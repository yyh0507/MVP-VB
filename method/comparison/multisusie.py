from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np
import MultiSuSiE


PIP_THRESHOLD = 0.9


def read_metadata(path: Path) -> dict[str, str]:
    with path.open("r", newline="", encoding="utf-8-sig") as handle:
        rows = list(csv.DictReader(handle))
    metadata = {row["field"]: row["value"] for row in rows}
    for field in ("K", "m", "q", "n_k"):
        if field not in metadata:
            raise ValueError(f"metadata.csv is missing {field}")
    return metadata


def read_snp_map(path: Path, p: int) -> list[dict[str, str]]:
    with path.open("r", newline="", encoding="utf-8-sig") as handle:
        rows = list(csv.DictReader(handle))
    if len(rows) != p:
        raise ValueError(f"Expected {p} SNP-map rows; observed {len(rows)}")
    if [int(row["SNP"]) for row in rows] != list(range(1, p + 1)):
        raise ValueError("snp_map.csv SNP indices are not consecutive")
    if len({row["variant_id"] for row in rows}) != p:
        raise ValueError("variant_id values must be unique")
    return rows


def load_and_validate(export_dir: Path):
    metadata = read_metadata(export_dir / "metadata.csv")
    k_count = int(metadata["K"])
    p = int(metadata["m"])
    q = int(metadata["q"])
    n_list = [int(value) for value in metadata["n_k"].split(",")]
    if k_count < 2 or p < 1 or q < 1 or len(n_list) != k_count:
        raise ValueError("Invalid K, m, q, or n_k in metadata.csv")

    snp_map = read_snp_map(export_dir / "snp_map.csv", p)
    x_list = []
    y_list = []
    truth_list = []
    for k, n in enumerate(n_list, start=1):
        x = np.loadtxt(export_dir / f"X_Group{k}.tsv", delimiter="\t")
        y = np.loadtxt(export_dir / f"Y_Group{k}.tsv", delimiter="\t")
        truth = np.loadtxt(
            export_dir / f"B_true_Group{k}.tsv", delimiter="\t"
        )
        x = np.asarray(x).reshape(n, p)
        y = np.asarray(y).reshape(n, q)
        truth = np.asarray(truth).reshape(p, q)
        if not (
            np.isfinite(x).all()
            and np.isfinite(y).all()
            and np.isfinite(truth).all()
        ):
            raise ValueError(f"Group{k} contains non-finite input values")
        if np.any((x < 0) | (x > 2) | (np.abs(x - np.rint(x)) > 1e-8)):
            raise ValueError(f"Group{k} X is not 0/1/2 allele dosage")
        if np.any(np.var(x, axis=0, ddof=1) <= 0):
            raise ValueError(f"Group{k} contains a monomorphic SNP")
        x_list.append(x)
        y_list.append(y)
        truth_list.append(truth)

    return metadata, snp_map, x_list, y_list, truth_list, n_list, p, q


def marginal_gwas(x: np.ndarray, y: np.ndarray):
    xc = x - x.mean(axis=0, keepdims=True)
    yc = y - y.mean()
    xtx = np.sum(xc * xc, axis=0)
    if np.any(xtx <= 0):
        raise ValueError("Monomorphic SNP encountered")
    beta_hat = (xc.T @ yc) / xtx
    residuals = yc[:, None] - xc * beta_hat[None, :]
    residual_sum_squares = np.sum(residuals**2, axis=0)
    se = np.sqrt(residual_sum_squares / ((x.shape[0] - 2) * xtx))
    if not np.isfinite(beta_hat).all() or not np.isfinite(se).all():
        raise ValueError("Non-finite marginal effect or standard error")
    return beta_hat, se


def write_csv(path: Path, rows: list[dict], fieldnames: list[str]):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def write_matrix_csv(path: Path, matrix: np.ndarray, trait_names: list[str]):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(trait_names)
        writer.writerows(matrix.tolist())


def finite_float_or_none(value):
    if value is None:
        return None
    value = float(value)
    return value if np.isfinite(value) else None


def indexed_value_or_none(values, index: int):
    if values is None or index >= len(values):
        return None
    return finite_float_or_none(values[index])


def extract_credible_sets(fit, variant_ids: list[str]):
    sets = getattr(fit, "sets", None)
    if sets is None or len(sets) == 0 or sets[0] is None:
        return [], np.zeros(len(variant_ids), dtype=bool)

    cs_indices = sets[0]
    purity = sets[1] if len(sets) > 1 else None
    claimed_coverage = sets[2] if len(sets) > 2 else None
    include_mask = sets[3] if len(sets) > 3 else None
    if include_mask is None:
        include_mask = np.ones(len(cs_indices), dtype=bool)
    include_mask = np.asarray(include_mask, dtype=bool)

    in_passed_cs = np.zeros(len(variant_ids), dtype=bool)
    rows = []
    for component, indices in enumerate(cs_indices):
        indices = np.asarray(indices, dtype=int)
        include = bool(include_mask[component]) if component < len(include_mask) else False
        if include:
            in_passed_cs[indices] = True
        rows.append(
            {
                "component": component + 1,
                "variant_indices_zero_based": indices.tolist(),
                "variant_ids": [variant_ids[i] for i in indices],
                "purity": indexed_value_or_none(purity, component),
                "claimed_coverage": indexed_value_or_none(
                    claimed_coverage, component
                ),
                "passes_purity_filter": include,
            }
        )
    return rows, in_passed_cs


def choose_mac_threshold(metadata: dict[str, str], n_list: list[int]) -> int:
    n_type = metadata.get("config.n_type", "").strip().lower()
    if n_type == "small":
        return 10
    if n_type == "big":
        return 20
    return 10 if min(n_list) < 100 else 20


def run_multisusie(export_dir: Path, result_dir: Path):
    result_dir.mkdir(parents=True, exist_ok=True)
    (
        metadata,
        snp_map,
        x_list,
        y_list,
        truth_list,
        n_list,
        p,
        q,
    ) = load_and_validate(export_dir)

    k_count = len(x_list)
    variant_ids = [row["variant_id"] for row in snp_map]
    trait_names = [f"Trait{r + 1}" for r in range(q)]
    r_list_master = [np.corrcoef(x, rowvar=False) for x in x_list]
    maf_list_master = [
        np.minimum(x.mean(axis=0) / 2, 1 - x.mean(axis=0) / 2)
        for x in x_list
    ]
    if not all(np.isfinite(r).all() for r in r_list_master):
        raise ValueError("Non-finite value in an LD matrix")

    rho = np.full((k_count, k_count), 0.75, dtype=np.float32)
    np.fill_diagonal(rho, 1.0)
    mac_threshold = choose_mac_threshold(metadata, n_list)
    shared_pip = np.zeros((p, q), dtype=np.float64)
    shared_zhat = np.zeros((p, q), dtype=np.int64)
    all_rows: list[dict] = []
    credible_sets: dict[str, list[dict]] = {}
    fit_diagnostics: list[dict] = []

    for trait_index in range(q):
        beta_hat_list = []
        se_list = []
        var_y_list = []
        for k in range(k_count):
            beta_hat, se = marginal_gwas(
                x_list[k], y_list[k][:, trait_index]
            )
            beta_hat_list.append(beta_hat.astype(np.float32, copy=True))
            se_list.append(se.astype(np.float32, copy=True))
            var_y_list.append(
                float(np.var(y_list[k][:, trait_index], ddof=1))
            )

        fit = MultiSuSiE.multisusie_rss(
            b_list=beta_hat_list,
            s_list=se_list,
            R_list=[r.astype(np.float32, copy=True) for r in r_list_master],
            varY_list=var_y_list,
            rho=rho,
            population_sizes=n_list,
            L=10,
            scaled_prior_variance=0.2,
            prior_weights=None,
            standardize=False,
            pop_spec_standardization=True,
            estimate_residual_variance=True,
            estimate_prior_variance=True,
            estimate_prior_method="early_EM",
            pop_spec_effect_priors=True,
            iter_before_zeroing_effects=5,
            prior_tol=1e-9,
            max_iter=100,
            tol=1e-3,
            verbose=False,
            coverage=0.95,
            min_abs_corr=0,
            float_type=np.float32,
            low_memory_mode=True,
            single_population_mac_thresh=mac_threshold,
            maf_list=[
                maf.astype(np.float32, copy=True) for maf in maf_list_master
            ],
            multi_population_maf_thresh=0.01,
            variant_ids=variant_ids,
        )

        if fit.pip.shape != (p,) or fit.coef.shape != (k_count, p):
            raise RuntimeError("Unexpected MultiSuSiE output dimensions")

        trait_name = trait_names[trait_index]
        cs_rows, in_passed_cs95 = extract_credible_sets(fit, variant_ids)
        credible_sets[trait_name] = cs_rows
        shared_pip[:, trait_index] = np.asarray(fit.pip, dtype=np.float64)
        shared_zhat[:, trait_index] = (
            shared_pip[:, trait_index] > PIP_THRESHOLD
        ).astype(np.int64)

        np.savez_compressed(
            result_dir / f"MultiSuSiE_{trait_name}_arrays.npz",
            pip=fit.pip,
            coef=fit.coef,
            coef_sd=fit.coef_sd,
            alpha=fit.alpha,
            V=fit.V,
            sigma2=fit.sigma2,
            elbo=fit.elbo,
        )
        fit_diagnostics.append(
            {
                "trait": trait_name,
                "converged": bool(fit.converged),
                "niter": int(fit.niter),
                "n_shared_PIP_gt_0_9": int(
                    np.sum(shared_zhat[:, trait_index])
                ),
                "n_passed_cs95_variants": int(np.sum(in_passed_cs95)),
            }
        )

        for k in range(k_count):
            for j in range(p):
                is_selected = bool(shared_zhat[j, trait_index])
                is_true = bool(truth_list[k][j, trait_index] != 0)
                all_rows.append(
                    {
                        "population": f"Group{k + 1}",
                        "SNP": j + 1,
                        "SNP_key": snp_map[j]["SNP_key"],
                        "variant_id": variant_ids[j],
                        "trait_index": trait_index + 1,
                        "trait": trait_name,
                        "overall_MultiSuSiE_PIP": float(fit.pip[j]),
                        "population_posterior_mean": float(fit.coef[k, j]),
                        "population_posterior_sd": float(fit.coef_sd[k, j]),
                        "in_passed_CS95": bool(in_passed_cs95[j]),
                        "shared_PIP_gt_0_9": is_selected,
                        "B_true": float(truth_list[k][j, trait_index]),
                        "is_true_population_pair": is_true,
                        "is_true_positive": bool(is_selected and is_true),
                    }
                )

    fields = [
        "population",
        "SNP",
        "SNP_key",
        "variant_id",
        "trait_index",
        "trait",
        "overall_MultiSuSiE_PIP",
        "population_posterior_mean",
        "population_posterior_sd",
        "in_passed_CS95",
        "shared_PIP_gt_0_9",
        "B_true",
        "is_true_population_pair",
        "is_true_positive",
    ]
    write_csv(
        result_dir / "MultiSuSiE_all_population_SNP_trait_pairs.csv",
        all_rows,
        fields,
    )
    selected_rows = [row for row in all_rows if row["shared_PIP_gt_0_9"]]
    write_csv(
        result_dir / "MultiSuSiE_shared_PIP_gt_0.9_population_pairs.csv",
        selected_rows,
        fields,
    )

    write_matrix_csv(
        result_dir / "MultiSuSiE_overall_PIP_matrix.csv",
        shared_pip,
        trait_names,
    )
    write_matrix_csv(
        result_dir / "MultiSuSiE_shared_Zhat_matrix.csv",
        shared_zhat,
        trait_names,
    )
    np.save(result_dir / "MultiSuSiE_overall_PIP_matrix.npy", shared_pip)
    np.save(result_dir / "MultiSuSiE_shared_Zhat_matrix.npy", shared_zhat)

    for k in range(k_count):
        write_matrix_csv(
            result_dir / f"MultiSuSiE_Zhat_Group{k + 1}.csv",
            shared_zhat,
            trait_names,
        )

    metric_rows = []
    for k in range(k_count):
        truth = truth_list[k] != 0
        pred = shared_zhat.astype(bool)
        tp = int(np.sum(pred & truth))
        fp = int(np.sum(pred & ~truth))
        fn = int(np.sum(~pred & truth))
        tn = int(np.sum(~pred & ~truth))
        precision = tp / (tp + fp) if tp + fp else 0.0
        recall = tp / (tp + fn) if tp + fn else 0.0
        f1 = (
            2 * precision * recall / (precision + recall)
            if precision + recall > 0
            else 0.0
        )
        metric_rows.append(
            {
                "method": "MultiSuSiE",
                "selection_rule": (
                    "overall_MultiSuSiE_PIP_gt_0.9_"
                    "replicated_across_populations"
                ),
                "population": f"Group{k + 1}",
                "TP": tp,
                "FP": fp,
                "FN": fn,
                "TN": tn,
                "precision": precision,
                "recall": recall,
                "F1": f1,
            }
        )
    write_csv(
        result_dir / "MultiSuSiE_population_pair_metrics.csv",
        metric_rows,
        list(metric_rows[0].keys()),
    )

    with (result_dir / "MultiSuSiE_credible_sets.json").open(
        "w", encoding="utf-8"
    ) as handle:
        json.dump(credible_sets, handle, indent=2, allow_nan=False)
    with (result_dir / "MultiSuSiE_fit_diagnostics.json").open(
        "w", encoding="utf-8"
    ) as handle:
        json.dump(fit_diagnostics, handle, indent=2, allow_nan=False)

    parameter_record = {
        "all_populations_fitted_jointly_for_each_trait": True,
        "population_specific_PIP": False,
        "selection_replication": "same overall-PIP selection in every population",
        "rho": 0.75,
        "L": 10,
        "scaled_prior_variance": 0.2,
        "prior_weights": None,
        "standardize": False,
        "pop_spec_standardization": True,
        "estimate_residual_variance": True,
        "estimate_prior_variance": True,
        "estimate_prior_method": "early_EM",
        "pop_spec_effect_priors": True,
        "iter_before_zeroing_effects": 5,
        "prior_tol": 1e-9,
        "max_iter": 100,
        "tol": 1e-3,
        "coverage": 0.95,
        "min_abs_corr": 0,
        "float_type": "numpy.float32",
        "low_memory_mode": True,
        "single_population_mac_thresh": mac_threshold,
        "single_population_mac_rule": (
            "10 for config.n_type=small or minimum n<100; otherwise 20"
        ),
        "multi_population_maf_thresh": 0.01,
        "fine_mapping_threshold": "overall PIP > 0.9",
    }
    with (result_dir / "MultiSuSiE_parameters.json").open(
        "w", encoding="utf-8"
    ) as handle:
        json.dump(parameter_record, handle, indent=2)

    print(f"MultiSuSiE completed: {result_dir}")
    print(f"K={k_count}, P={p}, Q={q}")
    print(f"Population-SNP-trait rows: {len(all_rows)}")
    print(f"Shared PIP > 0.9 rows per population: {int(shared_zhat.sum())}")
