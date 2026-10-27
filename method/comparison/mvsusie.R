library(mvsusieR)

## ============================================================
## Legacy manual pair-level PIP extractor retained only for comparison.
## The main reproduction path below uses the official exported function:
## mvsusieR::mvsusie_get_pip_per_outcome().
##
## Pair-level PIP:
## PIP_{j,r} = 1 - prod_l {1 - alpha_{l,j,r}}
## ============================================================

legacy_safe_mvsusie_get_pair_pip <- function(fit_k, q = NULL, verbose = TRUE) {
  if (is.null(fit_k$alpha)) {
    stop("fit_k$alpha is NULL. Invalid mvSuSiE fit.")
  }

  J <- ncol(fit_k$alpha)
  L <- nrow(fit_k$alpha)

  if (is.null(q)) {
    if (!is.null(fit_k$sigma2) && is.matrix(fit_k$sigma2)) {
      q <- ncol(fit_k$sigma2)
    } else if (!is.null(fit_k$lfsr) && is.matrix(fit_k$lfsr)) {
      q <- ncol(fit_k$lfsr)
    } else {
      stop("Cannot determine q in safe_mvsusie_get_pair_pip().")
    }
  }

  pmw <- fit_k$posterior_mixture_weights

  if (verbose) {
    cat("posterior_mixture_weights class:", class(pmw), "\n")
    cat(
      "posterior_mixture_weights dim:",
      ifelse(is.null(dim(pmw)), "NULL", paste(dim(pmw), collapse = " x ")),
      "\n"
    )
    if (!is.null(fit_k$V)) {
      cat("V summary:\n")
      print(summary(as.numeric(fit_k$V)))
    }
  }

  ## ------------------------------------------------------------
  ## 如果所有 single effects 都是 null，直接返回 0 矩阵
  ## ------------------------------------------------------------
  if (!is.null(fit_k$V)) {
    V_num <- suppressWarnings(as.numeric(fit_k$V))
    if (length(V_num) > 0 && all(is.finite(V_num)) && max(abs(V_num)) <= 1e-12) {
      if (verbose) {
        message("All mvSuSiE single-effect prior variances are zero. Return zero pair-level PIP matrix.")
      }
      return(matrix(0, nrow = J, ncol = q))
    }
  }

  ## ------------------------------------------------------------
  ## posterior_mixture_weights 必须是 L x J x K_prior 的 3D array
  ## 否则说明这个 fit 没有可用的 outcome-specific posterior mixture 信息
  ## ------------------------------------------------------------
  if (!is.array(pmw) || length(dim(pmw)) != 3) {
    if (verbose) {
      message("posterior_mixture_weights is not a valid 3D array. Return zero pair-level PIP matrix.")
    }
    return(matrix(0, nrow = J, ncol = q))
  }

  if (dim(pmw)[1] != L || dim(pmw)[2] != J) {
    if (verbose) {
      message("posterior_mixture_weights dimensions do not match alpha. Return zero pair-level PIP matrix.")
    }
    return(matrix(0, nrow = J, ncol = q))
  }

  active_matrices <- fit_k$V_structure

  if (is.null(active_matrices) || !is.list(active_matrices)) {
    if (verbose) {
      message("fit_k$V_structure is NULL or invalid. Return zero pair-level PIP matrix.")
    }
    return(matrix(0, nrow = J, ncol = q))
  }

  outcome_indicator <- do.call(
    rbind,
    lapply(active_matrices, function(U) {
      if (!is.matrix(U)) {
        return(rep(0L, q))
      }
      as.integer(diag(U) != 0)
    })
  )

  R <- ncol(outcome_indicator)
  K_prior <- nrow(outcome_indicator)

  if (R != q) {
    if (verbose) {
      message("Outcome dimension in V_structure does not match q. Return zero pair-level PIP matrix.")
    }
    return(matrix(0, nrow = J, ncol = q))
  }

  pmw_K <- dim(pmw)[3]

  ## 有些对象会保留 null component，有些不会。
  ## 如果 pmw_K = K_prior + 1，第一列是 null，非 null 从 2 开始。
  ## 如果 pmw_K = K_prior，说明 null 已经被去掉，非 null 从 1 开始。
  if (pmw_K == K_prior + 1) {
    offset <- 1L
  } else if (pmw_K == K_prior) {
    offset <- 0L
  } else {
    if (verbose) {
      message(
        "posterior_mixture_weights third dimension does not match V_structure. ",
        "pmw_K = ", pmw_K, ", K_prior = ", K_prior,
        ". Return zero pair-level PIP matrix."
      )
    }
    return(matrix(0, nrow = J, ncol = q))
  }

  alpha_out <- array(0, dim = c(L, J, R))

  for (r in seq_len(R)) {
    tmp <- matrix(0, nrow = L, ncol = J)

    for (kk in seq_len(K_prior)) {
      col_idx <- kk + offset
      tmp <- tmp + pmw[, , col_idx] * outcome_indicator[kk, r]
    }

    alpha_out[, , r] <- tmp * fit_k$alpha
  }

  pair_pip <- matrix(0, nrow = J, ncol = R)

  for (r in seq_len(R)) {
    A <- matrix(alpha_out[, , r], nrow = L, ncol = J)
    pair_pip[, r] <- 1 - apply(1 - A, 2, prod)
  }

  pair_pip[!is.finite(pair_pip)] <- 0
  pair_pip[pair_pip < 0] <- 0
  pair_pip[pair_pip > 1] <- 1

  colnames(pair_pip) <- paste0("Pheno", seq_len(R))
  rownames(pair_pip) <- paste0("SNP", seq_len(J))

  pair_pip
}


## ============================================================
## AUC by rank
## ============================================================

auc_rank <- function(score, truth) {
  score <- as.numeric(score)
  truth <- as.integer(as.numeric(truth) != 0)

  ok <- is.finite(score) & !is.na(truth)
  score <- score[ok]
  truth <- truth[ok]

  n1 <- sum(truth == 1)
  n0 <- sum(truth == 0)

  if (n1 == 0 || n0 == 0) return(NA_real_)

  r <- rank(score, ties.method = "average")
  auc <- (sum(r[truth == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
  auc
}


## ============================================================
## Pair-level metrics based on local false sign rate (lfsr)
## ============================================================

eval_pair_lfsr <- function(lfsr, truth, threshold = 0.01) {
  lfsr_vec <- as.numeric(lfsr)
  truth_vec <- as.integer(as.numeric(truth) != 0)

  if (length(lfsr_vec) != length(truth_vec)) {
    stop("lfsr and truth have different lengths.")
  }

  ## mvSuSiE source returns a J x R per-variable lfsr matrix.
  ## A SNP-phenotype pair is discovered only when lfsr < 0.01.
  pred_vec <- as.integer(lfsr_vec < threshold)

  TP <- sum(pred_vec == 1 & truth_vec == 1, na.rm = TRUE)
  FP <- sum(pred_vec == 1 & truth_vec == 0, na.rm = TRUE)
  FN <- sum(pred_vec == 0 & truth_vec == 1, na.rm = TRUE)
  TN <- sum(pred_vec == 0 & truth_vec == 0, na.rm = TRUE)

  precision <- ifelse(TP + FP > 0, TP / (TP + FP), 0)
  recall <- ifelse(TP + FN > 0, TP / (TP + FN), 0)
  specificity <- ifelse(TN + FP > 0, TN / (TN + FP), 0)
  youden <- recall + specificity - 1

  ## auc_rank expects larger scores to mean stronger evidence.
  auc <- auc_rank(1 - lfsr_vec, truth_vec)

  f1 <- ifelse(
    2 * TP + FP + FN > 0,
    2 * TP / (2 * TP + FP + FN),
    0
  )

  data.frame(
    TP = TP,
    FP = FP,
    FN = FN,
    TN = TN,
    Precision = precision,
    Recall = recall,
    Specificity = specificity,
    Youden = youden,
    F1 = f1,
    AUC = auc,
    threshold = threshold,
    n_selected = TP + FP,
    n_true = TP + FN
  )
}


## ============================================================
## Error inspection
## ============================================================

inspect_pair_lfsr_errors <- function(lfsr, truth_B, threshold = 0.01,
                                     population = 1, max_print = Inf) {
  truth <- truth_B != 0
  pred <- lfsr < threshold

  metric <- eval_pair_lfsr(lfsr, truth, threshold = threshold)

  cat(sprintf("\n--- Population %d mvSuSiE pair-level error inspection ---\n",
              population))
  cat(sprintf(
    "TP=%d, FP=%d, FN=%d, TN=%d\n",
    metric$TP, metric$FP, metric$FN, metric$TN
  ))
  cat(sprintf(
    "Precision=%.3f, Recall=%.3f, Specificity=%.3f, Youden=%.3f, AUC=%.3f, lfsr_threshold=%.6g\n",
    metric$Precision, metric$Recall, metric$Specificity,
    metric$Youden, metric$AUC, threshold
  ))

  err <- which(pred != truth, arr.ind = TRUE)
  if (nrow(err) == 0) {
    cat("No classification errors at this threshold.\n")
    return(invisible(metric))
  }

  n_show <- min(nrow(err), max_print)

  for (ii in seq_len(n_show)) {
    i <- err[ii, 1]
    r <- err[ii, 2]
    cat(sprintf(
      "Position (SNP=%d, Pheno=%d): B_true=%.4f, lfsr=%.6g, Predicted=%d, True=%d\n",
      i, r, truth_B[i, r], lfsr[i, r],
      as.integer(pred[i, r]), as.integer(truth[i, r])
    ))
  }

  if (n_show < nrow(err)) {
    cat(sprintf("... %d more errors not printed.\n", nrow(err) - n_show))
  }

  invisible(metric)
}


## ============================================================
## Author comparison setting: use the fixed upper bound L = 10.
## Do not tune L from target_causal_prop or any other truth-generating
## simulation parameter.
## ============================================================

choose_mvsusie_L <- function(common_data,
                             L_fixed = 10L) {
  m <- ncol(common_data$X_list[[1]])
  L <- as.integer(L_fixed)
  if (length(L) != 1L || is.na(L) || L < 1L || L > m) {
    stop("L_fixed must be one integer between 1 and the number of SNPs.")
  }
  L
}


## ============================================================
## Author-style trait-wise fine-mapped SNP-trait pairs
##
## A pair is reported when:
##   1. the SNP belongs to a retained 95% credible set;
##   2. that credible set has minimum absolute correlation >= 0.5; and
##   3. the corresponding single effect has average lfsr < 0.01 for
##      the trait.
## ============================================================

extract_author_snp_trait_pairs <- function(fit, population,
                                           lfsr_threshold = 0.01,
                                           min_abs_corr = 0.5) {
  if (is.null(fit$lfsr) || !is.matrix(fit$lfsr)) {
    stop("mvSuSiE fit does not contain a J x R lfsr matrix.")
  }
  if (is.null(fit$single_effect_lfsr) ||
      !is.matrix(fit$single_effect_lfsr)) {
    stop("mvSuSiE fit does not contain an L x R single_effect_lfsr matrix.")
  }

  J <- nrow(fit$lfsr)
  R <- ncol(fit$lfsr)
  if (ncol(fit$single_effect_lfsr) != R) {
    stop("single_effect_lfsr and lfsr have different trait dimensions.")
  }

  selected <- matrix(0L, nrow = J, ncol = R)
  rownames(selected) <- rownames(fit$lfsr)
  colnames(selected) <- colnames(fit$lfsr)

  empty_pairs <- data.frame(
    population = integer(0),
    snp_index = integer(0),
    snp = character(0),
    trait_index = integer(0),
    trait = character(0),
    credible_set = character(0),
    effect_index = integer(0),
    global_pip = numeric(0),
    snp_lfsr = numeric(0),
    average_lfsr = numeric(0),
    cs_min_abs_corr = numeric(0),
    stringsAsFactors = FALSE
  )

  if (is.null(fit$sets) || is.null(fit$sets$cs) ||
      length(fit$sets$cs) == 0L) {
    return(list(selected = selected, pairs = empty_pairs))
  }

  cs_list <- fit$sets$cs
  cs_names <- names(cs_list)
  if (is.null(cs_names)) cs_names <- paste0("L", seq_along(cs_list))

  cs_index <- fit$sets$cs_index
  if (is.null(cs_index) || length(cs_index) != length(cs_list)) {
    cs_index <- suppressWarnings(as.integer(sub("^L", "", cs_names)))
  }
  if (anyNA(cs_index)) {
    stop("Cannot map credible sets to mvSuSiE single-effect indices.")
  }

  purity <- fit$sets$purity
  if (is.null(purity) || !("min.abs.corr" %in% names(purity))) {
    stop("mvSuSiE fit does not contain credible-set purity information.")
  }
  if (!is.null(rownames(purity)) && all(cs_names %in% rownames(purity))) {
    cs_purity <- as.numeric(purity[cs_names, "min.abs.corr"])
  } else {
    if (nrow(purity) != length(cs_list)) {
      stop("Cannot align credible sets with purity rows.")
    }
    cs_purity <- as.numeric(purity$min.abs.corr)
  }

  snp_names <- rownames(fit$lfsr)
  if (is.null(snp_names)) snp_names <- paste0("SNP", seq_len(J))
  trait_names <- colnames(fit$lfsr)
  if (is.null(trait_names)) trait_names <- paste0("Pheno", seq_len(R))

  pair_rows <- vector("list", 0L)
  for (i in seq_along(cs_list)) {
    if (!is.finite(cs_purity[i]) || cs_purity[i] < min_abs_corr) next

    effect_index <- cs_index[i]
    if (effect_index < 1L || effect_index > nrow(fit$single_effect_lfsr)) {
      stop("Credible-set single-effect index is outside single_effect_lfsr.")
    }
    snp_index <- as.integer(cs_list[[i]])
    snp_index <- snp_index[snp_index >= 1L & snp_index <= J]
    trait_index <- which(
      fit$single_effect_lfsr[effect_index, ] < lfsr_threshold
    )
    if (length(snp_index) == 0L || length(trait_index) == 0L) next

    selected[snp_index, trait_index] <- 1L
    grid <- expand.grid(
      snp_index = snp_index,
      trait_index = trait_index,
      KEEP.OUT.ATTRS = FALSE,
      stringsAsFactors = FALSE
    )
    pair_rows[[length(pair_rows) + 1L]] <- data.frame(
      population = as.integer(population),
      snp_index = grid$snp_index,
      snp = snp_names[grid$snp_index],
      trait_index = grid$trait_index,
      trait = trait_names[grid$trait_index],
      credible_set = cs_names[i],
      effect_index = effect_index,
      global_pip = fit$pip[grid$snp_index],
      snp_lfsr = fit$lfsr[cbind(grid$snp_index, grid$trait_index)],
      average_lfsr = fit$single_effect_lfsr[
        effect_index, grid$trait_index
      ],
      cs_min_abs_corr = cs_purity[i],
      stringsAsFactors = FALSE
    )
  }

  pairs <- if (length(pair_rows) == 0L) {
    empty_pairs
  } else {
    pairs <- do.call(rbind, pair_rows)
    pairs[!duplicated(pairs[c(
      "population", "snp_index", "trait_index"
    )]), , drop = FALSE]
  }
  list(selected = selected, pairs = pairs)
}


eval_pair_selection <- function(selected, truth, score = NULL,
                                threshold = NA_real_) {
  pred_vec <- as.integer(as.logical(selected))
  truth_vec <- as.integer(as.numeric(truth) != 0)
  if (length(pred_vec) != length(truth_vec)) {
    stop("selected and truth have different lengths.")
  }

  TP <- sum(pred_vec == 1 & truth_vec == 1, na.rm = TRUE)
  FP <- sum(pred_vec == 1 & truth_vec == 0, na.rm = TRUE)
  FN <- sum(pred_vec == 0 & truth_vec == 1, na.rm = TRUE)
  TN <- sum(pred_vec == 0 & truth_vec == 0, na.rm = TRUE)
  precision <- ifelse(TP + FP > 0, TP / (TP + FP), 0)
  recall <- ifelse(TP + FN > 0, TP / (TP + FN), 0)
  specificity <- ifelse(TN + FP > 0, TN / (TN + FP), 0)
  youden <- recall + specificity - 1
  f1 <- ifelse(
    2 * TP + FP + FN > 0,
    2 * TP / (2 * TP + FP + FN), 0
  )
  auc <- if (is.null(score)) NA_real_ else auc_rank(score, truth_vec)

  data.frame(
    TP = TP, FP = FP, FN = FN, TN = TN,
    Precision = precision, Recall = recall,
    Specificity = specificity, Youden = youden,
    F1 = f1, AUC = auc, threshold = threshold,
    n_selected = TP + FP, n_true = TP + FN
  )
}


inspect_pair_selection_errors <- function(selected, truth_B,
                                          population = 1,
                                          max_print = Inf) {
  truth <- truth_B != 0
  metric <- eval_pair_selection(selected = selected, truth = truth)

  cat(sprintf(
    "\n--- Population %d author-style CS SNP-trait error inspection ---\n",
    population
  ))
  cat(sprintf(
    "TP=%d, FP=%d, FN=%d, TN=%d\n",
    metric$TP, metric$FP, metric$FN, metric$TN
  ))
  cat(sprintf(
    "Precision=%.3f, Recall=%.3f, Specificity=%.3f, Youden=%.3f\n",
    metric$Precision, metric$Recall, metric$Specificity, metric$Youden
  ))

  err <- which(selected != truth, arr.ind = TRUE)
  if (nrow(err) == 0L) {
    cat("No classification errors under the author-style CS rule.\n")
    return(invisible(metric))
  }

  n_show <- min(nrow(err), max_print)
  for (ii in seq_len(n_show)) {
    j <- err[ii, 1]
    r <- err[ii, 2]
    cat(sprintf(
      "Position (SNP=%d, Pheno=%d): B_true=%.4f, Predicted=%d, True=%d\n",
      j, r, truth_B[j, r], as.integer(selected[j, r]),
      as.integer(truth[j, r])
    ))
  }
  if (n_show < nrow(err)) {
    cat(sprintf("... %d more errors not printed.\n", nrow(err) - n_show))
  }
  invisible(metric)
}


## ============================================================
## Main function: run mvSuSiE on one common_data rds
##
## Fit stage:
##   only uses X_list and Y_list
##
## Evaluation stage:
##   uses B_true_list
## ============================================================

run_mvsusie_on_common_data <- function(common_data_file,
                                       out_file = NULL,
                                       csv_file = NULL,
                                       threshold = 0.01,
                                       eval_thresholds = c(0.01),
                                       L_fixed = 10L,
                                       seed = 2026,
                                       max_iter = 1000,
                                       tol = 1e-4,
                                       n_thread = 4,
                                       standardize = TRUE,
                                       intercept = TRUE,
                                       estimate_residual_variance = TRUE,
                                       estimate_prior_variance = TRUE,
                                       estimate_prior_method = "optim",
                                       estimate_prior_mixture_weights = TRUE,
                                       mixture_weight_method = "mixsqp",
                                       precompute_cache = TRUE,
                                       coverage = 0.95,
                                       min_abs_corr = 0.5,
                                       require_convergence = TRUE,
                                       verbose = TRUE,
                                       print_errors = TRUE,
                                       max_error_print = Inf) {
  if (
    length(threshold) != 1L ||
    !is.numeric(threshold) ||
    !is.finite(threshold) ||
    threshold < 0 ||
    threshold > 1
  ) {
    stop("threshold must be one finite numeric value between 0 and 1.")
  }
  eval_thresholds <- as.numeric(eval_thresholds)
  if (
    length(eval_thresholds) == 0L ||
    any(!is.finite(eval_thresholds)) ||
    any(eval_thresholds < 0 | eval_thresholds > 1)
  ) {
    stop("eval_thresholds must contain finite numeric values between 0 and 1.")
  }

  installed_version <- utils::packageVersion("mvsusieR")
  if (installed_version < "0.3.0") {
    stop(
      "This runner requires the current mvsusieR implementation (>= 0.3.0); ",
      "installed version is ", as.character(installed_version), "."
    )
  }

  common_data <- readRDS(common_data_file)

  scenario_name <- common_data$scenario_name
  if (is.null(scenario_name)) {
    scenario_name <- tools::file_path_sans_ext(basename(common_data_file))
    scenario_name <- sub("^common_data_", "", scenario_name)
  }

  X_list <- common_data$X_list
  Y_list <- common_data$Y_list

  ## B_true_list 只在评价阶段使用
  B_true_list <- common_data$B_true_list
  K <- length(X_list)
  q <- ncol(Y_list[[1]])
  m <- ncol(X_list[[1]])

  L <- choose_mvsusie_L(
    common_data = common_data,
    L_fixed = L_fixed
  )

  cat("\n============================================================\n")
  cat("Running mvSuSiE on common_data\n")
  cat("Scenario:", scenario_name, "\n")
  cat("K =", K, "| m =", m, "| q =", q, "| L =", L, "\n")
  cat("Input file:", common_data_file, "\n")
  cat("============================================================\n")

  fit_list <- vector("list", K)
  selected_pair_list <- vector("list", K)
  pair_table_list <- vector("list", K)
  lfsr_list <- vector("list", K)
  metric_list <- vector("list", K)
  error_metric_list <- vector("list", K)
  selection_status_list <- vector("list", K)

  prior <- mvsusieR::create_mixture_prior(R = q)

  for (k in seq_len(K)) {
    cat("\n------------------------------------------------------------\n")
    cat("mvSuSiE population", k, "of", K, "\n")
    cat("------------------------------------------------------------\n")

    set.seed(seed + k)

    ## ----------------------------
    ## Fit stage: only X and Y
    ## ----------------------------
    Xk <- as.matrix(X_list[[k]])
    Yk <- as.matrix(Y_list[[k]])

    if (!is.null(common_data$snp_key_vec) &&
        length(common_data$snp_key_vec) == ncol(Xk)) {
      colnames(Xk) <- as.character(common_data$snp_key_vec)
    } else if (is.null(colnames(Xk))) {
      colnames(Xk) <- paste0("SNP", seq_len(ncol(Xk)))
    }
    if (is.null(colnames(Yk))) {
      colnames(Yk) <- paste0("Pheno", seq_len(ncol(Yk)))
    }

    if (any(!is.finite(Xk)) || any(!is.finite(Yk))) {
      stop("Non-finite X or Y value in population ", k, ".")
    }

    fit_k <- mvsusieR::mvsusie(
      X = Xk,
      Y = Yk,
      L = L,
      prior_variance = prior,
      standardize = standardize,
      intercept = intercept,
      estimate_residual_variance = estimate_residual_variance,
      estimate_prior_variance = estimate_prior_variance,
      estimate_prior_method = estimate_prior_method,
      estimate_prior_mixture_weights = estimate_prior_mixture_weights,
      mixture_weight_method = mixture_weight_method,
      precompute_cache = precompute_cache,
      coverage = coverage,
      min_abs_corr = min_abs_corr,
      max_iter = max_iter,
      tol = tol,
      n_thread = n_thread,
      verbose = verbose
    )

    if (!isTRUE(fit_k$converged)) {
      msg <- paste0(
        "mvSuSiE did not converge for population ", k,
        " after ", fit_k$niter, " iterations."
      )
      if (require_convergence) stop(msg) else warning(msg)
    }

    ## ----------------------------
    ## Official mvSuSiE per-variable, per-outcome lfsr output.
    ## This is retained as an auxiliary SNP-wise result. The primary
    ## author-style result below additionally requires a retained CS and
    ## uses the corresponding single-effect (average) lfsr.
    ## ----------------------------
    expected_dim <- c(ncol(Xk), ncol(Yk))
    original_lfsr <- fit_k$lfsr
    original_lfsr_dim <- if (is.null(original_lfsr)) {
      NULL
    } else {
      dim(as.matrix(original_lfsr))
    }
    degenerate_lfsr <- is.null(original_lfsr) ||
      !identical(as.integer(original_lfsr_dim), as.integer(expected_dim))

    if (degenerate_lfsr) {
      status_message <- paste0(
        "mvSuSiE returned a degenerate lfsr object with dim ",
        if (is.null(original_lfsr_dim)) "NULL" else {
          paste(original_lfsr_dim, collapse = "x")
        },
        "; expected ", paste(expected_dim, collapse = "x"),
        ". This population is retained as a no-signal fit with an ",
        "all-zero SNP-trait selection matrix."
      )
      warning(status_message, call. = FALSE)

      ## Keep the original package return for diagnosis. The normalized
      ## fields below make all downstream thresholds produce zero selections.
      fit_k$original_degenerate_lfsr <- original_lfsr
      fit_k$lfsr <- matrix(
        1,
        nrow = expected_dim[1L],
        ncol = expected_dim[2L],
        dimnames = list(colnames(Xk), colnames(Yk))
      )
      fit_k$pip <- stats::setNames(rep(0, expected_dim[1L]), colnames(Xk))
      fit_k$single_effect_lfsr <- matrix(
        1,
        nrow = L,
        ncol = expected_dim[2L],
        dimnames = list(paste0("L", seq_len(L)), colnames(Yk))
      )
      fit_k$sets <- list(cs = list())
      fit_k$selection_status <- list(
        code = "degenerate_lfsr_no_signal",
        message = status_message,
        original_lfsr_dim = original_lfsr_dim,
        expected_lfsr_dim = expected_dim,
        forced_zero_selection = TRUE
      )
    } else {
      fit_k$selection_status <- list(
        code = "ok",
        message = "mvSuSiE returned an lfsr matrix with the expected dimensions.",
        original_lfsr_dim = original_lfsr_dim,
        expected_lfsr_dim = expected_dim,
        forced_zero_selection = FALSE
      )
    }
    selection_status_list[[k]] <- fit_k$selection_status

    lfsr_k <- as.matrix(fit_k$lfsr)
    lfsr_k[!is.finite(lfsr_k)] <- 1
    lfsr_k <- pmin(pmax(lfsr_k, 0), 1)

    if (!all(dim(lfsr_k) == dim(B_true_list[[k]]))) {
      stop(
        "Dimension mismatch in population ", k,
        ": lfsr has dim ",
        paste(dim(lfsr_k), collapse = "x"),
        ", B_true has dim ",
        paste(dim(B_true_list[[k]]), collapse = "x")
      )
    }

    author_pairs <- extract_author_snp_trait_pairs(
      fit = fit_k,
      population = k,
      lfsr_threshold = threshold,
      min_abs_corr = min_abs_corr
    )

    fit_list[[k]] <- fit_k
    lfsr_list[[k]] <- lfsr_k
    selected_pair_list[[k]] <- author_pairs$selected
    pair_table_list[[k]] <- author_pairs$pairs

    ## ----------------------------
    ## Evaluation stage: truth allowed
    ## ----------------------------
    truth_k <- B_true_list[[k]] != 0

    fixed_metrics <- do.call(rbind, lapply(eval_thresholds, function(th) {
      pairs_th <- extract_author_snp_trait_pairs(
        fit = fit_k,
        population = k,
        lfsr_threshold = th,
        min_abs_corr = min_abs_corr
      )
      ans <- eval_pair_selection(
        selected = pairs_th$selected,
        truth = truth_k,
        score = 1 - lfsr_k,
        threshold = th
      )

      cbind(
        method = "mvSuSiE_95CS_single_effect_lfsr01",
        scenario_name = scenario_name,
        population = k,
        L = L,
        threshold_type = "fixed",
        ans
      )
    }))

    metric_list[[k]] <- fixed_metrics

    if (print_errors) {
      error_metric_list[[k]] <- inspect_pair_selection_errors(
        selected = selected_pair_list[[k]],
        truth_B = B_true_list[[k]],
        population = k,
        max_print = max_error_print
      )
    }
  }

  metrics_by_population <- do.call(rbind, metric_list)

  ## ------------------------------------------------------------
  ## Overall metric across all populations
  ## 仅作补充，因为不同 population 的 PIP calibration 可能不同
  ## ------------------------------------------------------------
  lfsr_all <- unlist(lfsr_list)
  truth_all <- unlist(lapply(B_true_list, function(B) B != 0))

  overall_fixed_metrics <- do.call(rbind, lapply(eval_thresholds, function(th) {
    selected_th <- unlist(lapply(seq_len(K), function(k) {
      extract_author_snp_trait_pairs(
        fit = fit_list[[k]], population = k,
        lfsr_threshold = th, min_abs_corr = min_abs_corr
      )$selected
    }))
    ans <- eval_pair_selection(
      selected = selected_th,
      truth = truth_all,
      score = 1 - lfsr_all,
      threshold = th
    )

    cbind(
      method = "mvSuSiE_95CS_single_effect_lfsr01",
      scenario_name = scenario_name,
      population = "overall",
      L = L,
      threshold_type = "fixed",
      ans
    )
  }))

  metrics_overall <- overall_fixed_metrics

  out <- list(
    method = "mvSuSiE_95CS_single_effect_lfsr01",
    scenario_name = scenario_name,
    common_data_file = common_data_file,
    config = common_data$config,
    K = K,
    m = m,
    q = q,
    n_k = common_data$n_k,
    L = L,
    selection_rule = paste0(
      "SNP in retained ", coverage * 100,
      "% CS with min.abs.corr >= ", min_abs_corr,
      " and fit$single_effect_lfsr < threshold"
    ),
    selection_threshold = threshold,
    fit_settings = list(
      package_version = as.character(utils::packageVersion("mvsusieR")),
      L = L,
      max_iter = max_iter,
      tol = tol,
      n_thread = n_thread,
      standardize = standardize,
      intercept = intercept,
      residual_variance = "estimated from X and Y by mvsusieR",
      estimate_residual_variance = estimate_residual_variance,
      estimate_prior_variance = estimate_prior_variance,
      estimate_prior_method = estimate_prior_method,
      estimate_prior_mixture_weights = estimate_prior_mixture_weights,
      mixture_weight_method = mixture_weight_method,
      precompute_cache = precompute_cache,
      coverage = coverage,
      min_abs_corr = min_abs_corr,
      prior = "canonical mixture from create_mixture_prior(R = q)"
    ),
    fit_list = fit_list,
    Zhat_list = selected_pair_list,
    selected_pair_list = selected_pair_list,
    snp_trait_pairs = do.call(rbind, pair_table_list),
    pair_counts_by_population = data.frame(
      population = seq_len(K),
      n_snp_trait_pairs = vapply(
        selected_pair_list, function(z) sum(z == 1L), integer(1)
      ),
      stringsAsFactors = FALSE
    ),
    lfsr_list = lfsr_list,
    selection_status_list = selection_status_list,
    metrics_by_population = metrics_by_population,
    metrics_overall = metrics_overall,
    h2_table = common_data$h2_table,
    causal_info = common_data$causal_info
  )

  if (!is.null(out_file)) {
    dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
    saveRDS(out, out_file)
    cat("\nSaved mvSuSiE output:", out_file, "\n")
  }

  if (!is.null(csv_file)) {
    dir.create(dirname(csv_file), recursive = TRUE, showWarnings = FALSE)
    write.csv(metrics_by_population, csv_file, row.names = FALSE)
    cat("Saved mvSuSiE metrics csv:", csv_file, "\n")

    csv_overall_file <- sub("\\.csv$", "_overall.csv", csv_file)
    write.csv(metrics_overall, csv_overall_file, row.names = FALSE)
    cat("Saved mvSuSiE overall csv:", csv_overall_file, "\n")

    csv_pairs_file <- sub("\\.csv$", "_snp_trait_pairs.csv", csv_file)
    write.csv(out$snp_trait_pairs, csv_pairs_file, row.names = FALSE)
    cat("Saved mvSuSiE SNP-trait pairs csv:", csv_pairs_file, "\n")

    csv_pair_counts_file <- sub(
      "\\.csv$", "_snp_trait_pair_counts.csv", csv_file
    )
    write.csv(
      out$pair_counts_by_population,
      csv_pair_counts_file,
      row.names = FALSE
    )
    cat("Saved mvSuSiE pair counts csv:", csv_pair_counts_file, "\n")

    zhat_list_file <- sub("\\.csv$", "_Zhat_list.rds", csv_file)
    saveRDS(selected_pair_list, zhat_list_file)
    cat("Saved mvSuSiE Zhat list RDS:", zhat_list_file, "\n")

    for (k in seq_len(K)) {
      csv_pairs_k <- sub(
        "\\.csv$",
        paste0("_population", k, "_snp_trait_pairs.csv"),
        csv_file
      )
      write.csv(pair_table_list[[k]], csv_pairs_k, row.names = FALSE)
      cat("Saved population", k, "SNP-trait pairs csv:", csv_pairs_k, "\n")

      zhat_rds_k <- sub(
        "\\.csv$", paste0("_population", k, "_Zhat.rds"), csv_file
      )
      saveRDS(selected_pair_list[[k]], zhat_rds_k)
      cat("Saved population", k, "m x q Zhat RDS:", zhat_rds_k, "\n")

      zhat_csv_k <- sub(
        "\\.csv$", paste0("_population", k, "_Zhat.csv"), csv_file
      )
      write.csv(selected_pair_list[[k]], zhat_csv_k, row.names = FALSE)
      cat("Saved population", k, "m x q Zhat csv:", zhat_csv_k, "\n")
    }
  }

  cat("\n--- mvSuSiE metrics by population ---\n")
  print(metrics_by_population)

  cat("\n--- mvSuSiE overall metrics ---\n")
  print(metrics_overall)

  invisible(out)
}
