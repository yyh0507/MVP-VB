# MVP-VB

Code and result summaries for **Variational Bayes for Variable Selection with Adaptive Sharing across Populations and Responses**.

## Repository structure

```text
method/
  MVP-VB/
    initialization.R
    option_gate.R
    mvp_vb_core.R
    population_specific_variants.R
  comparison/
    mgflashfm.R
    multisusie_export.R
    multisusie.py
    mvsusie.R
    mvlasso.R
GLGC2021/
```

`method/MVP-VB` contains the initialization rules, adaptive sharing rules, the main variational Bayes routine, and the population-specific variants. Source the R files in the order shown above. The principal input is a list containing population-specific genotype and phenotype matrices (`X_list` and `Y_list`) together with their dimensions and sample sizes.

`method/comparison` contains the analysis functions used for MGflashfm, MultiSuSiE, mvSuSiE, and MVLasso. These files require the corresponding R or Python packages. MultiSuSiE inputs are prepared with `export_multisusie_inputs()` in `multisusie_export.R` and analysed with `run_multisusie()` in `multisusie.py`.

## GLGC 2021 results

The files in `GLGC2021` report variants with posterior inclusion probability `W > 0.5`.

- `single_Europe_*` contains the analysis restricted to the European population.
- `five_populations_joint_*` contains the joint analysis of AFR, EAS, EUR, HIS, and SAS.
- `*_selected_SNP_names_*` gives one row per selected variant.
- `*_population_trait_detail_*` gives population- and phenotype-specific posterior results, including `W`, posterior mean `M`, and the selected sharing option.

The result tables are provided as CSV files and retain the column order produced by the analysis.
