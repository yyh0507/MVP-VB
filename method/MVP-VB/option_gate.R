get_gate_rules <- function(
    n_k,
    q,
    threshold = 0.50,
    evidence_n_ratio_max = 5,
    option2_support_fraction = 0.50
) {
  ## ------------------------------------------------------------
  ## Basic checks
  ## ------------------------------------------------------------
  n_k <- as.numeric(n_k)

  if (length(n_k) < 1L ||
      any(!is.finite(n_k)) ||
      any(n_k <= 0)) {
    stop("n_k must contain positive finite sample sizes.")
  }

  K <- length(n_k)

  if (length(q) != 1L || !is.finite(q) || q < 1) {
    stop("q must be a positive scalar.")
  }

  q <- as.integer(q)

  ## ------------------------------------------------------------
  ## Evidence-only sample-size uplift
  ##
  ## Keep the largest population at its raw sample size and lift
  ## smaller populations to at least max(n_k) / evidence_n_ratio_max.
  ## ------------------------------------------------------------
  n_max <- max(n_k)

  n_k_evidence <- pmax(
    n_k,
    n_max / evidence_n_ratio_max
  )

  z_rescale <- sqrt(n_k_evidence / n_k)

  ## ------------------------------------------------------------
  ## Option 2:
  ## same SNP-same phenotype replication across populations
  ## ------------------------------------------------------------
  n_external <- max(K - 1L, 0L)

  if (K >= 2L) {
    option2_required_external <- max(
      2L,
      as.integer(ceiling(
        option2_support_fraction * n_external
      ))
    )
  } else {
    option2_required_external <- NA_integer_
  }

  option2_current_z_min <- 0.75

  ## ------------------------------------------------------------
  ## Return all gate rules
  ## ------------------------------------------------------------
  list(
    ## Dimensions
    K = K,
    q = q,

    ## Sample-size quantities
    n_k_raw = n_k,
    n_k_evidence = n_k_evidence,
    z_rescale = z_rescale,
    evidence_n_ratio_max = evidence_n_ratio_max,

    ## Final W-selection threshold / compatibility with current code
    c_W_gate = threshold,

    ## ----------------------------------------------------------
    ## Option 2
    ## ----------------------------------------------------------
    option2_enabled = (K >= 2L),

    ## One external population counts as weak support when p < 0.10
    option2_alpha = 0.10,

    ## Required breadth of replication
    option2_support_fraction = option2_support_fraction,
    option2_required_external = option2_required_external,

    ## If exactly one external population passes option2_alpha,
    ## require LOO-meta confirmation at p < 0.05
    option2_meta_alpha = 0.05,

    ## Target population must show at least weak marginal evidence
    option2_current_z_min = option2_current_z_min,
    option2_current_alpha =
      2 * pnorm(-option2_current_z_min),

    ## Strong opposite-direction evidence blocks borrowing
    option2_oppose_alpha = 0.05,

    ## ----------------------------------------------------------
    ## Option 3
    ## ----------------------------------------------------------
    option3_enabled = (q >= 2L),

    ## At least one other phenotype must show strong meta evidence
    option3_other_alpha = 1e-4,
    min_other_signal = 1L,

    ## Current phenotype must show some marginal/meta evidence
    option3_current_alpha = 0.10,

    ## ----------------------------------------------------------
    ## Priority
    ## ----------------------------------------------------------
    option_priority = c(3L, 2L, 1L),

    ## ----------------------------------------------------------
    ## Legacy diagnostic fields retained so the current core
    ## does not fail if it still reads these quantities.
    ## They do not determine the final hard switch.
    ## ----------------------------------------------------------
    gate_quantile = 0.90,
    local_ok_rule = "either",
    option3_margin = 0.05
  )
}
