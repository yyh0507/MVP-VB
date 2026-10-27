## ============================================================
## MVP-VB core iteration only
##
## This file deliberately contains NO initialization rule function and
## NO gate-rule generator function.  It consumes the outputs produced by
## those earlier pipeline stages:
##
##   initialization_rules <- get_initialization_rules(...)
##   prepared             <- prepare_fit_inputs(...)
##   gate_rules           <- get_gate_rules(...)
##
## The core then performs the active-set VB iteration.  The inclusion prior
## is deliberately global and identical for Options 1/2/3:
##   pi_fix = 1 / (m*q),
##   a       = initialization_rules$a (fixed),
##   b       = a * (1-pi_fix) / pi_fix.
## ============================================================

fit_MVP_VB_core <- function(
    common_data,
    prepared,
    initialization_rules,
    gate_rules,
    scenario_name = NULL,
    max_iter = 200L,
    tol_elbo = 1,
    active_controls = NULL,
    sigmaB2 = 5.0,
    nu0 = NULL,
    Psi0 = NULL,
    group_prior_n_c = NULL,       # deprecated; ignored
    pi_option1_by_pop = NULL,     # deprecated; ignored
    min_fully_supported_groups = 0L,
    min_fully_supported_option3_groups = 1L,
    option1_fallback_when_shared_mode_off = FALSE # deprecated; ignored
) {
  ## ============================================================
  ## 0. Validate and unpack inputs from earlier pipeline stages
  ## ============================================================

  required_common <- c("X_list", "Y_list", "n_k", "K", "m", "q")
  missing_common <- required_common[
    !vapply(required_common, function(nm) !is.null(common_data[[nm]]), logical(1))
  ]
  if (length(missing_common) > 0L) {
    stop("common_data is missing: ", paste(missing_common, collapse = ", "))
  }

  required_prepared <- c("XtX_list", "XtY_list", "W_list")
  missing_prepared <- required_prepared[
    !vapply(required_prepared, function(nm) !is.null(prepared[[nm]]), logical(1))
  ]
  if (length(missing_prepared) > 0L) {
    stop("prepared is missing: ", paste(missing_prepared, collapse = ", "))
  }

  required_init_rules <- c("a", "ld_threshold")
  missing_init_rules <- required_init_rules[
    !vapply(required_init_rules, function(nm) !is.null(initialization_rules[[nm]]), logical(1))
  ]
  if (length(missing_init_rules) > 0L) {
    stop(
      "initialization_rules is missing: ",
      paste(missing_init_rules, collapse = ", ")
    )
  }

  required_gate_rules <- c(
    "n_k_evidence", "z_rescale",
    "option2_enabled", "option2_alpha", "option2_required_external",
    "option2_meta_alpha", "option2_current_z_min", "option2_oppose_alpha",
    "option3_enabled", "option3_other_alpha", "option3_current_alpha",
    "min_other_signal"
  )
  missing_gate_rules <- required_gate_rules[
    !vapply(required_gate_rules, function(nm) !is.null(gate_rules[[nm]]), logical(1))
  ]
  if (length(missing_gate_rules) > 0L) {
    stop("gate_rules is missing: ", paste(missing_gate_rules, collapse = ", "))
  }

  X_list <- common_data$X_list
  Y_list <- common_data$Y_list
  n_k <- as.numeric(common_data$n_k)
  K <- as.integer(common_data$K)
  m <- as.integer(common_data$m)
  q <- as.integer(common_data$q)

  if (length(X_list) != K || length(Y_list) != K || length(n_k) != K) {
    stop("X_list, Y_list and n_k must all agree with K.")
  }
  if (length(prepared$XtX_list) != K ||
      length(prepared$XtY_list) != K ||
      length(prepared$W_list) != K) {
    stop("prepared$XtX_list, prepared$XtY_list and prepared$W_list must have length K.")
  }
  if (length(gate_rules$n_k_evidence) != K || length(gate_rules$z_rescale) != K) {
    stop("gate_rules$n_k_evidence and gate_rules$z_rescale must have length K.")
  }

  for (k in seq_len(K)) {
    if (!all(dim(X_list[[k]]) == c(n_k[k], m))) {
      stop("X_list[[", k, "]] has dimensions inconsistent with n_k and m.")
    }
    if (!all(dim(Y_list[[k]]) == c(n_k[k], q))) {
      stop("Y_list[[", k, "]] has dimensions inconsistent with n_k and q.")
    }
    if (!all(dim(prepared$XtX_list[[k]]) == c(m, m))) {
      stop("prepared$XtX_list[[", k, "]] must be m x m.")
    }
    if (!all(dim(prepared$XtY_list[[k]]) == c(m, q))) {
      stop("prepared$XtY_list[[", k, "]] must be m x q.")
    }
    if (!all(dim(prepared$W_list[[k]]) == c(m, q))) {
      stop("prepared$W_list[[", k, "]] must be m x q.")
    }
  }

  XtX_list <- prepared$XtX_list
  XtY_list <- prepared$XtY_list
  W_list <- lapply(prepared$W_list, function(x) as.matrix(x))

  a <- as.numeric(initialization_rules$a)
  pi_fix <- 1 / (m * q)
  b <- a * (1 - pi_fix) / pi_fix
  if (length(a) != 1L || !is.finite(a) || a <= 0) {
    stop("initialization_rules$a must be one positive finite scalar.")
  }
  if (!is.finite(pi_fix) || pi_fix <= 0 || pi_fix >= 1 ||
      !is.finite(b) || b <= 0) {
    stop("The global prior pi=1/(m*q), a, b is invalid.")
  }
  ## Enforce one prior irrespective of population, group, or sharing mode.
  ## Deprecated caller arguments remain in the signature only so an older
  ## runner cannot silently select a different code path.
  group_prior_n_c <- NULL
  pi_option1_by_pop <- rep(pi_fix, K)
  option1_fallback_when_shared_mode_off <- FALSE
  initialization_rules_used <- initialization_rules
  initialization_rules_used$a <- a
  initialization_rules_used$b <- b
  initialization_rules_used$pi <- pi_fix
  min_fully_supported_groups <- as.integer(min_fully_supported_groups)
  if (length(min_fully_supported_groups) != 1L ||
      is.na(min_fully_supported_groups) || min_fully_supported_groups < 0L) {
    stop("min_fully_supported_groups must be one non-negative integer.")
  }
  min_fully_supported_option3_groups <-
    as.integer(min_fully_supported_option3_groups)
  if (length(min_fully_supported_option3_groups) != 1L ||
      is.na(min_fully_supported_option3_groups) ||
      min_fully_supported_option3_groups < 0L) {
    stop("min_fully_supported_option3_groups must be one non-negative integer.")
  }
  if (!is.finite(sigmaB2) || sigmaB2 <= 0) {
    stop("sigmaB2 must be positive and finite.")
  }

  if (is.null(nu0)) nu0 <- q + 2
  if (is.null(Psi0)) Psi0 <- diag(q)
  if (!all(dim(Psi0) == c(q, q))) stop("Psi0 must be q x q.")

  n_k_evidence <- as.numeric(gate_rules$n_k_evidence)
  z_rescale <- as.numeric(gate_rules$z_rescale)

  block_id_vec <- if (!is.null(common_data$block_id_vec) &&
                      length(common_data$block_id_vec) == m) {
    common_data$block_id_vec
  } else {
    rep(1L, m)
  }

  ## ============================================================
  ## 1. Core helper functions
  ## ============================================================

  logit <- function(p) log(p) - log1p(-p)

  group_prior_ab <- function(member_populations) {
    ## member_populations is retained for a common interface, but every
    ## Option 2/3 group uses exactly the same global Beta(a,b) prior.
    c(a = a, b = b, n_eff = NA_real_, pi = pi_fix)
  }

  log_multivariate_gamma <- function(x, dimension) {
    dimension * (dimension - 1) * log(pi) / 4 +
      sum(lgamma(x + (1 - seq_len(dimension)) / 2))
  }

  top_n_idx <- function(score, n) {
    n <- min(length(score), max(0L, as.integer(n)))
    if (n <= 0L) return(integer(0))
    order(score, decreasing = TRUE, na.last = TRUE)[seq_len(n)]
  }

  precompute_top1_ld <- function(XtX_k, block_id_vec = NULL, ld_min = 0) {
    m_local <- nrow(XtX_k)
    out <- vector("list", m_local)
    snp_sd <- sqrt(pmax(diag(XtX_k), 0))

    if (is.null(block_id_vec) || length(block_id_vec) != m_local) {
      block_list <- list(seq_len(m_local))
    } else {
      block_list <- split(seq_len(m_local), block_id_vec)
    }

    for (idx_block in block_list) {
      idx_block <- as.integer(idx_block)
      if (length(idx_block) <= 1L) next

      denom <- outer(snp_sd[idx_block], snp_sd[idx_block], `*`)
      corr_abs <- abs(XtX_k[idx_block, idx_block, drop = FALSE] / denom)
      corr_abs[!is.finite(corr_abs)] <- 0
      diag(corr_abs) <- 0

      for (pos in seq_along(idx_block)) {
        jpos <- which.max(corr_abs[pos, ])
        if (length(jpos) == 1L && corr_abs[pos, jpos] > ld_min) {
          out[[idx_block[pos]]] <- idx_block[jpos]
        } else {
          out[[idx_block[pos]]] <- integer(0)
        }
      }
    }

    out
  }

  get_ld_top1_neighbors <- function(seed_idx, ld_top1_k) {
    if (length(seed_idx) == 0L) return(integer(0))
    unique(unlist(ld_top1_k[seed_idx], use.names = FALSE))
  }

  compute_ES_active <- function(M_C, W_C, XtX_C, V_C, Tau_k,
                                sigmaB2, m, q) {
    c_k <- nrow(M_C)
    E_mat <- matrix(0, q, q)
    S_mat <- matrix(0, q, q)

    for (r in seq_len(q)) {
      rows_r <- ((r - 1L) * c_k + 1L):(r * c_k)
      Mr <- M_C[, r]
      Wr <- W_C[, r]

      for (s in r:q) {
        rows_s <- ((s - 1L) * c_k + 1L):(s * c_k)
        Ms <- M_C[, s]
        Ws <- W_C[, s]

        V_rs <- V_C[rows_r, rows_s, drop = FALSE]
        Omega_rs <- outer(Wr, Ws)
        if (r == s) {
          diag(Omega_rs) <- diag(Omega_rs) + Wr * (1 - Wr)
        }

        MMt <- outer(Mr, Ms)
        val_E <- sum(XtX_C * ((MMt + V_rs) * Omega_rs))
        val_S <- sum(Mr * Ms) + sum(diag(V_rs))

        E_mat[r, s] <- val_E
        S_mat[r, s] <- val_S
        if (s != r) {
          E_mat[s, r] <- val_E
          S_mat[s, r] <- val_S
        }
      }
    }

    if (c_k < m) {
      S_mat <- S_mat + (m - c_k) * sigmaB2 * solve(Tau_k)
    }

    list(E_mat = E_mat, S_mat = S_mat)
  }

  compute_delta_coordinate_active <- function(
      pos_i, r, W_C, M_C, V_C, XtX_C, XtY_C, Tau_k,
      other_r_list) {
    c_k <- nrow(W_C)
    rows_r <- ((r - 1L) * c_k + 1L):(r * c_k)
    row_ir <- rows_r[pos_i]
    Mr <- M_C[, r]

    delta1 <- -0.5 * Tau_k[r, r] * XtX_C[pos_i, pos_i] *
      (M_C[pos_i, r]^2 + V_C[row_ir, row_ir])
    delta2 <- M_C[pos_i, r] * as.numeric(XtY_C[pos_i, ] %*% Tau_k[, r])

    t3_row <- M_C[pos_i, r] * Mr + V_C[row_ir, rows_r]
    t3_row[pos_i] <- 0
    delta3 <- -Tau_k[r, r] * sum(
      XtX_C[pos_i, ] * t3_row * W_C[, r]
    )

    delta4 <- 0
    for (s in other_r_list[[r]]) {
      rows_s <- ((s - 1L) * c_k + 1L):(s * c_k)
      t4_row <- M_C[pos_i, r] * M_C[, s] + V_C[row_ir, rows_s]
      delta4 <- delta4 - Tau_k[r, s] * sum(
        XtX_C[pos_i, ] * t4_row * W_C[, s]
      )
    }

    delta1 + delta2 + delta3 + delta4
  }

  compute_marginal_stats <- function(Xk, Yk) {
    Xc <- scale(Xk, center = TRUE, scale = FALSE)
    Yc <- scale(Yk, center = TRUE, scale = FALSE)

    Xc <- as.matrix(Xc)
    Yc <- as.matrix(Yc)
    Xc[!is.finite(Xc)] <- 0
    Yc[!is.finite(Yc)] <- 0

    n_obs <- nrow(Xc)
    xtx <- colSums(Xc^2)
    xty <- crossprod(Xc, Yc)
    yty <- colSums(Yc^2)

    beta_hat <- matrix(0, nrow = ncol(Xc), ncol = ncol(Yc))
    se_hat <- matrix(Inf, nrow = ncol(Xc), ncol = ncol(Yc))
    z_mat <- matrix(0, nrow = ncol(Xc), ncol = ncol(Yc))
    valid_snp <- is.finite(xtx) & xtx > 1e-12

    if (any(valid_snp)) {
      fit_ss <- matrix(0, nrow = ncol(Xc), ncol = ncol(Yc))

      beta_hat[valid_snp, ] <- sweep(
        xty[valid_snp, , drop = FALSE], 1, xtx[valid_snp], "/"
      )
      fit_ss[valid_snp, ] <- sweep(
        xty[valid_snp, , drop = FALSE]^2, 1, xtx[valid_snp], "/"
      )

      rss <- matrix(yty, nrow = ncol(Xc), ncol = ncol(Yc), byrow = TRUE) -
        fit_ss
      rss[!is.finite(rss)] <- 0
      rss <- pmax(rss, .Machine$double.eps)

      sigma2_hat <- pmax(rss / max(n_obs - 2L, 1L), .Machine$double.eps)
      se_hat[valid_snp, ] <- sqrt(
        sweep(
          sigma2_hat[valid_snp, , drop = FALSE],
          1,
          xtx[valid_snp],
          "/"
        )
      )

      z_mat <- beta_hat / se_hat
    }

    t_stat_mat <- z_mat
    z_mat[!is.finite(z_mat)] <- 0
    df <- n_obs - 2L
    p_t_mat <- matrix(NA_real_, nrow = ncol(Xc), ncol = ncol(Yc))
    valid_test <- is.finite(beta_hat) & is.finite(se_hat) & se_hat > 0 &
      is.finite(t_stat_mat) & df > 0L
    p_t_mat[valid_test] <- 2 * pt(-abs(t_stat_mat[valid_test]), df = df)

    list(
      beta_hat = beta_hat,
      se_hat = se_hat,
      z_mat = z_mat,
      p_t_mat = p_t_mat,
      df = as.integer(df)
    )
  }

  ## ============================================================
  ## 2. Fixed marginal evidence used by the gate and active set
  ## ============================================================

  marginal_stats_list <- vector("list", K)
  for (k in seq_len(K)) {
    marginal_stats_list[[k]] <- compute_marginal_stats(X_list[[k]], Y_list[[k]])
  }

  loo_k_list <- lapply(seq_len(K), function(k) setdiff(seq_len(K), k))
  other_r_list <- lapply(seq_len(q), function(r) setdiff(seq_len(q), r))

  marginal_z_raw_array <- array(NA_real_, dim = c(K, m, q))
  marginal_z_array <- array(NA_real_, dim = c(K, m, q))
  marginal_beta_hat_array <- array(NA_real_, dim = c(K, m, q))
  marginal_se_hat_array <- array(NA_real_, dim = c(K, m, q))
  donor_t_p_array <- array(NA_real_, dim = c(K, m, q))
  donor_bh_p_array <- array(NA_real_, dim = c(K, m, q))
  for (k in seq_len(K)) {
    marginal_z_raw_array[k, , ] <- marginal_stats_list[[k]]$z_mat
    marginal_z_array[k, , ] <- marginal_stats_list[[k]]$z_mat * z_rescale[k]
    marginal_beta_hat_array[k, , ] <- marginal_stats_list[[k]]$beta_hat
    marginal_se_hat_array[k, , ] <- marginal_stats_list[[k]]$se_hat
    donor_t_p_array[k, , ] <- marginal_stats_list[[k]]$p_t_mat
    donor_bh_p_array[k, , ] <- matrix(
      p.adjust(as.vector(marginal_stats_list[[k]]$p_t_mat), method = "BH"),
      nrow = m,
      ncol = q
    )
  }

  p_marg_fixed_array <- 2 * pnorm(-abs(marginal_z_array))
  z_abs_fixed_array <- abs(marginal_z_array)

  sqrt_n_evidence <- sqrt(n_k_evidence)
  p_meta_all_fixed_mat <- matrix(NA_real_, nrow = m, ncol = q)

  if (sum(n_k_evidence) <= 0) stop("gate_rules$n_k_evidence must be positive.")

  for (i in seq_len(m)) {
    for (r in seq_len(q)) {
      z_vec_all_ir <- marginal_z_array[, i, r]
      meta_z_all <- sum(sqrt_n_evidence * z_vec_all_ir) /
        sqrt(sum(n_k_evidence))
      p_meta_all_fixed_mat[i, r] <- 2 * pnorm(-abs(meta_z_all))
    }
  }

  ## Fixed group-level gate constants. These values are common to every
  ## scenario and are not tuned using simulation truth or fit metrics.
  tost_alpha <- 0.05
  tost_equiv_margin <- 0.10
  donor_bh_alpha <- 0.05
  final_option_priority <- c(3L, 2L, 1L)

  option2_enabled <- isTRUE(gate_rules$option2_enabled)
  option2_alpha <- as.numeric(gate_rules$option2_alpha)
  option2_required_external <- gate_rules$option2_required_external
  option2_meta_alpha <- as.numeric(gate_rules$option2_meta_alpha)
  option2_current_z_min <- as.numeric(gate_rules$option2_current_z_min)
  option2_oppose_alpha <- as.numeric(gate_rules$option2_oppose_alpha)
  option3_enabled <- isTRUE(gate_rules$option3_enabled)
  option3_other_alpha <- as.numeric(gate_rules$option3_other_alpha)
  option3_current_alpha <- as.numeric(gate_rules$option3_current_alpha)
  min_other_signal <- as.integer(gate_rules$min_other_signal)

  ## TOST is evaluated only for recipient-eligible populations in the
  ## unequal-n design and for every population in the equal-n fallback.
  n_mean <- mean(n_k)
  equal_n_fallback <- all(n_k == n_k[1L])
  if (equal_n_fallback) {
    donor_pool <- seq_len(K)
    recipient_pool <- seq_len(K)
  } else {
    donor_pool <- which(n_k > n_mean)
    recipient_pool <- which(n_k < n_mean)
  }

  tost_veto_array <- array(FALSE, dim = c(K, m, q))
  tost_evaluated_array <- array(FALSE, dim = c(K, m, q))
  tost_invalid_array <- array(FALSE, dim = c(K, m, q))
  for (k in recipient_pool) {
    beta_k <- marginal_beta_hat_array[k, , ]
    se_k <- marginal_se_hat_array[k, , ]
    df_k <- n_k[k] - 2
    valid_k <- is.finite(beta_k) & is.finite(se_k) & se_k > 0 & df_k > 0
    tost_evaluated_array[k, , ] <- TRUE
    tost_invalid_array[k, , ] <- !valid_k
    if (df_k > 0) {
      upper_abs_k <- abs(beta_k) + qt(1 - tost_alpha, df = df_k) * se_k
      veto_k <- matrix(FALSE, nrow = m, ncol = q)
      veto_k[valid_k] <- upper_abs_k[valid_k] < tost_equiv_margin
      tost_veto_array[k, , ] <- veto_k
    }
  }

  ## Preserve the original Option 2/3 evidence rules as candidate triggers.
  ## They continue to use the existing z-rescaled normal-approximation p-values.
  option2_candidate <- array(FALSE, dim = c(K, m, q))
  option3_candidate <- array(FALSE, dim = c(K, m, q))
  if (K > 1L) {
    for (k in seq_len(K)) {
      for (i in seq_len(m)) {
        for (r in seq_len(q)) {
          loo_idx <- loo_k_list[[k]]
          z_vec_loo <- marginal_z_array[loo_idx, i, r]
          p_vec_loo <- p_marg_fixed_array[loo_idx, i, r]
          sig_idx <- which(p_vec_loo < option2_alpha)
          N_p_loo <- length(sig_idx)
          n_vec_loo <- n_k_evidence[loo_idx]
          meta_z_loo <- if (length(loo_idx) > 0L && sum(n_vec_loo) > 0) {
            sum(sqrt(n_vec_loo) * z_vec_loo) / sqrt(sum(n_vec_loo))
          } else 0
          p_meta_loo <- 2 * pnorm(-abs(meta_z_loo))

          if (N_p_loo >= 2L) {
            sig_signs <- sign(z_vec_loo[sig_idx])
            sig_signs <- sig_signs[sig_signs != 0]
            opt2_sign_ok <- length(unique(sig_signs)) <= 1L
          } else if (N_p_loo == 1L) {
            opt2_sign_ok <- sign(z_vec_loo[sig_idx]) == sign(meta_z_loo)
          } else {
            opt2_sign_ok <- FALSE
          }

          opt2_replication_ok <- option2_enabled &&
            !is.na(option2_required_external) &&
            N_p_loo >= option2_required_external && opt2_sign_ok
          opt2_one_pop_meta_ok <- option2_enabled && N_p_loo == 1L &&
            opt2_sign_ok && p_meta_loo < option2_meta_alpha
          current_z <- marginal_z_array[k, i, r]
          current_p <- p_marg_fixed_array[k, i, r]
          current_support_weak <- abs(current_z) >= option2_current_z_min
          current_opposite_strong <- current_p < option2_oppose_alpha &&
            sign(current_z) != 0 && sign(meta_z_loo) != 0 &&
            sign(current_z) != sign(meta_z_loo)
          option2_candidate[k, i, r] <-
            (opt2_replication_ok || opt2_one_pop_meta_ok) &&
            current_support_weak && !current_opposite_strong

          other_pheno <- other_r_list[[r]]
          D_other_p <- sum(
            p_meta_all_fixed_mat[i, other_pheno] < option3_other_alpha
          )
          current_pheno_has_evidence <-
            p_meta_all_fixed_mat[i, r] < option3_current_alpha ||
            p_marg_fixed_array[k, i, r] < option3_current_alpha
          option3_candidate[k, i, r] <- option3_enabled &&
            D_other_p >= min_other_signal && current_pheno_has_evidence
        }
      }
    }
  }

  ## Build mutually exclusive fixed groups once, before the VB iterations.
  option_array <- array(1L, dim = c(K, m, q))
  candidate_member_list <- vector("list", m * q)
  donor_member_list <- vector("list", m * q)
  index_ir <- function(i, r) (i - 1L) * q + r
  for (i in seq_len(m)) {
    for (r in seq_len(q)) {
      if (equal_n_fallback) {
        eligible <- which(!tost_veto_array[, i, r])
        donors <- eligible[
          !is.na(donor_bh_p_array[eligible, i, r]) &
            donor_bh_p_array[eligible, i, r] <= donor_bh_alpha
        ]
        recipients <- eligible
      } else {
        donors <- donor_pool[
          !is.na(donor_bh_p_array[donor_pool, i, r]) &
            donor_bh_p_array[donor_pool, i, r] <= donor_bh_alpha
        ]
        recipients <- recipient_pool[!tost_veto_array[recipient_pool, i, r]]
      }
      candidate_member_list[[index_ir(i, r)]] <-
        sort(unique(c(donors, recipients)))
      donor_member_list[[index_ir(i, r)]] <- sort(unique(donors))
    }
  }

  option3_groups <- list()
  option2_groups <- list()
  for (i in seq_len(m)) {
    eligible_traits <- integer(0)
    for (r in seq_len(q)) {
      members <- candidate_member_list[[index_ir(i, r)]]
      donors <- donor_member_list[[index_ir(i, r)]]
      if (length(donors) >= 1L && length(members) >= 2L) {
        eligible_traits <- c(eligible_traits, r)
      }
    }
    any_option3_trigger <- any(vapply(eligible_traits, function(r) {
      members <- candidate_member_list[[index_ir(i, r)]]
      any(option3_candidate[members, i, r])
    }, logical(1)))
    if (length(eligible_traits) >= 2L && any_option3_trigger) {
      group_pairs <- matrix(integer(0), ncol = 2L,
                            dimnames = list(NULL, c("population", "trait")))
      for (r in eligible_traits) {
        members <- candidate_member_list[[index_ir(i, r)]]
        option_array[members, i, r] <- 3L
        group_pairs <- rbind(group_pairs, cbind(members, r))
      }
      option3_groups[[as.character(i)]] <- list(
        snp = i,
        traits = as.integer(eligible_traits),
        pairs = group_pairs
      )
    }
  }

  for (i in seq_len(m)) {
    for (r in seq_len(q)) {
      members <- candidate_member_list[[index_ir(i, r)]]
      donors <- donor_member_list[[index_ir(i, r)]]
      remaining <- members[option_array[members, i, r] != 3L]
      remaining_donors <- intersect(remaining, donors)
      if (length(remaining_donors) >= 1L && length(remaining) >= 2L &&
          any(option2_candidate[remaining, i, r])) {
        option_array[remaining, i, r] <- 2L
        option2_groups[[paste(i, r, sep = ":")]] <- list(
          snp = i,
          trait = r,
          populations = as.integer(remaining)
        )
      }
    }
  }

  if (length(option_array) != K * m * q ||
      any(!option_array %in% c(1L, 2L, 3L))) {
    stop("Fixed group construction did not assign exactly one valid option per indicator.")
  }
  listed_option2 <- array(FALSE, dim = c(K, m, q))
  for (group in option2_groups) {
    listed_option2[group$populations, group$snp, group$trait] <- TRUE
  }
  listed_option3 <- array(FALSE, dim = c(K, m, q))
  for (group in option3_groups) {
    pairs <- group$pairs
    listed_option3[cbind(pairs[, "population"], group$snp, pairs[, "trait"])] <- TRUE
  }
  if (!identical(listed_option2, option_array == 2L) ||
      !identical(listed_option3, option_array == 3L) ||
      any(listed_option2 & listed_option3)) {
    stop("Fixed group membership is inconsistent or overlapping.")
  }
  fully_supported_option2 <- vapply(option2_groups, function(group) {
    all(option2_candidate[group$populations, group$snp, group$trait])
  }, logical(1))
  fully_supported_option3 <- vapply(option3_groups, function(group) {
    pairs <- group$pairs
    idx3 <- cbind(pairs[, "population"], group$snp, pairs[, "trait"])
    all(option3_candidate[idx3])
  }, logical(1))
  n_fully_supported_option2_groups <- sum(fully_supported_option2)
  n_fully_supported_option3_groups <- sum(fully_supported_option3)
  n_fully_supported_groups <- n_fully_supported_option2_groups +
    n_fully_supported_option3_groups
  shared_mode_enabled <- min_fully_supported_groups == 0L ||
    n_fully_supported_option2_groups >= min_fully_supported_groups ||
    (min_fully_supported_option3_groups > 0L &&
       n_fully_supported_option3_groups >= min_fully_supported_option3_groups)
  if (!shared_mode_enabled) {
    option_array[] <- 1L
    option2_groups <- list()
    option3_groups <- list()
  }
  protected_snp_list <- lapply(seq_len(K), function(k) {
    which(apply(option_array[k, , , drop = FALSE] != 1L, 2L, any))
  })

  ## ============================================================
  ## 3. Computational active set
  ## ============================================================

  active_controls_default <- list(
    use_active = TRUE,
    open_iter = 10L,
    drop_eps = 5e-4,
    c_min = 50L,
    c_frac_n = 1,
    z_frac = 0.20,
    c_max = if (m <= 1000L) 300L else 500L,
    use_ld_neighbor = TRUE,
    ld_min = initialization_rules$ld_threshold,
    store_active_history = FALSE
  )

  active_controls_used <- active_controls_default
  if (!is.null(active_controls)) {
    active_controls_used[names(active_controls)] <- active_controls
  }

  active_use <- isTRUE(active_controls_used$use_active)
  active_open_iter <- as.integer(active_controls_used$open_iter)
  active_drop_eps <- as.numeric(active_controls_used$drop_eps)
  active_c_min <- as.integer(active_controls_used$c_min)
  active_c_frac_n <- as.numeric(active_controls_used$c_frac_n)
  active_z_frac <- as.numeric(active_controls_used$z_frac)
  active_use_ld_neighbor <- isTRUE(active_controls_used$use_ld_neighbor)
  active_ld_min <- as.numeric(active_controls_used$ld_min)
  active_c_max <- as.integer(active_controls_used$c_max)

  if (!is.finite(active_c_max) || active_c_max < 1L) active_c_max <- m

  if (!active_use) {
    c_target_k <- rep(m, K)
  } else {
    c_target_k <- pmin(
      m,
      active_c_max,
      pmax(active_c_min, floor(n_k * active_c_frac_n))
    )
  }
  c_target_k <- as.integer(pmax(1L, c_target_k))

  c_z_k <- as.integer(
    pmin(c_target_k, pmax(5L, ceiling(active_z_frac * c_target_k)))
  )
  c_w_k <- as.integer(pmax(1L, c_target_k - c_z_k))
  active_min_size <- as.integer(min(m, max(1L, ceiling(0.05 * m))))

  z_score_snp_list <- vector("list", K)
  Z_watch_list <- vector("list", K)
  C_list <- vector("list", K)
  drop_count_list <- vector("list", K)
  V_active_list <- vector("list", K)
  active_idx_list <- vector("list", K)
  active_frozen <- rep(FALSE, K)

  ld_top1_list <- vector("list", K)
  if (active_use && active_use_ld_neighbor) {
    for (k in seq_len(K)) {
      ld_top1_list[[k]] <- precompute_top1_ld(
        XtX_k = XtX_list[[k]],
        block_id_vec = block_id_vec,
        ld_min = active_ld_min
      )
    }
  } else {
    for (k in seq_len(K)) ld_top1_list[[k]] <- vector("list", m)
  }

  for (k in seq_len(K)) {
    z_score_snp_list[[k]] <- apply(
      abs(marginal_z_raw_array[k, , , drop = FALSE]),
      2,
      max
    )

    if (!active_use) {
      C_init <- seq_len(m)
      Z_watch_list[[k]] <- seq_len(m)
    } else {
      Z_watch_list[[k]] <- top_n_idx(z_score_snp_list[[k]], c_z_k[k])

      W_score <- apply(W_list[[k]], 1, max)
      W_top <- top_n_idx(W_score, c_w_k[k])
      seed_idx <- sort(unique(c(W_top, Z_watch_list[[k]])))
      LD_top <- if (active_use_ld_neighbor) {
        get_ld_top1_neighbors(seed_idx, ld_top1_list[[k]])
      } else {
        integer(0)
      }

      C_init <- sort(unique(c(seed_idx, LD_top, protected_snp_list[[k]])))

      if (length(C_init) < c_target_k[k]) {
        filler_score <- W_score + 1e-6 * z_score_snp_list[[k]]
        filler <- setdiff(top_n_idx(filler_score, m), C_init)
        C_init <- sort(unique(c(
          C_init,
          head(filler, c_target_k[k] - length(C_init))
        )))
      }
    }

    C_list[[k]] <- C_init
    active_idx_list[[k]] <- C_init
    drop_count_list[[k]] <- integer(m)

    inactive_idx <- setdiff(seq_len(m), C_init)
    if (length(inactive_idx) > 0L) {
      W_list[[k]][inactive_idx, ] <- 0
    }
  }

  ## ============================================================
  ## 4. Initialize variational quantities used by the core
  ## ============================================================

  Tau_list <- vector("list", K)
  M_list <- vector("list", K)
  Tau0 <- nu0 * chol2inv(chol(Psi0))
  for (k in seq_len(K)) {
    ## E_q[Sigma_k^{-1}] under the initial IW(nu0, Psi0) factor.
    Tau_list[[k]] <- Tau0
    M_list[[k]] <- matrix(0, m, q)
  }

  YtY_list <- lapply(Y_list, crossprod)

  lambda_opt2_group_array <- array(NA_real_, dim = c(K, m, q))
  lambda_opt3_group_array <- array(NA_real_, dim = c(K, m, q))
  lambda_update_array <- array(NA_real_, dim = c(K, m, q))
  W_after_update_array <- array(NA_real_, dim = c(K, m, q))

  iter_summary <- data.frame(
    iter = seq_len(max_iter),
    ELBO = NA_real_,
    ELBO_change = NA_real_,
    ELBO_continuous = NA_real_,
    ELBO_entropy = NA_real_,
    ELBO_option1 = NA_real_,
    ELBO_option2 = NA_real_,
    ELBO_option3 = NA_real_,
    ELBO_trace_correction = NA_real_,
    beta_sync_max_error = NA_real_,
    option1_count = NA_integer_,
    option2_count = NA_integer_,
    option3_count = NA_integer_,
    option_change = NA_integer_,
    n_option2_ok = NA_integer_,
    n_option3_ok = NA_integer_,
    mean_delta_gate_opt2 = NA_real_,
    mean_lambda_opt2 = NA_real_,
    mean_W_after_opt2 = NA_real_,
    stringsAsFactors = FALSE
  )

  active_size_history <- matrix(
    NA_integer_,
    nrow = max_iter,
    ncol = K,
    dimnames = list(NULL, paste0("pop", seq_len(K)))
  )
  protected_active_count_history <- matrix(
    NA_integer_, nrow = max_iter, ncol = K,
    dimnames = list(NULL, paste0("pop", seq_len(K)))
  )
  protected_active_all_history <- matrix(
    NA, nrow = max_iter, ncol = K,
    dimnames = list(NULL, paste0("pop", seq_len(K)))
  )

  store_active_history <- isTRUE(active_controls_used$store_active_history)
  active_idx_history <- if (store_active_history) vector("list", max_iter) else NULL

  lambda_opt1 <- logit(pi_fix)
  elbo_prev <- -Inf
  converged <- FALSE
  convergence_reason <- "max_iter"

  ## ============================================================
  ## 5. Main VB iteration
  ## ============================================================

  runtime_sec <- system.time({
    for (it in seq_len(max_iter)) {
      logdet_V_list <- vector("list", K)
      logdet_Psi_list <- vector("list", K)
      T_tilde_list <- vector("list", K)
      T_current_list <- vector("list", K)
      Tau_gaussian_list <- vector("list", K)
      trace_correction_list <- numeric(K)
      delta_gate_array <- array(-Inf, dim = c(K, m, q))

      ## ----------------------------------------------------------
      ## 5.1 Dynamic active-set update during the opening iterations
      ## ----------------------------------------------------------
      if (active_use && it > 1L && it <= active_open_iter) {
        for (k in seq_len(K)) {
          if (active_frozen[k]) next

          current_active <- C_list[[k]]
          W_score <- apply(W_list[[k]], 1, max)
          W_top <- top_n_idx(W_score, c_w_k[k])

          seed_idx <- sort(unique(c(W_top, Z_watch_list[[k]])))
          LD_top <- if (active_use_ld_neighbor) {
            get_ld_top1_neighbors(seed_idx, ld_top1_list[[k]])
          } else {
            integer(0)
          }
          C_new <- sort(unique(c(seed_idx, LD_top)))

          low_now <- rep(FALSE, m)
          if (length(current_active) > 0L) {
            low_now[current_active] <- W_score[current_active] < active_drop_eps
          }
          drop_count_list[[k]][low_now] <- drop_count_list[[k]][low_now] + 1L
          drop_count_list[[k]][!low_now] <- 0L

          drop_idx <- which(drop_count_list[[k]] >= 2L)
          C_new <- setdiff(C_new, drop_idx)
          C_new <- sort(unique(c(C_new, protected_snp_list[[k]])))

          if (length(C_new) < active_min_size) {
            C_list[[k]] <- sort(unique(c(
              top_n_idx(z_score_snp_list[[k]], active_min_size),
              protected_snp_list[[k]]
            )))
            drop_count_list[[k]] <- integer(m)
            active_frozen[k] <- TRUE
            next
          }

          if (length(C_new) < c_target_k[k]) {
            filler_score <- W_score + 1e-6 * z_score_snp_list[[k]]
            filler <- setdiff(
              top_n_idx(filler_score, m),
              c(C_new, drop_idx)
            )
            C_new <- sort(unique(c(
              C_new,
              head(filler, c_target_k[k] - length(C_new))
            )))
          }

          C_list[[k]] <- sort(C_new)
        }
      }

      for (k in seq_len(K)) {
        active_size_history[it, k] <- length(C_list[[k]])
        protected_active_count_history[it, k] <-
          sum(protected_snp_list[[k]] %in% C_list[[k]])
        protected_active_all_history[it, k] <-
          all(protected_snp_list[[k]] %in% C_list[[k]])
      }
      if (store_active_history) {
        active_idx_history[[it]] <- lapply(C_list, identity)
      }

      ## ----------------------------------------------------------
      ## 5.2 Gaussian coefficient update and covariance update
      ## ----------------------------------------------------------
      for (k in seq_len(K)) {
        idx <- C_list[[k]]
        c_k <- length(idx)
        active_idx_list[[k]] <- idx

        W_full <- matrix(0, m, q)
        W_full[idx, ] <- W_list[[k]][idx, , drop = FALSE]
        W_list[[k]] <- W_full

        W_C <- W_full[idx, , drop = FALSE]
        wk <- as.vector(W_C)
        Tau_k <- Tau_list[[k]]
        Tau_gaussian_list[[k]] <- Tau_k
        XtX_C <- XtX_list[[k]][idx, idx, drop = FALSE]
        XtY_C <- XtY_list[[k]][idx, , drop = FALSE]
        X_k <- X_list[[k]]
        Y_k <- Y_list[[k]]

        K_prior_k <- kronecker(Tau_k, diag(c_k) / sigmaB2)
        K_like_k <- kronecker(Tau_k, XtX_C)

        Vinv_k <- K_prior_k +
          sweep(sweep(K_like_k, 1, wk, `*`), 2, wk, `*`)

        diag_like_k <- rep(diag(Tau_k), each = c_k) *
          rep(diag(XtX_C), times = q)
        diag(Vinv_k) <- diag(Vinv_k) +
          diag_like_k * wk * (1 - wk)

        chol_Vinv_k <- chol(Vinv_k)
        logdet_active <- -2 * sum(log(diag(chol_Vinv_k)))
        logdet_tau <- 2 * sum(log(diag(chol(Tau_k))))
        logdet_inactive <- (m - c_k) *
          (q * log(sigmaB2) - logdet_tau)
        logdet_V_list[[k]] <- logdet_active + logdet_inactive

        V_C <- chol2inv(chol_Vinv_k)
        V_active_list[[k]] <- V_C

        RHS_mat_C <- XtY_C %*% Tau_k
        rhs_C <- as.vector(W_C * RHS_mat_C)
        tmp_k <- forwardsolve(t(chol_Vinv_k), rhs_C)
        m_new_C <- as.numeric(backsolve(chol_Vinv_k, tmp_k))

        M_C <- matrix(m_new_C, c_k, q)
        M_full <- matrix(0, m, q)
        M_full[idx, ] <- M_C
        M_list[[k]] <- M_full

        ES_k <- compute_ES_active(
          M_C = M_C,
          W_C = W_C,
          XtX_C = XtX_C,
          V_C = V_C,
          Tau_k = Tau_k,
          sigmaB2 = sigmaB2,
          m = m,
          q = q
        )

        Bhat_C <- W_C * M_C
        XBhat_C <- X_k[, idx, drop = FALSE] %*% Bhat_C
        YtX_Bhat_k <- crossprod(Y_k, XBhat_C)
        Bt_XtY_k <- t(Bhat_C) %*% XtY_C
        YtY_k <- YtY_list[[k]]

        Psi_scale_k <- Psi0 + YtY_k -
          Bt_XtY_k -
          YtX_Bhat_k +
          ES_k$E_mat +
          (1 / sigmaB2) * ES_k$S_mat

        Psi_scale_k[!is.finite(Psi_scale_k)] <- 0
        ## Numerical covariance ridge: added only to the stored IW scale.
        ## It is not folded into Psi0 or the theoretical T_current below.
        Psi_scale_k <- Psi_scale_k + 1e-6 * diag(q)
        T_tilde_list[[k]] <- Psi_scale_k

        chol_Psi <- chol(Psi_scale_k)
        logdet_Psi_list[[k]] <- 2 * sum(log(diag(chol_Psi)))

        nu_star <- nu0 + n_k[k] + m
        Tau_list[[k]] <- nu_star * chol2inv(chol_Psi)

      }

      ## ----------------------------------------------------------
      ## 5.3 Initialize fixed-group Beta-factor records
      ## ----------------------------------------------------------
      W_old_array <- array(0, dim = c(K, m, q))
      for (k in seq_len(K)) W_old_array[k, , ] <- W_list[[k]]

      lambda_opt2_group_array[] <- NA_real_
      lambda_opt3_group_array[] <- NA_real_
      lambda_update_array[] <- NA_real_
      W_after_update_array[] <- NA_real_
      beta_factor_option2 <- list()
      beta_factor_option3 <- list()

      for (group_name in names(option2_groups)) {
        group <- option2_groups[[group_name]]
        i <- group$snp
        r <- group$trait
        members <- group$populations
        J2 <- length(members)
        u2_old <- sum(W_old_array[members, i, r])
        prior2 <- group_prior_ab(members)
        a2_prior <- unname(prior2["a"])
        b2_prior <- unname(prior2["b"])
        a2_old <- a2_prior + u2_old
        b2_old <- b2_prior + J2 - u2_old
        beta_factor_option2[[group_name]] <- list(
          J = J2, u_old = u2_old, a_old = a2_old, b_old = b2_old,
          prior_a = a2_prior, prior_b = b2_prior,
          n_eff = unname(prior2["n_eff"]), pi_group = unname(prior2["pi"])
        )
      }

      for (group_name in names(option3_groups)) {
        group <- option3_groups[[group_name]]
        pairs <- group$pairs
        idx3 <- cbind(pairs[, "population"], group$snp, pairs[, "trait"])
        J3 <- nrow(idx3)
        u3_old <- sum(W_old_array[idx3])
        prior3 <- group_prior_ab(pairs[, "population"])
        a3_prior <- unname(prior3["a"])
        b3_prior <- unname(prior3["b"])
        a3_old <- a3_prior + u3_old
        b3_old <- b3_prior + J3 - u3_old
        beta_factor_option3[[group_name]] <- list(
          J = J3, u_old = u3_old, a_old = a3_old, b_old = b3_old,
          prior_a = a3_prior, prior_b = b3_prior,
          n_eff = unname(prior3["n_eff"]), pi_group = unname(prior3["pi"])
        )
      }

      ## ----------------------------------------------------------
      ## 5.4 Sequential coordinate update in k, then i, then r order
      ## ----------------------------------------------------------
      ## M, V and Tau remain fixed during this block.  Every delta, shared
      ## Beta factor and lambda is refreshed from the latest W_temp, so an
      ## updated coordinate is available to all later coordinates in the
      ## same VB iteration.
      W_temp <- lapply(W_list, function(W) W)

      for (k in seq_len(K)) {
        active_idx <- active_idx_list[[k]]
        XtX_C <- XtX_list[[k]][active_idx, active_idx, drop = FALSE]
        XtY_C <- XtY_list[[k]][active_idx, , drop = FALSE]
        M_C <- M_list[[k]][active_idx, , drop = FALSE]
        V_C <- V_active_list[[k]]
        Tau_k <- Tau_list[[k]]
        W_C <- W_temp[[k]][active_idx, , drop = FALSE]

        for (pos_i in seq_along(active_idx)) {
          i <- active_idx[pos_i]
          for (r in seq_len(q)) {
            delta_ir <- compute_delta_coordinate_active(
              pos_i = pos_i,
              r = r,
              W_C = W_C,
              M_C = M_C,
              V_C = V_C,
              XtX_C = XtX_C,
              XtY_C = XtY_C,
              Tau_k = Tau_k,
              other_r_list = other_r_list
            )
            delta_gate_array[k, i, r] <- delta_ir
            option_kir <- option_array[k, i, r]

            if (option_kir == 1L) {
              lambda_ir <- lambda_opt1
            } else if (option_kir == 2L) {
              group_name <- paste(i, r, sep = ":")
              group <- option2_groups[[group_name]]
              if (is.null(group) || !k %in% group$populations) {
                stop("Missing Option 2 group for coordinate ", k, ":", i, ":", r, ".")
              }
              members <- group$populations
              u_current <- sum(vapply(
                members,
                function(member_k) W_temp[[member_k]][i, r],
                numeric(1)
              ))
              prior <- group_prior_ab(members)
              a_current <- unname(prior["a"]) + u_current
              b_current <- unname(prior["b"]) + length(members) - u_current
              lambda_ir <- digamma(a_current) - digamma(b_current)
              lambda_opt2_group_array[k, i, r] <- lambda_ir
              beta_factor_option2[[group_name]]$u_last_coordinate <- u_current
              beta_factor_option2[[group_name]]$a_last_coordinate <- a_current
              beta_factor_option2[[group_name]]$b_last_coordinate <- b_current
            } else if (option_kir == 3L) {
              group_name <- as.character(i)
              group <- option3_groups[[group_name]]
              if (is.null(group)) {
                stop("Missing Option 3 group for coordinate ", k, ":", i, ":", r, ".")
              }
              pairs <- group$pairs
              pair_match <- pairs[, "population"] == k & pairs[, "trait"] == r
              if (!any(pair_match)) {
                stop("Coordinate is absent from its Option 3 group: ", k, ":", i, ":", r, ".")
              }
              u_current <- sum(vapply(
                seq_len(nrow(pairs)),
                function(j) {
                  W_temp[[pairs[j, "population"]]][i, pairs[j, "trait"]]
                },
                numeric(1)
              ))
              prior <- group_prior_ab(pairs[, "population"])
              a_current <- unname(prior["a"]) + u_current
              b_current <- unname(prior["b"]) + nrow(pairs) - u_current
              lambda_ir <- digamma(a_current) - digamma(b_current)
              lambda_opt3_group_array[k, i, r] <- lambda_ir
              beta_factor_option3[[group_name]]$u_last_coordinate <- u_current
              beta_factor_option3[[group_name]]$a_last_coordinate <- a_current
              beta_factor_option3[[group_name]]$b_last_coordinate <- b_current
            } else {
              stop("Unknown option_kir.")
            }

            W_new <- plogis(lambda_ir + delta_ir)
            lambda_update_array[k, i, r] <- lambda_ir
            W_after_update_array[k, i, r] <- W_new
            W_temp[[k]][i, r] <- W_new
            W_C[pos_i, r] <- W_new
          }
        }
      }
      W_list <- W_temp

      ## ----------------------------------------------------------
      ## 5.5 ELBO
      ## ----------------------------------------------------------
      elbo_continuous <- 0
      covariance_ridge <- 1e-6
      logdet_Psi0 <- 2 * sum(log(diag(chol(Psi0))))
      for (k in seq_len(K)) {
        idx <- active_idx_list[[k]]
        W_C_new <- W_list[[k]][idx, , drop = FALSE]
        M_C <- M_list[[k]][idx, , drop = FALSE]
        V_C <- V_active_list[[k]]
        XtX_C <- XtX_list[[k]][idx, idx, drop = FALSE]
        XtY_C <- XtY_list[[k]][idx, , drop = FALSE]

        ES_current <- compute_ES_active(
          M_C = M_C,
          W_C = W_C_new,
          XtX_C = XtX_C,
          V_C = V_C,
          Tau_k = Tau_gaussian_list[[k]],
          sigmaB2 = sigmaB2,
          m = m,
          q = q
        )
        Bhat_C_new <- W_C_new * M_C
        XBhat_C_new <- X_list[[k]][, idx, drop = FALSE] %*% Bhat_C_new
        T_current_k <- Psi0 + YtY_list[[k]] -
          t(Bhat_C_new) %*% XtY_C -
          crossprod(Y_list[[k]], XBhat_C_new) +
          ES_current$E_mat + (1 / sigmaB2) * ES_current$S_mat
        if (any(!is.finite(T_current_k))) {
          stop("Non-finite theoretical covariance scale in population ", k, ".")
        }
        T_current_list[[k]] <- T_current_k

        d_k <- nu0 + n_k[k] + m
        Psi_star_k <- Tau_list[[k]]
        trace_correction_list[k] <- 0.5 * sum(diag(
          (T_tilde_list[[k]] - T_current_k) %*% Psi_star_k
        ))
        elbo_continuous <- elbo_continuous +
          -(n_k[k] * q / 2) * log(2 * pi) -
          (m * q / 2) * log(sigmaB2) +
          (m * q / 2) +
          0.5 * logdet_V_list[[k]] +
          (nu0 / 2) * logdet_Psi0 +
          ((n_k[k] + m) * q / 2) * log(2) +
          log_multivariate_gamma(d_k / 2, q) -
          log_multivariate_gamma(nu0 / 2, q) -
          (d_k / 2) * logdet_Psi_list[[k]] +
          trace_correction_list[k]
      }

      ## Bernoulli entropy.
      W_new_array <- array(0, dim = c(K, m, q))
      for (k in seq_len(K)) W_new_array[k, , ] <- W_list[[k]]
      entropy_mask <- W_new_array > 0 & W_new_array < 1
      elbo_entropy <- -sum(
        W_new_array[entropy_mask] * log(W_new_array[entropy_mask]) +
          (1 - W_new_array[entropy_mask]) *
          log1p(-W_new_array[entropy_mask])
      )

      ## Option 1 prior applies only to final Option 1 indicators.
      option1_mask <- option_array == 1L
      w_option1 <- W_new_array[option1_mask]
      elbo_option1 <- sum(
        w_option1 * log(pi_fix) +
          (1 - w_option1) * log1p(-pi_fix)
      )

      ## Synchronize every shared Beta factor to the new W, then use the
      ## matching actual group for both the simplified and general formulas.
      elbo_option2 <- 0
      elbo_option3 <- 0
      beta_sync_errors <- numeric(0)
      for (group_name in names(option2_groups)) {
        group <- option2_groups[[group_name]]
        J2 <- length(group$populations)
        u2_new <- sum(W_new_array[group$populations, group$snp, group$trait])
        a2_prior <- beta_factor_option2[[group_name]]$prior_a
        b2_prior <- beta_factor_option2[[group_name]]$prior_b
        a2_new <- a2_prior + u2_new
        b2_new <- b2_prior + J2 - u2_new
        simple2 <- lbeta(a2_new, b2_new) - lbeta(a2_prior, b2_prior)
        general2 <- lbeta(a2_new, b2_new) - lbeta(a2_prior, b2_prior) +
          (a2_prior + u2_new - a2_new) * (digamma(a2_new) - digamma(a2_new + b2_new)) +
          (b2_prior + J2 - u2_new - b2_new) *
          (digamma(b2_new) - digamma(a2_new + b2_new))
        beta_sync_errors <- c(beta_sync_errors, abs(simple2 - general2))
        elbo_option2 <- elbo_option2 + simple2
        beta_factor_option2[[group_name]]$u_new <- u2_new
        beta_factor_option2[[group_name]]$a_new <- a2_new
        beta_factor_option2[[group_name]]$b_new <- b2_new
      }
      for (group_name in names(option3_groups)) {
        group <- option3_groups[[group_name]]
        pairs <- group$pairs
        idx3 <- cbind(pairs[, "population"], group$snp, pairs[, "trait"])
        J3 <- nrow(idx3)
        u3_new <- sum(W_new_array[idx3])
        a3_prior <- beta_factor_option3[[group_name]]$prior_a
        b3_prior <- beta_factor_option3[[group_name]]$prior_b
        a3_new <- a3_prior + u3_new
        b3_new <- b3_prior + J3 - u3_new
        simple3 <- lbeta(a3_new, b3_new) - lbeta(a3_prior, b3_prior)
        general3 <- lbeta(a3_new, b3_new) - lbeta(a3_prior, b3_prior) +
          (a3_prior + u3_new - a3_new) * (digamma(a3_new) - digamma(a3_new + b3_new)) +
          (b3_prior + J3 - u3_new - b3_new) *
          (digamma(b3_new) - digamma(a3_new + b3_new))
        beta_sync_errors <- c(beta_sync_errors, abs(simple3 - general3))
        elbo_option3 <- elbo_option3 + simple3
        beta_factor_option3[[group_name]]$u_new <- u3_new
        beta_factor_option3[[group_name]]$a_new <- a3_new
        beta_factor_option3[[group_name]]$b_new <- b3_new
      }
      beta_sync_max_error <- if (length(beta_sync_errors)) {
        max(beta_sync_errors)
      } else 0
      if (!is.finite(beta_sync_max_error) || beta_sync_max_error > 1e-10) {
        stop("Synchronized Beta-factor ELBO check failed.")
      }

      elbo_curr <- elbo_continuous + elbo_entropy + elbo_option1 +
        elbo_option2 + elbo_option3

      ## ----------------------------------------------------------
      ## 5.6 Diagnostics and convergence
      ## ----------------------------------------------------------
      option1_count_now <- sum(option_array == 1L)
      option2_count_now <- sum(option_array == 2L)
      option3_count_now <- sum(option_array == 3L)

      option_change_now <- 0L

      opt2_mask <- option_array == 2L
      mean_delta_gate_opt2 <- if (any(opt2_mask)) {
        mean(delta_gate_array[opt2_mask], na.rm = TRUE)
      } else {
        NA_real_
      }
      mean_lambda_opt2 <- if (any(opt2_mask)) {
        mean(lambda_update_array[opt2_mask], na.rm = TRUE)
      } else {
        NA_real_
      }
      mean_W_after_opt2 <- if (any(opt2_mask)) {
        mean(W_after_update_array[opt2_mask], na.rm = TRUE)
      } else {
        NA_real_
      }

      iter_summary$ELBO[it] <- elbo_curr
      iter_summary$ELBO_change[it] <- if (it > 1L) elbo_curr - elbo_prev else NA_real_
      iter_summary$ELBO_continuous[it] <- elbo_continuous
      iter_summary$ELBO_entropy[it] <- elbo_entropy
      iter_summary$ELBO_option1[it] <- elbo_option1
      iter_summary$ELBO_option2[it] <- elbo_option2
      iter_summary$ELBO_option3[it] <- elbo_option3
      iter_summary$ELBO_trace_correction[it] <- sum(trace_correction_list)
      iter_summary$beta_sync_max_error[it] <- beta_sync_max_error
      iter_summary$option1_count[it] <- option1_count_now
      iter_summary$option2_count[it] <- option2_count_now
      iter_summary$option3_count[it] <- option3_count_now
      iter_summary$option_change[it] <- option_change_now
      iter_summary$n_option2_ok[it] <- sum(option2_candidate, na.rm = TRUE)
      iter_summary$n_option3_ok[it] <- sum(option3_candidate, na.rm = TRUE)
      iter_summary$mean_delta_gate_opt2[it] <- mean_delta_gate_opt2
      iter_summary$mean_lambda_opt2[it] <- mean_lambda_opt2
      iter_summary$mean_W_after_opt2[it] <- mean_W_after_opt2

      cat(sprintf(
        paste0(
          "Iter %d: ELBO=%.7f | opt1=%d | opt2=%d | opt3=%d | ",
          "option_change=%s | mean_delta_opt2=%.3f | ",
          "mean_lambda_opt2=%.3f | mean_W_after_opt2=%.3f\n"
        ),
        it,
        elbo_curr,
        option1_count_now,
        option2_count_now,
        option3_count_now,
        ifelse(is.na(option_change_now), "NA", as.character(option_change_now)),
        mean_delta_gate_opt2,
        mean_lambda_opt2,
        mean_W_after_opt2
      ))

      if (it > 1L) {
        elbo_diff <- abs(elbo_curr - elbo_prev)
        if (is.finite(elbo_diff) && elbo_diff < tol_elbo) {
          converged <- TRUE
          convergence_reason <- "elbo"
          cat(sprintf(
            "Converged at iter %d | Delta ELBO = %.3e\n",
            it,
            elbo_diff
          ))
          break
        }
      }

      elbo_prev <- elbo_curr
    }

    n_iter_done <- it
  })[["elapsed"]]

  ## ============================================================
  ## 6. Return core outputs
  ## ============================================================

  iter_summary <- iter_summary[seq_len(n_iter_done), , drop = FALSE]
  active_size_history <- active_size_history[seq_len(n_iter_done), , drop = FALSE]
  protected_active_count_history <-
    protected_active_count_history[seq_len(n_iter_done), , drop = FALSE]
  protected_active_all_history <-
    protected_active_all_history[seq_len(n_iter_done), , drop = FALSE]
  if (store_active_history) {
    active_idx_history <- active_idx_history[seq_len(n_iter_done)]
  }

  performance_by_population <- NULL
  if (!is.null(common_data$B_true_list) &&
      length(common_data$B_true_list) == K) {
    performance_by_population <- do.call(rbind, lapply(seq_len(K), function(k) {
      truth <- as.matrix(common_data$B_true_list[[k]]) != 0
      pred <- W_list[[k]] > 0.5
      if (!identical(dim(truth), c(m, q))) {
        stop("B_true_list dimensions do not agree with m and q in population ", k, ".")
      }
      TP <- sum(pred & truth)
      FP <- sum(pred & !truth)
      FN <- sum(!pred & truth)
      precision <- if (TP + FP > 0) TP / (TP + FP) else 0
      recall <- if (TP + FN > 0) TP / (TP + FN) else 0
      F1 <- if (2 * TP + FP + FN > 0) 2 * TP / (2 * TP + FP + FN) else 0
      data.frame(
        population = k, TP = TP, FP = FP, FN = FN,
        precision = precision, recall = recall, F1 = F1
      )
    }))
  }

  group_diagnostics <- list(
    rules = list(
      tost_alpha = tost_alpha,
      tost_equiv_margin = tost_equiv_margin,
      donor_bh_alpha = donor_bh_alpha,
      final_option_priority = final_option_priority,
      equal_n_fallback = equal_n_fallback,
      n_mean = n_mean,
      donor_pool = donor_pool,
      recipient_pool = recipient_pool,
      covariance_ridge = covariance_ridge,
      prior_mode = "global_uniform_1_over_mq",
      global_pi = pi_fix,
      global_a = a,
      global_b = b,
      group_prior_n_c = group_prior_n_c,
      min_fully_supported_groups = min_fully_supported_groups,
      min_fully_supported_option3_groups = min_fully_supported_option3_groups,
      n_fully_supported_option2_groups = n_fully_supported_option2_groups,
      n_fully_supported_option3_groups = n_fully_supported_option3_groups,
      n_fully_supported_groups = n_fully_supported_groups,
      shared_mode_enabled = shared_mode_enabled,
      option1_fallback_when_shared_mode_off = option1_fallback_when_shared_mode_off
    ),
    marginal_beta_hat_array = marginal_beta_hat_array,
    marginal_se_hat_array = marginal_se_hat_array,
    marginal_z_raw_array = marginal_z_raw_array,
    p_marg_fixed_array = p_marg_fixed_array,
    donor_t_p_array = donor_t_p_array,
    donor_bh_p_array = donor_bh_p_array,
    tost_veto_array = tost_veto_array,
    tost_evaluated_array = tost_evaluated_array,
    tost_invalid_array = tost_invalid_array,
    tost_veto_count = sum(tost_veto_array),
    tost_invalid_count = sum(tost_invalid_array),
    option2_candidate = option2_candidate,
    option3_candidate = option3_candidate,
    option2_groups = option2_groups,
    option3_groups = option3_groups,
    protected_snp_list = protected_snp_list,
    protected_active_count_history = protected_active_count_history,
    protected_active_all_history = protected_active_all_history,
    beta_factor_option2 = beta_factor_option2,
    beta_factor_option3 = beta_factor_option3,
    T_tilde_list = T_tilde_list,
    T_current_list = T_current_list,
    Tau_gaussian_list = Tau_gaussian_list,
    trace_correction_list = trace_correction_list,
    delta_gate_array = delta_gate_array,
    lambda_update_array = lambda_update_array,
    W_after_update_array = W_after_update_array,
    performance_by_population = performance_by_population
  )

  list(
    method = "MVP-VB",
    scenario_name = scenario_name,
    config = if (!is.null(common_data$config)) common_data$config else NULL,
    W_list = W_list,
    M_list = M_list,
    Tau_list = Tau_list,
    V_active_list = V_active_list,
    option_array = option_array,
    iter_summary = iter_summary,
    n_iter_done = n_iter_done,
    converged = converged,
    convergence_reason = convergence_reason,
    runtime_sec = runtime_sec,

    ## Record the actual prior used by this core, including the enforced
    ## global pi/a/b values.
    initialization_rules = initialization_rules_used,
    gate_rules = gate_rules,
    prior_parameters = list(
      a = a, b = b, pi = pi_fix,
      pi_option1_by_pop = pi_option1_by_pop,
      group_prior_n_c = group_prior_n_c,
      prior_mode = "global_uniform_1_over_mq"
    ),
    n_k_evidence = n_k_evidence,

    active_idx_list = active_idx_list,
    C_list = C_list,
    Z_watch_list = Z_watch_list,
    active_controls = active_controls_used,
    c_target_k = c_target_k,
    c_z_k = c_z_k,
    c_w_k = c_w_k,
    active_size_history = active_size_history,
    active_idx_history = active_idx_history,
    group_diagnostics = group_diagnostics
  )
}
