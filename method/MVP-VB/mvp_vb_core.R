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

  compute_delta_one_k_active <- function(k, Wcur, M_list, V_C, active_idx,
                                         Tau_list, XtX_list, XtY_list,
                                         m, q, other_r_list) {
    delta_mat <- matrix(-Inf, nrow = m, ncol = q)
    c_k <- length(active_idx)
    if (c_k == 0L) return(delta_mat)

    XtX_C <- XtX_list[[k]][active_idx, active_idx, drop = FALSE]
    XtY_C <- XtY_list[[k]][active_idx, , drop = FALSE]
    M_C <- M_list[[k]][active_idx, , drop = FALSE]
    Tau_k <- Tau_list[[k]]
    W_C <- Wcur[[k]][active_idx, , drop = FALSE]
    s_diag_C <- diag(XtX_C)
    delta_C <- matrix(0, nrow = c_k, ncol = q)

    for (r in seq_len(q)) {
      rows_r <- ((r - 1L) * c_k + 1L):(r * c_k)
      Mr <- M_C[, r]
      Wr <- W_C[, r]
      Vrr <- V_C[rows_r, rows_r, drop = FALSE]

      Delta1 <- -0.5 * Tau_k[r, r] * s_diag_C * (Mr^2 + diag(Vrr))
      Delta2 <- Mr * as.numeric(XtY_C %*% Tau_k[, r])

      T3 <- outer(Mr, Mr) + Vrr
      diag(T3) <- 0
      Delta3 <- -Tau_k[r, r] *
        rowSums(XtX_C * sweep(T3, 2, Wr, `*`))

      Delta4 <- numeric(c_k)
      for (s in other_r_list[[r]]) {
        rows_s <- ((s - 1L) * c_k + 1L):(s * c_k)
        Ms <- M_C[, s]
        Ws <- W_C[, s]
        Vrs <- V_C[rows_r, rows_s, drop = FALSE]
        T4 <- outer(Mr, Ms) + Vrs

        Delta4 <- Delta4 -
          Tau_k[r, s] * rowSums(XtX_C * sweep(T4, 2, Ws, `*`))
      }

      delta_C[, r] <- Delta1 + Delta2 + Delta3 + Delta4
    }

    delta_mat[active_idx, ] <- delta_C
    delta_mat
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
      z_vec_all_ir <- marïÎø¶‰žËkºwµçAÐÍ}É½ÕÁ}…ÉÉ…åm¥‘àÍt€ð´±…µ‰‘„Ì4(€€€€€€€‰•Ñ…}™…Ñ½É}½ÁÑ¥½¸ÍmmÉ½ÕÁ}¹…µ•ut€ð´±¥ÍÐ 4(€€€€€€€€€(€ô(Ì°Õ}½±€ôÔÍ}½±°…}½±€ô„Í}½±°‰}½±€ôˆÍ}½±°4(€€€€€€€€€ÁÉ¥½É}„€ô„Í}ÁÉ¥½È°ÁÉ¥½É}ˆ€ôˆÍ}ÁÉ¥½È°4(€€€€€€€€€¹}•™˜€ôÕ¹¹…µ”¡ÁÉ¥½ÈÍl‰¹}•™˜‰t¤°Á¥}É½ÕÀ€ôÕ¹¹…µ”¡ÁÉ¥½ÈÍl‰Á¤‰t¤4(€€€€€€€€¤4(€€€€€ô4(4(€€€€€€ŒŒ€´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´4(€€€€€€ŒŒ€Ô¸Ð±•µ•¹ÐµÝ¥Í”\ÕÁ‘…Ñ”4(€€€€€€ŒŒ€´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´4(€€€€€]}Ñ•µÀ€ð´±…ÁÁ±ä¡Í•Å}±•¸¡,¤°™Õ¹Ñ¥½¸¡¬¤µ…ÑÉ¥à À°´°Ä¤¤4(4(€€€€€™½È€¡¬¥¸Í•Å}±•¸¡,¤¤ì4(€€€€€€€™½È€¡¤¥¸…Ñ¥Ù•}¥‘á}±¥ÍÑmm­ut¤ì4(€€€€€€€€€™½È€¡È¥¸Í•Å}±•¸¡Ä¤¤ì4(€€€€€€€€€€€‘•±Ñ…}¥È€ð´‘•±Ñ…}…Ñ•}…ÉÉ…åm¬°¤°Ét4(€€€€€€€€€€€½ÁÑ¥½¹}­¥È€ð´½ÁÑ¥½¹}…ÉÉ…åm¬°¤°Ét4(4(€€€€€€€€€€€¥˜€¡½ÁÑ¥½¹}­¥È€ôô€Å0¤ì4(€€€€€€€€€€€€€±…µ‰‘…}¥È€ð´±…µ‰‘…}½ÁÐÄ4(€€€€€€€€€€€ô•±Í”¥˜€¡½ÁÑ¥½¹}­¥È€ôô€É0¤ì4(€€€€€€€€€€€€€±…µ‰‘…}¥È€ð´±…µ‰‘…}½ÁÐÉ}É½ÕÁ}…ÉÉ…åm¬°¤°Ét4(€€€€€€€€€€€ô•±Í”¥˜€¡½ÁÑ¥½¹}­¥È€ôô€Í0¤ì4(€€€€€€€€€€€€€±…µ‰‘…}¥È€ð´±…µ‰‘…}½ÁÐÍ}É½ÕÁ}…ÉÉ…åm¬°¤°Ét4(€€€€€€€€€€€ô•±Í”ì4(€€€€€€€€€€€€€ÍÑ½À ‰U¹­¹½Ý¸½ÁÑ¥½¹}­¥È¸ˆ¤4(€€€€€€€€€€€ô4(4(€€€€€€€€€€€]}¹•Ü€ð´Á±½¥Ì¡±…µ‰‘…}¥È€¬‘•±Ñ…}¥È¤4(€€€€€€€€€€€±…µ‰‘…}ÕÁ‘…Ñ•}…ÉÉ…åm¬°¤°Ét€ð´±…µ‰‘…}¥È4(€€€€€€€€€€€]}…™Ñ•É}ÕÁ‘…Ñ•}…ÉÉ…åm¬°¤°Ét€ð´]}¹•Ü4(€€€€€€€€€€€]}Ñ•µÁmm­uum¤°Ét€ð´]}¹•Ü4(€€€€€€€€€ô4(€€€€€€€ô4(€€€€€ô4(€€€€€]}±¥ÍÐ€ð´]}Ñ•µÀ4(4(€€€€€€ŒŒ€´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´4(€€€€€€ŒŒ€Ô¸Ô1	<4(€€€€€€ŒŒ€´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´4(€€€€€•±‰½}½¹Ñ¥¹Õ½ÕÌ€ð´€À4(€€€€€½Ù…É¥…¹•}É¥‘”€ð´€Å”´Ø4(€€€€€±½‘•Ñ}AÍ¤À€ð´€È€¨ÍÕ´¡±½œ¡‘¥…œ¡¡½°¡AÍ¤À¤¤¤¤4(€€€€€™½È€¡¬¥¸Í•Å}±•¸¡,¤¤ì4(€€€€€€€¥‘à€ð´…Ñ¥Ù•}¥‘á}±¥ÍÑmm­ut4(€€€€€€€]}}¹•Ü€ð´]}±¥ÍÑmm­uum¥‘à°€°‘É½À€ô1Mt4(€€€€€€€5}€ð´5}±¥ÍÑmm­uum¥‘à°€°‘É½À€ô1Mt4(€€€€€€€Y}€ð´Y}…Ñ¥Ù•}±¥ÍÑmm­ut4(€€€€€€€aÑa}€ð´aÑa}±¥ÍÑmm­uum¥‘à°¥‘à°‘É½À€ô1Mt4(€€€€€€€aÑe}€ð´aÑe}±¥ÍÑmm­uum¥‘à°€°‘É½À€ô1Mt4(4(€€€€€€€M}ÕÉÉ•¹Ð€ð´½µÁÕÑ•}M}…Ñ¥Ù” 4(€€€€€€€€€5}€ô5}°4(€€€€€€€€€]}€ô]}}¹•Ü°4(€€€€€€€€€aÑa}€ôaÑa}°4(€€€€€€€€€Y}€ôY}°4(€€€€€€€€€Q…Õ}¬€ôQ…Õ}…ÕÍÍ¥…¹}±¥ÍÑmm­ut°4(€€€€€€€€€Í¥µ…È€ôÍ¥µ…È°4(€€€€€€€€€´€ô´°4(€€€€€€€€€Ä€ôÄ4(€€€€€€€€¤4(€€€€€€€	¡…Ñ}}¹•Ü€ð´]}}¹•Ü€¨5}4(€€€€€€€a	¡…Ñ}}¹•Ü€ð´a}±¥ÍÑmm­uul°¥‘à°‘É½À€ô1Mt€”¨”	¡…Ñ}}¹•Ü4(€€€€€€€Q}ÕÉÉ•¹Ñ}¬€ð´AÍ¤À€¬eÑe}±¥ÍÑmm­ut€´4(€€€€€€€€€Ð¡	¡…Ñ}}¹•Ü¤€”¨”aÑe}€´4(€€€€€€€€€É½ÍÍÁÉ½¡e}±¥ÍÑmm­ut°a	¡…Ñ}}¹•Ü¤€¬4(€€€€€€€€€M}ÕÉÉ•¹Ð‘}µ…Ð€¬€ Ä€¼Í¥µ…È¤€¨M}ÕÉÉ•¹Ð‘M}µ…Ð4(€€€€€€€¥˜€¡…¹ä …¥Ì¹™¥¹¥Ñ”¡Q}ÕÉÉ•¹Ñ}¬¤¤¤ì4(€€€€€€€€€ÍÑ½À ‰9½¸µ™¥¹¥Ñ”Ñ¡•½É•Ñ¥…°½Ù…É¥…¹”Í…±”¥¸Á½ÁÕ±…Ñ¥½¸€ˆ°¬°€ˆ¸ˆ¤4(€€€€€€€ô4(€€€€€€€Q}ÕÉÉ•¹Ñ}±¥ÍÑmm­ut€ð´Q}ÕÉÉ•¹Ñ}¬4(4(€€€€€€€‘}¬€ð´¹ÔÀ€¬¹}­m­t€¬´4(€€€€€€€AÍ¥}ÍÑ…É}¬€ð´Q…Õ}±¥ÍÑmm­ut4(€€€€€€€ÑÉ…•}½ÉÉ•Ñ¥½¹}±¥ÍÑm­t€ð´€À¸Ô€¨ÍÕ´¡‘¥…œ 4(€€€€€€€€€€¡Q}Ñ¥±‘•}±¥ÍÑmm­ut€´Q}ÕÉÉ•¹Ñ}¬¤€”¨”AÍ¥}ÍÑ…É}¬4(€€€€€€€€¤¤4(€€€€€€€•±‰½}½¹Ñ¥¹Õ½ÕÌ€ð´•±‰½}½¹Ñ¥¹Õ½ÕÌ€¬4(€€€€€€€€€€´¡¹}­m­t€¨Ä€¼€È¤€¨±½œ È€¨Á¤¤€´4(€€€€€€€€€€¡´€¨Ä€¼€È¤€¨±½œ¡Í¥µ…È¤€¬4(€€€€€€€€€€¡´€¨Ä€¼€È¤€¬4(€€€€€€€€€€À¸Ô€¨±½‘•Ñ}Y}±¥ÍÑmm­ut€¬4(€€€€€€€€€€¡¹ÔÀ€¼€È¤€¨±½‘•Ñ}AÍ¤À€¬4(€€€€€€€€€€ ¡¹}­m­t€¬´¤€¨Ä€¼€È¤€¨±½œ È¤€¬4(€€€€€€€€€±½}µÕ±Ñ¥Ù…É¥…Ñ•}…µµ„¡‘}¬€¼€È°Ä¤€´4(€€€€€€€€€±½}µÕ±Ñ¥Ù…É¥…Ñ•}…µµ„¡¹ÔÀ€¼€È°Ä¤€´4(€€€€€€€€€€¡‘}¬€¼€È¤€¨±½‘•Ñ}AÍ¥}±¥ÍÑmm­ut€¬4(€€€€€€€€€ÑÉ…•}½ÉÉ•Ñ¥½¹}±¥ÍÑm­t4(€€€€€ô4(4(€€€€€€ŒŒ	•É¹½Õ±±¤•¹ÑÉ½Áä¸4(€€€€€]}¹•Ý}…ÉÉ…ä€ð´…ÉÉ…ä À°‘¥´€ôŒ¡,°´°Ä¤¤4(€€€€€™½È€¡¬¥¸Í•Å}±•¸¡,¤¤]}¹•Ý}…ÉÉ…åm¬°€°t€ð´]}±¥ÍÑmm­ut4(€€€€€•¹ÑÉ½Áå}µ…Í¬€ð´]}¹•Ý}…ÉÉ…ä€ø€À€˜]}¹•Ý}…ÉÉ…ä€ð€Ä4(€€€€€•±‰½}•¹ÑÉ½Áä€ð´€µÍÕ´ 4(€€€€€€€]}¹•Ý}…ÉÉ…åm•¹ÑÉ½Áå}µ…Í­t€¨±½œ¡]}¹•Ý}…ÉÉ…åm•¹ÑÉ½Áå}µ…Í­t¤€¬4(€€€€€€€€€€ Ä€´]}¹•Ý}…ÉÉ…åm•¹ÑÉ½Áå}µ…Í­t¤€¨4(€€€€€€€€€±½œÅÀ µ]}¹•Ý}…ÉÉ…åm•¹ÑÉ½Áå}µ…Í­t¤4(€€€€€€¤4(4(€€€€€€ŒŒ=ÁÑ¥½¸€ÄÁÉ¥½È…ÁÁ±¥•Ì½¹±äÑ¼™¥¹…°=ÁÑ¥½¸€Ä¥¹‘¥…Ñ½ÉÌ¸4(€€€€€½ÁÑ¥½¸Å}µ…Í¬€ð´½ÁÑ¥½¹}…ÉÉ…ä€ôô€Å04(€€€€€Ý}½ÁÑ¥½¸Ä€ð´]}¹•Ý}…ÉÉ…åm½ÁÑ¥½¸Å}µ…Í­t4(€€€€€•±‰½}½ÁÑ¥½¸Ä€ð´ÍÕ´ 4(€€€€€€€Ý}½ÁÑ¥½¸Ä€¨±½œ¡Á¥}™¥à¤€¬4(€€€€€€€€€€ Ä€´Ý}½ÁÑ¥½¸Ä¤€¨±½œÅÀ µÁ¥}™¥à¤4(€€€€€€¤4(4(€€€€€€ŒŒMå¹¡É½¹¥é”•Ù•ÉäÍ¡…É•	•Ñ„™…Ñ½ÈÑ¼Ñ¡”¹•Ü\°Ñ¡•¸ÕÍ”Ñ¡”4(€€€€€€ŒŒµ…Ñ¡¥¹œ…ÑÕ…°É½ÕÀ™½È‰½Ñ Ñ¡”Í¥µÁ±¥™¥•…¹•¹•É…°™½ÉµÕ±…Ì¸4(€€€€€•±‰½}½ÁÑ¥½¸È€ð´€À4(€€€€€•±‰½}½ÁÑ¥½¸Ì€ð´€À4(€€€€€‰•Ñ…}Íå¹}•ÉÉ½ÉÌ€ð´¹Õµ•É¥Œ À¤4(€€€€€™½È€¡É½ÕÁ}¹…µ”¥¸¹…µ•Ì¡½ÁÑ¥½¸É}É½ÕÁÌ¤¤ì4(€€€€€€€É½ÕÀ€ð´½ÁÑ¥½¸É}É½ÕÁÍmmÉ½ÕÁ}¹…µ•ut4(€€€€€€€(È€ð´±•¹Ñ ¡É½ÕÀ‘Á½ÁÕ±…Ñ¥½¹Ì¤4(€€€€€€€ÔÉ}¹•Ü€ð´ÍÕ´¡]}¹•Ý}…ÉÉ…åmÉ½ÕÀ‘Á½ÁÕ±…Ñ¥½¹Ì°É½ÕÀ‘Í¹À°É½ÕÀ‘ÑÉ…¥Ñt¤4(€€€€€€€„É}ÁÉ¥½È€ð´‰•Ñ…}™…Ñ½É}½ÁÑ¥½¸ÉmmÉ½ÕÁ}¹…µ•ut‘ÁÉ¥½É}„4(€€€€€€€ˆÉ}ÁÉ¥½È€ð´‰•Ñ…}™…Ñ½É}½ÁÑ¥½¸ÉmmÉ½ÕÁ}¹…µ•ut‘ÁÉ¥½É}ˆ4(€€€€€€€„É}¹•Ü€ð´„É}ÁÉ¥½È€¬ÔÉ}¹•Ü4(€€€€€€€ˆÉ}¹•Ü€ð´ˆÉ}ÁÉ¥½È€¬(È€´ÔÉ}¹•Ü4(€€€€€€€Í¥µÁ±”È€ð´±‰•Ñ„¡„É}¹•Ü°ˆÉ}¹•Ü¤€´±‰•Ñ„¡„É}ÁÉ¥½È°ˆÉ}ÁÉ¥½È¤4(€€€€€€€•¹•É…°È€ð´±‰•Ñ„¡„É}¹•Ü°ˆÉ}¹•Ü¤€´±‰•Ñ„¡„É}ÁÉ¥½È°ˆÉ}ÁÉ¥½È¤€¬4(€€€€€€€€€€¡„É}ÁÉ¥½È€¬ÔÉ}¹•Ü€´„É}¹•Ü¤€¨€¡‘¥…µµ„¡„É}¹•Ü¤€´‘¥…µµ„¡„É}¹•Ü€¬ˆÉ}¹•Ü¤¤€¬4(€€€€€€€€€€¡ˆÉ}ÁÉ¥½È€¬(È€´ÔÉ}¹•Ü€´ˆÉ}¹•Ü¤€¨4(€€€€€€€€€€¡‘¥…µµ„¡ˆÉ}¹•Ü¤€´‘¥…µµ„¡„É}¹•Ü€¬ˆÉ}¹•Ü¤¤4(€€€€€€€‰•Ñ…}Íå¹}•ÉÉ½ÉÌ€ð´Œ¡‰•Ñ…}Íå¹}•ÉÉ½ÉÌ°…‰Ì¡Í¥µÁ±”È€´•¹•É…°È¤¤4(€€€€€€€•±‰½}½ÁÑ¥½¸È€ð´•±‰½}½ÁÑ¥½¸È€¬Í¥µÁ±”È4(€€€€€€€‰•Ñ…}™…Ñ½É}½ÁÑ¥½¸ÉmmÉ½ÕÁ}¹…µ•ut‘Õ}¹•Ü€ð´ÔÉ}¹•Ü4(€€€€€€€‰•Ñ…}™…Ñ½É}½ÁÑ¥½¸ÉmmÉ½ÕÁ}¹…µ•ut‘…}¹•Ü€ð´„É}¹•Ü4(€€€€€€€‰•Ñ…}™…Ñ½É}½ÁÑ¥½¸ÉmmÉ½ÕÁ}¹…µ•ut‘‰}¹•Ü€ð´ˆÉ}¹•Ü4(€€€€€ô4(€€€€€™½È€¡É½ÕÁ}¹…µ”¥¸¹…µ•Ì¡½ÁÑ¥½¸Í}É½ÕÁÌ¤¤ì4(€€€€€€€É½ÕÀ€ð´½ÁÑ¥½¸Í}É½ÕÁÍmmÉ½ÕÁ}¹…µ•ut4(€€€€€€€Á…¥ÉÌ€ð´É½ÕÀ‘Á…¥ÉÌ4(€€€€€€€¥‘àÌ€ð´‰¥¹¡Á…¥ÉÍl°€‰Á½ÁÕ±…Ñ¥½¸‰t°É½ÕÀ‘Í¹À°Á…¥ÉÍl°€‰ÑÉ…¥Ð‰t¤4(€€€€€€€(Ì€ð´¹É½Ü¡¥‘àÌ¤4(€€€€€€€ÔÍ}¹•Ü€ð´ÍÕ´¡]}¹•Ý}…ÉÉ…åm¥‘àÍt¤4(€€€€€€€„Í}ÁÉ¥½È€ð´‰•Ñ…}™…Ñ½É}½ÁÑ¥½¸ÍmmÉ½ÕÁ}¹…µ•ut‘ÁÉ¥½É}„4(€€€€€€€ˆÍ}ÁÉ¥½È€ð´‰•Ñ…}™…Ñ½É}½ÁÑ¥½¸ÍmmÉ½ÕÁ}¹…µ•ut‘ÁÉ¥½É}ˆ4(€€€€€€€„Í}¹•Ü€ð´„Í}ÁÉ¥½È€¬ÔÍ}¹•Ü4(€€€€€€€ˆÍ}¹•Ü€ð´ˆÍ}ÁÉ¥½È€¬(Ì€´ÔÍ}¹•Ü4(€€€€€€€Í¥µÁ±”Ì€ð´±‰•Ñ„¡„Í}¹•Ü°ˆÍ}¹•Ü¤€´±‰•Ñ„¡„Í}ÁÉ¥½È°ˆÍ}ÁÉ¥½È¤4(€€€€€€€•¹•É…°Ì€ð´±‰•Ñ„¡„Í}¹•Ü°ˆÍ}¹•Ü¤€´±‰•Ñ„¡„Í}ÁÉ¥½È°ˆÍ}ÁÉ¥½È¤€¬4(€€€€€€€€€€¡„Í}ÁÉ¥½È€¬ÔÍ}¹•Ü€´„Í}¹•Ü¤€¨€¡‘¥…µµ„¡„Í}¹•Ü¤€´‘¥…µµ„¡„Í}¹•Ü€¬ˆÍ}¹•Ü¤¤€¬4(€€€€€€€€€€¡ˆÍ}ÁÉ¥½È€¬(Ì€´ÔÍ}¹•Ü€´ˆÍ}¹•Ü¤€¨4(€€€€€€€€€€¡‘¥…µµ„¡ˆÍ}¹•Ü¤€´‘¥…µµ„¡„Í}¹•Ü€¬ˆÍ}¹•Ü¤¤4(€€€€€€€‰•Ñ…}Íå¹}•ÉÉ½ÉÌ€ð´Œ¡‰•Ñ…}Íå¹}•ÉÉ½ÉÌ°…‰Ì¡Í¥µÁ±”Ì€´•¹•É…°Ì¤¤4(€€€€€€€•±‰½}½ÁÑ¥½¸Ì€ð´•±‰½}½ÁÑ¥½¸Ì€¬Í¥µÁ±”Ì4(€€€€€€€‰•Ñ…}™…Ñ½É}½ÁÑ¥½¸ÍmmÉ½ÕÁ}¹…µ•ut‘Õ}¹•Ü€ð´ÔÍ}¹•Ü4(€€€€€€€‰•Ñ…}™…Ñ½É}½ÁÑ¥½¸ÍmmÉ½ÕÁ}¹…µ•ut‘…}¹•Ü€ð´„Í}¹•Ü4(€€€€€€€‰•Ñ…}™…Ñ½É}½ÁÑ¥½¸ÍmmÉ½ÕÁ}¹…µ•ut‘‰}¹•Ü€ð´ˆÍ}¹•Ü4(€€€€€ô4(€€€€€‰•Ñ…}Íå¹}µ…á}•ÉÉ½È€ð´¥˜€¡±•¹Ñ ¡‰•Ñ…}Íå¹}•ÉÉ½ÉÌ¤¤ì4(€€€€€€€µ…à¡‰•Ñ…}Íå¹}•ÉÉ½ÉÌ¤4(€€€€€ô•±Í”€À4(€€€€€¥˜€ …¥Ì¹™¥¹¥Ñ”¡‰•Ñ…}Íå¹}µ…á}•ÉÉ½È¤ñð‰•Ñ…}Íå¹}µ…á}•ÉÉ½È€ø€Å”´ÄÀ¤ì4(€€€€€€€ÍÑ½À ‰Må¹¡É½¹¥é•	•Ñ„µ™…Ñ½È1	<¡•¬™…¥±•¸ˆ¤4(€€€€€ô4(4(€€€€€•±‰½}ÕÉÈ€ð´•±‰½}½¹Ñ¥¹Õ½ÕÌ€¬•±‰½}•¹ÑÉ½Áä€¬•±‰½}½ÁÑ¥½¸Ä€¬4(€€€€€€€•±‰½}½ÁÑ¥½¸È€¬•±‰½}½ÁÑ¥½¸Ì4(4(€€€€€€ŒŒ€´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´4(€€€€€€ŒŒ€Ô¸Ø¥…¹½ÍÑ¥Ì…¹½¹Ù•É•¹”4(€€€€€€ŒŒ€´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´4(€€€€€½ÁÑ¥½¸Å}½Õ¹Ñ}¹½Ü€ð´ÍÕ´¡½ÁÑ¥½¹}…ÉÉ…ä€ôô€Å0¤4(€€€€€½ÁÑ¥½¸É}½Õ¹Ñ}¹½Ü€ð´ÍÕ´¡½ÁÑ¥½¹}…ÉÉ…ä€ôô€É0¤4(€€€€€½ÁÑ¥½¸Í}½Õ¹Ñ}¹½Ü€ð´ÍÕ´¡½ÁÑ¥½¹}…ÉÉ…ä€ôô€Í0¤4(4(€€€€€½ÁÑ¥½¹}¡…¹•}¹½Ü€ð´€Á04(4(€€€€€½ÁÐÉ}µ…Í¬€ð´½ÁÑ¥½¹}…ÉÉ…ä€ôô€É04(€€€€€µ•…¹}‘•±Ñ…}…Ñ•}½ÁÐÈ€ð´¥˜€¡…¹ä¡½ÁÐÉ}µ…Í¬¤¤ì4(€€€€€€€µ•…¸¡‘•±Ñ…}…Ñ•}…ÉÉ…åm½ÁÐÉ}µ…Í­t°¹„¹É´€ôQIU¤4(€€€€€ô•±Í”ì4(€€€€€€€9}É•…±|4(€€€€€ô4(€€€€€µ•…¹}±…µ‰‘…}½ÁÐÈ€ð´¥˜€¡…¹ä¡½ÁÐÉ}µ…Í¬¤¤ì4(€€€€€€€µ•…¸¡±…µ‰‘…}ÕÁ‘…Ñ•}…ÉÉ…åm½ÁÐÉ}µ…Í­t°¹„¹É´€ôQIU¤4(€€€€€ô•±Í”ì4(€€€€€€€9}É•…±|4(€€€€€ô4(€€€€€µ•…¹}]}…™Ñ•É}½ÁÐÈ€ð´¥˜€¡…¹ä¡½ÁÐÉ}µ…Í¬¤¤ì4(€€€€€€€µ•…¸¡]}…™Ñ•É}ÕÁ‘…Ñ•}…ÉÉ…åm½ÁÐÉ}µ…Í­t°¹„¹É´€ôQIU¤4(€€€€€ô•±Í”ì4(€€€€€€€9}É•…±|4(€€€€€ô4(4(€€€€€¥Ñ•É}ÍÕµµ…Éä‘1	=m¥Ñt€ð´•±‰½}ÕÉÈ4(€€€€€¥Ñ•É}ÍÕµµ…Éä‘1	=}¡…¹•m¥Ñt€ð´¥˜€¡¥Ð€ø€Å0¤•±‰½}ÕÉÈ€´•±‰½}ÁÉ•Ø•±Í”9}É•…±|4(€€€€€¥Ñ•É}ÍÕµµ…Éä‘1	=}½¹Ñ¥¹Õ½ÕÍm¥Ñt€ð´•±‰½}½¹Ñ¥¹Õ½ÕÌ4(€€€€€¥Ñ•É}ÍÕµµ…Éä‘1	=}•¹ÑÉ½Áåm¥Ñt€ð´•±‰½}•¹ÑÉ½Áä4(€€€€€¥Ñ•É}ÍÕµµ…Éä‘1	=}½ÁÑ¥½¸Åm¥Ñt€ð´•±‰½}½ÁÑ¥½¸Ä4(€€€€€¥Ñ•É}ÍÕµµ…Éä‘1	=}½ÁÑ¥½¸Ém¥Ñt€ð´•±‰½}½ÁÑ¥½¸È4(€€€€€¥Ñ•É}ÍÕµµ…Éä‘1	=}½ÁÑ¥½¸Ím¥Ñt€ð´•±‰½}½ÁÑ¥½¸Ì4(€€€€€¥Ñ•É}ÍÕµµ…Éä‘1	=}ÑÉ…•}½ÉÉ•Ñ¥½¹m¥Ñt€ð´ÍÕ´¡ÑÉ…•}½ÉÉ•Ñ¥½¹}±¥ÍÐ¤4(€€€€€¥Ñ•É}ÍÕµµ…Éä‘‰•Ñ…}Íå¹}µ…á}•ÉÉ½Ém¥Ñt€ð´‰•Ñ…}Íå¹}µ…á}•ÉÉ½È4(€€€€€¥Ñ•É}ÍÕµµ…Éä‘½ÁÑ¥½¸Å}½Õ¹Ñm¥Ñt€ð´½ÁÑ¥½¸Å}½Õ¹Ñ}¹½Ü4(€€€€€¥Ñ•É}ÍÕµµ…Éä‘½ÁÑ¥½¸É}½Õ¹Ñm¥Ñt€ð´½ÁÑ¥½¸É}½Õ¹Ñ}¹½Ü4(€€€€€¥Ñ•É}ÍÕµµ…Éä‘½ÁÑ¥½¸Í}½Õ¹Ñm¥Ñt€ð´½ÁÑ¥½¸Í}½Õ¹Ñ}¹½Ü4(€€€€€¥Ñ•É}ÍÕµµ…Éä‘½ÁÑ¥½¹}¡…¹•m¥Ñt€ð´½ÁÑ¥½¹}¡…¹•}¹½Ü4(€€€€€¥Ñ•É}ÍÕµµ…Éä‘¹}½ÁÑ¥½¸É}½­m¥Ñt€ð´ÍÕ´¡½ÁÑ¥½¸É}…¹‘¥‘…Ñ”°¹„¹É´€ôQIU¤4(€€€€€¥Ñ•É}ÍÕµµ…Éä‘¹}½ÁÑ¥½¸Í}½­m¥Ñt€ð´ÍÕ´¡½ÁÑ¥½¸Í}…¹‘¥‘…Ñ”°¹„¹É´€ôQIU¤4(€€€€€¥Ñ•É}ÍÕµµ…Éä‘µ•…¹}‘•±Ñ…}…Ñ•}½ÁÐÉm¥Ñt€ð´µ•…¹}‘•±Ñ…}…Ñ•}½ÁÐÈ4(€€€€€¥Ñ•É}ÍÕµµ…Éä‘µ•…¹}±…µ‰‘…}½ÁÐÉm¥Ñt€ð´µ•…¹}±…µ‰‘…}½ÁÐÈ4(€€€€€¥Ñ•É}ÍÕµµ…Éä‘µ•…¹}]}…™Ñ•É}½ÁÐÉm¥Ñt€ð´µ•…¹}]}…™Ñ•É}½ÁÐÈ4(4(€€€€€…Ð¡ÍÁÉ¥¹Ñ˜ 4(€€€€€€€Á…ÍÑ”À 4(€€€€€€€€€€‰%Ñ•È€•è1	<ô”¸Ý˜ð½ÁÐÄô•ð½ÁÐÈô•ð½ÁÐÌô•ð€ˆ°4(€€€€€€€€€€‰½ÁÑ¥½¹}¡…¹”ô•Ìðµ•…¹}‘•±Ñ…}½ÁÐÈô”¸Í˜ð€ˆ°4(€€€€€€€€€€‰µ•…¹}±…µ‰‘…}½ÁÐÈô”¸Í˜ðµ•…¹}]}…™Ñ•É}½ÁÐÈô”¸Í™q¸ˆ4(€€€€€€€€¤°4(€€€€€€€¥Ð°4(€€€€€€€•±‰½}ÕÉÈ°4(€€€€€€€½ÁÑ¥½¸Å}½Õ¹Ñ}¹½Ü°4(€€€€€€€½ÁÑ¥½¸É}½Õ¹Ñ}¹½Ü°4(€€€€€€€½ÁÑ¥½¸Í}½Õ¹Ñ}¹½Ü°4(€€€€€€€¥™•±Í”¡¥Ì¹¹„¡½ÁÑ¥½¹}¡…¹•}¹½Ü¤°€‰9ˆ°…Ì¹¡…É…Ñ•È¡½ÁÑ¥½¹}¡…¹•}¹½Ü¤¤°4(€€€€€€€µ•…¹}‘•±Ñ…}…Ñ•}½ÁÐÈ°4(€€€€€€€µ•…¹}±…µ‰‘…}½ÁÐÈ°4(€€€€€€€µ•…¹}]}…™Ñ•É}½ÁÐÈ4(€€€€€€¤¤4(4(€€€€€¥˜€¡¥Ð€ø€Å0¤ì4(€€€€€€€•±‰½}‘¥™˜€ð´…‰Ì¡•±‰½}ÕÉÈ€´•±‰½}ÁÉ•Ø¤4(€€€€€€€¥˜€¡¥Ì¹™¥¹¥Ñ”¡•±‰½}‘¥™˜¤€˜˜•±‰½}‘¥™˜€ðÑ½±}•±‰¼¤ì4(€€€€€€€€€½¹Ù•É•€ð´QIU4(€€€€€€€€€½¹Ù•É•¹•}É•…Í½¸€ð´€‰•±‰¼ˆ4(€€€€€€€€€…Ð¡ÍÁÉ¥¹Ñ˜ 4(€€€€€€€€€€€€‰½¹Ù•É•…Ð¥Ñ•È€•ð•±Ñ„1	<€ô€”¸Í•q¸ˆ°4(€€€€€€€€€€€¥Ð°4(€€€€€€€€€€€•±‰½}‘¥™˜4(€€€€€€€€€€¤¤4(€€€€€€€€€‰É•…¬4(€€€€€€€ô4(€€€€€ô4(4(€€€€€•±‰½}ÁÉ•Ø€ð´•±‰½}ÕÉÈ4(€€€ô4(4(€€€¹}¥Ñ•É}‘½¹”€ð´¥Ð4(€ô¥ml‰•±…ÁÍ•‰ut4(4(€€ŒŒ€ôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôô4(€€ŒŒ€Ø¸I•ÑÕÉ¸½É”½ÕÑÁÕÑÌ4(€€ŒŒ€ôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôôô4(4(€¥Ñ•É}ÍÕµµ…Éä€ð´¥Ñ•É}ÍÕµµ…ÉåmÍ•Å}±•¸¡¹}¥Ñ•É}‘½¹”¤°€°‘É½À€ô1Mt4(€…Ñ¥Ù•}Í¥é•}¡¥ÍÑ½Éä€ð´…Ñ¥Ù•}Í¥é•}¡¥ÍÑ½ÉåmÍ•Å}±•¸¡¹}¥Ñ•É}‘½¹”¤°€°‘É½À€ô1Mt4(€ÁÉ½Ñ•Ñ•‘}…Ñ¥Ù•}½Õ¹Ñ}¡¥ÍÑ½Éä€ð´4(€€€ÁÉ½Ñ•Ñ•‘}…Ñ¥Ù•}½Õ¹Ñ}¡¥ÍÑ½ÉåmÍ•Å}±•¸¡¹}¥Ñ•É}‘½¹”¤°€°‘É½À€ô1Mt4(€ÁÉ½Ñ•Ñ•‘}…Ñ¥Ù•}…±±}¡¥ÍÑ½Éä€ð´4(€€€ÁÉ½Ñ•Ñ•‘}…Ñ¥Ù•}…±±}¡¥ÍÑ½ÉåmÍ•Å}±•¸¡¹}¥Ñ•É}‘½¹”¤°€°‘É½À€ô1Mt4(€¥˜€¡ÍÑ½É•}…Ñ¥Ù•}¡¥ÍÑ½Éä¤ì4(€€€…Ñ¥Ù•}¥‘á}¡¥ÍÑ½Éä€ð´…Ñ¥Ù•}¥‘á}¡¥ÍÑ½ÉåmÍ•Å}±•¸¡¹}¥Ñ•É}‘½¹”¥t4(€ô4(4(€Á•É™½Éµ…¹•}‰å}Á½ÁÕ±…Ñ¥½¸€ð´9U104(€¥˜€ …¥Ì¹¹Õ±°¡½µµ½¹}‘…Ñ„‘	}ÑÉÕ•}±¥ÍÐ¤€˜˜4(€€€€€±•¹Ñ ¡½µµ½¹}‘…Ñ„‘	}ÑÉÕ•}±¥ÍÐ¤€ôô,¤ì4(€€€Á•É™½Éµ…¹•}‰å}Á½ÁÕ±…Ñ¥½¸€ð´‘¼¹…±°¡É‰¥¹°±…ÁÁ±ä¡Í•Å}±•¸¡,¤°™Õ¹Ñ¥½¸¡¬¤ì4(€€€€€ÑÉÕÑ €ð´…Ì¹µ…ÑÉ¥à¡½µµ½¹}‘…Ñ„‘	}ÑÉÕ•}±¥ÍÑmm­ut¤€„ô€À4(€€€€€ÁÉ•€ð´]}±¥ÍÑmm­ut€ø€À¸Ô4(€€€€€¥˜€ …¥‘•¹Ñ¥…°¡‘¥´¡ÑÉÕÑ ¤°Œ¡´°Ä¤¤¤ì4(€€€€€€€ÍÑ½À ‰	}ÑÉÕ•}±¥ÍÐ‘¥µ•¹Í¥½¹Ì‘¼¹½Ð…É•”Ý¥Ñ ´…¹Ä¥¸Á½ÁÕ±…Ñ¥½¸€ˆ°¬°€ˆ¸ˆ¤4(€€€€€ô4(€€€€€Q@€ð´ÍÕ´¡ÁÉ•€˜ÑÉÕÑ ¤4(€€€€€@€ð´ÍÕ´¡ÁÉ•€˜€…ÑÉÕÑ ¤4(€€€€€8€ð´ÍÕ´ …ÁÉ•€˜ÑÉÕÑ ¤4(€€€€€ÁÉ•¥Í¥½¸€ð´¥˜€¡Q@€¬@€ø€À¤Q@€¼€¡Q@€¬@¤•±Í”€À4(€€€€€É•…±°€ð´¥˜€¡Q@€¬8€ø€À¤Q@€¼€¡Q@€¬8¤•±Í”€À4(€€€€€Ä€ð´¥˜€ È€¨Q@€¬@€¬8€ø€À¤€È€¨Q@€¼€ È€¨Q@€¬@€¬8¤•±Í”€À4(€€€€€‘…Ñ„¹™É…µ” 4(€€€€€€€Á½ÁÕ±…Ñ¥½¸€ô¬°Q@€ôQ@°@€ô@°8€ô8°4(€€€€€€€ÁÉ•¥Í¥½¸€ôÁÉ•¥Í¥½¸°É•…±°€ôÉ•…±°°Ä€ôÄ4(€€€€€€¤4(€€€ô¤¤4(€ô4(4(€É½ÕÁ}‘¥…¹½ÍÑ¥Ì€ð´±¥ÍÐ 4(€€€ÉÕ±•Ì€ô±¥ÍÐ 4(€€€€€Ñ½ÍÑ}…±Á¡„€ôÑ½ÍÑ}…±Á¡„°4(€€€€€Ñ½ÍÑ}•ÅÕ¥Ù}µ…É¥¸€ôÑ½ÍÑ}•ÅÕ¥Ù}µ…É¥¸°4(€€€€€‘½¹½É}‰¡}…±Á¡„€ô‘½¹½É}‰¡}…±Á¡„°4(€€€€€™¥¹…±}½ÁÑ¥½¹}ÁÉ¥½É¥Ñä€ô™¥¹…±}½ÁÑ¥½¹}ÁÉ¥½É¥Ñä°4(€€€€€•ÅÕ…±}¹}™…±±‰…¬€ô•ÅÕ…±}¹}™…±±‰…¬°4(€€€€€¹}µ•…¸€ô¹}µ•…¸°4(€€€€€‘½¹½É}Á½½°€ô‘½¹½É}Á½½°°4(€€€€€É•¥Á¥•¹Ñ}Á½½°€ôÉ•¥Á¥•¹Ñ}Á½½°°4(€€€€€½Ù…É¥…¹•}É¥‘”€ô½Ù…É¥…¹•}É¥‘”°4(€€€€€ÁÉ¥½É}µ½‘”€ô€‰±½‰…±}Õ¹¥™½Éµ|Å}½Ù•É}µÄˆ°4(€€€€€±½‰…±}Á¤€ôÁ¥}™¥à°4(€€€€€±½‰…±}„€ô„°4(€€€€€±½‰…±}ˆ€ôˆ°4(€€€€€É½ÕÁ}ÁÉ¥½É}¹}Œ€ôÉ½ÕÁ}ÁÉ¥½É}¹}Œ°4(€€€€€µ¥¹}™Õ±±å}ÍÕÁÁ½ÉÑ•‘}É½ÕÁÌ€ôµ¥¹}™Õ±±å}ÍÕÁÁ½ÉÑ•‘}É½ÕÁÌ°4(€€€€€µ¥¹}™Õ±±å}ÍÕÁÁ½ÉÑ•‘}½ÁÑ¥½¸Í}É½ÕÁÌ€ôµ¥¹}™Õ±±å}ÍÕÁÁ½ÉÑ•‘}½ÁÑ¥½¸Í}É½ÕÁÌ°4(€€€€€¹}™Õ±±å}ÍÕÁÁ½ÉÑ•‘}½ÁÑ¥½¸É}É½ÕÁÌ€ô¹}™Õ±±å}ÍÕÁÁ½ÉÑ•‘}½ÁÑ¥½¸É}É½ÕÁÌ°4(€€€€€¹}™Õ±±å}ÍÕÁÁ½ÉÑ•‘}½ÁÑ¥½¸Í}É½ÕÁÌ€ô¹}™Õ±±å}ÍÕÁÁ½ÉÑ•‘}½ÁÑ¥½¸Í}É½ÕÁÌ°4(€€€€€¹}™Õ±±å}ÍÕÁÁ½ÉÑ•‘}É½ÕÁÌ€ô¹}™Õ±±å}ÍÕÁÁ½ÉÑ•‘}É½ÕÁÌ°4(€€€€€Í¡…É•‘}µ½‘•}•¹…‰±•€ôÍ¡…É•‘}µ½‘•}•¹…‰±•°4(€€€€€½ÁÑ¥½¸Å}™…±±‰…­}Ý¡•¹}Í¡…É•‘}µ½‘•}½™˜€ô½ÁÑ¥½¸Å}™…±±‰…­}Ý¡•¹}Í¡…É•‘}µ½‘•}½™˜4(€€€€¤°4(€€€µ…É¥¹…±}‰•Ñ…}¡…Ñ}…ÉÉ…ä€ôµ…É¥¹…±}‰•Ñ…}¡…Ñ}…ÉÉ…ä°4(€€€µ…É¥¹…±}Í•}¡…Ñ}…ÉÉ…ä€ôµ…É¥¹…±}Í•}¡…Ñ}…ÉÉ…ä°4(€€€µ…É¥¹…±}é}É…Ý}…ÉÉ…ä€ôµ…É¥¹…±}é}É…Ý}…ÉÉ…ä°4(€€€Á}µ…É}™¥á•‘}…ÉÉ…ä€ôÁ}µ…É}™¥á•‘}…ÉÉ…ä°4(€€€‘½¹½É}Ñ}Á}…ÉÉ…ä€ô‘½¹½É}Ñ}Á}…ÉÉ…ä°4(€€€‘½¹½É}‰¡}Á}…ÉÉ…ä€ô‘½¹½É}‰¡}Á}…ÉÉ…ä°4(€€€Ñ½ÍÑ}Ù•Ñ½}…ÉÉ…ä€ôÑ½ÍÑ}Ù•Ñ½}…ÉÉ…ä°4(€€€Ñ½ÍÑ}•Ù…±Õ…Ñ•‘}…ÉÉ…ä€ôÑ½ÍÑ}•Ù…±Õ…Ñ•‘}…ÉÉ…ä°4(€€€Ñ½ÍÑ}¥¹Ù…±¥‘}…ÉÉ…ä€ôÑ½ÍÑ}¥¹Ù…±¥‘}…ÉÉ…ä°4(€€€Ñ½ÍÑ}Ù•Ñ½}½Õ¹Ð€ôÍÕ´¡Ñ½ÍÑ}Ù•Ñ½}…ÉÉ…ä¤°4(€€€Ñ½ÍÑ}¥¹Ù…±¥‘}½Õ¹Ð€ôÍÕ´¡Ñ½ÍÑ}¥¹Ù…±¥‘}…ÉÉ…ä¤°4(€€€½ÁÑ¥½¸É}…¹‘¥‘…Ñ”€ô½ÁÑ¥½¸É}…¹‘¥‘…Ñ”°4(€€€½ÁÑ¥½¸Í}…¹‘¥‘…Ñ”€ô½ÁÑ¥½¸Í}…¹‘¥‘…Ñ”°4(€€€½ÁÑ¥½¸É}É½ÕÁÌ€ô½ÁÑ¥½¸É}É½ÕÁÌ°4(€€€½ÁÑ¥½¸Í}É½ÕÁÌ€ô½ÁÑ¥½¸Í}É½ÕÁÌ°4(€€€ÁÉ½Ñ•Ñ•‘}Í¹Á}±¥ÍÐ€ôÁÉ½Ñ•Ñ•‘}Í¹Á}±¥ÍÐ°4(€€€ÁÉ½Ñ•Ñ•‘}…Ñ¥Ù•}½Õ¹Ñ}¡¥ÍÑ½Éä€ôÁÉ½Ñ•Ñ•‘}…Ñ¥Ù•}½Õ¹Ñ}¡¥ÍÑ½Éä°4(€€€ÁÉ½Ñ•Ñ•‘}…Ñ¥Ù•}…±±}¡¥ÍÑ½Éä€ôÁÉ½Ñ•Ñ•‘}…Ñ¥Ù•}…±±}¡¥ÍÑ½Éä°4(€€€‰•Ñ…}™…Ñ½É}½ÁÑ¥½¸È€ô‰•Ñ…}™…Ñ½É}½ÁÑ¥½¸È°4(€€€‰•Ñ…}™…Ñ½É}½ÁÑ¥½¸Ì€ô‰•Ñ…}™…Ñ½É}½ÁÑ¥½¸Ì°4(€€€Q}Ñ¥±‘•}±¥ÍÐ€ôQ}Ñ¥±‘•}±¥ÍÐ°4(€€€Q}ÕÉÉ•¹Ñ}±¥ÍÐ€ôQ}ÕÉÉ•¹Ñ}±¥ÍÐ°4(€€€Q…Õ}…ÕÍÍ¥…¹}±¥ÍÐ€ôQ…Õ}…ÕÍÍ¥…¹}±¥ÍÐ°4(€€€ÑÉ…•}½ÉÉ•Ñ¥½¹}±¥ÍÐ€ôÑÉ…•}½ÉÉ•Ñ¥½¹}±¥ÍÐ°4(€€€‘•±Ñ…}…Ñ•}…ÉÉ…ä€ô‘•±Ñ…}…Ñ•}…ÉÉ…ä°4(€€€±…µ‰‘…}ÕÁ‘…Ñ•}…ÉÉ…ä€ô±…µ‰‘…}ÕÁ‘…Ñ•}…ÉÉ…ä°4(€€€]}…™Ñ•É}ÕÁ‘…Ñ•}…ÉÉ…ä€ô]}…™Ñ•É}ÕÁ‘…Ñ•}…ÉÉ…ä°4(€€€Á•É™½Éµ…¹•}‰å}Á½ÁÕ±…Ñ¥½¸€ôÁ•É™½Éµ…¹•}‰å}Á½ÁÕ±…Ñ¥½¸4(€€¤4(4(€±¥ÍÐ 4(€€€µ•Ñ¡½€ô€‰5Y@µYˆ°4(€€€Í•¹…É¥½}¹…µ”€ôÍ•¹…É¥½}¹…µ”°4(€€€½¹™¥œ€ô¥˜€ …¥Ì¹¹Õ±°¡½µµ½¹}‘…Ñ„‘½¹™¥œ¤¤½µµ½¹}‘…Ñ„‘½¹™¥œ•±Í”9U10°4(€€€]}±¥ÍÐ€ô]}±¥ÍÐ°4(€€€5}±¥ÍÐ€ô5}±¥ÍÐ°4(€€€Q…Õ}±¥ÍÐ€ôQ…Õ}±¥ÍÐ°4(€€€Y}…Ñ¥Ù•}±¥ÍÐ€ôY}…Ñ¥Ù•}±¥ÍÐ°4(€€€½ÁÑ¥½¹}…ÉÉ…ä€ô½ÁÑ¥½¹}…ÉÉ…ä°4(€€€¥Ñ•É}ÍÕµµ…Éä€ô¥Ñ•É}ÍÕµµ…Éä°4(€€€¹}¥Ñ•É}‘½¹”€ô¹}¥Ñ•É}‘½¹”°4(€€€½¹Ù•É•€ô½¹Ù•É•°4(€€€½¹Ù•É•¹•}É•…Í½¸€ô½¹Ù•É•¹•}É•…Í½¸°4(€€€ÉÕ¹Ñ¥µ•}Í•Œ€ôÉÕ¹Ñ¥µ•}Í•Œ°4(4(€€€€ŒŒI•½ÉÑ¡”…ÑÕ…°ÁÉ¥½ÈÕÍ•‰äÑ¡¥Ì½É”°¥¹±Õ‘¥¹œÑ¡”•¹™½É•4(€€€€ŒŒ±½‰…°Á¤½„½ˆÙ…±Õ•Ì¸4(€€€¥¹¥Ñ¥…±¥é…Ñ¥½¹}ÉÕ±•Ì€ô¥¹¥Ñ¥…±¥é…Ñ¥½¹}ÉÕ±•Í}ÕÍ•°4(€€€…Ñ•}ÉÕ±•Ì€ô…Ñ•}ÉÕ±•Ì°4(€€€ÁÉ¥½É}Á…É…µ•Ñ•ÉÌ€ô±¥ÍÐ 4(€€€€€„€ô„°ˆ€ôˆ°Á¤€ôÁ¥}™¥à°4(€€€€€Á¥}½ÁÑ¥½¸Å}‰å}Á½À€ôÁ¥}½ÁÑ¥½¸Å}‰å}Á½À°4(€€€€€É½ÕÁ}ÁÉ¥½É}¹}Œ€ôÉ½ÕÁ}ÁÉ¥½É}¹}Œ°4(€€€€€ÁÉ¥½É}µ½‘”€ô€‰±½‰…±}Õ¹¥™½Éµ|Å}½Ù•É}µÄˆ4(€€€€¤°4(€€€¹}­}•Ù¥‘•¹”€ô¹}­}•Ù¥‘•¹”°4(4(€€€…Ñ¥Ù•}¥‘á}±¥ÍÐ€ô…Ñ¥Ù•}¥‘á}±¥ÍÐ°4(€€€}±¥ÍÐ€ô}±¥ÍÐ°4(€€€i}Ý…Ñ¡}±¥ÍÐ€ôi}Ý…Ñ¡}±¥ÍÐ°4(€€€…Ñ¥Ù•}½¹ÑÉ½±Ì€ô…Ñ¥Ù•}½¹ÑÉ½±Í}ÕÍ•°4(€€€}Ñ…É•Ñ}¬€ô}Ñ…É•Ñ}¬°4(€€€}é}¬€ô}é}¬°4(€€€}Ý}¬€ô}Ý}¬°4(€€€…Ñ¥Ù•}Í¥é•}¡¥ÍÑ½Éä€ô…Ñ¥Ù•}Í¥é•}¡¥ÍÑ½Éä°4(€€€…Ñ¥Ù•}¥‘á}¡¥ÍÑ½Éä€ô…Ñ¥Ù•}¥‘á}¡¥ÍÑ½Éä°4(€€€É½ÕÁ}‘¥…¹½ÍÑ¥Ì€ôÉ½ÕÁ}‘¥…¹½ÍÑ¥Ì4(€€¤4)ô4(