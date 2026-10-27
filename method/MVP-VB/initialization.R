get_initialization_rules <- function(
    n_k,
    m,
    q,
    m_global = m,
    seed_fraction = 0.10,
    max_steps = 10L,
    evidence_n_ratio_max = 5,
    direct_weight = 0.9,
    halo_weight = 0.10,
    cross_trait_weight = 0.10,
    background_weight = 0.001,
    ld_threshold = 0.80,
    prior_mass = 2500
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

  if (length(m) != 1L || !is.finite(m) || m < 1) {
    stop("m must be a positive scalar.")
  }

  if (length(m_global) != 1L ||
      !is.finite(m_global) ||
      m_global < m) {
    stop("m_global must be a positive scalar and at least m.")
  }

  if (length(q) != 1L || !is.finite(q) || q < 1) {
    stop("q must be a positive scalar.")
  }

  m <- as.integer(m)
  m_global <- as.integer(m_global)
  q <- as.integer(q)

  ## ------------------------------------------------------------
  ## Evidence-only sample-size balancing
  ##
  ## Evidence-only sample-size uplift.
  ## The largest population remains at its raw n; smaller populations
  ## are lifted to at least max(n_k) / evidence_n_ratio_max.
  ## ------------------------------------------------------------
  n_max <- max(n_k)

  n_k_evidence <- pmax(
    n_k,
    n_max / evidence_n_ratio_max
  )

  z_rescale <- sqrt(n_k_evidence / n_k)

  ## ------------------------------------------------------------
  ## Population-specific, phenotype-specific screening
  ##
  ## Rank marginal p-values and define the top 5% as the eligible
  ## pool. The original stepwise residual screening is then run
  ## only inside this pool.
  ## ------------------------------------------------------------
  n_screen_pool <- max(
    1L,
    as.integer(ceiling(seed_fraction * m))
  )

  n_stepwise_max <- min(
    as.integer(max_steps),
    n_screen_pool
  )

  ## ------------------------------------------------------------
  ## Inclusion-prior hyperparameters
  ##
  ## m_global, not window-specific m, controls multiplicity.
  ## ------------------------------------------------------------
  a <- 1

  pi <- prior_mass / (m_global * q)

  if (pi <= 0 || pi >= 1) {
    stop("The induced prior inclusion probability pi is invalid.")
  }

  b <- a * (1 - pi) / pi

  ## ------------------------------------------------------------
  ## Return all initialization rules
  ## ------------------------------------------------------------
  list(
    ## Dimensions
    K = K,
    q = q,
    m = m,
    m_global = m_global,

    ## Sample-size quantities
    n_k_raw = n_k,
    n_k_evidence = n_k_evidence,
    z_rescale = z_rescale,
    evidence_n_ratio_max = evidence_n_ratio_max,

    ## Screening
    screening_method = "population_response_specific_top_fraction_then_stepwise",
    seed_fraction = seed_fraction,
    n_screen_pool = n_screen_pool,
    max_steps = as.integer(max_steps),
    n_stepwise_max = n_stepwise_max,

    ## Initialization W
    W_background = background_weight,
    W_direct = direct_weight,
    W_ld_halo = halo_weight,
    W_cross_trait = cross_trait_weight,

    ## LD
    ld_threshold = ld_threshold,

    ## Prior
    a = a,
    pi = pi,
    b = b
  )
}

precompute_XtX_XtY <- function(X_list, Y_list) {
  K <- length(X_list)
  XtX_list <- vector("list", K)
  XtY_list <- vector("list", K)

  for (k in seq_len(K)) {
    XtX_list[[k]] <- crossprod(X_list[[k]])
    XtY_list[[k]] <- crossprod(X_list[[k]], Y_list[[k]])
  }

  list(XtX_list = XtX_list, XtY_list = XtY_list)
}

.initialize_population_W <- function(X, Y, rules) {
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  m <- ncol(X)
  q <- ncol(Y)

  Xc <- scale(X, center = TRUE, scale = FALSE)
  Yc <- scale(Y, center = TRUE, scale = FALSE)
  Xc <- as.matrix(Xc)
  Yc <- as.matrix(Yc)

  W <- matrix(rules$W_background, nrow = m, ncol = q)
  W_init_type <- matrix("background", nrow = m, ncol = q)

  x_ss <- colSums(Xc^2)
  for (r in seq_len(q)) {
    y <- Yc[, r]
    y_ss <- sum(y^2)
    score <- abs(as.vector(crossprod(Xc, y))) /
      sqrt(pmax(x_ss * y_ss, .Machine$double.eps))
    score[!is.finite(score)] <- -Inf

    pool_n <- min(m, as.integer(rules$n_screen_pool))
    pool <- order(score, decreasing = TRUE)[seq_len(pool_n)]
    selected <- integer(0)
    residual <- y

    for (step in seq_len(min(as.integer(rules$n_stepwise_max), length(pool)))) {
      remaining <- setdiff(pool, selected)
      if (length(remaining) == 0L) break

      candidate_score <- vapply(
        remaining,
        function(i) abs(sum(Xc[, i] * residual)) /
          sqrt(max(sum(Xc[, i]^2) * sum(residual^2), .Machine$double.eps)),
        numeric(1)
      )
      best <- remaining[which.max(candidate_score)]
      if (!is.finite(candidate_score[which.max(candidate_score)]) ||
          candidate_score[which.max(candidate_score)] <= 0) break

      selected <- c(selected, best)
      fit <- lm.fit(
        x = cbind(1, Xc[, selected, drop = FALSE]),
        y = y
      )
      residual <- fit$residuals
    }

    if (length(selected) > 0L) {
      W[selected, r] <- rules$W_direct
      W_init_type[selected, r] <- "direct"
    }
  }

  ## LD halo and cross-trait initialization are evidence-only starting values.
  x_sd <- sqrt(pmax(x_ss, 0))
  Cor <- crossprod(Xc) / outer(x_sd, x_sd)
  Cor[!is.finite(Cor)] <- 0
  diag(Cor) <- 1

  for (r in seq_len(q)) {
    direct <- which(W[, r] >= rules$W_direct)
    if (length(direct) == 0L) next

    halo <- unique(unlist(lapply(direct, function(i) {
      which(abs(Cor[i, ]) >= rules$ld_threshold)
    }), use.names = FALSE))
    halo <- setdiff(halo, direct)
    if (length(halo) > 0L) {
      use <- halo[W[halo, r] < rules$W_ld_halo]
      W[use, r] <- rules$W_ld_halo
      W_init_type[use, r] <- "ld_halo"
    }

    other_traits <- setdiff(seq_len(q), r)
    for (i in direct) {
      use <- other_traits[W[i, other_traits] < rules$W_cross_trait]
      if (length(use) > 0L) {
        W[i, use] <- rules$W_cross_trait
        W_init_type[i, use] <- "cross_trait"
      }
    }
  }

  list(W = W, W_init_type = W_init_type)
}

prepare_fit_inputs <- function(
    X_list,
    Y_list,
    q,
    K,
    initialization_rules
) {
  if (length(X_list) != K || length(Y_list) != K) {
    stop("X_list and Y_list must both have length K.")
  }
  if (initialization_rules$K != K || initialization_rules$q != q) {
    stop("initialization_rules dimensions do not match K and q.")
  }

  for (k in seq_len(K)) {
    X_list[[k]] <- as.matrix(X_list[[k]])
    Y_list[[k]] <- as.matrix(Y_list[[k]])
    if (nrow(X_list[[k]]) != nrow(Y_list[[k]]) ||
        ncol(X_list[[k]]) != initialization_rules$m ||
        ncol(Y_list[[k]]) != q) {
      stop("X_list[[", k, "]] and Y_list[[", k, "]] have incompatible dimensions.")
    }
  }

  pre <- precompute_XtX_XtY(X_list, Y_list)
  init <- lapply(seq_len(K), function(k) {
    .initialize_population_W(X_list[[k]], Y_list[[k]], initialization_rules)
  })

  c(
    pre,
    list(
      W_list = lapply(init, function(z) z$W),
      W_base_list = lapply(init, function(z) z$W),
      W_init_type_list = lapply(init, function(z) z$W_init_type)
    )
  )
}
