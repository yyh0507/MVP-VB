suppressPackageStartupMessages({
  if (!requireNamespace("MSGLasso", quietly = TRUE)) {
    stop("Package 'MSGLasso' is required. Install the official MSGLasso 2.1 package first.")
  }
})

MSGLasso_fixed <- function(X.m, Y.m, grp.WTs, Pen.L, Pen.G, PQ.grps, GR.grps,
                           grp_Norm0, lam1, lam.G, Beta0 = NULL) {
  X.m <- as.matrix(X.m); Y.m <- as.matrix(Y.m); grp.WTs <- as.matrix(grp.WTs)
  Pen.L <- as.matrix(Pen.L); Pen.G <- as.matrix(Pen.G)
  PQ.grps <- as.matrix(PQ.grps); GR.grps <- as.matrix(GR.grps)
  grp_Norm0 <- as.matrix(grp_Norm0); lam.G <- as.matrix(lam.G)
  N <- nrow(X.m); P <- ncol(X.m); Q <- ncol(Y.m)
  G <- nrow(grp.WTs); R <- ncol(grp.WTs)
  gmax <- ncol(PQ.grps); cmax <- ncol(GR.grps)
  if (gmax < 2L || cmax < 2L) stop("Invalid MSGLasso group metadata.")
  args <- list(as.integer(N), as.integer(P), as.integer(Q), as.integer(G), as.integer(R),
               as.double(as.vector(t(X.m))), as.double(as.vector(t(Y.m))),
               as.double(as.vector(t(grp.WTs))), as.integer(as.vector(t(Pen.L))),
               as.integer(as.vector(t(Pen.G))), as.integer(gmax),
               as.integer(as.vector(t(PQ.grps))), as.integer(cmax),
               as.integer(as.vector(t(GR.grps))), as.double(lam1),
               as.double(as.vector(t(lam.G))))
  if (is.null(Beta0)) {
    z <- .C("MSGLasso", args[[1]], args[[2]], args[[3]], args[[4]], args[[5]],
            args[[6]], args[[7]], args[[8]], args[[9]], args[[10]], args[[11]],
            args[[12]], args[[13]], args[[14]], args[[15]], args[[16]],
            grp_Norm = as.double(as.vector(t(grp_Norm0))),
            Beta.out = as.double(rep(0, P * Q)), n.iter = as.integer(0),
            RSS = as.double(0), Edebug = as.double(rep(0, N * Q)), PACKAGE = "MSGLasso")
  } else {
    if (!identical(dim(Beta0), c(P, Q))) stop("Beta0 has the wrong dimension.")
    z <- .C("MSGLasso_Ini", args[[1]], args[[2]], args[[3]], args[[4]], args[[5]],
            args[[6]], args[[7]], args[[8]], args[[9]], args[[10]], args[[11]],
            args[[12]], args[[13]], args[[14]], args[[15]], args[[16]],
            Beta.ini = as.double(as.vector(t(Beta0))),
            grp_Norm = as.double(as.vector(t(grp_Norm0))),
            Beta.out = as.double(rep(0, P * Q)), n.iter = as.integer(0),
            RSS = as.double(0), Edebug = as.double(rep(0, N * Q)), PACKAGE = "MSGLasso")
  }
  B <- matrix(z$Beta.out, nrow = P, byrow = TRUE)
  E <- matrix(z$Edebug, nrow = N, ncol = Q, byrow = TRUE)
  list(Beta = B, grpNorm = matrix(z$grp_Norm, nrow = G, byrow = TRUE),
       E = E, rss.v = colSums(E^2), rss = sum(E^2), iter = as.integer(z$n.iter))
}

build_groups <- function(p, q) {
  G <- p; R <- q; gmax <- 1L; cmax <- max(p, q)
  gs <- 0:(p - 1L); ge <- gs; rs <- 0:(q - 1L); re <- rs
  pq <- MSGLasso::FindingPQGrps(p, q, G, R, gmax, gs, ge, rs, re)$PQgrps
  gr <- MSGLasso::FindingGRGrps(p, q, G, R, cmax, gs, ge, rs, re)$GRgrps
  wt <- MSGLasso::Cal_grpWTs(p, q, G, R, gmax, pq)$grpWTs
  list(G = G, R = R, PQgrps = pq, GRgrps = gr, grpWTs = wt)
}

standardize_xy <- function(X, Y, tol = 1e-12) {
  X <- as.matrix(X); Y <- as.matrix(Y)
  if (nrow(X) != nrow(Y) || any(!is.finite(X)) || any(!is.finite(Y))) stop("Invalid X/Y.")
  xm <- colMeans(X); ym <- colMeans(Y)
  Xc <- sweep(X, 2L, xm, "-")
  Yc <- sweep(Y, 2L, ym, "-")
  x_norm <- sqrt(colSums(Xc^2)); y_norm <- sqrt(colSums(Yc^2))
  keep <- is.finite(x_norm) & x_norm > tol
  if (!any(keep)) stop("No polymorphic SNP remains.")
  if (any(!is.finite(y_norm)) || any(y_norm <= tol)) stop("Every trait must have positive variance.")
  Xs <- sweep(Xc[, keep, drop = FALSE], 2L, x_norm[keep], "/")
  Ys <- sweep(Yc, 2L, y_norm, "/")
  list(X = Xs, Y = Ys, keep = keep, x_mean = xm[keep],
       x_norm = x_norm[keep], y_mean = ym, y_norm = y_norm)
}

fit_one_population <- function(X, Y, gamma = 0.5, n_lambda = 50L,
                               log10_lambda_min = -15,
                               log10_lambda_max = 5,
                               verbose = FALSE) {
  st <- standardize_xy(X, Y); Xs <- st$X; Ys <- st$Y
  n <- nrow(Xs); p <- ncol(Xs); q <- ncol(Ys)
  g <- build_groups(p, q); penL <- matrix(1L, p, q); penG <- matrix(0L, g$G, g$R)
  lamG <- matrix(0, g$G, g$R); gn0 <- matrix(0, g$G, g$R)
  lambda_max <- max(abs(crossprod(Xs, Ys) / n))
  if (!is.finite(lambda_max) || lambda_max <= 0) stop("Invalid lambda_max.")
  ## EBIC search requested for the comparison: 50 uniformly spaced points
  ## on the absolute log10(lambda) interval [-15, 5].
  if (length(n_lambda) != 1L || n_lambda != 50L ||
      !isTRUE(all.equal(log10_lambda_min, -15)) ||
      !isTRUE(all.equal(log10_lambda_max, 5))) {
    stop("MVLasso EBIC grid must use 50 points on log10(lambda) [-15, 5].")
  }
  lambda_grid <- 10^seq(log10_lambda_min, log10_lambda_max,
                        length.out = n_lambda)
  prev <- NULL; gn <- gn0; best <- NULL; best_ebic <- Inf; best_i <- NA_integer_; rows <- vector("list", n_lambda)
  for (i in seq_along(lambda_grid)) {
    z <- MSGLasso_fixed(Xs, Ys, g$grpWTs, penL, penG, g$PQgrps, g$GRgrps, gn0,
                         lambda_grid[i], lamG, prev)
    B <- z$Beta; d <- sum(B != 0); rss <- pmax(z$rss.v, .Machine$double.xmin)
    ebic <- n * sum(log(rss / n)) + d * log(n) + 2 * gamma * lchoose(p * q, d)
    rows[[i]] <- data.frame(lambda_index = i, lambda = lambda_grid[i], RSS = sum(z$rss.v),
                             d_lambda = d, EBIC = ebic, iterations = z$iter)
    if (ebic < best_ebic) { best_ebic <- ebic; best <- z; best_i <- i }
    prev <- B; gn <- z$grpNorm
  }
  B <- best$Beta; full <- matrix(0, ncol(X), q); full[st$keep, ] <- B
  beta_original <- matrix(0, ncol(X), q)
  beta_original[st$keep, ] <- sweep(sweep(B, 1L, st$x_norm, "/"), 2L, st$y_norm, "*")
  list(beta_std = full, beta_original = beta_original, selected = full != 0,
       best_lambda = lambda_grid[best_i], best_ebic = best_ebic, best_index = best_i,
       gamma = gamma, lambda_max = lambda_max, lambda_grid = lambda_grid,
       path = do.call(rbind, rows), standardization = st)
}

pair_metrics <- function(selected, truth) {
  selected <- as.logical(selected); truth <- as.logical(truth)
  tp <- sum(selected & truth); fp <- sum(selected & !truth); fn <- sum(!selected & truth)
  data.frame(TP = tp, FP = fp, FN = fn,
             Precision = if ((tp + fp) > 0) tp / (tp + fp) else 0,
             Recall = if ((tp + fn) > 0) tp / (tp + fn) else 0,
             F1 = if ((2 * tp + fp + fn) > 0) 2 * tp / (2 * tp + fp + fn) else 0)
}

run_MVLasso_any_scenario <- function(common_data, output_dir, scenario_name = "") {
  req <- c("X_list", "Y_list", "B_true_list")
  if (!all(req %in% names(common_data))) stop("common_data must contain X_list, Y_list, B_true_list.")
  K <- length(common_data$X_list); if (K < 1L || length(common_data$Y_list) != K || length(common_data$B_true_list) != K) stop("Invalid population lists.")
  m <- ncol(as.matrix(common_data$X_list[[1]])); q <- ncol(as.matrix(common_data$Y_list[[1]]))
  fits <- vector("list", K); pair_rows <- vector("list", K); metric_rows <- vector("list", K)
  for (k in seq_len(K)) {
    X <- as.matrix(common_data$X_list[[k]]); Y <- as.matrix(common_data$Y_list[[k]]); truth <- as.matrix(common_data$B_true_list[[k]]) != 0
    if (!identical(dim(X), c(nrow(X), m)) || !identical(dim(Y), c(nrow(X), q)) || !identical(dim(truth), c(m, q))) stop("Dimension mismatch in population ", k)
    fit <- fit_one_population(X, Y); fits[[k]] <- fit
    metric <- pair_metrics(fit$selected, truth); metric$population <- k; metric$scenario_name <- scenario_name; metric_rows[[k]] <- metric
    snp_key <- if (!is.null(common_data$snp_key_vec)) as.character(common_data$snp_key_vec) else paste0("snp", seq_len(m))
    pair_rows[[k]] <- data.frame(population = k, SNP = rep(seq_len(m), q), SNP_key = rep(snp_key, q), trait = rep(seq_len(q), each = m),
                                 beta_hat = as.vector(fit$beta_original), is_selected = as.vector(fit$selected),
                                 B_true = as.vector(common_data$B_true_list[[k]]), is_true = as.vector(truth))
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  pairs <- do.call(rbind, pair_rows); metrics <- do.call(rbind, metric_rows)
  saveRDS(fits, file.path(output_dir, "MVLasso_fits.rds")); write.csv(pairs, file.path(output_dir, "MVLasso_population_SNP_trait_pairs.csv"), row.names = FALSE)
  write.csv(pairs[pairs$is_selected, , drop = FALSE], file.path(output_dir, "MVLasso_selected_SNP_trait_pairs.csv"), row.names = FALSE)
  write.csv(metrics, file.path(output_dir, "MVLasso_population_metrics.csv"), row.names = FALSE)
  invisible(list(fits = fits, pairs = pairs, metrics = metrics))
}
