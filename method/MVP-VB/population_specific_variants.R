## ============================================================
## Population-wise Option-1-only MVP-VB wrappers
##
## Required functions (source before using this file):
##   get_initialization_rules()
##   prepare_fit_inputs()
##   get_gate_rules()
##   fit_MVP_VB_core()
##
## Preprocessing contract:
##   X is already column-standardized by the fixed-X generator.
##   Y is already column-centered by data_generate.R.
##   These wrappers do not scale X or Y and do not scale Y variance.
## ============================================================

.MVP_VB_required_functions <- c(
  "get_initialization_rules",
  "prepare_fit_inputs",
  "get_gate_rules",
  "fit_MVP_VB_core"
)

.check_MVP_VB_dependencies <- function() {
  missing <- .MVP_VB_required_functions[
    !vapply(.MVP_VB_required_functions, exists, logical(1), mode = "function")
  ]
  if (length(missing) > 0L) {
    stop(
      "Source initialization.R, option_gate.R, and mvp_vb_core.R first. ",
      "Missing functions: ", paste(missing, collapse = ", ")
    )
  }
}

.subset_common_data_MVP_VB <- function(common_data, population, phenotype = NULL) {
  population <- as.integer(population)
  if (length(population) != 1L || population < 1L || population > common_data$K) {
    stop("population is out of range.")
  }

  out <- common_data
  out$K <- 1L
  out$n_k <- as.numeric(common_data$n_k[population])
  out$X_list <- list(as.matrix(common_data$X_list[[population]]))
  out$Y_list <- list(as.matrix(common_data$Y_list[[population]]))
  out$B_true_list <- list(as.matrix(common_data$B_true_list[[population]]))

  if (!is.null(common_data$E_list)) {
    out$E_list <- list(as.matrix(common_data$E_list[[population]]))
  }
  if (!is.null(common_data$Sigma_true_list)) {
    out$Sigma_true_list <- list(common_data$Sigma_true_list[[population]])
  }
  if (!is.null(common_data$target_h2)) {
    out$target_h2 <- common_data$target_h2[population]
  }
  if (!is.null(common_data$h2_population)) {
    out$h2_population <- common_data$h2_population[population]
  }

  if (!is.null(phenotype)) {
    phenotype <- as.integer(phenotype)
    if (length(phenotype) != 1L || phenotype < 1L || phenotype > common_data$q) {
      stop("phenotype is out of range.")
    }
    out$q <- 1L
    out$Y_list <- list(out$Y_list[[1L]][, phenotype, drop = FALSE])
    out$B_true_list <- list(out$B_true_list[[1L]][, phenotype, drop = FALSE])
    if (!is.null(out$E_list)) {
      out$E_list <- list(out$E_list[[1L]][, phenotype, drop = FALSE])
    }
    if (!is.null(out$Sigma_true_list)) {
      out$Sigma_true_list <- list(
        matrix(out$Sigma_true_list[[1L]][phenotype, phenotype], 1L, 1L)
      )
    }
  } else {
    out$q <- as.integer(common_data$q)
  }

  if (!is.null(out$config)) {
    out$config$K <- 1L
    out$config$q <- out$q
    out$config$n_k <- out$n_k
  }
  out$scenario_name <- paste0(
    common_data$scenario_name,
    "_population", population,
    if (is.null(phenotype)) "" else paste0("_phenotype", phenotype)
  )
  out
}

.make_MVP_VB_option1_inputs <- function(common_one, initialization_args = list()) {
  base_args <- list(
    n_k = common_one$n_k,
    m = common_one$m,
    q = common_one$q,
    m_global = common_one$m
  )
  rules <- do.call(
    get_initialization_rules,
    utils::modifyList(base_args, initialization_args)
  )
  prepared <- prepare_fit_inputs(
    X_list = common_one$X_list,
    Y_list = common_one$Y_list,
    q = common_one$q,
    K = common_one$K,
    initialization_rules = rules
  )

  gate <- get_gate_rules(n_k = common_one$n_k, q = common_one$q)
  gate$option2_enabled <- FALSE
  gate$option3_enabled <- FALSE
  gate$option_priority <- 1L

  list(rules = rules, prepared = prepared, gate = gate)
}

.fit_MVP_VB_option1_one <- function(
    common_one,
    initialization_args = list(),
    max_iter = 200L,
    tol_elbo = 1,
    active_controls = NULL,
    sigmaB2 = 5,
    verbose = TRUE
) {
  inputs <- .make_MVP_VB_option1_inputs(common_one, initialization_args)
  if (isTRUE(verbose)) {
    return(fit_MVP_VB_core(
      common_data = common_one,
      prepared = inputs$prepared,
      initialization_rules = inputs$rules,
      gate_rules = inputs$gate,
      scenario_name = common_one$scenario_name,
      max_iter = max_iter,
      tol_elbo = tol_elbo,
      active_controls = active_controls,
      sigmaB2 = sigmaB2
    ))
  }

  fit <- NULL
  invisible(capture.output({
    fit <- fit_MVP_VB_core(
      common_data = common_one,
      prepared = inputs$prepared,
      initialization_rules = inputs$rules,
      gate_rules = inputs$gate,
      scenario_name = common_one$scenario_name,
      max_iter = max_iter,
      tol_elbo = tol_elbo,
      active_controls = active_controls,
      sigmaB2 = sigmaB2
    )
  }))
  fit
}

fit_MVP_VB_pop <- function(
    common_data,
    initialization_args = list(),
    max_iter = 200L,
    tol_elbo = 1,
    active_controls = NULL,
    sigmaB2 = 5,
    verbose = TRUE
) {
  .check_MVP_VB_dependencies()
  K <- as.integer(common_data$K)
  q <- as.integer(common_data$q)
  m <- as.integer(common_data$m)

  fits <- vector("list", K)
  elapsed <- system.time({
    for (k in seq_len(K)) {
      one <- .subset_common_data_MVP_VB(common_data, population = k)
      fits[[k]] <- .fit_MVP_VB_option1_one(
        one,
        initialization_args = initialization_args,
        max_iter = max_iter,
        tol_elbo = tol_elbo,
        active_controls = active_controls,
        sigmaB2 = sigmaB2,
        verbose = verbose
      )
      if (any(fits[[k]]$option_array != 1L)) {
        stop("MVP-VB-pop produced a non-Option-1 entry in population ", k, ".")
      }
    }
  })[["elapsed"]]

  list(
    method = "MVP-VB-pop",
    scenario_name = common_data$scenario_name,
    W_list = lapply(fits, function(x) x$W_list[[1L]]),
    Zhat_list = lapply(fits, function(x) x$W_list[[1L]] > 0.5),
    fits_by_population = fits,
    option_array = array(1L, dim = c(K, m, q)),
    runtime_sec = elapsed
  )
}

fit_MVP_VB_pop_pheno <- function(
    common_data,
    initialization_args = list(),
    max_iter = 200L,
    tol_elbo = 1,
    active_controls = NULL,
    sigmaB2 = 5,
    verbose = TRUE
) {
  .check_MVP_VB_dependencies()
  K <- as.integer(common_data$K)
  q <- as.integer(common_data$q)
  m <- as.integer(common_data$m)

  fits <- lapply(seq_len(K), function(k) vector("list", q))
  W_list <- lapply(seq_len(K), function(k) matrix(0, m, q))

  elapsed <- system.time({
    for (k in seq_len(K)) {
      for (r in seq_len(q)) {
        one <- .subset_common_data_MVP_VB(
          common_data,
          population = k,
          phenotype = r
        )
        fit_kr <- .fit_MVP_VB_option1_one(
          one,
          initialization_args = initialization_args,
          max_iter = max_iter,
          tol_elbo = tol_elbo,
          active_controls = active_controls,
          sigmaB2 = sigmaB2,
          verbose = verbose
        )
        if (any(fit_kr$option_array != 1L)) {
          stop(
            "MVP-VB-pop-pheno produced a non-Option-1 entry in population ",
            k, ", phenotype ", r, "."
          )
        }
        fits[[k]][[r]] <- fit_kr
        W_list[[k]][, r] <- fit_kr$W_list[[1L]][, 1L]
      }
    }
  })[["elapsed"]]

  list(
    method = "MVP-VB-pop-pheno",
    scenario_name = common_data$scenario_name,
    W_list = W_list,
    Zhat_list = lapply(W_list, function(W) W > 0.5),
    fits_by_population_and_phenotype = fits,
    option_array = array(1L, dim = c(K, m, q)),
    runtime_sec = elapsed
  )
}
