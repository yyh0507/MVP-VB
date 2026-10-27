## ============================================================
## MGflashfm wrapper for simulation data with arbitrary positive K, m and q
##
## Important implementation details:
##   1. MGflashfm receives raw 0/1/2 genotype dosages.
##   2. Marginal betas, LD, EAF and cov(Y) are computed population-wise.
##   3. Expensive population-level FLASHFM/JAM results are checkpointed.
##   4. The official MGFLASHFMwithJAM post-processing bug for q > K is
##      avoided by using a K x q sample-size matrix in the final MAF step.
##   5. The returned object keeps the same main components used by the
##      existing pipeline: result, mgMPP, Zhat, selected_pairs, metrics.
## ============================================================

make_mgflashfm_inputs <- function(common_data) {

  required <- c(
    "Y_list",
    "B_true_list",
    "snp_key_vec",
    "K",
    "m",
    "q",
    "n_k",
    "config"
  )

  missing_fields <- setdiff(required, names(common_data))
  if (length(missing_fields) > 0L) {
    stop(
      "common_data is missing: ",
      paste(missing_fields, collapse = ", ")
    )
  }

  ## Raw 0/1/2 genotype must be used by MGflashfm.
  ## In the current data generator X_list is standardized.
  if (!is.null(common_data$X_raw_list)) {
    X_raw_list <- common_data$X_raw_list
  } else if (!is.null(common_data$X_dosage_list)) {
    X_raw_list <- common_data$X_dosage_list
  } else {
    stop(
      "MGflashfm requires raw 0/1/2 genotype matrices in ",
      "common_data$X_raw_list or common_data$X_dosage_list."
    )
  }

  K <- as.integer(common_data$K)
  m <- as.integer(common_data$m)
  q <- as.integer(common_data$q)
  n_k <- as.integer(common_data$n_k)

  if (length(K) != 1L || is.na(K) || K < 1L) {
    stop("K must be one positive integer; observed K=", K)
  }
  if (length(m) != 1L || is.na(m) || m < 1L) {
    stop("m must be one positive integer; observed m=", m)
  }
  if (length(q) != 1L || is.na(q) || q < 1L) {
    stop("q must be one positive integer; observed q=", q)
  }
  if (length(n_k) != K) {
    stop("length(n_k) must equal K.")
  }
  if (length(common_data$snp_key_vec) != m) {
    stop("length(snp_key_vec) must equal m.")
  }

  if (
    length(X_raw_list) != K ||
    length(common_data$Y_list) != K ||
    length(common_data$B_true_list) != K
  ) {
    stop(
      "X_raw_list, Y_list, and B_true_list must each contain K populations."
    )
  }

  truth_masks <- lapply(
    common_data$B_true_list,
    function(B) as.matrix(B != 0)
  )

  for (k in seq_len(K)) {
    if (!identical(dim(truth_masks[[k]]), c(m, q))) {
      stop("B_true dimension mismatch in population ", k, ".")
    }
    if (any(colSums(truth_masks[[k]]) < 1L)) {
      stop(
        "Population ", k,
        " contains a phenotype with no causal SNP."
      )
    }
  }

  ## Fixed-width IDs prevent substring collisions in older flashfm code.
  snp_ids <- sprintf("snp%04d", seq_len(m))
  trait_names <- paste0("Trait", seq_len(q))
  group_names <- paste0("Group", seq_len(K))

  gwas.list <- vector("list", K)
  LD.list <- vector("list", K)
  covY.list <- vector("list", K)
  Nall <- vector("list", K)

  for (k in seq_len(K)) {

    X <- as.matrix(X_raw_list[[k]])
    Y <- as.matrix(common_data$Y_list[[k]])
    B <- as.matrix(common_data$B_true_list[[k]])

    if (
      !identical(dim(X), c(n_k[k], m)) ||
      !identical(dim(Y), c(n_k[k], q)) ||
      !identical(dim(B), c(m, q))
    ) {
      stop("Dimension mismatch in Group", k, ".")
    }

    if (any(!is.finite(X)) || any(!is.finite(Y))) {
      stop("Non-finite X or Y value in Group", k, ".")
    }

    if (
      any(X < 0 | X > 2) ||
      any(abs(X - round(X)) > 1e-8)
    ) {
      stop(
        "Group", k,
        " raw X is not coded as HapGen2 0/1/2 allele dosage."
      )
    }

    ## Marginal beta from centered raw dosage and centered Y.
    Xc <- scale(X, center = TRUE, scale = FALSE)
    Yc <- scale(Y, center = TRUE, scale = FALSE)
    xss <- colSums(Xc^2)

    if (any(!is.finite(xss)) || any(xss <= 0)) {
      stop("Group", k, " contains a monomorphic SNP.")
    }

    marginal_beta <- sweep(
      crossprod(Xc, Yc),
      1L,
      xss,
      "/"
    )

    EAF <- colMeans(X) / 2
    if (any(!is.finite(EAF)) || any(EAF <= 0 | EAF >= 1)) {
      stop("Group", k, " has an invalid effect-allele frequency.")
    }

    trait_gwas <- lapply(
      seq_len(q),
      function(r) {
        out <- data.frame(
          rsID = snp_ids,
          beta = as.numeric(marginal_beta[, r]),
          EAF = as.numeric(EAF),
          stringsAsFactors = FALSE
        )
        rownames(out) <- out$rsID
        out
      }
    )
    names(trait_gwas) <- trait_names
    gwas.list[[k]] <- trait_gwas

    LD <- cor(X)
    if (any(!is.finite(LD))) {
      stop("Non-finite LD in Group", k, ".")
    }
    rownames(LD) <- colnames(LD) <- snp_ids
    LD.list[[k]] <- LD

    covY <- cov(Y)
    if (any(!is.finite(covY))) {
      stop("Non-finite phenotype covariance in Group", k, ".")
    }
    rownames(covY) <- colnames(covY) <- trait_names
    covY.list[[k]] <- covY

    Nall[[k]] <- setNames(rep(nrow(X), q), trait_names)
  }

  names(gwas.list) <- group_names
  names(LD.list) <- group_names
  names(covY.list) <- group_names
  names(Nall) <- group_names

  list(
    gwas.list = gwas.list,
    LD.list = LD.list,
    covY.list = covY.list,
    Nall = Nall,
    snp_map = data.frame(
      SNP = seq_len(m),
      mgflashfm_rsID = snp_ids,
      SNP_key = as.character(common_data$snp_key_vec),
      stringsAsFactors = FALSE
    ),
    trait_names = trait_names,
    group_names = group_names,
    truth_masks = truth_masks
  )
}

## ------------------------------------------------------------
## Preserve the legal null model "1" in MGflashfm's correlation pruning.
##
## MGflashfm::JAMmulti2() uses "1" as the no-SNP model and its downstream
## tagexpand.mod()/multibeta() functions explicitly support it.  The package
## version used here mistakenly removes "1" because it is not a row name of
## corX; an all-null JAM result then leaves expmods with zero columns and
## fails at expmods[, 2].  This session-local namespace patch changes only
## that case.  All non-null pruning operations reproduce the package code.
## ------------------------------------------------------------
patch_mgflashfm_null_model <- function() {

  if (!requireNamespace("MGflashfm", quietly = TRUE)) {
    stop("Package MGflashfm is not installed.")
  }

  ns <- asNamespace("MGflashfm")
  current <- get("prune_models_by_corr", envir = ns, inherits = FALSE)

  required_formals <- c("snpmods", "corX", "beta", "thr")
  if (!identical(names(formals(current)), required_formals)) {
    stop(
      "Unsupported MGflashfm::prune_models_by_corr() signature: ",
      paste(names(formals(current)), collapse = ", "),
      "."
    )
  }

  probe_cor <- matrix(
    1,
    nrow = 1L,
    ncol = 1L,
    dimnames = list("snp0001", "snp0001")
  )
  probe_beta <- c(snp0001 = 1)
  probe <- current(
    snpmods = "1",
    corX = probe_cor,
    beta = probe_beta,
    thr = 0.8
  )

  probe_empty <- current(
    snpmods = character(0),
    corX = probe_cor,
    beta = probe_beta,
    thr = 0.8
  )

  ## A newer package version may already contain both corrections.
  if (
    identical(as.character(probe), "1") &&
    identical(as.character(probe_empty), "1")
  ) {
    return(invisible(FALSE))
  }

  fixed <- function(snpmods, corX, beta, thr) {
    if (is.null(names(beta))) stop("beta must be a named vector")

    ## JAM.tries.maxcv() can legitimately return no retained non-null
    ## model for a trait with essentially no signal.  JAMmulti2() then
    ## constructs a zero-row binout matrix and apply(..., 1, mod.fn)
    ## returns character(0).  MGflashfm's legal representation of that
    ## null-only result is the model marker "1".
    if (length(snpmods) == 0L) {
      return("1")
    }

    pruned <- character(length(snpmods))

    for (i in seq_along(snpmods)) {
      ## "1" is flashfm/MGflashfm's legal null-model marker.
      if (identical(as.character(snpmods[i]), "1")) {
        pruned[i] <- "1"
        next
      }

      snps <- unlist(strsplit(snpmods[i], "%", fixed = TRUE))
      snps <- snps[snps %in% rownames(corX)]

      if (length(snps) <= 1L) {
        pruned[i] <- paste(snps, collapse = "%")
        next
      }

      repeat {
        subcor <- abs(corX[snps, snps, drop = FALSE])
        diag(subcor) <- 0
        max_r <- max(subcor, na.rm = TRUE)
        if (is.na(max_r) || max_r < thr) break

        pair_idx <- which(subcor == max_r, arr.ind = TRUE)[1L, , drop = TRUE]
        s1 <- snps[pair_idx[1L]]
        s2 <- snps[pair_idx[2L]]
        b1 <- if (!is.na(beta[s1])) abs(beta[s1]) else 0
        b2 <- if (!is.na(beta[s2])) abs(beta[s2]) else 0
        drop_snp <- if (b1 < b2) s1 else s2
        snps <- setdiff(snps, drop_snp)
        if (length(snps) <= 1L) break
      }

      pruned[i] <- paste(snps, collapse = "%")
    }

    pruned_models <- unique(pruned)
    pruned_models <- pruned_models[nzchar(pruned_models)]
    pruned_models
  }

  environment(fixed) <- ns
  assignInNamespace(
    x = "prune_models_by_corr",
    value = fixed,
    ns = "MGflashfm"
  )

  patched <- get(
    "prune_models_by_corr",
    envir = asNamespace("MGflashfm"),
    inherits = FALSE
  )
  check <- patched(
    snpmods = "1",
    corX = probe_cor,
    beta = probe_beta,
    thr = 0.8
  )
  if (!identical(as.character(check), "1")) {
    stop("Failed to preserve MGflashfm's explicit null model.")
  }

  check_empty <- patched(
    snpmods = character(0),
    corX = probe_cor,
    beta = probe_beta,
    thr = 0.8
  )
  if (!identical(as.character(check_empty), "1")) {
    stop("Failed to convert an empty JAM model set to the null model.")
  }

  check_non_null <- patched(
    snpmods = "snp0001",
    corX = probe_cor,
    beta = probe_beta,
    thr = 0.8
  )
  if (!identical(as.character(check_non_null), "snp0001")) {
    stop("MGflashfm null-model patch changed an ordinary non-null model.")
  }

  message(
    "MGflashfm null-model patch applied: explicit and empty null models supported."
  )
  invisible(TRUE)
}


## ------------------------------------------------------------
## Remove exact duplicate canonical models returned by JAMmulti2().
##
## tagexpand.mod()/multibeta() can map two explored tag-model paths to the
## same canonical model. JAMmulti2() otherwise keeps both copies. Downstream
## flashfmU() then loses the duplicated model names when it constructs a data
## frame; makeSNPgroups2U()/PP2snpmod() consequently interprets automatic row
## names such as "2" and "3" as SNP IDs and groupmultiU() fails while indexing
## the LD matrix. Intervene only for exact duplicate copies: their posterior
## probabilities and mbeta objects must agree. Ambiguous duplicates stop with
## an explicit error instead of being merged heuristically.
## ------------------------------------------------------------
patch_mgflashfm_duplicate_canonical_models <- function() {
  if (!requireNamespace("MGflashfm", quietly = TRUE)) {
    stop("Package MGflashfm is not installed.")
  }

  ns <- asNamespace("MGflashfm")
  current <- get("JAMmulti2", envir = ns, inherits = FALSE)
  required_formals <- c(
    "gwas.list", "corX", "ybar", "Vy", "N", "r2", "save.path",
    "maxcv", "maxcv_stop", "maxcv_autocheck", "jam.nM.iter",
    "extra.java.arguments"
  )
  if (!identical(names(formals(current)), required_formals)) {
    stop(
      "Unsupported MGflashfm::JAMmulti2() signature: ",
      paste(names(formals(current)), collapse = ", "),
      "."
    )
  }
  if (isTRUE(attr(current, "mvpvb_duplicate_canonical_model_patch"))) {
    return(invisible(FALSE))
  }

  original <- current
  PP2snpmod_function <- get("PP2snpmod", envir = ns, inherits = TRUE)

  repair_output <- function(out) {
    if (
      !is.list(out) ||
      !is.list(out$SM) ||
      !is.list(out$mbeta) ||
      length(out$SM) != length(out$mbeta)
    ) {
      stop("MGflashfm JAMmulti2() returned an invalid model structure.")
    }

    repaired_traits <- integer()

    for (i in seq_along(out$SM)) {
      sm <- out$SM[[i]]
      if (!methods::is(sm, "snpmod")) {
        stop("MGflashfm JAMmulti2() returned a non-snpmod object for trait ", i, ".")
      }
      models <- sm@models
      if (
        !is.data.frame(models) ||
        !all(c("str", "PP") %in% names(models))
      ) {
        stop("MGflashfm JAMmulti2() returned an invalid model table for trait ", i, ".")
      }

      model_names <- as.character(models$str)
      duplicated_names <- unique(model_names[
        duplicated(model_names) | duplicated(model_names, fromLast = TRUE)
      ])
      if (!length(duplicated_names)) next

      mbeta <- out$mbeta[[i]]
      mbeta_names <- names(mbeta)
      if (is.null(mbeta_names)) {
        stop("MGflashfm JAMmulti2() returned unnamed mbeta objects for trait ", i, ".")
      }

      for (model_name in duplicated_names) {
        model_idx <- which(model_names == model_name)
        model_pp <- as.numeric(models$PP[model_idx])
        if (
          any(!is.finite(model_pp)) ||
          any(model_pp < 0) ||
          !isTRUE(all.equal(
            model_pp,
            rep(model_pp[[1L]], length(model_pp)),
            tolerance = 100 * .Machine$double.eps,
            check.attributes = FALSE
          ))
        ) {
          stop(
            "MGflashfm JAMmulti2() returned non-identical duplicate posterior ",
            "copies for trait ", i, ", model ", dQuote(model_name), "."
          )
        }

        beta_idx <- which(mbeta_names == model_name)
        if (length(beta_idx) != length(model_idx)) {
          stop(
            "MGflashfm JAMmulti2() duplicate model/mbeta counts disagree for ",
            "trait ", i, ", model ", dQuote(model_name), "."
          )
        }
        beta_reference <- mbeta[[beta_idx[[1L]]]]
        beta_equal <- vapply(
          beta_idx,
          function(j) isTRUE(all.equal(
            mbeta[[j]],
            beta_reference,
            tolerance = 0,
            check.attributes = TRUE
          )),
          logical(1)
        )
        if (!all(beta_equal)) {
          stop(
            "MGflashfm JAMmulti2() returned non-identical duplicate mbeta ",
            "copies for trait ", i, ", model ", dQuote(model_name), "."
          )
        }
      }

      keep_model <- !duplicated(model_names)
      repaired_pp <- as.numeric(models$PP[keep_model])
      pp_total <- sum(repaired_pp)
      if (!is.finite(pp_total) || pp_total <= 0) {
        stop("MGflashfm duplicate-model repair found invalid posterior mass for trait ", i, ".")
      }
      repaired_pp <- repaired_pp / pp_total

      ppdf <- data.frame(
        str = model_names[keep_model],
        PP = repaired_pp,
        stringsAsFactors = FALSE
      )
      out$SM[[i]] <- PP2snpmod_function(ppdf)
      out$mbeta[[i]] <- mbeta[!duplicated(mbeta_names)]
      repaired_traits <- c(repaired_traits, i)
    }

    if (length(repaired_traits)) {
      attr(out, "mvpvb_duplicate_model_traits") <- repaired_traits
    }
    out
  }

  fixed <- function(
    gwas.list,
    corX,
    ybar,
    Vy,
    N,
    r2 = 0.99,
    save.path,
    maxcv = 2,
    maxcv_stop = 20,
    maxcv_autocheck = TRUE,
    jam.nM.iter = 1,
    extra.java.arguments = NULL
  ) {
    out <- original(
      gwas.list = gwas.list,
      corX = corX,
      ybar = ybar,
      Vy = Vy,
      N = N,
      r2 = r2,
      save.path = save.path,
      maxcv = maxcv,
      maxcv_stop = maxcv_stop,
      maxcv_autocheck = maxcv_autocheck,
      jam.nM.iter = jam.nM.iter,
      extra.java.arguments = extra.java.arguments
    )
    repair_output(out)
  }

  probe_sm <- PP2snpmod_function(data.frame(
    str = c("snp0001", "snp0001", "snp0002"),
    PP = c(0.4, 0.4, 0.2),
    stringsAsFactors = FALSE
  ))
  probe <- repair_output(list(
    SM = list(probe_sm),
    mbeta = list(stats::setNames(
      list(1, 1, 2),
      c("snp0001", "snp0001", "snp0002")
    ))
  ))
  ordinary_probe <- list(
    SM = list(PP2snpmod_function(data.frame(
      str = c("snp0001", "snp0002"),
      PP = c(0.6, 0.4),
      stringsAsFactors = FALSE
    ))),
    mbeta = list(stats::setNames(
      list(1, 2),
      c("snp0001", "snp0002")
    ))
  )
  ordinary_check <- repair_output(ordinary_probe)
  if (
    !identical(as.character(probe$SM[[1L]]@models$str), c("snp0001", "snp0002")) ||
    !isTRUE(all.equal(
      as.numeric(probe$SM[[1L]]@models$PP),
      c(2 / 3, 1 / 3),
      tolerance = 100 * .Machine$double.eps
    )) ||
    anyDuplicated(names(probe$mbeta[[1L]])) ||
    !identical(ordinary_check, ordinary_probe)
  ) {
    stop("Failed to validate the MGflashfm duplicate canonical-model patch.")
  }

  attr(fixed, "mvpvb_duplicate_canonical_model_patch") <- TRUE
  assignInNamespace(
    x = "JAMmulti2",
    value = fixed,
    ns = "MGflashfm"
  )

  patched <- get("JAMmulti2", envir = ns, inherits = FALSE)
  if (!isTRUE(attr(patched, "mvpvb_duplicate_canonical_model_patch"))) {
    stop("Failed to install the MGflashfm duplicate canonical-model patch.")
  }

  message(
    "MGflashfm duplicate-model patch applied: exact canonical copies are removed."
  )
  invisible(TRUE)
}


## ------------------------------------------------------------
## Bound the exponential six-trait joint-model calculation.
##
## flashfm::ppadjT6fast() explicitly enumerates the Cartesian product of
## the retained models for all six traits, and flashfmU() repeats that work
## once per trait.  Diffuse JAM posteriors can therefore turn a modest list
## of per-trait models into billions of joint evaluations.  flashfmU()
## already exposes the official maxmod argument for this purpose, but
## FLASHFMwithJAMd() always calls it with maxmod = NULL.
##
## This patch leaves the ordinary cpp-based path unchanged when its predicted
## work is within budget.  Only an over-budget six-trait analysis receives
## the largest official maxmod cap whose predicted work stays within budget.
## ------------------------------------------------------------
patch_mgflashfm_six_trait_runtime_guard <- function(
    max_joint_evaluations = 1200000
) {

  if (!requireNamespace("MGflashfm", quietly = TRUE)) {
    stop("Package MGflashfm is not installed.")
  }

  if (
    length(max_joint_evaluations) != 1L ||
    !is.finite(max_joint_evaluations) ||
    max_joint_evaluations < 6
  ) {
    stop("max_joint_evaluations must be one finite number >= 6.")
  }
  max_joint_evaluations <- as.numeric(max_joint_evaluations)

  ns <- asNamespace("MGflashfm")
  current <- get("flashfmU", envir = ns, inherits = FALSE)
  required_formals <- c(
    "main.input", "TOdds", "covY", "ss.stats", "cpp", "maxmod",
    "fastapprox", "NCORES"
  )
  if (!identical(names(formals(current)), required_formals)) {
    stop(
      "Unsupported MGflashfm::flashfmU() signature: ",
      paste(names(formals(current)), collapse = ", "),
      "."
    )
  }

  if (isTRUE(attr(current, "mvpvb_six_trait_runtime_guard"))) {
    installed_budget <- attr(current, "mvpvb_max_joint_evaluations")
    if (!identical(as.numeric(installed_budget), max_joint_evaluations)) {
      stop(
        "MGflashfm six-trait runtime guard is already installed with budget ",
        installed_budget,
        "; requested ",
        max_joint_evaluations,
        "."
      )
    }
    return(invisible(TRUE))
  }

  original <- current

  fixed <- function(
    main.input,
    TOdds,
    covY,
    ss.stats,
    cpp = 0.99,
    maxmod = NULL,
    fastapprox = FALSE,
    NCORES
  ) {
    effective_maxmod <- maxmod

    if (is.null(maxmod) && length(main.input$SM) == 6L) {
      model_tables <- lapply(main.input$SM, function(x) x@models)
      raw_counts <- vapply(model_tables, nrow, integer(1L))

      if (
        any(raw_counts < 1L) ||
        any(vapply(
          model_tables,
          function(x) {
            !"PP" %in% names(x) ||
              any(!is.finite(x$PP)) ||
              any(x$PP < 0) ||
              sum(x$PP) <= 0
          },
          logical(1L)
        ))
      ) {
        stop("Invalid six-trait model posterior supplied to flashfmU().")
      }

      cpp_counts <- vapply(
        model_tables,
        function(x) {
          pp <- sort(as.numeric(x$PP), decreasing = TRUE)
          wh <- which(cumsum(pp) <= cpp)
          if (!length(wh)) wh <- 1L
          length(unique(c(wh, min(max(wh) + 1L, length(pp)))))
        },
        integer(1L)
      )

      uncapped_states <- prod(as.numeric(cpp_counts))
      uncapped_evaluations <- 6 * uncapped_states

      if (uncapped_evaluations > max_joint_evaluations) {
        candidate_caps <- seq_len(max(raw_counts))
        candidate_evaluations <- vapply(
          candidate_caps,
          function(cap) {
            6 * prod(as.numeric(pmin(raw_counts, cap)))
          },
          numeric(1L)
        )
        feasible <- candidate_caps[
          candidate_evaluations <= max_joint_evaluations
        ]
        if (!length(feasible)) feasible <- 1L
        effective_maxmod <- max(feasible)
        capped_counts <- pmin(raw_counts, effective_maxmod)
        capped_evaluations <- 6 * prod(as.numeric(capped_counts))
        retained_mass <- vapply(
          model_tables,
          function(x) {
            pp <- sort(as.numeric(x$PP), decreasing = TRUE)
            sum(pp[seq_len(min(effective_maxmod, length(pp)))]) / sum(pp)
          },
          numeric(1L)
        )

        message(
          "MGflashfm six-trait runtime guard triggered: cpp model counts = ",
          paste(cpp_counts, collapse = "/"),
          "; predicted evaluations = ",
          format(uncapped_evaluations, scientific = FALSE, trim = TRUE),
          "."
        )
        message(
          "Using official flashfmU maxmod = ",
          effective_maxmod,
          "; capped model counts = ",
          paste(capped_counts, collapse = "/"),
          "; predicted evaluations = ",
          format(capped_evaluations, scientific = FALSE, trim = TRUE),
          "; retained posterior mass = ",
          paste(formatC(retained_mass, digits = 6L, format = "f"), collapse = "/"),
          "."
        )
      }
    }

    original(
      main.input = main.input,
      TOdds = TOdds,
      covY = covY,
      ss.stats = ss.stats,
      cpp = cpp,
      maxmod = effective_maxmod,
      fastapprox = fastapprox,
      NCORES = NCORES
    )
  }

  attr(fixed, "mvpvb_six_trait_runtime_guard") <- TRUE
  attr(fixed, "mvpvb_max_joint_evaluations") <- max_joint_evaluations
  assignInNamespace(x = "flashfmU", value = fixed, ns = "MGflashfm")

  patched <- get("flashfmU", envir = ns, inherits = FALSE)
  if (
    !isTRUE(attr(patched, "mvpvb_six_trait_runtime_guard")) ||
    !identical(
      as.numeric(attr(patched, "mvpvb_max_joint_evaluations")),
      max_joint_evaluations
    )
  ) {
    stop("Failed to install the MGflashfm six-trait runtime guard.")
  }

  message(
    "MGflashfm six-trait runtime guard applied: joint-evaluation budget = ",
    format(max_joint_evaluations, scientific = FALSE, trim = TRUE),
    "."
  )
  invisible(TRUE)
}


## ------------------------------------------------------------
## Prevent best.models.cpp() from appending an out-of-range row.
##
## In the MGflashfm version used on the server, cpp selection always appends
## max(wh) + 1.  When JAM returns only the legal null model ("1", PP = 1),
## that index is 2 although the model table has one row.  R then creates an
## all-NA model row; flashfmU() later reaches FALSE | NA inside calcQ12() and
## stops with "missing value where TRUE/FALSE needed".
##
## The wrapper below delegates every ordinary case to the installed package
## unchanged.  It intervenes only when the package's next-row index would be
## outside the model table, retaining the valid in-range rows and applying the
## same posterior-probability normalization/output construction.
## ------------------------------------------------------------
patch_mgflashfm_best_models_cpp <- function() {

  if (!requireNamespace("MGflashfm", quietly = TRUE)) {
    stop("Package MGflashfm is not installed.")
  }

  ns <- asNamespace("MGflashfm")
  current <- get("best.models.cpp", envir = ns, inherits = FALSE)

  required_formals <- c("d", "cpp.thr", "maxmod")
  if (!identical(names(formals(current)), required_formals)) {
    stop(
      "Unsupported MGflashfm::best.models.cpp() signature: ",
      paste(names(formals(current)), collapse = ", "),
      "."
    )
  }

  ## Make repeated calls within the same R process harmless.
  if (isTRUE(attr(current, "mvpvb_cpp_boundary_patch"))) {
    return(invisible(FALSE))
  }

  original <- current
  fixed <- function(d, cpp.thr = 0.99, maxmod = NULL) {
    if (is.null(maxmod)) {
      models <- d@models[
        order(d@models$PP, decreasing = TRUE),
        ,
        drop = FALSE
      ]

      if (nrow(models) < 1L) {
        stop("MGflashfm best.models.cpp() received an empty model table.")
      }
      if (
        any(!is.finite(models$PP)) ||
        any(models$PP < 0)
      ) {
        stop("MGflashfm best.models.cpp() received invalid model probabilities.")
      }

      cumulative_pp <- cumsum(models$PP)
      keep <- which(cumulative_pp <= cpp.thr)
      if (!length(keep)) keep <- 1L
      next_index <- max(keep) + 1L

      ## Only replace the package path that would create an all-NA row.
      if (next_index > nrow(models)) {
        keep <- unique(keep[
          keep >= 1L & keep <= nrow(models)
        ])
        if (!length(keep)) {
          stop("MGflashfm CPP selection retained no valid model row.")
        }

        d2 <- d
        d2@models <- models[keep, , drop = FALSE]
        pp_total <- sum(d2@models$PP)
        if (!is.finite(pp_total) || pp_total <= 0) {
          stop("MGflashfm CPP-selected probabilities have invalid total mass.")
        }
        d2@models$PP <- d2@models$PP / pp_total

        makestr <- get("makestr", envir = ns, inherits = TRUE)
        snps <- unlist(
          lapply(
            strsplit(as.character(d2@models$str), "%", fixed = TRUE),
            makestr
          ),
          use.names = FALSE
        )

        return(list(
          models = cbind(d2@models, snps = snps),
          old.cpp = cpp.thr
        ))
      }
    }

    original(d = d, cpp.thr = cpp.thr, maxmod = maxmod)
  }

  attr(fixed, "mvpvb_cpp_boundary_patch") <- TRUE
  assignInNamespace(
    x = "best.models.cpp",
    value = fixed,
    ns = "MGflashfm"
  )

  patched <- get(
    "best.models.cpp",
    envir = asNamespace("MGflashfm"),
    inherits = FALSE
  )
  if (!isTRUE(attr(patched, "mvpvb_cpp_boundary_patch"))) {
    stop("Failed to install the MGflashfm CPP boundary patch.")
  }

  message(
    "MGflashfm CPP boundary patch applied: out-of-range model rows are blocked."
  )
  invisible(TRUE)
}


## ------------------------------------------------------------
## Preserve matrix dimensions in MGflashfm::MPP.fn() when a trait has
## exactly one retained model.
##
## The installed MGflashfm version uses apply() twice.  For a one-model
## trait, the first apply() simplifies its result to a scalar/vector; the
## second apply(tmp1, 1, sum) then fails because tmp1 has no dim attribute.
## This patch delegates every multi-model input to the package unchanged.
## For a one-model input it reproduces the package function, restoring only
## the dimension that apply() dropped before performing the same row sums.
## ------------------------------------------------------------
patch_mgflashfm_mpp_single_model <- function() {

  if (!requireNamespace("MGflashfm", quietly = TRUE)) {
    stop("Package MGflashfm is not installed.")
  }

  ns <- asNamespace("MGflashfm")
  current <- get("MPP.fn", envir = ns, inherits = FALSE)

  if (!identical(names(formals(current)), "PP1")) {
    stop(
      "Unsupported MGflashfm::MPP.fn() signature: ",
      paste(names(formals(current)), collapse = ", "),
      "."
    )
  }

  ## Make repeated calls within the same R process harmless.
  if (isTRUE(attr(current, "mvpvb_mpp_single_model_patch"))) {
    return(invisible(FALSE))
  }

  probe <- matrix(
    c(1, 1),
    nrow = 1L,
    dimnames = list("1", c("pp.null", "pp.1"))
  )
  probe_result <- tryCatch(current(probe), error = identity)

  ## Do not replace a future package version that already handles this case.
  if (
    !inherits(probe_result, "error") &&
    is.matrix(probe_result) &&
    identical(dim(probe_result), c(1L, 2L)) &&
    identical(rownames(probe_result), "1") &&
    identical(colnames(probe_result), c("pp.null", "pp.1")) &&
    !anyNA(probe_result) &&
    all(is.finite(probe_result))
  ) {
    return(invisible(FALSE))
  }

  original <- current
  original_environment <- environment(current)
  if (!is.environment(original_environment)) {
    stop("MGflashfm::MPP.fn() has no valid function environment.")
  }

  ## Resolve the helpers exactly where the installed MPP.fn() resolves its
  ## free variables.  In MGflashfm 0.0.0.9000 these helpers are not exported
  ## bindings in asNamespace("MGflashfm"), although the official function can
  ## reach them through its own lexical environment/import chain.
  sep_function <- get(
    "sep.fn",
    envir = original_environment,
    inherits = TRUE
  )
  check_function <- get(
    "check.fn",
    envir = original_environment,
    inherits = TRUE
  )

  fixed <- function(PP1) {
    ## Preserve the installed package behavior for every ordinary input.
    if (!is.matrix(PP1) || nrow(PP1) != 1L) {
      return(original(PP1))
    }

    if (
      ncol(PP1) < 1L ||
      is.null(rownames(PP1)) ||
      is.null(colnames(PP1))
    ) {
      stop("MGflashfm MPP.fn() received an invalid one-model matrix.")
    }

    mnames <- rownames(PP1)
    model_index <- matrix(seq_along(mnames), ncol = 1L)
    msep <- apply(model_index, 1L, sep_function, mnames)
    gnames <- unique(unlist(msep, use.names = FALSE))

    if (length(gnames) < 1L || anyNA(gnames)) {
      stop("MGflashfm MPP.fn() produced invalid model names.")
    }

    mpp1 <- NULL
    for (k in seq_len(ncol(PP1))) {
      tmp1 <- apply(
        model_index,
        1L,
        check_function,
        msep,
        PP1[, k],
        gnames
      )

      ## apply() drops the matrix dimension when there is one model.
      if (is.null(dim(tmp1))) {
        expected_length <- length(gnames) * length(mnames)
        if (length(tmp1) != expected_length) {
          stop(
            "MGflashfm MPP.fn() returned an unexpected one-model result: ",
            "length ", length(tmp1), ", expected ", expected_length, "."
          )
        }
        tmp1 <- matrix(
          tmp1,
          nrow = length(gnames),
          ncol = length(mnames)
        )
      }

      mpp1 <- rbind(mpp1, rowSums(tmp1))
    }

    mpp1 <- data.frame(mpp1, row.names = colnames(PP1))
    names(mpp1) <- gnames
    t(mpp1)
  }

  attr(fixed, "mvpvb_mpp_single_model_patch") <- TRUE
  assignInNamespace(
    x = "MPP.fn",
    value = fixed,
    ns = "MGflashfm"
  )

  patched <- get("MPP.fn", envir = ns, inherits = FALSE)
  check <- patched(probe)
  if (
    !is.matrix(check) ||
    !identical(dim(check), c(1L, 2L)) ||
    !identical(rownames(check), "1") ||
    !identical(colnames(check), c("pp.null", "pp.1")) ||
    anyNA(check) ||
    any(!is.finite(check)) ||
    !isTRUE(all.equal(unname(check), unname(probe), tolerance = 0))
  ) {
    stop("Failed to install the MGflashfm one-model MPP patch.")
  }

  message(
    "MGflashfm MPP patch applied: one-model apply() dimensions preserved."
  )
  invisible(TRUE)
}


## ------------------------------------------------------------
## Preserve the model name passed from flashfm::PPsummarise() to
## flashfm::PPmodGroups() when a trait has exactly one retained model.
##
## The installed flashfm version extracts a posterior column with
## fm.multi$PP[[i]][, j].  For a 1 x 2 matrix, R simplifies this to an
## unnamed scalar, so names(modPP) is NULL inside PPmodGroups() and
## as.matrix(NULL) fails.  Every input without a one-model trait is delegated
## to the package unchanged.  On the boundary path, the package function is
## reproduced with only the two column extractions made explicitly named.
## ------------------------------------------------------------
patch_mgflashfm_ppsummarise_single_model <- function() {
  if (
    !requireNamespace("MGflashfm", quietly = TRUE) ||
    !requireNamespace("flashfm", quietly = TRUE)
  ) {
    stop("Packages MGflashfm and flashfm must both be installed.")
  }

  mg_ns <- asNamespace("MGflashfm")
  flash_ns <- asNamespace("flashfm")
  mg_binding_env <- mg_ns
  while (
    !identical(mg_binding_env, emptyenv()) &&
    !exists("PPsummarise", envir = mg_binding_env, inherits = FALSE)
  ) {
    mg_binding_env <- parent.env(mg_binding_env)
  }
  if (
    identical(mg_binding_env, emptyenv()) ||
    !exists("PPsummarise", envir = mg_binding_env, inherits = FALSE)
  ) {
    stop("MGflashfm cannot resolve its imported PPsummarise binding.")
  }
  current <- get("PPsummarise", envir = mg_binding_env, inherits = FALSE)
  required_formals <- c("fm.multi", "snpGroups", "minPP")

  if (!identical(names(formals(current)), required_formals)) {
    stop(
      "Unsupported flashfm::PPsummarise() signature: ",
      paste(names(formals(current)), collapse = ", "),
      "."
    )
  }

  if (isTRUE(attr(current, "mvpvb_ppsummarise_single_model_patch"))) {
    return(invisible(FALSE))
  }

  original <- current
  original_environment <- environment(current)
  if (!is.environment(original_environment)) {
    stop("flashfm::PPsummarise() has no valid function environment.")
  }
  PPmodGroups_function <- get(
    "PPmodGroups",
    envir = original_environment,
    inherits = TRUE
  )
  MPPmodGroups_function <- get(
    "MPPmodGroups",
    envir = original_environment,
    inherits = TRUE
  )
  smartbind_function <- get(
    "smartbind",
    envir = original_environment,
    inherits = TRUE
  )

  fixed <- function(fm.multi, snpGroups, minPP = 0.05) {
    one_model <- vapply(
      fm.multi$PP,
      function(x) is.matrix(x) && nrow(x) == 1L,
      logical(1)
    )

    ## Preserve the installed package behavior for every ordinary input.
    if (!any(one_model)) {
      return(original(fm.multi, snpGroups, minPP = minPP))
    }

    MPP <- fm.multi$MPP
    pp <- fm.multi$PP
    M <- length(pp)
    snpgroups <- snpGroups[[1]]
    snpgroups2 <- snpGroups[[2]]
    tnames <- names(fm.multi$PP)
    PPout <- PPout2 <- outPPg <- outPP <- vector("list", M)
    MPPout <- MPPout2 <- outMPPg <- outMPP <- vector("list", M)

    extract_named_pp <- function(model_pp, column_index) {
      if (
        !is.matrix(model_pp) ||
        ncol(model_pp) < column_index ||
        is.null(rownames(model_pp)) ||
        length(rownames(model_pp)) != nrow(model_pp)
      ) {
        stop("flashfm::PPsummarise() received an invalid posterior matrix.")
      }
      stats::setNames(
        as.numeric(model_pp[, column_index]),
        rownames(model_pp)
      )
    }

    for (i in seq_len(M)) {
      PPout[[i]] <- PPmodGroups_function(
        extract_named_pp(fm.multi$PP[[i]], 1L),
        snpgroups,
        minPP = 0
      )
      PPout2[[i]] <- PPmodGroups_function(
        extract_named_pp(fm.multi$PP[[i]], 2L),
        snpgroups2,
        minPP = 0
      )
      pplist <- list(
        t(data.frame(
          PPout[[i]]$group[, 2],
          row.names = rownames(PPout[[i]]$group)
        )),
        t(data.frame(
          PPout2[[i]]$group[, 2],
          row.names = rownames(PPout2[[i]]$group)
        ))
      )
      outPPg[[i]] <- t(do.call(smartbind_function, c(pplist, fill = 0)))
      colnames(outPPg[[i]]) <- fm.multi$sharing
      pplist <- list(
        t(data.frame(
          PPout[[i]]$snp[, 3],
          row.names = rownames(PPout[[i]]$snp)
        )),
        t(data.frame(
          PPout2[[i]]$snp[, 3],
          row.names = rownames(PPout2[[i]]$snp)
        ))
      )
      outPP[[i]] <- t(do.call(smartbind_function, c(pplist, fill = 0)))
      colnames(outPP[[i]]) <- fm.multi$sharing
      MPPout[[i]] <- MPPmodGroups_function(PPout[[i]]$snp)
      MPPout2[[i]] <- MPPmodGroups_function(PPout2[[i]]$snp)
      pplist <- list(
        t(data.frame(
          MPPout[[i]][, 1],
          row.names = rownames(MPPout[[i]])
        )),
        t(data.frame(
          MPPout2[[i]][, 1],
          row.names = rownames(MPPout2[[i]])
        ))
      )
      outMPP[[i]] <- t(do.call(smartbind_function, c(pplist, fill = 0)))
      colnames(outMPP[[i]]) <- fm.multi$sharing
      MPPout[[i]] <- MPPmodGroups_function(PPout[[i]]$group)
      MPPout2[[i]] <- MPPmodGroups_function(PPout2[[i]]$group)
      pplist <- list(
        t(data.frame(
          MPPout[[i]][, 1],
          row.names = rownames(MPPout[[i]])
        )),
        t(data.frame(
          MPPout2[[i]][, 1],
          row.names = rownames(MPPout2[[i]])
        ))
      )
      outMPPg[[i]] <- t(do.call(smartbind_function, c(pplist, fill = 0)))
      colnames(outMPPg[[i]]) <- fm.multi$sharing
      ind <- which(outPPg[[i]][, 1] > minPP | outPPg[[i]][, 2] > minPP)
      rnames <- rownames(outPPg[[i]])
      if (length(ind) > 1) {
        tmp <- outPPg[[i]][ind, ]
      } else {
        tmp <- t(as.matrix(outPPg[[i]][ind, ]))
      }
      rownames(tmp) <- rnames[ind]
      outPPg[[i]] <- tmp
      ind <- which(outMPPg[[i]][, 1] > minPP | outMPPg[[i]][, 2] > minPP)
      rnames <- rownames(outMPPg[[i]])
      if (length(ind) > 1) {
        tmp <- outMPPg[[i]][ind, ]
      } else {
        tmp <- t(as.matrix(outMPPg[[i]][ind, ]))
      }
      rownames(tmp) <- rnames[ind]
      outMPPg[[i]] <- tmp
    }

    names(outMPPg) <- names(outMPP) <- names(outPP) <- names(outPPg) <- tnames
    list(MPP = outMPP, MPPg = outMPPg, PP = outPP, PPg = outPPg)
  }

  attr(fixed, "mvpvb_ppsummarise_single_model_patch") <- TRUE

  replace_binding <- function(ns, name, value) {
    was_locked <- bindingIsLocked(name, ns)
    if (was_locked) unlockBinding(name, ns)
    assign(name, value, envir = ns)
    if (was_locked) lockBinding(name, ns)
  }

  ## FLASHFMwithJAMd() resolves the imported MGflashfm binding, while direct
  ## flashfm calls resolve the defining flashfm binding.  Keep both identical.
  replace_binding(flash_ns, "PPsummarise", fixed)
  replace_binding(mg_binding_env, "PPsummarise", fixed)

  patched_mg <- get("PPsummarise", envir = mg_binding_env, inherits = FALSE)
  patched_flash <- get("PPsummarise", envir = flash_ns, inherits = FALSE)
  if (
    !isTRUE(attr(patched_mg, "mvpvb_ppsummarise_single_model_patch")) ||
    !identical(patched_mg, patched_flash)
  ) {
    stop("Failed to install the flashfm one-model PPsummarise patch.")
  }

  message(
    "MGflashfm PPsummarise patch applied: one-model posterior names preserved."
  )
  invisible(TRUE)
}


## ------------------------------------------------------------
## Resolve SNP IDs returned by MGflashfm without fuzzy matching.
##
## Only harmless textual variants that map unambiguously to one supplied
## SNP ID are accepted. A genuinely foreign ID still stops with the
## offending value in the error message.
## ------------------------------------------------------------
resolve_mgflashfm_snp_ids <- function(ids, inputs, context) {

  known_ids <- as.character(inputs$snp_map$mgflashfm_rsID)
  raw_ids <- as.character(ids)

  if (
    length(raw_ids) == 0L ||
    length(known_ids) == 0L ||
    anyDuplicated(known_ids)
  ) {
    stop("Invalid SNP-ID map while processing ", context, ".")
  }

  clean_ids <- trimws(raw_ids)
  wrapped <- grepl("^([\"'`])(.*)\\1$", clean_ids, perl = TRUE)
  clean_ids[wrapped] <- sub(
    "^([\"'`])(.*)\\1$",
    "\\2",
    clean_ids[wrapped],
    perl = TRUE
  )

  snp_index <- match(clean_ids, known_ids)

  unresolved <- which(is.na(snp_index))
  if (length(unresolved) > 0L) {
    lower_known <- tolower(known_ids)
    if (anyDuplicated(lower_known)) {
      stop("Case-insensitive SNP IDs are not unique in the input map.")
    }
    snp_index[unresolved] <- match(
      tolower(clean_ids[unresolved]),
      lower_known
    )
  }

  ## Accept an anchored SNP-number form, optionally followed by the suffix
  ## added by make.unique(). Do not extract SNP substrings from model labels.
  unresolved <- which(is.na(snp_index))
  if (length(unresolved) > 0L) {
    candidate <- clean_ids[unresolved]
    valid_shape <- grepl(
      "(?i)^snp0*[0-9]+(?:\\.[0-9]+)?$",
      candidate,
      perl = TRUE
    )
    known_number <- suppressWarnings(as.integer(sub(
      "(?i)^snp0*",
      "",
      known_ids,
      perl = TRUE
    )))

    if (anyNA(known_number) || anyDuplicated(known_number)) {
      stop("Numeric SNP IDs are not unique in the input map.")
    }

    if (any(valid_shape)) {
      candidate_base <- sub(
        "\\.[0-9]+$",
        "",
        candidate[valid_shape]
      )
      candidate_number <- suppressWarnings(as.integer(sub(
        "(?i)^snp0*",
        "",
        candidate_base,
        perl = TRUE
      )))
      snp_index[unresolved[valid_shape]] <- match(
        candidate_number,
        known_number
      )
    }
  }

  if (anyNA(snp_index)) {
    bad <- unique(raw_ids[is.na(snp_index)])
    bad_shown <- paste(utils::head(dQuote(bad), 10L), collapse = ", ")
    if (length(bad) > 10L) bad_shown <- paste0(bad_shown, ", ...")
    stop(
      "MGflashfm returned unknown SNP ID(s) while processing ",
      context,
      ": ",
      bad_shown,
      "."
    )
  }

  as.integer(snp_index)
}


## ------------------------------------------------------------
## Final MGflashfm post-processing.
##
## This reproduces the output layout of MGFLASHFMwithJAM but uses
## a K x q sample-size matrix. The upstream package constructs
## data.frame(Nall), which is q x K and fails when q > K.
## ------------------------------------------------------------
postprocess_mgflashfm_core <- function(
    csta,
    flashfm.out,
    inputs
) {

  K <- length(inputs$gwas.list)
  q <- length(inputs$trait_names)

  if (!is.list(csta) || !is.list(csta$summary) || !is.list(csta$details)) {
    stop("Invalid MGflashfm core output.")
  }
  if (length(csta$summary) != q || length(csta$details) != q) {
    stop("MGflashfm core output does not contain q trait results.")
  }

  maf_by_trait <- vector("list", q)
  names(maf_by_trait) <- inputs$trait_names

  for (r in seq_len(q)) {

    cs <- csta$summary[[r]]

    if (is.data.frame(cs) && nrow(cs) > 0L) {

      tsnps <- rownames(cs)
      if (is.null(tsnps) || length(tsnps) != nrow(cs)) {
        stop("MGflashfm credible-set SNP names are invalid for Trait", r, ".")
      }

      snp_index <- resolve_mgflashfm_snp_ids(
        tsnps,
        inputs,
        paste0("credible-set summary for Trait", r)
      )
      canonical_ids <- inputs$snp_map$mgflashfm_rsID[snp_index]

      mafdf <- matrix(
        NA_real_,
        nrow = nrow(cs),
        ncol = K,
        dimnames = list(tsnps, inputs$group_names)
      )

      for (k in seq_len(K)) {
        gwas_r <- inputs$gwas.list[[k]][[r]]
        idx <- match(canonical_ids, gwas_r$rsID)
        if (anyNA(idx)) {
          stop(
            "Resolved SNP ID is absent from Group", k,
            " GWAS input for Trait", r, "."
          )
        }
        eaf <- as.numeric(gwas_r$EAF[idx])
        mafdf[, k] <- ifelse(
          is.na(eaf),
          NA_real_,
          ifelse(eaf < 0.5, eaf, 1 - eaf)
        )
      }

      maf_by_trait[[r]] <- mafdf

    } else {
      ## Single-bracket assignment preserves the trait slot when cs is NULL.
      maf_by_trait[r] <- list(NULL)
    }
  }

  ## Correct orientation: populations x traits = K x q.
  Nmat <- do.call(
    rbind,
    lapply(inputs$Nall, as.numeric)
  )
  rownames(Nmat) <- inputs$group_names
  colnames(Nmat) <- inputs$trait_names

  if (!identical(dim(Nmat), c(K, q))) {
    stop(
      "Internal Nall dimension error: expected ", K, "x", q,
      ", observed ", paste(dim(Nmat), collapse = "x"), "."
    )
  }

  cstaS <- vector("list", q)
  names(cstaS) <- inputs$trait_names

  for (r in seq_len(q)) {

    cs <- csta$summary[[r]]

    if (is.data.frame(cs) && nrow(cs) > 0L) {

      ## Keep MAF separate from cs. If cs already contains a Group* column,
      ## selecting names from data.frame(cs, mafdf) can raise
      ## "undefined columns selected" or select the wrong duplicate column.
      maf_mat <- as.matrix(maf_by_trait[[r]])

      if (!identical(dim(maf_mat), c(nrow(cs), K))) {
        stop("Internal MAF dimension error for Trait", r, ".")
      }

      Nw <- as.numeric(Nmat[, r])
      if (any(!is.finite(Nw)) || any(Nw <= 0)) {
        stop("Invalid sample size in Nall for Trait", r, ".")
      }

      pmaf <- apply(
        maf_mat,
        1L,
        function(x) {
          keep <- is.finite(x) & is.finite(Nw)
          if (!any(keep)) return(NA_real_)
          weighted.mean(x[keep], w = Nw[keep])
        }
      )

      Aprop <- apply(
        maf_mat,
        1L,
        function(x) mean(is.finite(x))
      )

      Aanc <- apply(
        maf_mat,
        1L,
        function(x) {
          paste(inputs$group_names[is.finite(x)], collapse = ",")
        }
      )

      cstaS[[r]] <- data.frame(
        cs,
        maf_pooled = pmaf,
        prop_groups = Aprop,
        groups = Aanc,
        check.names = FALSE
      )

    } else {
      ## Do not use [[r]] <- NULL: that deletes the list element and shifts
      ## every later trait index.
      cstaS[r] <- list(cs)
    }
  }

  list(
    CSsummary = cstaS,
    CSdetail = csta,
    flashfm.out = flashfm.out
  )
}


## ------------------------------------------------------------
## Checkpointed implementation of the expensive official workflow.
##
## The statistical core is unchanged:
##   FLASHFMwithJAMd is run once per population;
##   MGflashfm combines the population-specific posterior model output.
##
## Checkpoints are stored inside result_dir, not tempdir().
## ------------------------------------------------------------
run_mgflashfm_core_checkpointed <- function(
    inputs,
    result_dir,
    jam.nM.iter = 1,
    jam.seed = 20260929L,
    cpp = 0.99,
    cred = 0.99,
    TOdds = 1,
    maxcv = 1,
    maxcv_stop = 20,
    r2 = 0.80,
    NCORES = 1,
    extra.java.arguments = NULL,
    resume = TRUE,
    scenario_name = NULL
) {

  ## Must run inside every worker process before FLASHFMwithJAMd().
  patch_mgflashfm_null_model()
  patch_mgflashfm_duplicate_canonical_models()
  patch_mgflashfm_six_trait_runtime_guard()
  patch_mgflashfm_best_models_cpp()
  patch_mgflashfm_mpp_single_model()
  patch_mgflashfm_ppsummarise_single_model()

  K <- length(inputs$gwas.list)
  q <- length(inputs$trait_names)
  m <- nrow(inputs$snp_map)

  iter_tag <- format(jam.nM.iter, scientific = FALSE, trim = TRUE)
  iter_tag <- gsub("\\.", "p", iter_tag)
  seed_tag <- as.character(as.integer(jam.seed))
  run_tag <- paste0("jam_", iter_tag, "M_seed_", seed_tag)

  checkpoint_dir <- file.path(result_dir, "checkpoints", run_tag)
  jam_dir <- file.path(result_dir, "tmpJAM", run_tag)
  dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(jam_dir, recursive = TRUE, showWarnings = FALSE)

  settings <- list(
    scenario_name = if (is.null(scenario_name)) NA_character_ else as.character(scenario_name),
    K = K,
    m = m,
    q = q,
    jam.nM.iter = as.numeric(jam.nM.iter),
    jam.seed = as.integer(jam.seed),
    cpp = cpp,
    cred = cred,
    TOdds = TOdds,
    maxcv = maxcv,
    maxcv_stop = maxcv_stop,
    r2 = r2,
    NCORES = NCORES
  )

  same_settings <- function(x) {
    is.list(x) && isTRUE(all.equal(x, settings, check.attributes = FALSE))
  }

  set.seed(as.integer(jam.seed))

  flashfm.out <- vector("list", K)
  names(flashfm.out) <- inputs$group_names

  ssnps <- vector("list", K)
  allflashfmPP <- vector("list", K)
  names(allflashfmPP) <- inputs$group_names

  ybar <- setNames(rep(0, q), inputs$trait_names)

  for (k in seq_len(K)) {

    group_checkpoint <- file.path(
      checkpoint_dir,
      sprintf("flashfm_group_%d.rds", k)
    )

    loaded <- FALSE

    if (isTRUE(resume) && file.exists(group_checkpoint)) {
      chk <- tryCatch(readRDS(group_checkpoint), error = function(e) NULL)

      if (
        is.list(chk) &&
        same_settings(chk$settings) &&
        identical(as.integer(chk$group_index), as.integer(k)) &&
        !is.null(chk$flashfm.out)
      ) {
        flashfm.out[[k]] <- chk$flashfm.out
        if (!is.null(chk$rng_state_after)) {
          assign(".Random.seed", chk$rng_state_after, envir = .GlobalEnv)
        }
        loaded <- TRUE
        cat(
          "MGflashfm: loaded population ", k,
          " checkpoint (", run_tag, ").\n",
          sep = ""
        )
      }
    }

    if (!loaded) {
      cat(
        "MGflashfm: running population ", k, " of ", K,
        " | q=", q,
        " | jam.nM.iter=", jam.nM.iter, "M\n",
        sep = ""
      )

      flashfm.out[[k]] <- MGflashfm::FLASHFMwithJAMd(
        gwas.list = inputs$gwas.list[[k]],
        corX = inputs$LD.list[[k]],
        ybar = ybar,
        N = inputs$Nall[[k]],
        save.path = jam_dir,
        TOdds = TOdds,
        covY = inputs$covY.list[[k]],
        cpp = cpp,
        NCORES = NCORES,
        maxcv = maxcv,
        maxcv_stop = maxcv_stop,
        jam.nM.iter = jam.nM.iter,
        r2 = r2,
        extra.java.arguments = extra.java.arguments
      )

      saveRDS(
        list(
          settings = settings,
          group_index = as.integer(k),
          group_name = inputs$group_names[[k]],
          flashfm.out = flashfm.out[[k]],
          rng_state_after = if (exists(".Random.seed", envir = .GlobalEnv)) {
            get(".Random.seed", envir = .GlobalEnv)
          } else {
            NULL
          }
        ),
        group_checkpoint
      )

      cat(
        "MGflashfm: saved population ", k,
        " checkpoint.\n",
        sep = ""
      )
    }

    ssnps[[k]] <- Reduce(
      intersect,
      lapply(inputs$gwas.list[[k]], function(x) x$rsID)
    )

    allflashfmPP[[k]] <- flashfm.out[[k]]$mpp.pp

    if (is.null(allflashfmPP[[k]])) {
      stop("Population ", k, " flashfm output does not contain mpp.pp.")
    }
  }

  all_groups_checkpoint <- file.path(
    checkpoint_dir,
    "flashfm_all_groups.rds"
  )

  saveRDS(
    list(
      settings = settings,
      flashfm.out = flashfm.out,
      allflashfmPP = allflashfmPP,
      ssnps = ssnps
    ),
    all_groups_checkpoint
  )

  allsnps <- Reduce(union, ssnps)

  csta_checkpoint <- file.path(
    checkpoint_dir,
    "MGflashfm_core_csta.rds"
  )

  csta <- NULL

  if (isTRUE(resume) && file.exists(csta_checkpoint)) {
    chk <- tryCatch(readRDS(csta_checkpoint), error = function(e) NULL)
    if (is.list(chk) && same_settings(chk$settings) && !is.null(chk$csta)) {
      csta <- chk$csta
      cat("MGflashfm: loaded cross-population core checkpoint.\n")
    }
  }

  if (is.null(csta)) {
    cat("MGflashfm: combining population-specific flashfm results.\n")

    csta <- MGflashfm::MGflashfm(
      allflashfmPP = allflashfmPP,
      Nall = inputs$Nall,
      snps = allsnps,
      cred = cred,
      multi = TRUE,
      cpp = cpp,
      NCORES = NCORES
    )

    saveRDS(
      list(
        settings = settings,
        csta = csta
      ),
      csta_checkpoint
    )

    cat("MGflashfm: saved cross-population core checkpoint.\n")
  }

  mgmtCS <- postprocess_mgflashfm_core(
    csta = csta,
    flashfm.out = flashfm.out,
    inputs = inputs
  )

  list(
    mgmtCS = mgmtCS,
    checkpoint_dir = checkpoint_dir,
    jam_dir = jam_dir,
    settings = settings
  )
}


credible_set_pairs <- function(mgmtCS, inputs, common_data) {
  rows <- list()
  out_i <- 1L

  for (r in seq_along(inputs$trait_names)) {
    trait <- inputs$trait_names[[r]]
    cs <- mgmtCS$CSsummary[[trait]]
    if (is.null(cs)) cs <- mgmtCS$CSsummary[[r]]
    if (!is.data.frame(cs) || nrow(cs) == 0L) next

    ids <- rownames(cs)
    snp_index <- resolve_mgflashfm_snp_ids(
      ids,
      inputs,
      paste0("credible-set pairs for ", trait)
    )

    one <- data.frame(
      SNP = snp_index,
      SNP_key = inputs$snp_map$SNP_key[snp_index],
      trait_index = r,
      trait = trait,
      mgMPP = as.numeric(cs$mgMPP),
      maf_pooled = as.numeric(cs$maf_pooled),
      prop_groups = as.numeric(cs$prop_groups),
      groups = as.character(cs$groups),
      in_CS99 = TRUE,
      is_true_in_population1 = inputs$truth_masks[[1L]][cbind(snp_index, r)],
      stringsAsFactors = FALSE
    )

    ## Backward-compatible name retained for old downstream code.
    one$is_true_shared_pair <- one$is_true_in_population1

    for (k in seq_len(common_data$K)) {
      one[[paste0("B_true_Group", k)]] <-
        common_data$B_true_list[[k]][cbind(snp_index, r)]
      one[[paste0("is_true_Group", k)]] <-
        inputs$truth_masks[[k]][cbind(snp_index, r)]
    }

    rows[[out_i]] <- one
    out_i <- out_i + 1L
  }

  if (length(rows) == 0L) {
    out <- data.frame(
      SNP = integer(),
      SNP_key = character(),
      trait_index = integer(),
      trait = character(),
      mgMPP = numeric(),
      maf_pooled = numeric(),
      prop_groups = numeric(),
      groups = character(),
      in_CS99 = logical(),
      is_true_in_population1 = logical(),
      is_true_shared_pair = logical(),
      stringsAsFactors = FALSE
    )
    return(out)
  }

  do.call(rbind, rows)
}


extract_mgmpp_selection <- function(
    mgmtCS,
    inputs,
    common_data,
    threshold = 0.90
) {

  m <- as.integer(common_data$m)
  q <- as.integer(common_data$q)

  mgMPP <- matrix(
    0,
    nrow = m,
    ncol = q,
    dimnames = list(inputs$snp_map$mgflashfm_rsID, inputs$trait_names)
  )

  details <- mgmtCS$CSdetail$details
  if (!is.list(details)) {
    stop("MGflashfm output does not contain CSdetail$details.")
  }
  if (length(details) != q) {
    stop("MGflashfm CSdetail$details does not contain q traits.")
  }

  for (r in seq_len(q)) {

    trait <- inputs$trait_names[[r]]
    detail <- details[[trait]]
    if (is.null(detail)) detail <- details[[r]]

    ## If the core method returns no mgMPP object for a trait, preserve
    ## the original all-zero row/column contribution. Non-finite numeric
    ## mgMPP values are still treated as an error below.
    if (is.null(detail) || is.null(detail$mgMPP)) next

    values <- as.numeric(detail$mgMPP)
    ids <- names(detail$mgMPP)

    if (is.null(ids) || length(ids) != length(values)) {
      stop("MGflashfm mgMPP values have no valid SNP names for ", trait, ".")
    }
    ## Numerical boundary tolerance only. Values outside [0, 1] merely by
    ## floating-point rounding are projected back to the probability range.
    mgmpp_tol <- 100 * .Machine$double.eps

    if (any(!is.finite(values))) {
      stop("MGflashfm returned non-finite mgMPP values for ", trait, ".")
    }

    large_bad <- values < -mgmpp_tol | values > 1 + mgmpp_tol
    if (any(large_bad)) {
      stop(
        "MGflashfm returned materially out-of-range mgMPP values for ",
        trait,
        ": range = [",
        format(min(values), digits = 17),
        ", ",
        format(max(values), digits = 17),
        "]."
      )
    }

    values <- pmin(1, pmax(0, values))

    ## MGflashfm can emit a zero-probability entry whose name is an actual
    ## NA_character_ sentinel. It is not an SNP and cannot affect the
    ## mgMPP > threshold selection. Drop only zero-valued missing-name
    ## sentinels; a missing name carrying nonzero mass remains an error.
    missing_id <- is.na(ids) | !nzchar(trimws(ids))
    if (any(missing_id & values > mgmpp_tol)) {
      stop(
        "MGflashfm returned a missing SNP ID with nonzero mgMPP for ",
        trait,
        "."
      )
    }
    if (any(missing_id)) {
      values <- values[!missing_id]
      ids <- ids[!missing_id]
    }
    if (length(values) == 0L) next

    snp_index <- resolve_mgflashfm_snp_ids(
      ids,
      inputs,
      paste0("mgMPP for ", trait)
    )

    if (anyDuplicated(snp_index)) {
      values <- vapply(
        split(values, snp_index),
        max,
        numeric(1),
        na.rm = TRUE
      )
      snp_index <- as.integer(names(values))
    }

    mgMPP[snp_index, r] <- values
  }

  Zhat <- matrix(
    as.integer(mgMPP > threshold),
    nrow = m,
    ncol = q,
    dimnames = dimnames(mgMPP)
  )

  all_grid <- expand.grid(
    SNP = seq_len(m),
    trait_index = seq_len(q),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )

  all_evidence <- data.frame(
    SNP = all_grid$SNP,
    SNP_key = inputs$snp_map$SNP_key[all_grid$SNP],
    trait_index = all_grid$trait_index,
    trait = inputs$trait_names[all_grid$trait_index],
    mgMPP = mgMPP[cbind(all_grid$SNP, all_grid$trait_index)],
    selected = Zhat[cbind(all_grid$SNP, all_grid$trait_index)] == 1L,
    stringsAsFactors = FALSE
  )

  selected_pairs <- all_evidence[
    all_evidence$selected,
    ,
    drop = FALSE
  ]

  list(
    mgMPP = mgMPP,
    Zhat = Zhat,
    all_evidence = all_evidence,
    selected_pairs = selected_pairs,
    threshold = threshold
  )
}


evaluate_pair_selection <- function(Zhat, truth, population = NA_integer_) {

  pred <- as.matrix(Zhat != 0)
  truth <- as.matrix(truth != 0)

  if (!identical(dim(pred), dim(truth))) {
    stop("Prediction and truth matrices have different dimensions.")
  }

  TP <- sum(pred & truth)
  FP <- sum(pred & !truth)
  FN <- sum(!pred & truth)
  TN <- sum(!pred & !truth)

  precision <- if (TP + FP == 0L) 0 else TP / (TP + FP)
  recall <- if (TP + FN == 0L) 0 else TP / (TP + FN)
  specificity <- if (TN + FP == 0L) 0 else TN / (TN + FP)
  F1 <- if (precision + recall == 0) {
    0
  } else {
    2 * precision * recall / (precision + recall)
  }

  data.frame(
    method = "MGflashfm",
    population = as.integer(population),
    selection_rule = "mgMPP > 0.90",
    TP = TP,
    FP = FP,
    FN = FN,
    TN = TN,
    precision = precision,
    recall = recall,
    specificity = specificity,
    F1 = F1,
    stringsAsFactors = FALSE
  )
}


## Backward-compatible helper retained for old callers.
evaluate_shared_pairs <- function(Zhat, truth) {
  out <- evaluate_pair_selection(Zhat, truth, population = 1L)
  out$population <- NULL
  out
}


build_truth_pairs_by_population <- function(
    common_data,
    inputs,
    selected_pairs
) {

  K <- as.integer(common_data$K)
  rows <- vector("list", K)

  predicted_key <- paste(
    selected_pairs$SNP,
    selected_pairs$trait_index,
    sep = ":"
  )

  for (k in seq_len(K)) {

    truth_idx <- which(
      inputs$truth_masks[[k]],
      arr.ind = TRUE
    )

    truth_pairs <- data.frame(
      population = k,
      SNP = truth_idx[, 1L],
      SNP_key = common_data$snp_key_vec[truth_idx[, 1L]],
      trait_index = truth_idx[, 2L],
      trait = inputs$trait_names[truth_idx[, 2L]],
      B_true = common_data$B_true_list[[k]][truth_idx],
      stringsAsFactors = FALSE
    )

    truth_key <- paste(
      truth_pairs$SNP,
      truth_pairs$trait_index,
      sep = ":"
    )

    truth_pairs$detected_mgMPP90 <- truth_key %in% predicted_key
    rows[[k]] <- truth_pairs
  }

  do.call(rbind, rows)
}


run_mgflashfm <- function(
    common_file,
    result_dir,
    jam.nM.iter = getOption("MGflashfm.jam.nM.iter", 1),
    jam.seed = getOption("MGflashfm.jam.seed", 20260929L),
    resume = TRUE
) {

  if (!requireNamespace("MGflashfm", quietly = TRUE)) {
    stop("Package MGflashfm is not installed.")
  }
  if (!requireNamespace("R2BGLiMS", quietly = TRUE)) {
    stop("Package R2BGLiMS is not installed; JAM also requires Java.")
  }

  if (
    length(jam.nM.iter) != 1L ||
    !is.numeric(jam.nM.iter) ||
    !is.finite(jam.nM.iter) ||
    jam.nM.iter <= 0
  ) {
    stop("jam.nM.iter must be one positive finite number.")
  }

  if (
    length(jam.seed) != 1L ||
    !is.numeric(jam.seed) ||
    !is.finite(jam.seed)
  ) {
    stop("jam.seed must be one finite numeric value.")
  }

  suppressPackageStartupMessages(library(R2BGLiMS))
  suppressPackageStartupMessages(library(MGflashfm))

  common_file <- normalizePath(common_file, mustWork = TRUE)
  dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
  result_dir <- normalizePath(result_dir, mustWork = TRUE)

  common_data <- readRDS(common_file)
  inputs <- make_mgflashfm_inputs(common_data)

  K <- as.integer(common_data$K)
  m <- as.integer(common_data$m)
  q <- as.integer(common_data$q)

  cat("\n============================================================\n")
  cat("Running MGflashfm\n")
  cat("============================================================\n")
  cat("Scenario:", common_data$scenario_name, "\n")
  cat("K =", K, "| m =", m, "| q =", q, "\n")
  cat(
    "jam.nM.iter =", jam.nM.iter,
    "million =",
    format(jam.nM.iter * 1e6, scientific = FALSE),
    "iterations per JAM run\n"
  )
  cat("JAM R seed =", as.integer(jam.seed), "\n")
  cat("resume =", isTRUE(resume), "\n")

  write.csv(
    inputs$snp_map,
    file.path(result_dir, "snp_id_map.csv"),
    row.names = FALSE
  )

  core_time <- system.time({
    core <- run_mgflashfm_core_checkpointed(
      inputs = inputs,
      result_dir = result_dir,
      jam.nM.iter = jam.nM.iter,
      jam.seed = as.integer(jam.seed),
      cpp = 0.99,
      cred = 0.99,
      TOdds = 1,
      maxcv = 1,
      maxcv_stop = 20,
      r2 = 0.80,
      NCORES = 1,
      extra.java.arguments = NULL,
      resume = resume,
      scenario_name = common_data$scenario_name
    )
  })

  mgmtCS <- core$mgmtCS

  ## Save the complete official-layout result before any later diagnostic
  ## or CSV code. This prevents losing an expensive fit to post-processing.
  saveRDS(
    mgmtCS,
    file.path(result_dir, "MGflashfm_official_result.rds")
  )

  selection <- extract_mgmpp_selection(
    mgmtCS = mgmtCS,
    inputs = inputs,
    common_data = common_data,
    threshold = 0.90
  )

  ## Save selection immediately as another recovery point.
  saveRDS(
    selection$mgMPP,
    file.path(result_dir, "MGflashfm_mgMPP_matrix.rds")
  )
  saveRDS(
    selection$Zhat,
    file.path(result_dir, "MGflashfm_Zhat_mgMPP90.rds")
  )

  cs_pairs <- credible_set_pairs(
    mgmtCS = mgmtCS,
    inputs = inputs,
    common_data = common_data
  )

  metrics_by_population <- do.call(
    rbind,
    lapply(
      seq_len(K),
      function(k) {
        evaluate_pair_selection(
          Zhat = selection$Zhat,
          truth = inputs$truth_masks[[k]],
          population = k
        )
      }
    )
  )

  truth_pairs <- build_truth_pairs_by_population(
    common_data = common_data,
    inputs = inputs,
    selected_pairs = selection$selected_pairs
  )

  ## Credible-set coverage by population. This remains diagnostic only;
  ## the simulation selection rule is mgMPP > 0.90.
  cs_key <- paste(cs_pairs$SNP, cs_pairs$trait_index, sep = ":")

  cs_metrics_by_population <- do.call(
    rbind,
    lapply(
      seq_len(K),
      function(k) {
        truth_k <- truth_pairs[truth_pairs$population == k, , drop = FALSE]
        truth_key <- paste(truth_k$SNP, truth_k$trait_index, sep = ":")

        data.frame(
          method = "MGflashfm_CS99",
          population = k,
          purpose = "credible_set_coverage_resolution_only",
          n_CS99_SNP_trait_pairs = nrow(cs_pairs),
          n_true_pairs_covered = sum(truth_key %in% cs_key),
          truth_pair_coverage = if (length(truth_key) == 0L) {
            NA_real_
          } else {
            mean(truth_key %in% cs_key)
          },
          mean_CS99_size = if (nrow(cs_pairs) == 0L) {
            0
          } else {
            mean(as.numeric(table(cs_pairs$trait_index)))
          },
          stringsAsFactors = FALSE
        )
      }
    )
  )

  ## Nonessential human-readable outputs are written only after all RDS
  ## recovery points above exist.
  write.csv(
    cs_pairs,
    file.path(result_dir, "MGflashfm_CS99_SNP_trait_pairs.csv"),
    row.names = FALSE
  )
  write.csv(
    selection$mgMPP,
    file.path(result_dir, "MGflashfm_mgMPP_matrix.csv"),
    row.names = FALSE
  )
  write.csv(
    selection$Zhat,
    file.path(result_dir, "MGflashfm_Zhat_mgMPP90.csv"),
    row.names = FALSE
  )
  write.csv(
    selection$all_evidence,
    file.path(result_dir, "MGflashfm_all_SNP_trait_mgMPP.csv"),
    row.names = FALSE
  )
  write.csv(
    selection$selected_pairs,
    file.path(result_dir, "MGflashfm_mgMPP90_SNP_trait_pairs.csv"),
    row.names = FALSE
  )
  write.csv(
    truth_pairs,
    file.path(result_dir, "truth_SNP_trait_pairs_by_population.csv"),
    row.names = FALSE
  )
  write.csv(
    metrics_by_population,
    file.path(result_dir, "MGflashfm_pair_metrics.csv"),
    row.names = FALSE
  )
  write.csv(
    cs_metrics_by_population,
    file.path(result_dir, "MGflashfm_CS99_coverage_resolution.csv"),
    row.names = FALSE
  )

  runtime <- list(
    core_elapsed_sec = unname(core_time["elapsed"]),
    core_user_sec = unname(core_time["user.self"]),
    core_system_sec = unname(core_time["sys.self"])
  )

  saveRDS(
    list(
      scenario_name = common_data$scenario_name,
      K = K,
      m = m,
      q = q,
      jam.nM.iter = as.numeric(jam.nM.iter),
      jam.iterations = as.numeric(jam.nM.iter) * 1e6,
      jam.seed = as.integer(jam.seed),
      runtime = runtime,
      checkpoint_dir = core$checkpoint_dir,
      jam_dir = core$jam_dir,
      metrics_by_population = metrics_by_population
    ),
    file.path(result_dir, "MGflashfm_run_metadata.rds")
  )

  cat("\nMGflashfm completed.\n")
  cat("mgMPP > 0.90 SNP-trait pairs:", nrow(selection$selected_pairs), "\n")
  cat("CS99 SNP-trait pairs:", nrow(cs_pairs), "\n")
  cat("Core elapsed seconds:", runtime$core_elapsed_sec, "\n")
  cat("Results:", result_dir, "\n")

  invisible(
    list(
      result = mgmtCS,
      mgMPP = selection$mgMPP,
      Zhat = selection$Zhat,
      selected_pairs = selection$selected_pairs,
      CS99_pairs = cs_pairs,
      truth = truth_pairs,
      metrics = metrics_by_population,
      metrics_by_population = metrics_by_population,
      CS99_metrics = cs_metrics_by_population,
      jam.nM.iter = as.numeric(jam.nM.iter),
      jam.iterations = as.numeric(jam.nM.iter) * 1e6,
      runtime = runtime,
      checkpoint_dir = core$checkpoint_dir,
      jam_dir = core$jam_dir
    )
  )
}
