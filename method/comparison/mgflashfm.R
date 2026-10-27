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
      repaired_traits <- c(repaired_trai×mzÖÚ$z{-®éÜj×7&÷72×÷VÆF–öâ6÷&R6†V6·ö–çBåÆâ"¢Ð ¢Öv×D52ÂÒ÷7G&ö6W75öÖvfÆ6†fÕö6÷&R€¢77FÒ77FÀ¢fÆ6†fÒæ÷WBÒfÆ6†fÒæ÷WBÀ¢–çWG2Ò–çWG0¢ ¢Æ—7B€¢Öv×D52ÒÖv×D52À¢6†V6·ö–çEöF—"Ò6†V6·ö–çEöF—"À¢¦ÕöF—"Ò¦ÕöF—"À¢6WGF–æw2Ò6WGF–æw0¢§Ð  ¦7&VF–&ÆU÷6WE÷—'2ÂÒgVæ7F–öâ†Öv×D52Â–çWG2Â6öÖÖöåöFF’°¢&÷w2ÂÒÆ—7B‚¢÷WEö’ÂÒÀ ¢f÷"‡"–â6WöÆöær†–çWG2GG&—EöæÖW2’’°¢G&—BÂÒ–çWG2GG&—EöæÖW5µ·%ÕÐ¢72ÂÒÖv×D52D577VÖÖ'•µ·G&—EÕÐ¢–b†—2æçVÆÂ†72’’72ÂÒÖv×D52D577VÖÖ'•µ·%ÕÐ¢–b‚—2æFFæg&ÖR†72’ÇÂç&÷r†72’ÓÒÂ’æW‡@ ¢–G2ÂÒ&÷væÖW2†72¢6çö–æFW‚ÂÒ&W6öÇfUöÖvfÆ6†fÕ÷6çö–G2€¢–G2À¢–çWG2À¢7FS‚&7&VF–&ÆR×6WB—'2f÷""ÂG&—B¢ ¢öæRÂÒFFæg&ÖR€¢4åÒ6çö–æFW‚À¢4åö¶W’Ò–çWG2G6çöÖE4åö¶W•·6çö–æFW…ÒÀ¢G&—Eö–æFW‚Ò"À¢G&—BÒG&—BÀ¢ÖtÕÒ2æçVÖW&–2†72FÖtÕ’À¢Öe÷ööÆVBÒ2æçVÖW&–2†72FÖe÷ööÆVB’À¢&÷öw&÷W2Ò2æçVÖW&–2†72G&÷öw&÷W2’À¢w&÷W2Ò2æ6†&7FW"†72Fw&÷W2’À¢–åô53“’ÒE%TRÀ¢—5÷G'VUö–å÷÷VÆF–öãÒ–çWG2GG'WF…öÖ6·5µ³ÅÕÕ¶6&–æB‡6çö–æFW‚Â"•ÒÀ¢7G&–æw44f7F÷'2ÒdÅ4P¢ ¢22&6·v&BÖ6ö×F–&ÆRæÖR&WF–æVBf÷"öÆBF÷vç7G&VÒ6öFRà¢öæRF—5÷G'VU÷6†&VE÷—"ÂÒöæRF—5÷G'VUö–å÷÷VÆF–öã ¢f÷"†²–â6WöÆVâ†6öÖÖöåöFFD²’’°¢öæUµ·7FS‚$%÷G'VUôw&÷W"Â²•ÕÒÂÐ¢6öÖÖöåöFFD%÷G'VUöÆ—7Eµ¶µÕÕ¶6&–æB‡6çö–æFW‚Â"•Ð¢öæUµ·7FS‚&—5÷G'VUôw&÷W"Â²•ÕÒÂÐ¢–çWG2GG'WF…öÖ6·5µ¶µÕÕ¶6&–æB‡6çö–æFW‚Â"•Ð¢Ð ¢&÷w5µ¶÷WEö•ÕÒÂÒöæP¢÷WEö’ÂÒ÷WEö’²À¢Ð ¢–b†ÆVæwF‚‡&÷w2’ÓÒÂ’°¢÷WBÂÒFFæg&ÖR€¢4åÒ–çFVvW"‚’À¢4åö¶W’Ò6†&7FW"‚’À¢G&—Eö–æFW‚Ò–çFVvW"‚’À¢G&—BÒ6†&7FW"‚’À¢ÖtÕÒçVÖW&–2‚’À¢Öe÷ööÆVBÒçVÖW&–2‚’À¢&÷öw&÷W2ÒçVÖW&–2‚’À¢w&÷W2Ò6†&7FW"‚’À¢–åô53“’ÒÆöv–6Â‚’À¢—5÷G'VUö–å÷÷VÆF–öãÒÆöv–6Â‚’À¢—5÷G'VU÷6†&VE÷—"ÒÆöv–6Â‚’À¢7G&–æw44f7F÷'2ÒdÅ4P¢¢&WGW&â†÷WB¢Ð ¢Fòæ6ÆÂ‡&&–æBÂ&÷w2§Ð  ¦W‡G&7EöÖv×÷6VÆV7F–öâÂÒgVæ7F–öâ€¢Öv×D52À¢–çWG2À¢6öÖÖöåöFFÀ¢F‡&W6†öÆBÒã“ ¢’° ¢ÒÂÒ2æ–çFVvW"†6öÖÖöåöFFFÒ¢ÂÒ2æ–çFVvW"†6öÖÖöåöFFG ¢ÖtÕÂÒÖG&—‚€¢À¢ç&÷rÒÒÀ¢æ6öÂÒÀ¢F–ÖæÖW2ÒÆ—7B†–çWG2G6çöÖFÖvfÆ6†fÕ÷'4”BÂ–çWG2GG&—EöæÖW2¢ ¢FWF–Ç2ÂÒÖv×D52D56FWF–ÂFFWF–Ç0¢–b‚—2æÆ—7B†FWF–Ç2’’°¢7F÷‚$ÔvfÆ6†fÒ÷WGWBFöW2æ÷B6öçF–â56FWF–ÂFFWF–Ç2â"¢Ð¢–b†ÆVæwF‚†FWF–Ç2’Ò’°¢7F÷‚$ÔvfÆ6†fÒ56FWF–ÂFFWF–Ç2FöW2æ÷B6öçF–âG&—G2â"¢Ð ¢f÷"‡"–â6WöÆVâ‡’’° ¢G&—BÂÒ–çWG2GG&—EöæÖW5µ·%ÕÐ¢FWF–ÂÂÒFWF–Ç5µ·G&—EÕÐ¢–b†—2æçVÆÂ†FWF–Â’’FWF–ÂÂÒFWF–Ç5µ·%ÕÐ ¢22–bF†R6÷&RÖWF†öB&WGW&ç2æòÖtÕö&¦V7Bf÷"G&—BÂ&W6W'fP¢22F†R÷&–v–æÂÆÂ×¦W&ò&÷rö6öÇVÖâ6öçG&–'WF–öââæöâÖf–æ—FRçVÖW&–0¢22ÖtÕfÇVW2&R7F–ÆÂG&VFVB2âW'&÷"&VÆ÷rà¢–b†—2æçVÆÂ†FWF–Â’ÇÂ—2æçVÆÂ†FWF–ÂFÖtÕ’’æW‡@ ¢fÇVW2ÂÒ2æçVÖW&–2†FWF–ÂFÖtÕ¢–G2ÂÒæÖW2†FWF–ÂFÖtÕ ¢–b†—2æçVÆÂ†–G2’ÇÂÆVæwF‚†–G2’ÒÆVæwF‚‡fÇVW2’’°¢7F÷‚$ÔvfÆ6†fÒÖtÕfÇVW2†fRæòfÆ–B4åæÖW2f÷""ÂG&—BÂ"â"¢Ð¢22çVÖW&–6Â&÷VæF'’FöÆW&æ6RöæÇ’âfÇVW2÷WG6–FR³ÂÒÖW&VÇ’'¢22fÆöF–ær×ö–çB&÷VæF–ær&R&ö¦V7FVB&6²FòF†R&ö&&–Æ—G’&ævRà¢Öv×÷FöÂÂÒ¢äÖ6†–æRFF÷V&ÆRæW0 ¢–b†ç’‚—2æf–æ—FR‡fÇVW2’’’°¢7F÷‚$ÔvfÆ6†fÒ&WGW&æVBæöâÖf–æ—FRÖtÕfÇVW2f÷""ÂG&—BÂ"â"¢Ð ¢Æ&vUö&BÂÒfÇVW2ÂÖÖv×÷FöÂÂfÇVW2â²Öv×÷FöÀ¢–b†ç’†Æ&vUö&B’’°¢7F÷€¢$ÔvfÆ6†fÒ&WGW&æVBÖFW&–ÆÇ’÷WBÖöb×&ævRÖtÕfÇVW2f÷""À¢G&—BÀ¢#¢&ævRÒ²"À¢f÷&ÖB†Ö–â‡fÇVW2’ÂF–v—G2Òr’À¢"Â"À¢f÷&ÖB†Ö‚‡fÇVW2’ÂF–v—G2Òr’À¢%Òâ ¢¢Ð ¢fÇVW2ÂÒÖ–âƒÂÖ‚ƒÂfÇVW2’ ¢22ÔvfÆ6†fÒ6âVÖ—B¦W&ò×&ö&&–Æ—G’VçG'’v†÷6RæÖR—2â7GVÀ¢22äö6†&7FW%ò6VçF–æVÂâ—B—2æ÷Bâ4åæB6ææ÷BffV7BF†P¢22ÖtÕâF‡&W6†öÆB6VÆV7F–öââG&÷öæÇ’¦W&ò×fÇVVBÖ—76–ærÖæÖP¢226VçF–æVÇ3²Ö—76–æræÖR6''––æræöç¦W&òÖ72&VÖ–ç2âW'&÷"à¢Ö—76–æuö–BÂÒ—2ææ†–G2’Âç¦6†"‡G&–×w2†–G2’¢–b†ç’†Ö—76–æuö–BbfÇVW2âÖv×÷FöÂ’’°¢7F÷€¢$ÔvfÆ6†fÒ&WGW&æVBÖ—76–ær4å”Bv—F‚æöç¦W&òÖtÕf÷""À¢G&—BÀ¢"â ¢¢Ð¢–b†ç’†Ö—76–æuö–B’’°¢fÇVW2ÂÒfÇVW5²Ö—76–æuö–EÐ¢–G2ÂÒ–G5²Ö—76–æuö–EÐ¢Ð¢–b†ÆVæwF‚‡fÇVW2’ÓÒÂ’æW‡@ ¢6çö–æFW‚ÂÒ&W6öÇfUöÖvfÆ6†fÕ÷6çö–G2€¢–G2À¢–çWG2À¢7FS‚&ÖtÕf÷""ÂG&—B¢ ¢–b†ç”GWÆ–6FVB‡6çö–æFW‚’’°¢fÇVW2ÂÒfÇ’€¢7Æ—B‡fÇVW2Â6çö–æFW‚’À¢Ö‚À¢çVÖW&–2ƒ’À¢æç&ÒÒE%TP¢¢6çö–æFW‚ÂÒ2æ–çFVvW"†æÖW2‡fÇVW2’¢Ð ¢ÖtÕ·6çö–æFW‚Â%ÒÂÒfÇVW0¢Ð ¢¦†BÂÒÖG&—‚€¢2æ–çFVvW"†ÖtÕâF‡&W6†öÆB’À¢ç&÷rÒÒÀ¢æ6öÂÒÀ¢F–ÖæÖW2ÒF–ÖæÖW2†ÖtÕ¢ ¢ÆÅöw&–BÂÒW‡æBæw&–B€¢4åÒ6WöÆVâ†Ò’À¢G&—Eö–æFW‚Ò6WöÆVâ‡’À¢´TUäõUBäEE%2ÒdÅ4RÀ¢7G&–æw44f7F÷'2ÒdÅ4P¢ ¢ÆÅöWf–FVæ6RÂÒFFæg&ÖR€¢4åÒÆÅöw&–BE4åÀ¢4åö¶W’Ò–çWG2G6çöÖE4åö¶W•¶ÆÅöw&–BE4åÒÀ¢G&—Eö–æFW‚ÒÆÅöw&–BGG&—Eö–æFW‚À¢G&—BÒ–çWG2GG&—EöæÖW5¶ÆÅöw&–BGG&—Eö–æFW…ÒÀ¢ÖtÕÒÖtÕ¶6&–æB†ÆÅöw&–BE4åÂÆÅöw&–BGG&—Eö–æFW‚•ÒÀ¢6VÆV7FVBÒ¦†E¶6&–æB†ÆÅöw&–BE4åÂÆÅöw&–BGG&—Eö–æFW‚•ÒÓÒÂÀ¢7G&–æw44f7F÷'2ÒdÅ4P¢ ¢6VÆV7FVE÷—'2ÂÒÆÅöWf–FVæ6U°¢ÆÅöWf–FVæ6RG6VÆV7FVBÀ¢À¢G&÷ÒdÅ4P¢Ð ¢Æ—7B€¢ÖtÕÒÖtÕÀ¢¦†BÒ¦†BÀ¢ÆÅöWf–FVæ6RÒÆÅöWf–FVæ6RÀ¢6VÆV7FVE÷—'2Ò6VÆV7FVE÷—'2À¢F‡&W6†öÆBÒF‡&W6†öÆ@¢§Ð  ¦WfÇVFU÷—%÷6VÆV7F–öâÂÒgVæ7F–öâ…¦†BÂG'WF‚Â÷VÆF–öâÒäö–çFVvW%ò’° ¢&VBÂÒ2æÖG&—‚…¦†BÒ¢G'WF‚ÂÒ2æÖG&—‚‡G'WF‚Ò ¢–b‚–FVçF–6Â†F–Ò‡&VB’ÂF–Ò‡G'WF‚’’’°¢7F÷‚%&VF–7F–öâæBG'WF‚ÖG&–6W2†fRF–ffW&VçBF–ÖVç6–öç2â"¢Ð ¢EÂÒ7VÒ‡&VBbG'WF‚¢eÂÒ7VÒ‡&VBbG'WF‚¢dâÂÒ7VÒ‚&VBbG'WF‚¢DâÂÒ7VÒ‚&VBbG'WF‚ ¢&V6—6–öâÂÒ–b…E²eÓÒÂ’VÇ6REò…E²e¢&V6ÆÂÂÒ–b…E²dâÓÒÂ’VÇ6REò…E²dâ¢7V6–f–6—G’ÂÒ–b…Dâ²eÓÒÂ’VÇ6RDâò…Dâ²e¢cÂÒ–b‡&V6—6–öâ²&V6ÆÂÓÒ’°¢ ¢ÒVÇ6R°¢"¢&V6—6–öâ¢&V6ÆÂò‡&V6—6–öâ²&V6ÆÂ¢Ð ¢FFæg&ÖR€¢ÖWF†öBÒ$ÔvfÆ6†fÒ"À¢÷VÆF–öâÒ2æ–çFVvW"‡÷VÆF–öâ’À¢6VÆV7F–öå÷'VÆRÒ&ÖtÕâã“"À¢EÒEÀ¢eÒeÀ¢dâÒdâÀ¢DâÒDâÀ¢&V6—6–öâÒ&V6—6–öâÀ¢&V6ÆÂÒ&V6ÆÂÀ¢7V6–f–6—G’Ò7V6–f–6—G’À¢cÒcÀ¢7G&–æw44f7F÷'2ÒdÅ4P¢§Ð  ¢22&6·v&BÖ6ö×F–&ÆR†VÇW"&WF–æVBf÷"öÆB6ÆÆW'2à¦WfÇVFU÷6†&VE÷—'2ÂÒgVæ7F–öâ…¦†BÂG'WF‚’°¢÷WBÂÒWfÇVFU÷—%÷6VÆV7F–öâ…¦†BÂG'WF‚Â÷VÆF–öâÒÂ¢÷WBG÷VÆF–öâÂÒåTÄÀ¢÷W@§Ð  ¦'V–ÆE÷G'WF…÷—'5ö'•÷÷VÆF–öâÂÒgVæ7F–öâ€¢6öÖÖöåöFFÀ¢–çWG2À¢6VÆV7FVE÷—'0¢’° ¢²ÂÒ2æ–çFVvW"†6öÖÖöåöFFD²¢&÷w2ÂÒfV7F÷"‚&Æ—7B"Â² ¢&VF–7FVEö¶W’ÂÒ7FR€¢6VÆV7FVE÷—'2E4åÀ¢6VÆV7FVE÷—'2GG&—Eö–æFW‚À¢6WÒ#¢ ¢ ¢f÷"†²–â6WöÆVâ„²’’° ¢G'WF…ö–G‚ÂÒv†–6‚€¢–çWG2GG'WF…öÖ6·5µ¶µÕÒÀ¢'"æ–æBÒE%TP¢ ¢G'WF…÷—'2ÂÒFFæg&ÖR€¢÷VÆF–öâÒ²À¢4åÒG'WF…ö–G…²ÂÅÒÀ¢4åö¶W’Ò6öÖÖöåöFFG6çö¶W•÷fV5·G'WF…ö–G…²ÂÅÕÒÀ¢G&—Eö–æFW‚ÒG'WF…ö–G…²Â$ÅÒÀ¢G&—BÒ–çWG2GG&—EöæÖW5·G'WF…ö–G…²Â$ÅÕÒÀ¢%÷G'VRÒ6öÖÖöåöFFD%÷G'VUöÆ—7Eµ¶µÕÕ·G'WF…ö–G…ÒÀ¢7G&–æw44f7F÷'2ÒdÅ4P¢ ¢G'WF…ö¶W’ÂÒ7FR€¢G'WF…÷—'2E4åÀ¢G'WF…÷—'2GG&—Eö–æFW‚À¢6WÒ#¢ ¢ ¢G'WF…÷—'2FFWFV7FVEöÖtÕ“ÂÒG'WF…ö¶W’V–âR&VF–7FVEö¶W¢&÷w5µ¶µÕÒÂÒG'WF…÷—'0¢Ð ¢Fòæ6ÆÂ‡&&–æBÂ&÷w2§Ð  §'VåöÖvfÆ6†fÒÂÒgVæ7F–öâ€¢6öÖÖöåöf–ÆRÀ¢&W7VÇEöF—"À¢¦ÒæäÒæ—FW"ÒvWD÷F–öâ‚$ÔvfÆ6†fÒæ¦ÒæäÒæ—FW""Â’À¢¦Òç6VVBÒvWD÷F–öâ‚$ÔvfÆ6†fÒæ¦Òç6VVB"Â##c“#”Â’À¢&W7VÖRÒE%TP¢’° ¢–b‚&WV—&TæÖW76R‚$ÔvfÆ6†fÒ"ÂV–WFÇ’ÒE%TR’’°¢7F÷‚%6¶vRÔvfÆ6†fÒ—2æ÷B–ç7FÆÆVBâ"¢Ð¢–b‚&WV—&TæÖW76R‚%#$$tÆ”Õ2"ÂV–WFÇ’ÒE%TR’’°¢7F÷‚%6¶vR#$$tÆ”Õ2—2æ÷B–ç7FÆÆVC²¤ÒÇ6ò&WV—&W2¦fâ"¢Ð ¢–b€¢ÆVæwF‚†¦ÒæäÒæ—FW"’ÒÂÇÀ¢—2æçVÖW&–2†¦ÒæäÒæ—FW"’ÇÀ¢—2æf–æ—FR†¦ÒæäÒæ—FW"’ÇÀ¢¦ÒæäÒæ—FW"ÃÒ ¢’°¢7F÷‚&¦ÒæäÒæ—FW"×W7B&RöæR÷6—F—fRf–æ—FRçVÖ&W"â"¢Ð ¢–b€¢ÆVæwF‚†¦Òç6VVB’ÒÂÇÀ¢—2æçVÖW&–2†¦Òç6VVB’ÇÀ¢—2æf–æ—FR†¦Òç6VVB¢’°¢7F÷‚&¦Òç6VVB×W7B&RöæRf–æ—FRçVÖW&–2fÇVRâ"¢Ð ¢7W&W756¶vU7F'GWÖW76vW2†Æ–'&'’…#$$tÆ”Õ2’¢7W&W756¶vU7F'GWÖW76vW2†Æ–'&'’„ÔvfÆ6†fÒ’ ¢6öÖÖöåöf–ÆRÂÒæ÷&ÖÆ—¦UF‚†6öÖÖöåöf–ÆRÂ×W7Ev÷&²ÒE%TR¢F—"æ7&VFR‡&W7VÇEöF—"Â&V7W'6—fRÒE%TRÂ6†÷uv&æ–æw2ÒdÅ4R¢&W7VÇEöF—"ÂÒæ÷&ÖÆ—¦UF‚‡&W7VÇEöF—"Â×W7Ev÷&²ÒE%TR ¢6öÖÖöåöFFÂÒ&VE$E2†6öÖÖöåöf–ÆR¢–çWG2ÂÒÖ¶UöÖvfÆ6†fÕö–çWG2†6öÖÖöåöFF ¢²ÂÒ2æ–çFVvW"†6öÖÖöåöFFD²¢ÒÂÒ2æ–çFVvW"†6öÖÖöåöFFFÒ¢ÂÒ2æ–çFVvW"†6öÖÖöåöFFG ¢6B‚%ÆãÓÓÓÓÓÓÓÓÓÓÓÓÓÓÓÓÓÓÓÓÓÓÓÓÓÓÓÓÓÓÓÓÓÓÓÓÓÓÓÓÓÓÓÓÓÓÓÓÓÓÓÓÓÓÓÓÓÓÓÕÆâ"¢6B‚%'Vææ–ærÔvfÆ6†fÕÆâ"¢6B‚#ÓÓÓÓÓÓÓÓÓÓÓÓÓÓÓÓÓÓÓÓÓÓÓÓÓÓÓÓÓÓÓÓÓÓÓÓÓÓÓÓÓÓÓÓÓÓÓÓÓÓÓÓÓÓÓÓÓÓÓÕÆâ"¢6B‚%66Væ&–ó¢"Â6öÖÖöåöFFG66Væ&–õöæÖRÂ%Æâ"¢6B‚$²Ò"Â²Â'ÂÒÒ"ÂÒÂ'ÂÒ"ÂÂ%Æâ"¢6B€¢&¦ÒæäÒæ—FW"Ò"Â¦ÒæäÒæ—FW"À¢&Ö–ÆÆ–öâÒ"À¢f÷&ÖB†¦ÒæäÒæ—FW"¢SbÂ66–VçF–f–2ÒdÅ4R’À¢&—FW&F–öç2W"¤Ò'VåÆâ ¢¢6B‚$¤Ò"6VVBÒ"Â2æ–çFVvW"†¦Òç6VVB’Â%Æâ"¢6B‚'&W7VÖRÒ"Â—5E%TR‡&W7VÖR’Â%Æâ" ¢w&—FRæ77b€¢–çWG2G6çöÖÀ¢f–ÆRçF‚‡&W7VÇEöF—"Â'6çö–EöÖæ77b"’À¢&÷rææÖW2ÒdÅ4P¢ ¢6÷&U÷F–ÖRÂÒ7—7FVÒçF–ÖR‡°¢6÷&RÂÒ'VåöÖvfÆ6†fÕö6÷&Uö6†V6·ö–çFVB€¢–çWG2Ò–çWG2À¢&W7VÇEöF—"Ò&W7VÇEöF—"À¢¦ÒæäÒæ—FW"Ò¦ÒæäÒæ—FW"À¢¦Òç6VVBÒ2æ–çFVvW"†¦Òç6VVB’À¢7Òã“’À¢7&VBÒã“’À¢DöFG2ÒÀ¢Ö†7bÒÀ¢Ö†7e÷7F÷Ò#À¢#"ÒãƒÀ¢ä4õ$U2ÒÀ¢W‡G&æ¦fæ&wVÖVçG2ÒåTÄÂÀ¢&W7VÖRÒ&W7VÖRÀ¢66Væ&–õöæÖRÒ6öÖÖöåöFFG66Væ&–õöæÖP¢¢Ò ¢Öv×D52ÂÒ6÷&RFÖv×D50 ¢226fRF†R6ö×ÆWFRöff–6–ÂÖÆ–÷WB&W7VÇB&Vf÷&Rç’ÆFW"F–væ÷7F–0¢22÷"55b6öFRâF†—2&WfVçG2Æ÷6–ærâW‡Vç6—fRf—BFò÷7B×&ö6W76–ærà¢6fU$E2€¢Öv×D52À¢f–ÆRçF‚‡&W7VÇEöF—"Â$ÔvfÆ6†fÕööff–6–Å÷&W7VÇBç&G2"¢ ¢6VÆV7F–öâÂÒW‡G&7EöÖv×÷6VÆV7F–öâ€¢Öv×D52ÒÖv×D52À¢–çWG2Ò–çWG2À¢6öÖÖöåöFFÒ6öÖÖöåöFFÀ¢F‡&W6†öÆBÒã“ ¢ ¢226fR6VÆV7F–öâ–ÖÖVF–FVÇ’2æ÷F†W"&V6÷fW'’ö–çBà¢6fU$E2€¢6VÆV7F–öâFÖtÕÀ¢f–ÆRçF‚‡&W7VÇEöF—"Â$ÔvfÆ6†fÕöÖtÕöÖG&—‚ç&G2"¢¢6fU$E2€¢6VÆV7F–öâE¦†BÀ¢f–ÆRçF‚‡&W7VÇEöF—"Â$ÔvfÆ6†fÕõ¦†EöÖtÕ“ç&G2"¢ ¢75÷—'2ÂÒ7&VF–&ÆU÷6WE÷—'2€¢Öv×D52ÒÖv×D52À¢–çWG2Ò–çWG2À¢6öÖÖöåöFFÒ6öÖÖöåöFF¢ ¢ÖWG&–75ö'•÷÷VÆF–öâÂÒFòæ6ÆÂ€¢&&–æBÀ¢ÆÇ’€¢6WöÆVâ„²’À¢gVæ7F–öâ†²’°¢WfÇVFU÷—%÷6VÆV7F–öâ€¢¦†BÒ6VÆV7F–öâE¦†BÀ¢G'WF‚Ò–çWG2GG'WF…öÖ6·5µ¶µÕÒÀ¢÷VÆF–öâÒ°¢¢Ð¢¢ ¢G'WF…÷—'2ÂÒ'V–ÆE÷G'WF…÷—'5ö'•÷÷VÆF–öâ€¢6öÖÖöåöFFÒ6öÖÖöåöFFÀ¢–çWG2Ò–çWG2À¢6VÆV7FVE÷—'2Ò6VÆV7F–öâG6VÆV7FVE÷—'0¢ ¢227&VF–&ÆR×6WB6÷fW&vR'’÷VÆF–öââF†—2&VÖ–ç2F–væ÷7F–2öæÇ“°¢22F†R6–×VÆF–öâ6VÆV7F–öâ'VÆR—2ÖtÕâã“à¢75ö¶W’ÂÒ7FR†75÷—'2E4åÂ75÷—'2GG&—Eö–æFW‚Â6WÒ#¢" ¢75öÖWG&–75ö'•÷÷VÆF–öâÂÒFòæ6ÆÂ€¢&&–æBÀ¢ÆÇ’€¢6WöÆVâ„²’À¢gVæ7F–öâ†²’°¢G'WF…ö²ÂÒG'WF…÷—'5·G'WF…÷—'2G÷VÆF–öâÓÒ²ÂÂG&÷ÒdÅ4UÐ¢G'WF…ö¶W’ÂÒ7FR‡G'WF…ö²E4åÂG'WF…ö²GG&—Eö–æFW‚Â6WÒ#¢" ¢FFæg&ÖR€¢ÖWF†öBÒ$ÔvfÆ6†fÕô53“’"À¢÷VÆF–öâÒ²À¢W'÷6RÒ&7&VF–&ÆU÷6WEö6÷fW&vU÷&W6öÇWF–öåööæÇ’"À¢åô53“•õ4å÷G&—E÷—'2Òç&÷r†75÷—'2’À¢å÷G'VU÷—'5ö6÷fW&VBÒ7VÒ‡G'WF…ö¶W’V–âR75ö¶W’’À¢G'WF…÷—%ö6÷fW&vRÒ–b†ÆVæwF‚‡G'WF…ö¶W’’ÓÒÂ’°¢ä÷&VÅð¢ÒVÇ6R°¢ÖVâ‡G'WF…ö¶W’V–âR75ö¶W’¢ÒÀ¢ÖVåô53“•÷6—¦RÒ–b†ç&÷r†75÷—'2’ÓÒÂ’°¢ ¢ÒVÇ6R°¢ÖVâ†2æçVÖW&–2‡F&ÆR†75÷—'2GG&—Eö–æFW‚’’¢ÒÀ¢7G&–æw44f7F÷'2ÒdÅ4P¢¢Ð¢¢ ¢22æöæW76VçF–Â‡VÖâ×&VF&ÆR÷WGWG2&Rw&—GFVâöæÇ’gFW"ÆÂ$E0¢22&V6÷fW'’ö–çG2&÷fRW†—7Bà¢w&—FRæ77b€¢75÷—'2À¢f–ÆRçF‚‡&W7VÇEöF—"Â$ÔvfÆ6†fÕô53“•õ4å÷G&—E÷—'2æ77b"’À¢&÷rææÖW2ÒdÅ4P¢¢w&—FRæ77b€¢6VÆV7F–öâFÖtÕÀ¢f–ÆRçF‚‡&W7VÇEöF—"Â$ÔvfÆ6†fÕöÖtÕöÖG&—‚æ77b"’À¢&÷rææÖW2ÒdÅ4P¢¢w&—FRæ77b€¢6VÆV7F–öâE¦†BÀ¢f–ÆRçF‚‡&W7VÇEöF—"Â$ÔvfÆ6†fÕõ¦†EöÖtÕ“æ77b"’À¢&÷rææÖW2ÒdÅ4P¢¢w&—FRæ77b€¢6VÆV7F–öâFÆÅöWf–FVæ6RÀ¢f–ÆRçF‚‡&W7VÇEöF—"Â$ÔvfÆ6†fÕöÆÅõ4å÷G&—EöÖtÕæ77b"’À¢&÷rææÖW2ÒdÅ4P¢¢w&—FRæ77b€¢6VÆV7F–öâG6VÆV7FVE÷—'2À¢f–ÆRçF‚‡&W7VÇEöF—"Â$ÔvfÆ6†fÕöÖtÕ“õ4å÷G&—E÷—'2æ77b"’À¢&÷rææÖW2ÒdÅ4P¢¢w&—FRæ77b€¢G'WF…÷—'2À¢f–ÆRçF‚‡&W7VÇEöF—"Â'G'WF…õ4å÷G&—E÷—'5ö'•÷÷VÆF–öâæ77b"’À¢&÷rææÖW2ÒdÅ4P¢¢w&—FRæ77b€¢ÖWG&–75ö'•÷÷VÆF–öâÀ¢f–ÆRçF‚‡&W7VÇEöF—"Â$ÔvfÆ6†fÕ÷—%öÖWG&–72æ77b"’À¢&÷rææÖW2ÒdÅ4P¢¢w&—FRæ77b€¢75öÖWG&–75ö'•÷÷VÆF–öâÀ¢f–ÆRçF‚‡&W7VÇEöF—"Â$ÔvfÆ6†fÕô53“•ö6÷fW&vU÷&W6öÇWF–öâæ77b"’À¢&÷rææÖW2ÒdÅ4P¢ ¢'VçF–ÖRÂÒÆ—7B€¢6÷&UöVÆ6VE÷6V2ÒVææÖR†6÷&U÷F–ÖU²&VÆ6VB%Ò’À¢6÷&U÷W6W%÷6V2ÒVææÖR†6÷&U÷F–ÖU²'W6W"ç6VÆb%Ò’À¢6÷&U÷7—7FVÕ÷6V2ÒVææÖR†6÷&U÷F–ÖU²'7—2ç6VÆb%Ò¢ ¢6fU$E2€¢Æ—7B€¢66Væ&–õöæÖRÒ6öÖÖöåöFFG66Væ&–õöæÖRÀ¢²Ò²À¢ÒÒÒÀ¢ÒÀ¢¦ÒæäÒæ—FW"Ò2æçVÖW&–2†¦ÒæäÒæ—FW"’À¢¦Òæ—FW&F–öç2Ò2æçVÖW&–2†¦ÒæäÒæ—FW"’¢SbÀ¢¦Òç6VVBÒ2æ–çFVvW"†¦Òç6VVB’À¢'VçF–ÖRÒ'VçF–ÖRÀ¢6†V6·ö–çEöF—"Ò6÷&RF6†V6·ö–çEöF—"À¢¦ÕöF—"Ò6÷&RF¦ÕöF—"À¢ÖWG&–75ö'•÷÷VÆF–öâÒÖWG&–75ö'•÷÷VÆF–öà¢’À¢f–ÆRçF‚‡&W7VÇEöF—"Â$ÔvfÆ6†fÕ÷'VåöÖWFFFç&G2"¢ ¢6B‚%ÆäÔvfÆ6†fÒ6ö×ÆWFVBåÆâ"¢6B‚&ÖtÕâã“4å×G&—B—'3¢"Âç&÷r‡6VÆV7F–öâG6VÆV7FVE÷—'2’Â%Æâ"¢6B‚$53“’4å×G&—B—'3¢"Âç&÷r†75÷—'2’Â%Æâ"¢6B‚$6÷&RVÆ6VB6V6öæG3¢"Â'VçF–ÖRF6÷&UöVÆ6VE÷6V2Â%Æâ"¢6B‚%&W7VÇG3¢"Â&W7VÇEöF—"Â%Æâ" ¢–çf—6–&ÆR€¢Æ—7B€¢&W7VÇBÒÖv×D52À¢ÖtÕÒ6VÆV7F–öâFÖtÕÀ¢¦†BÒ6VÆV7F–öâE¦†BÀ¢6VÆV7FVE÷—'2Ò6VÆV7F–öâG6VÆV7FVE÷—'2À¢53“•÷—'2Ò75÷—'2À¢G'WF‚ÒG'WF…÷—'2À¢ÖWG&–72ÒÖWG&–75ö'•÷÷VÆF–öâÀ¢ÖWG&–75ö'•÷÷VÆF–öâÒÖWG&–75ö'•÷÷VÆF–öâÀ¢53“•öÖWG&–72Ò75öÖWG&–75ö'•÷÷VÆF–öâÀ¢¦ÒæäÒæ—FW"Ò2æçVÖW&–2†¦ÒæäÒæ—FW"’À¢¦Òæ—FW&F–öç2Ò2æçVÖW&–2†¦ÒæäÒæ—FW"’¢SbÀ¢'VçF–ÖRÒ'VçF–ÖRÀ¢6†V6·ö–çEöF—"Ò6÷&RF6†V6·ö–çEöF—"À¢¦ÕöF—"Ò6÷&RF¦ÕöF— ¢¢§Ð