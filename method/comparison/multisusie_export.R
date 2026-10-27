export_multisusie_inputs <- function(common_data, export_dir) {
dir.create(export_dir, recursive = TRUE, showWarnings = FALSE)
export_dir <- normalizePath(export_dir, mustWork = TRUE)

required <- c(
  "X_list", "Y_list", "B_true_list", "snp_key_vec",
  "K", "m", "q", "n_k"
)
missing_fields <- setdiff(required, names(common_data))
if (length(missing_fields) > 0L) {
  stop("common_data is missing: ", paste(missing_fields, collapse = ", "))
}

K <- as.integer(common_data$K)
m <- as.integer(common_data$m)
q <- as.integer(common_data$q)
n_k <- as.integer(common_data$n_k)

if (length(K) != 1L || K < 2L || length(m) != 1L || m < 1L ||
    length(q) != 1L || q < 1L || length(n_k) != K || any(n_k < 3L)) {
  stop("Invalid K, m, q, or n_k in common_data.")
}
if (length(common_data$X_list) != K ||
    length(common_data$Y_list) != K ||
    length(common_data$B_true_list) != K ||
    length(common_data$snp_key_vec) != m) {
  stop("common_data list lengths do not agree with K or m.")
}

snp_map <- data.frame(
  SNP = seq_len(m),
  variant_id = paste0("snp", seq_len(m)),
  SNP_key = as.character(common_data$snp_key_vec),
  stringsAsFactors = FALSE
)
write.csv(snp_map, file.path(export_dir, "snp_map.csv"), row.names = FALSE)

for (k in seq_len(K)) {
  X <- as.matrix(common_data$X_list[[k]])
  Y <- as.matrix(common_data$Y_list[[k]])
  B <- as.matrix(common_data$B_true_list[[k]])

  if (!identical(dim(X), c(n_k[k], m)) ||
      !identical(dim(Y), c(n_k[k], q)) ||
      !identical(dim(B), c(m, q))) {
    stop("Dimension mismatch in Group", k, ".")
  }
  if (any(!is.finite(X)) || any(!is.finite(Y)) || any(!is.finite(B))) {
    stop("Non-finite X, Y, or B_true value in Group", k, ".")
  }
  if (any(X < 0 | X > 2) || any(abs(X - round(X)) > 1e-8)) {
    stop("Group", k, " X is not coded as 0/1/2 allele dosage.")
  }
  if (any(apply(X, 2L, var) <= 0)) {
    stop("Group", k, " contains a monomorphic SNP.")
  }

  write.table(
    X,
    file.path(export_dir, paste0("X_Group", k, ".tsv")),
    sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE
  )
  write.table(
    Y,
    file.path(export_dir, paste0("Y_Group", k, ".tsv")),
    sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE
  )
  write.table(
    B,
    file.path(export_dir, paste0("B_true_Group", k, ".tsv")),
    sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE
  )
}

scenario_name <- if (is.null(common_data$scenario_name)) {
  "unspecified"
} else {
  as.character(common_data$scenario_name)
}

metadata <- data.frame(
  field = c("scenario_name", "K", "m", "q", "n_k"),
  value = c(scenario_name, K, m, q, paste(n_k, collapse = ",")),
  stringsAsFactors = FALSE
)

if (!is.null(common_data$config) && is.list(common_data$config)) {
  scalar_config <- vapply(
    common_data$config,
    function(x) length(x) == 1L && is.atomic(x) && !is.na(x),
    logical(1)
  )
  if (any(scalar_config)) {
    config_values <- common_data$config[scalar_config]
    metadata <- rbind(
      metadata,
      data.frame(
        field = paste0("config.", names(config_values)),
        value = vapply(config_values, as.character, character(1)),
        stringsAsFactors = FALSE
      )
    )
  }
}

write.csv(metadata, file.path(export_dir, "metadata.csv"), row.names = FALSE)
invisible(export_dir)
}
