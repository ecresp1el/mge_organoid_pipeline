#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
})

timestamp <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")
log_msg <- function(...) message("[R ", timestamp(), "] ", paste0(..., collapse = ""))

parse_args <- function(args) {
  out <- list(
    reference = NULL,
    query = NULL,
    labels = NULL,
    outdir = NULL,
    reference_assay = "RNA",
    query_assay = "RNA",
    transfer_name = "shi_label",
    normalization_method = "LogNormalize",
    dims = "50",
    nfeatures = "3000",
    npcs = "50",
    seed = "0",
    output_prefix = "div30"
  )
  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, "--")) stop("Unknown argument: ", key, call. = FALSE)
    name <- substring(key, 3L)
    if (!(name %in% names(out))) stop("Unknown argument: ", key, call. = FALSE)
    if (i == length(args)) stop("Missing value for argument: ", key, call. = FALSE)
    out[[name]] <- args[[i + 1L]]
    i <- i + 2L
  }
  out
}

read_tsv <- function(path) {
  utils::read.delim(path, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
}

write_tsv <- function(x, path) {
  utils::write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
}

write_tsv_gz <- function(x, path) {
  con <- gzfile(path, open = "wt")
  on.exit(close(con), add = TRUE)
  utils::write.table(x, con, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
}

join_layers_if_needed <- function(obj, assay) {
  if (
    inherits(obj[[assay]], "Assay5") &&
    exists("JoinLayers", where = asNamespace("SeuratObject"), inherits = FALSE)
  ) {
    log_msg("Joining Assay5 layers for assay: ", assay)
    obj <- SeuratObject::JoinLayers(obj, assay = assay)
  }
  obj
}

get_data_like_matrix <- function(obj, assay) {
  last_error <- NULL
  for (layer in c("data", "logcounts", "counts")) {
    mat <- tryCatch(
      SeuratObject::GetAssayData(obj, assay = assay, layer = layer),
      error = function(e) {
        last_error <<- e
        NULL
      }
    )
    if (!is.null(mat)) {
      return(list(matrix = mat, layer = layer))
    }
  }
  stop("Could not extract data-like matrix for assay ", assay, ": ", conditionMessage(last_error), call. = FALSE)
}

key_map <- function(labels, key_col, value_col = "shi_label") {
  if (!(key_col %in% colnames(labels))) return(setNames(character(0), character(0)))
  if (!(value_col %in% colnames(labels))) return(setNames(character(0), character(0)))
  keys <- as.character(labels[[key_col]])
  values <- as.character(labels[[value_col]])
  keep <- !is.na(keys) & nzchar(keys) & !is.na(values) & nzchar(values)
  key_counts <- table(keys[keep])
  unique_keys <- names(key_counts)[key_counts == 1L]
  keep <- keep & keys %in% unique_keys
  stats::setNames(values[keep], keys[keep])
}

attach_reference_labels <- function(reference, labels) {
  maps <- list(
    colname = key_map(labels, "reference_obs_name")
  )
  week_maps <- list(
    colname = key_map(labels, "reference_obs_name", "shi_week_label")
  )
  week_numeric_maps <- list(
    colname = key_map(labels, "reference_obs_name", "shi_week_numeric")
  )
  metadata <- reference@meta.data
  result <- rep(NA_character_, ncol(reference))
  source <- rep(NA_character_, ncol(reference))
  week_result <- rep(NA_character_, ncol(reference))
  week_numeric_result <- rep(NA_character_, ncol(reference))
  cells <- colnames(reference)

  for (i in seq_along(cells)) {
    candidates <- c(cells[[i]])
    names(candidates) <- "colname"
    for (j in seq_along(candidates)) {
      key <- candidates[[j]]
      map_name <- names(candidates)[[j]]
      if (!is.na(key) && nzchar(key) && key %in% names(maps[[map_name]])) {
        result[[i]] <- unname(maps[[map_name]][[key]])
        source[[i]] <- map_name
        if (key %in% names(week_maps[[map_name]])) {
          week_result[[i]] <- unname(week_maps[[map_name]][[key]])
        }
        if (key %in% names(week_numeric_maps[[map_name]])) {
          week_numeric_result[[i]] <- unname(week_numeric_maps[[map_name]][[key]])
        }
        break
      }
    }
  }
  reference$shi_transfer_label <- result
  reference$shi_transfer_label_source <- source
  reference$shi_transfer_week_label <- week_result
  reference$shi_transfer_week_numeric <- suppressWarnings(as.numeric(week_numeric_result))
  reference
}

ensure_log_normalized <- function(obj, assay, nfeatures, verbose = FALSE) {
  Seurat::DefaultAssay(obj) <- assay
  data_res <- get_data_like_matrix(obj, assay)
  if (data_res$layer == "counts") {
    log_msg("No data/logcounts layer found for assay ", assay, "; running NormalizeData")
    obj <- Seurat::NormalizeData(obj, assay = assay, normalization.method = "LogNormalize", verbose = verbose)
  } else {
    log_msg("Using existing ", data_res$layer, " layer for assay ", assay)
  }
  if (length(Seurat::VariableFeatures(obj, assay = assay)) == 0) {
    log_msg("No variable features found for assay ", assay, "; running FindVariableFeatures")
    obj <- Seurat::FindVariableFeatures(obj, assay = assay, nfeatures = nfeatures, verbose = verbose)
  }
  obj
}

sanitize_token <- function(x) {
  y <- gsub("[^A-Za-z0-9]+", "_", trimws(as.character(x)))
  y <- gsub("_+", "_", y)
  y <- gsub("^_|_$", "", y)
  ifelse(nzchar(y), y, "value")
}

prediction_score_cols <- function(df) {
  setdiff(grep("^prediction\\.score\\.", colnames(df), value = TRUE), "prediction.score.max")
}

rename_score_columns <- function(df, prefix, score_cols) {
  out <- df[, c("cell_id", score_cols), drop = FALSE]
  new_names <- colnames(out)
  for (col in score_cols) {
    label <- sub("^prediction\\.score\\.", "", col)
    new_names[new_names == col] <- paste0(prefix, sanitize_token(label))
  }
  colnames(out) <- new_names
  out
}

anchor_query_cell_count <- function(anchors) {
  anchor_df <- as.data.frame(anchors@anchors)
  if ("cell2" %in% colnames(anchor_df)) {
    return(length(unique(anchor_df$cell2)))
  }
  NA_integer_
}

choose_transfer_k_weight <- function(n_anchors, requested = 50L) {
  if (n_anchors <= requested) {
    used <- max(1L, n_anchors - 1L)
    reason <- paste0("lowered_from_", requested, "_because_anchor_count_was_", n_anchors)
  } else {
    used <- requested
    reason <- "requested_k_weight_used"
  }
  list(used = used, reason = reason, requested = requested)
}

validate_transfer_predictions <- function(predictions, transfer_label) {
  required <- c("cell_id", "predicted.id", "prediction.score.max")
  missing <- setdiff(required, colnames(predictions))
  if (length(missing) > 0) {
    stop(transfer_label, " TransferData output is missing: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  score_cols <- prediction_score_cols(predictions)
  if (length(score_cols) == 0) {
    stop(transfer_label, " TransferData returned no per-label score columns", call. = FALSE)
  }
  max_score <- as.numeric(predictions$prediction.score.max)
  if (any(!is.finite(max_score))) {
    stop(transfer_label, " prediction.score.max contains non-finite values", call. = FALSE)
  }
  if (any(max_score < -1e-8 | max_score > 1 + 1e-8)) {
    stop(transfer_label, " prediction.score.max contains values outside [0,1]", call. = FALSE)
  }
  score_mat <- as.matrix(predictions[, score_cols, drop = FALSE])
  storage.mode(score_mat) <- "double"
  if (any(!is.finite(score_mat))) {
    stop(transfer_label, " per-label prediction scores contain non-finite values", call. = FALSE)
  }
  if (any(abs(max_score - apply(score_mat, 1, max)) > 1e-8)) {
    stop(transfer_label, " prediction.score.max is not the row-wise max of per-label scores", call. = FALSE)
  }
  score_cols
}

validate_canonical_scores <- function(scores, max_score, uncertainty, score_label) {
  score_cols <- setdiff(colnames(scores), "cell_id")
  if (length(score_cols) == 0) stop("No canonical ", score_label, " score columns to validate", call. = FALSE)
  score_mat <- as.matrix(scores[, score_cols, drop = FALSE])
  storage.mode(score_mat) <- "double"
  max_score <- as.numeric(max_score)
  uncertainty <- as.numeric(uncertainty)
  if (any(!is.finite(score_mat)) || any(!is.finite(max_score)) || any(!is.finite(uncertainty))) {
    stop("Non-finite canonical ", score_label, " scores", call. = FALSE)
  }
  if (any(abs(max_score - apply(score_mat, 1, max)) > 1e-8)) {
    stop("Canonical ", score_label, " max score is not the row-wise max of per-label scores", call. = FALSE)
  }
  if (any(abs(uncertainty - (1 - max_score)) > 1e-8)) {
    stop("Canonical ", score_label, " uncertainty is not 1 - max score", call. = FALSE)
  }
}

opt <- parse_args(commandArgs(trailingOnly = TRUE))
required <- c("reference", "query", "labels", "outdir")
for (name in required) {
  if (is.null(opt[[name]]) || !nzchar(opt[[name]])) stop("--", name, " is required", call. = FALSE)
}
for (path in c(opt$reference, opt$query, opt$labels)) {
  if (!file.exists(path)) stop("Missing input: ", path, call. = FALSE)
}
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

dims_n <- as.integer(opt$dims)
nfeatures <- as.integer(opt$nfeatures)
npcs <- as.integer(opt$npcs)
seed <- as.integer(opt$seed)
set.seed(seed)

log_msg("R version: ", R.version.string)
log_msg("Seurat version: ", as.character(utils::packageVersion("Seurat")))
log_msg("SeuratObject version: ", as.character(utils::packageVersion("SeuratObject")))
log_msg("Reference RDS: ", opt$reference)
log_msg("Query RDS: ", opt$query)
log_msg("Labels TSV: ", opt$labels)
log_msg("Output dir: ", opt$outdir)
log_msg("Output prefix: ", opt$output_prefix)

labels <- read_tsv(opt$labels)
if (!("shi_label" %in% colnames(labels))) stop("Labels TSV lacks shi_label column", call. = FALSE)
labels <- labels[!is.na(labels$shi_label) & nzchar(labels$shi_label), , drop = FALSE]
log_msg("Usable label rows: ", nrow(labels))

log_msg("Reading reference Seurat object")
reference <- readRDS(opt$reference)
if (!inherits(reference, "Seurat")) stop("Reference is not a Seurat object", call. = FALSE)
if (!(opt$reference_assay %in% Seurat::Assays(reference))) {
  stop("Reference assay missing: ", opt$reference_assay, call. = FALSE)
}
reference <- join_layers_if_needed(reference, opt$reference_assay)
reference <- attach_reference_labels(reference, labels)
labelled_cells <- colnames(reference)[
  !is.na(reference$shi_transfer_label) &
    nzchar(reference$shi_transfer_label) &
    !is.na(reference$shi_transfer_week_label) &
    nzchar(reference$shi_transfer_week_label)
]
if (length(labelled_cells) < 10) {
  diagnostics <- data.frame(metric = c("reference_cells", "labelled_reference_cells"), value = c(ncol(reference), length(labelled_cells)))
  write_tsv(diagnostics, file.path(opt$outdir, paste0(opt$output_prefix, "_shi_seurat_full_transfer_diagnostics.tsv")))
  stop("Fewer than 10 reference cells received Shi labels in Seurat object", call. = FALSE)
}
reference <- subset(reference, cells = labelled_cells)
reference$shi_transfer_label <- factor(reference$shi_transfer_label)
reference$shi_transfer_week_label <- factor(reference$shi_transfer_week_label)
log_msg("Reference cells after label subset: ", ncol(reference))
log_msg("Reference labels: ", paste(levels(reference$shi_transfer_label), collapse = ", "))
log_msg("Reference week labels: ", paste(levels(reference$shi_transfer_week_label), collapse = ", "))

log_msg("Reading query Seurat object")
query <- readRDS(opt$query)
if (!inherits(query, "Seurat")) stop("Query is not a Seurat object", call. = FALSE)
if (!(opt$query_assay %in% Seurat::Assays(query))) {
  stop("Query assay missing: ", opt$query_assay, call. = FALSE)
}
query <- join_layers_if_needed(query, opt$query_assay)
log_msg("Query cells: ", ncol(query))

reference <- ensure_log_normalized(reference, opt$reference_assay, nfeatures)
query <- ensure_log_normalized(query, opt$query_assay, nfeatures)

common_features <- intersect(rownames(reference[[opt$reference_assay]]), rownames(query[[opt$query_assay]]))
variable_features <- intersect(Seurat::VariableFeatures(reference, assay = opt$reference_assay), common_features)
if (length(variable_features) < 500) {
  log_msg("Only ", length(variable_features), " shared reference variable features; using common features capped at nfeatures")
  variable_features <- common_features[seq_len(min(length(common_features), nfeatures))]
}
if (length(variable_features) < 50) stop("Too few shared features for transfer: ", length(variable_features), call. = FALSE)
log_msg("Shared transfer features: ", length(variable_features))

npcs_use <- min(npcs, dims_n, length(variable_features) - 1L, ncol(reference) - 1L)
if (npcs_use < 2L) stop("Too few PCs available for transfer", call. = FALSE)
dims_use <- seq_len(npcs_use)

log_msg("Scaling reference and running PCA with ", npcs_use, " components")
reference <- Seurat::ScaleData(reference, assay = opt$reference_assay, features = variable_features, verbose = FALSE)
reference <- Seurat::RunPCA(
  reference,
  assay = opt$reference_assay,
  features = variable_features,
  npcs = npcs_use,
  verbose = FALSE
)

log_msg("Finding transfer anchors")
anchors <- Seurat::FindTransferAnchors(
  reference = reference,
  query = query,
  normalization.method = opt$normalization_method,
  reference.assay = opt$reference_assay,
  query.assay = opt$query_assay,
  features = variable_features,
  reference.reduction = "pca",
  reduction = "pcaproject",
  dims = dims_use,
  verbose = TRUE
)
n_anchors <- nrow(anchors@anchors)
if (n_anchors < 2L) stop("Produced too few transfer anchors: ", n_anchors, call. = FALSE)
n_unique_query_anchor_cells <- anchor_query_cell_count(anchors)
transfer_k <- choose_transfer_k_weight(n_anchors, requested = 50L)
log_msg(
  "TransferData k.weight: ", transfer_k$used,
  " (requested=", transfer_k$requested,
  ", anchors=", n_anchors,
  ", unique_query_anchor_cells=", n_unique_query_anchor_cells,
  ", reason=", transfer_k$reason, ")"
)

log_msg("Running TransferData for Shi labels")
predictions <- Seurat::TransferData(
  anchorset = anchors,
  refdata = reference$shi_transfer_label,
  dims = dims_use,
  k.weight = transfer_k$used,
  verbose = TRUE
)
predictions <- as.data.frame(predictions, stringsAsFactors = FALSE)
predictions <- data.frame(cell_id = rownames(predictions), predictions, check.names = FALSE)
score_cols <- validate_transfer_predictions(predictions, "major-label")

log_msg("Running TransferData for Shi gestational week labels")
week_predictions <- Seurat::TransferData(
  anchorset = anchors,
  refdata = reference$shi_transfer_week_label,
  dims = dims_use,
  k.weight = transfer_k$used,
  verbose = TRUE
)
week_predictions <- as.data.frame(week_predictions, stringsAsFactors = FALSE)
week_predictions <- data.frame(cell_id = rownames(week_predictions), week_predictions, check.names = FALSE)
week_score_cols <- validate_transfer_predictions(week_predictions, "week")

prediction_path <- file.path(opt$outdir, paste0(opt$output_prefix, "_shi_seurat_full_predictions.tsv.gz"))
score_path <- file.path(opt$outdir, paste0(opt$output_prefix, "_shi_seurat_full_prediction_scores.tsv.gz"))
week_prediction_path <- file.path(opt$outdir, paste0(opt$output_prefix, "_shi_seurat_full_week_predictions.tsv.gz"))
week_score_path <- file.path(opt$outdir, paste0(opt$output_prefix, "_shi_seurat_full_week_prediction_scores.tsv.gz"))
diagnostics_path <- file.path(opt$outdir, paste0(opt$output_prefix, "_shi_seurat_full_transfer_diagnostics.tsv"))

scores <- rename_score_columns(predictions, "shi_seurat_full_prediction_score_", score_cols)
week_scores <- rename_score_columns(week_predictions, "shi_seurat_full_week_prediction_score_", week_score_cols)
validate_canonical_scores(
  scores,
  predictions$prediction.score.max,
  1 - as.numeric(predictions$prediction.score.max),
  "major-label"
)
validate_canonical_scores(
  week_scores,
  week_predictions$prediction.score.max,
  1 - as.numeric(week_predictions$prediction.score.max),
  "week"
)

log_msg("Writing predictions: ", prediction_path)
write_tsv_gz(predictions, prediction_path)
log_msg("Writing score matrix: ", score_path)
write_tsv_gz(scores, score_path)
log_msg("Writing week predictions: ", week_prediction_path)
write_tsv_gz(week_predictions, week_prediction_path)
log_msg("Writing week score matrix: ", week_score_path)
write_tsv_gz(week_scores, week_score_path)

label_counts <- as.data.frame(table(reference$shi_transfer_label), stringsAsFactors = FALSE)
colnames(label_counts) <- c("shi_label", "n_reference_cells")
label_counts_path <- file.path(opt$outdir, "shi_reference_labels_used_by_seurat.tsv")
write_tsv(label_counts, label_counts_path)

week_counts <- as.data.frame(table(reference$shi_transfer_week_label), stringsAsFactors = FALSE)
colnames(week_counts) <- c("shi_week_label", "n_reference_cells")
week_counts_path <- file.path(opt$outdir, "shi_reference_weeks_used_by_seurat.tsv")
write_tsv(week_counts, week_counts_path)

diagnostics <- data.frame(
  metric = c(
    "reference_cells_labelled",
    "reference_cells_with_week",
    "query_cells",
    "shared_features_used",
    "pca_components_used",
    "anchors",
    "unique_query_anchor_cells",
    "k_weight_requested",
    "k_weight_used",
    "k_weight_reason",
    "median_prediction_score_max",
    "fraction_prediction_score_max_ge_0_75_qc_only",
    "reference_assay",
    "query_assay",
    "normalization_method",
    "transfer_name"
  ),
  value = c(
    as.character(ncol(reference)),
    as.character(sum(!is.na(reference$shi_transfer_week_label) & nzchar(as.character(reference$shi_transfer_week_label)))),
    as.character(ncol(query)),
    as.character(length(variable_features)),
    as.character(length(dims_use)),
    as.character(n_anchors),
    as.character(n_unique_query_anchor_cells),
    as.character(transfer_k$requested),
    as.character(transfer_k$used),
    transfer_k$reason,
    as.character(stats::median(as.numeric(predictions$prediction.score.max), na.rm = TRUE)),
    as.character(mean(as.numeric(predictions$prediction.score.max) >= 0.75, na.rm = TRUE)),
    opt$reference_assay,
    opt$query_assay,
    opt$normalization_method,
    opt$transfer_name
  ),
  stringsAsFactors = FALSE
)
write_tsv(diagnostics, diagnostics_path)

log_msg("Finished Seurat Shi label transfer export")
