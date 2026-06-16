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
    bridge_dir = NULL,
    outdir = NULL,
    label_column = "candidate_jia_group",
    reduction = "rpca",
    dims = 20L,
    nfeatures = 3000L,
    npcs = 50L,
    seed = 0L,
    k_weight = 50L,
    min_shared_features = 500L
  )
  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, "--")) stop("Unknown argument: ", key, call. = FALSE)
    name <- gsub("-", "_", substring(key, 3L))
    if (!(name %in% names(out))) stop("Unknown argument: ", key, call. = FALSE)
    if (i == length(args)) stop("Missing value for argument: ", key, call. = FALSE)
    out[[name]] <- args[[i + 1L]]
    i <- i + 2L
  }
  for (name in c("dims", "nfeatures", "npcs", "seed", "k_weight", "min_shared_features")) {
    out[[name]] <- as.integer(out[[name]])
  }
  out
}
opt <- parse_args(commandArgs(trailingOnly = TRUE))

required <- c("bridge_dir", "outdir")
for (name in required) {
  if (is.null(opt[[name]]) || !nzchar(opt[[name]])) stop("--", gsub("_", "-", name), " is required", call. = FALSE)
}

read_tsv <- function(path) {
  utils::read.delim(path, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
}

write_tsv <- function(x, path) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  utils::write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
}

write_tsv_gz <- function(x, path) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  con <- gzfile(path, open = "wt")
  on.exit(close(con), add = TRUE)
  utils::write.table(x, con, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
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

anchor_query_cell_count <- function(anchors) {
  anchor_df <- as.data.frame(anchors@anchors)
  if ("cell2" %in% colnames(anchor_df)) length(unique(anchor_df$cell2)) else NA_integer_
}

choose_transfer_k_weight <- function(n_anchors, requested = 50L) {
  if (n_anchors <= requested) {
    list(used = max(1L, n_anchors - 1L), reason = paste0("lowered_from_", requested, "_because_anchor_count_was_", n_anchors))
  } else {
    list(used = requested, reason = "requested_k_weight_used")
  }
}

load_bridge_object <- function(prefix, bridge_dir) {
  mat <- Matrix::readMM(file.path(bridge_dir, paste0(prefix, "_counts.mtx")))
  genes <- read_tsv(file.path(bridge_dir, paste0(prefix, "_genes.tsv")))$gene
  barcodes <- read_tsv(file.path(bridge_dir, paste0(prefix, "_barcodes.tsv")))$cell_id
  meta <- read_tsv(file.path(bridge_dir, paste0(prefix, "_metadata.tsv.gz")))
  if (nrow(mat) != length(genes)) stop(prefix, ": matrix rows do not match genes", call. = FALSE)
  if (ncol(mat) != length(barcodes)) stop(prefix, ": matrix columns do not match barcodes", call. = FALSE)
  rownames(mat) <- make.unique(as.character(genes))
  colnames(mat) <- as.character(barcodes)
  rownames(meta) <- as.character(meta$seurat_cell_id)
  obj <- CreateSeuratObject(counts = mat, meta.data = meta)
  obj
}

set.seed(opt$seed)
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(opt$outdir, "seurat"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(opt$outdir, "tables"), showWarnings = FALSE, recursive = TRUE)

log_msg("Seurat version: ", as.character(packageVersion("Seurat")))
log_msg("SeuratObject version: ", as.character(packageVersion("SeuratObject")))
log_msg("Bridge dir: ", opt$bridge_dir)
log_msg("Label column: ", opt$label_column)
log_msg("Reduction: ", opt$reduction)
log_msg("Dims: 1:", opt$dims)
log_msg("nfeatures: ", opt$nfeatures)
log_msg("npcs: ", opt$npcs)
log_msg("seed: ", opt$seed)

reference <- load_bridge_object("reference", opt$bridge_dir)
query <- load_bridge_object("query", opt$bridge_dir)

if (!(opt$label_column %in% colnames(reference@meta.data))) {
  stop("Reference metadata is missing label column: ", opt$label_column, call. = FALSE)
}
reference[[opt$label_column]][is.na(reference[[opt$label_column]][, 1]) | reference[[opt$label_column]][, 1] == ""] <- "unlabeled_or_na"

shared_features <- intersect(rownames(reference), rownames(query))
if (length(shared_features) < opt$min_shared_features) {
  stop("Only ", length(shared_features), " shared features; minimum is ", opt$min_shared_features, call. = FALSE)
}

reference <- NormalizeData(reference, normalization.method = "LogNormalize", verbose = FALSE)
query <- NormalizeData(query, normalization.method = "LogNormalize", verbose = FALSE)

reference <- FindVariableFeatures(reference, selection.method = "vst", nfeatures = opt$nfeatures, verbose = FALSE)
query <- FindVariableFeatures(query, selection.method = "vst", nfeatures = opt$nfeatures, verbose = FALSE)
features <- unique(c(VariableFeatures(reference), VariableFeatures(query)))
features <- intersect(features, shared_features)
features <- features[seq_len(min(length(features), opt$nfeatures))]
if (length(features) < opt$min_shared_features) {
  stop("Only ", length(features), " selected variable shared features; minimum is ", opt$min_shared_features, call. = FALSE)
}
write_tsv(data.frame(feature = features), file.path(opt$outdir, "seurat", "selected_transfer_features.tsv"))

VariableFeatures(reference) <- features
VariableFeatures(query) <- features
reference <- ScaleData(reference, features = features, verbose = FALSE)
query <- ScaleData(query, features = features, verbose = FALSE)
reference <- RunPCA(reference, features = features, npcs = opt$npcs, verbose = FALSE)
query <- RunPCA(query, features = features, npcs = opt$npcs, verbose = FALSE)

dims_use <- seq_len(min(opt$dims, opt$npcs))
anchor_args <- list(
  reference = reference,
  query = query,
  normalization.method = "LogNormalize",
  reference.assay = "RNA",
  query.assay = "RNA",
  reduction = opt$reduction,
  features = features,
  dims = dims_use,
  npcs = opt$npcs,
  verbose = FALSE
)
if (opt$reduction == "rpca") {
  anchor_args$reference.reduction <- "pca"
}
anchors <- do.call(FindTransferAnchors, anchor_args)

anchor_count <- nrow(as.data.frame(anchors@anchors))
k_choice <- choose_transfer_k_weight(anchor_count, requested = opt$k_weight)
predictions <- as.data.frame(TransferData(
  anchorset = anchors,
  refdata = reference@meta.data[[opt$label_column]],
  dims = dims_use,
  k.weight = k_choice$used,
  verbose = FALSE
))
predictions$cell_id <- colnames(query)
predictions <- predictions[, c("cell_id", setdiff(colnames(predictions), "cell_id")), drop = FALSE]

score_cols <- prediction_score_cols(predictions)
scores <- predictions[, c("cell_id", score_cols), drop = FALSE]
colnames(scores) <- c("cell_id", paste0("score_", sanitize_token(sub("^prediction\\.score\\.", "", score_cols))))

diag <- data.frame(
  label_column = opt$label_column,
  reduction = opt$reduction,
  dims = paste0("1:", max(dims_use)),
  nfeatures_requested = opt$nfeatures,
  nfeatures_used = length(features),
  npcs = opt$npcs,
  seed = opt$seed,
  n_reference_cells = ncol(reference),
  n_query_cells = ncol(query),
  n_anchors = anchor_count,
  n_unique_query_anchor_cells = anchor_query_cell_count(anchors),
  transfer_k_weight_requested = opt$k_weight,
  transfer_k_weight_used = k_choice$used,
  transfer_k_weight_reason = k_choice$reason,
  n_prediction_score_columns = length(score_cols)
)

prefix <- paste("siletti_div90", opt$label_column, opt$reduction, paste0("dims", opt$dims), sep = "__")
prefix <- sanitize_token(prefix)
write_tsv_gz(predictions, file.path(opt$outdir, "seurat", paste0(prefix, "_predictions.tsv.gz")))
write_tsv_gz(scores, file.path(opt$outdir, "seurat", paste0(prefix, "_prediction_scores.tsv.gz")))
write_tsv(diag, file.path(opt$outdir, "seurat", paste0(prefix, "_transfer_diagnostics.tsv")))

query_meta <- query@meta.data
query_meta$cell_id <- colnames(query)
obs <- cbind(query_meta, predictions[, c("predicted.id", "prediction.score.max"), drop = FALSE])
write_tsv_gz(obs, file.path(opt$outdir, "tables", paste0(prefix, "_query_obs_with_predictions.tsv.gz")))

cluster_col <- if ("cluster_id" %in% colnames(obs)) "cluster_id" else if ("cluster_id_manifest" %in% colnames(obs)) "cluster_id_manifest" else NA_character_
if (!is.na(cluster_col)) {
  summary <- as.data.frame(table(obs[[cluster_col]], obs$predicted.id), stringsAsFactors = FALSE)
  colnames(summary) <- c("cluster_id", "predicted_label", "n_cells")
  summary <- summary[summary$n_cells > 0, , drop = FALSE]
  write_tsv(summary, file.path(opt$outdir, "tables", paste0(prefix, "_cluster_label_counts.tsv")))
}

log_msg("Finished transfer: ", prefix)
print(diag)
