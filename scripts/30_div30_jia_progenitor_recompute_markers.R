#!/usr/bin/env Rscript

# Recompute DIV30 progenitor reclustering markers from a saved progenitor-only
# Seurat object. This is a post hoc repair path for Seurat/future environments
# where FindAllMarkers was blocked by future.globals.maxSize.

suppressPackageStartupMessages({
  library(Seurat)
})

log_msg <- function(...) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste0(..., collapse = "")))
  flush.console()
}

parse_args <- function(args) {
  out <- list(
    `seurat-rds` = NULL,
    outdir = NULL,
    assay = "RNA",
    `cluster-col` = "div30_progenitor_cluster",
    `min-pct` = "0.1",
    `logfc-threshold` = "0.25",
    `top-n-markers` = "50",
    `future-max-gb` = "16",
    help = FALSE
  )
  i <- 1L
  while (i <= length(args)) {
    a <- args[[i]]
    if (a %in% c("--help", "-h")) {
      out$help <- TRUE
      i <- i + 1L
      next
    }
    if (!startsWith(a, "--")) stop("Unknown argument: ", a, call. = FALSE)
    key <- substring(a, 3L)
    if (!(key %in% names(out))) stop("Unknown argument: ", a, call. = FALSE)
    if (i == length(args)) stop("Missing value for argument: ", a, call. = FALSE)
    out[[key]] <- args[[i + 1L]]
    i <- i + 2L
  }
  out
}

as_num <- function(x, name) {
  value <- suppressWarnings(as.numeric(x))
  if (is.na(value)) stop(name, " must be numeric; got ", x, call. = FALSE)
  value
}

as_int <- function(x, name) {
  value <- suppressWarnings(as.integer(x))
  if (is.na(value)) stop(name, " must be an integer; got ", x, call. = FALSE)
  value
}

write_tsv <- function(x, path) {
  con <- if (grepl("\\.gz$", path)) gzfile(path, open = "wt") else file(path, open = "wt")
  on.exit(close(con), add = TRUE)
  write.table(x, con, sep = "\t", row.names = FALSE, quote = FALSE, na = "")
}

top_markers <- function(markers, top_n) {
  if (nrow(markers) == 0) return(markers)
  pieces <- split(markers, markers$cluster)
  out <- do.call(rbind, lapply(pieces, function(df) head(df, top_n)))
  rownames(out) <- NULL
  out
}

opt <- parse_args(commandArgs(trailingOnly = TRUE))
if (isTRUE(opt$help)) {
  cat("Usage: Rscript scripts/30_div30_jia_progenitor_recompute_markers.R --seurat-rds <rds> --outdir <run_root>\n")
  quit(save = "no", status = 0)
}
if (is.null(opt$`seurat-rds`) || !nzchar(opt$`seurat-rds`)) stop("--seurat-rds is required", call. = FALSE)
if (is.null(opt$outdir) || !nzchar(opt$outdir)) stop("--outdir is required", call. = FALSE)

cfg <- list(
  seurat_rds = opt$`seurat-rds`,
  outdir = opt$outdir,
  table_dir = file.path(opt$outdir, "tables"),
  assay = opt$assay,
  cluster_col = opt$`cluster-col`,
  min_pct = as_num(opt$`min-pct`, "min-pct"),
  logfc_threshold = as_num(opt$`logfc-threshold`, "logfc-threshold"),
  top_n_markers = as_int(opt$`top-n-markers`, "top-n-markers"),
  future_max_gb = as_num(opt$`future-max-gb`, "future-max-gb")
)
dir.create(cfg$table_dir, recursive = TRUE, showWarnings = FALSE)

log_msg("Reading Seurat object: ", cfg$seurat_rds)
obj <- readRDS(cfg$seurat_rds)
if (!(cfg$cluster_col %in% colnames(obj@meta.data))) stop("Cluster column not found: ", cfg$cluster_col, call. = FALSE)
DefaultAssay(obj) <- cfg$assay
Idents(obj) <- obj@meta.data[[cfg$cluster_col]]

options(future.globals.maxSize = cfg$future_max_gb * 1024^3)
if (requireNamespace("future", quietly = TRUE)) {
  future::plan("sequential")
}

log_msg("Running FindAllMarkers with future.globals.maxSize=", cfg$future_max_gb, " GiB")
markers <- FindAllMarkers(
  obj,
  assay = cfg$assay,
  only.pos = TRUE,
  min.pct = cfg$min_pct,
  logfc.threshold = cfg$logfc_threshold,
  verbose = TRUE
)
if (nrow(markers) > 0 && "avg_logFC" %in% colnames(markers) && !("avg_log2FC" %in% colnames(markers))) {
  markers$avg_log2FC <- markers$avg_logFC
}
if (nrow(markers) > 0 && all(c("cluster", "avg_log2FC", "p_val_adj") %in% colnames(markers))) {
  markers <- markers[order(markers$cluster, -markers$avg_log2FC, markers$p_val_adj), , drop = FALSE]
}
rownames(markers) <- NULL

write_tsv(markers, file.path(cfg$table_dir, "div30_jia_progenitor_phase1_all_markers.tsv.gz"))
marker_cols <- intersect(c("cluster", "gene", "avg_log2FC", "pct.1", "pct.2", "p_val_adj"), colnames(markers))
write_tsv(top_markers(markers, cfg$top_n_markers)[, marker_cols, drop = FALSE], file.path(cfg$table_dir, "div30_jia_progenitor_phase1_top_markers.tsv"))
write_tsv(
  data.frame(
    completed_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    seurat_rds = cfg$seurat_rds,
    n_markers = nrow(markers),
    n_clusters = length(unique(as.character(obj@meta.data[[cfg$cluster_col]]))),
    future_max_gb = cfg$future_max_gb,
    stringsAsFactors = FALSE
  ),
  file.path(cfg$table_dir, "div30_jia_progenitor_phase1_marker_recompute_complete.tsv")
)
log_msg("Marker recompute complete. n_markers=", nrow(markers))
