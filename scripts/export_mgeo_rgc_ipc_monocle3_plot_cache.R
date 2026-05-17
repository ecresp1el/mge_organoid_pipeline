#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
options(expressions = 500000)

suppressPackageStartupMessages({
  library(Matrix)
  library(monocle3)
})

PROJECT_ROOT <- Sys.getenv(
  "PROJECT_ROOT",
  unset = "/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder"
)

RUN_ROOT <- file.path(PROJECT_ROOT, "results/mgeo_rgc_ipc_monocle3")
INTERACTIVE_DIR <- file.path(RUN_ROOT, "interactive")
CACHE_DIR <- file.path(INTERACTIVE_DIR, "plot_cache")
dir.create(CACHE_DIR, recursive = TRUE, showWarnings = FALSE)

CDS_PREPROCESS_RDS <- file.path(INTERACTIVE_DIR, "cds_after_preprocess.rds")

marker_genes <- c(
  "NES",
  "VIM",
  "HES1",
  "ASCL1",
  "DLX2",
  "DCX",
  "STMN2",
  "FABP7",
  "HES5",
  "SOX2",
  "DLX1",
  "NKX2-1"
)

log_step <- function(...) {
  message(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste0(..., collapse = ""))
  flush.console()
}

partition_path <- function(prefix, partition_label, ext = ".csv") {
  file.path(INTERACTIVE_DIR, paste0(prefix, "_partition_", partition_label, ext))
}

cache_path <- function(prefix, partition_label = NULL, ext = ".csv") {
  suffix <- if (is.null(partition_label)) "" else paste0("_partition_", partition_label)
  file.path(CACHE_DIR, paste0(prefix, suffix, ext))
}

write_plot_ready_cache <- function(partition_label) {
  umap_path <- partition_path("mgeo_rgc_ipc_monocle3_umap_clusters", partition_label)
  pseudotime_path <- partition_path("mgeo_rgc_ipc_monocle3_pseudotime_interactive", partition_label)

  stopifnot(file.exists(umap_path))
  stopifnot(file.exists(pseudotime_path))

  log_step("Reading UMAP CSV for partition_", partition_label)
  umap_df <- read.csv(umap_path, check.names = FALSE)
  log_step("Reading pseudotime CSV for partition_", partition_label)
  pseudotime_df <- read.csv(pseudotime_path, check.names = FALSE)

  pseudotime_df <- pseudotime_df[match(umap_df$cell_id, pseudotime_df$cell_id), , drop = FALSE]
  stopifnot(identical(umap_df$cell_id, pseudotime_df$cell_id))

  duplicated_columns <- intersect(
    setdiff(colnames(pseudotime_df), "cell_id"),
    setdiff(colnames(umap_df), "cell_id")
  )
  pseudotime_keep <- setdiff(colnames(pseudotime_df), c("cell_id", duplicated_columns))

  plot_ready_df <- cbind(
    umap_df,
    pseudotime_df[, pseudotime_keep, drop = FALSE]
  )

  out_path <- cache_path("mgeo_rgc_ipc_monocle3_plot_ready", partition_label)
  write.csv(plot_ready_df, out_path, row.names = FALSE)
  log_step("Wrote plot-ready cache: ", out_path)

  out_path
}

log_step("PROJECT_ROOT=", PROJECT_ROOT)
log_step("CACHE_DIR=", CACHE_DIR)

plot_ready_true <- write_plot_ready_cache("true")
plot_ready_false <- write_plot_ready_cache("false")

stopifnot(file.exists(CDS_PREPROCESS_RDS))
log_step("Reading preprocessed CDS only to export marker expression cache")
cds <- readRDS(CDS_PREPROCESS_RDS)

available_marker_genes <- intersect(marker_genes, rownames(cds))
missing_marker_genes <- setdiff(marker_genes, available_marker_genes)
if (length(missing_marker_genes) > 0) {
  warning("Missing marker genes: ", paste(missing_marker_genes, collapse = ", "))
}
stopifnot(length(available_marker_genes) > 0)

plot_ready_df <- read.csv(plot_ready_true, check.names = FALSE)
cell_ids <- plot_ready_df$cell_id
stopifnot(all(cell_ids %in% colnames(cds)))

log_step("Extracting marker expression for ", length(available_marker_genes), " genes x ", length(cell_ids), " cells")
expr_mat <- exprs(cds)[available_marker_genes, cell_ids, drop = FALSE]
if (inherits(expr_mat, "sparseMatrix")) {
  expr_mat <- as.matrix(expr_mat)
}

expr_df <- data.frame(cell_id = cell_ids, t(expr_mat), check.names = FALSE)
expression_path <- cache_path("mgeo_rgc_ipc_marker_expression_log_normalized")
write.csv(expr_df, expression_path, row.names = FALSE)
log_step("Wrote marker expression cache: ", expression_path)

manifest <- data.frame(
  key = c(
    "project_root",
    "cache_dir",
    "plot_ready_partition_true",
    "plot_ready_partition_false",
    "marker_expression_log_normalized",
    "marker_genes",
    "source_expression",
    "note"
  ),
  value = c(
    PROJECT_ROOT,
    CACHE_DIR,
    plot_ready_true,
    plot_ready_false,
    expression_path,
    paste(available_marker_genes, collapse = ","),
    CDS_PREPROCESS_RDS,
    "Marker values are the same log-normalized expression values used by Monocle; UMAP/pseudotime are cached flat files, so plotting does not rerun Monocle."
  )
)

manifest_path <- cache_path("mgeo_rgc_ipc_monocle3_plot_cache_manifest", ext = ".tsv")
write.table(manifest, manifest_path, sep = "\t", quote = FALSE, row.names = FALSE)
log_step("Wrote manifest: ", manifest_path)
log_step("Done")
