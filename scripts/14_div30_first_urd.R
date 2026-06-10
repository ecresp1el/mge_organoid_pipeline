#!/usr/bin/env Rscript

# First-pass URD pseudotime for DIV30 cells.
#
# Inputs are plain Matrix Market files exported by
# python_notebooks/scripts/export_div30_first_urd_inputs.py. This keeps the URD
# run independent of the Seurat object and records the exact cells used for the
# pilot lineage reconstruction.

log_msg <- function(...) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste0(..., collapse = " ")))
  flush.console()
}

parse_args <- function(args) {
  out <- list(
    `project-root` = Sys.getenv("PROJECT_ROOT"),
    `input-dir` = NULL,
    outdir = NULL,
    `run-label` = Sys.getenv("RUN_LABEL", "div30_first_urd_paper_radial_glia_v1"),
    `root-label` = Sys.getenv("ROOT_LABEL", "Radial glia"),
    knn = Sys.getenv("URD_KNN", "100"),
    sigma = Sys.getenv("URD_SIGMA", "local"),
    `n-floods` = Sys.getenv("URD_N_FLOODS", "20"),
    `num-variable-genes` = Sys.getenv("URD_NUM_VARIABLE_GENES", "3000"),
    `min-genes` = Sys.getenv("URD_MIN_GENES", "500"),
    `min-cells` = Sys.getenv("URD_MIN_CELLS", "3"),
    `min-counts` = Sys.getenv("URD_MIN_COUNTS", "10"),
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

print_usage <- function() {
  cat(paste(
    "Usage:",
    "  Rscript scripts/14_div30_first_urd.R --project-root <PROJECT_ROOT> [--input-dir <dir>] [--outdir <dir>]",
    "",
    "Inputs in --input-dir:",
    "  div30_first_urd_counts.mtx",
    "  div30_first_urd_features.tsv",
    "  div30_first_urd_cell_metadata.tsv",
    "",
    "Outputs:",
    "  div30_first_urd_pseudotime.tsv",
    "  div30_first_urd_summary.tsv",
    "  div30_first_urd_marker_correlations.tsv",
    "  div30_first_urd_object.rds",
    sep = "\n"
  ))
}

opt <- parse_args(commandArgs(trailingOnly = TRUE))
if (isTRUE(opt$help)) {
  print_usage()
  quit(save = "no", status = 0)
}

required <- c("Matrix", "URD", "ggplot2")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) {
  stop("Missing required R package(s): ", paste(missing, collapse = ", "), call. = FALSE)
}

suppressPackageStartupMessages({
  library(Matrix)
  library(URD)
  library(ggplot2)
})

trim_trailing_slash <- function(x) sub("/+$", "", x)
if (!nzchar(opt$`project-root`)) stop("PROJECT_ROOT or --project-root is required", call. = FALSE)
project_root <- trim_trailing_slash(opt$`project-root`)
run_label <- opt$`run-label`
input_dir <- if (is.null(opt$`input-dir`) || !nzchar(opt$`input-dir`)) {
  file.path(project_root, "results/div30_first_urd", run_label, "inputs")
} else {
  opt$`input-dir`
}
outdir <- if (is.null(opt$outdir) || !nzchar(opt$outdir)) {
  file.path(project_root, "results/div30_first_urd", run_label)
} else {
  opt$outdir
}
plot_dir <- file.path(outdir, "plots")
table_dir <- file.path(outdir, "tables")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

knn <- as.integer(opt$knn)
n_floods <- as.integer(opt$`n-floods`)
num_variable_genes <- as.integer(opt$`num-variable-genes`)
min_genes <- as.integer(opt$`min-genes`)
min_cells <- as.integer(opt$`min-cells`)
min_counts <- as.integer(opt$`min-counts`)
sigma_text <- opt$sigma
sigma <- suppressWarnings(as.numeric(sigma_text))
if (is.na(sigma)) {
  sigma <- if (toupper(sigma_text) == "NULL") NULL else sigma_text
}

counts_path <- file.path(input_dir, "div30_first_urd_counts.mtx")
features_path <- file.path(input_dir, "div30_first_urd_features.tsv")
metadata_path <- file.path(input_dir, "div30_first_urd_cell_metadata.tsv")
needed <- c(counts_path, features_path, metadata_path)
missing_files <- needed[!file.exists(needed)]
if (length(missing_files) > 0) {
  stop("Missing URD input file(s): ", paste(missing_files, collapse = ", "), call. = FALSE)
}

log_msg("Reading counts: ", counts_path)
counts <- Matrix::readMM(counts_path)
if (!inherits(counts, "dgCMatrix")) counts <- as(counts, "dgCMatrix")

features <- read.delim(features_path, stringsAsFactors = FALSE, check.names = FALSE)
metadata <- read.delim(metadata_path, stringsAsFactors = FALSE, check.names = FALSE)
if (!("feature_id" %in% colnames(features))) stop("Feature table needs feature_id column", call. = FALSE)
if (!("cell_id" %in% colnames(metadata))) stop("Metadata table needs cell_id column", call. = FALSE)
if (nrow(counts) != nrow(features)) stop("Counts rows do not match features", call. = FALSE)
if (ncol(counts) != nrow(metadata)) stop("Counts columns do not match metadata rows", call. = FALSE)

rownames(counts) <- make.unique(as.character(features$feature_id))
colnames(counts) <- as.character(metadata$cell_id)
rownames(metadata) <- metadata$cell_id
metadata$urd_root_candidate <- as.logical(metadata$urd_root_candidate)

log_msg("Creating URD object")
urd <- createURD(
  count.data = counts,
  meta = metadata,
  min.cells = min_cells,
  min.genes = min_genes,
  min.counts = min_counts,
  verbose = TRUE
)
rm(counts)
invisible(gc())

cells <- rownames(urd@meta)
root_cells <- cells[as.logical(urd@meta[cells, "urd_root_candidate"])]
if (length(root_cells) == 0) {
  stop("No root candidate cells remained after URD filtering", call. = FALSE)
}
log_msg("Root cells retained after filtering: ", length(root_cells))

log_msg("Selecting top variable genes from URD logupx.data")
log_data <- urd@logupx.data
gene_means <- Matrix::rowMeans(log_data)
log_data_sq <- log_data
log_data_sq@x <- log_data_sq@x^2
gene_vars <- Matrix::rowMeans(log_data_sq) - gene_means^2
gene_vars <- gene_vars[is.finite(gene_vars)]
gene_vars <- sort(gene_vars, decreasing = TRUE)
urd@var.genes <- names(gene_vars)[seq_len(min(num_variable_genes, length(gene_vars)))]
log_msg("Variable genes stored in URD object: ", length(urd@var.genes))

set.seed(7)
log_msg("Running calcPCA")
urd <- calcPCA(urd, mp.factor = 2)
log_msg("Running calcDM with knn=", knn, " sigma=", sigma_text)
urd <- calcDM(urd, knn = knn, sigma = sigma)

pseudotime_name <- "paper_radial_glia_root"
log_msg("Running floodPseudotime n=", n_floods)
floods <- floodPseudotime(
  urd,
  root.cells = root_cells,
  n = n_floods,
  minimum.cells.flooded = 2,
  verbose = TRUE
)
log_msg("Processing flood pseudotime")
urd <- floodPseudotimeProcess(
  urd,
  floods,
  floods.name = pseudotime_name,
  max.frac.NA = 0.4,
  pseudotime.fun = mean,
  stability.div = min(10, n_floods)
)

extract_pseudotime <- function(object, name) {
  pt <- object@pseudotime
  if (is.data.frame(pt) || is.matrix(pt)) return(as.numeric(pt[rownames(object@meta), name]))
  if (is.list(pt) && name %in% names(pt)) return(as.numeric(pt[[name]][rownames(object@meta)]))
  stop("Could not extract pseudotime named ", name, call. = FALSE)
}

pt_values <- extract_pseudotime(urd, pseudotime_name)
out_meta <- urd@meta[rownames(urd@meta), , drop = FALSE]
pt_df <- data.frame(
  cell_id = rownames(out_meta),
  urd_pseudotime = pt_values,
  urd_pseudotime_name = pseudotime_name,
  urd_root_label = opt$`root-label`,
  stringsAsFactors = FALSE
)
pt_df <- cbind(pt_df, out_meta)

marker_genes <- c("DLX2", "ASCL1", "DCX", "STMN2", "TUBB3", "MAP2", "RBFOX3", "SYT1", "SNAP25")
find_gene <- function(gene, available) {
  hit <- which(toupper(available) == toupper(gene))
  if (length(hit) == 0) return(NA_integer_)
  hit[[1]]
}
available_genes <- rownames(urd@logupx.data)
for (gene in marker_genes) {
  idx <- find_gene(gene, available_genes)
  col <- paste0("logupx_", gene)
  if (is.na(idx)) {
    pt_df[[col]] <- NA_real_
  } else {
    pt_df[[col]] <- as.numeric(urd@logupx.data[idx, pt_df$cell_id])
  }
}
neuron_cols <- paste0("logupx_", c("DCX", "STMN2", "TUBB3", "MAP2", "RBFOX3", "SYT1", "SNAP25"))
neuron_cols <- neuron_cols[neuron_cols %in% colnames(pt_df)]
pt_df$neuronal_maturation_score <- rowMeans(pt_df[, neuron_cols, drop = FALSE], na.rm = TRUE)

pseudotime_path <- file.path(table_dir, "div30_first_urd_pseudotime.tsv")
write.table(pt_df, pseudotime_path, sep = "\t", row.names = FALSE, quote = FALSE)

safe_cor <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3) return(NA_real_)
  suppressWarnings(cor(x[ok], y[ok], method = "spearman"))
}
score_cols <- intersect(
  c(
    "jia_score_RGC1",
    "jia_score_RGC2",
    "jia_score_IPC",
    "logupx_DLX2",
    "logupx_ASCL1",
    "logupx_DCX",
    "neuronal_maturation_score",
    "shi_seurat_full_prediction_score_MGE",
    "shi_seurat_full_prediction_score_progenitor",
    "shi_seurat_full_expected_shi_week_numeric"
  ),
  colnames(pt_df)
)
cor_df <- data.frame(
  feature = score_cols,
  spearman_with_urd_pseudotime = vapply(score_cols, function(col) safe_cor(pt_df$urd_pseudotime, as.numeric(pt_df[[col]])), numeric(1)),
  stringsAsFactors = FALSE
)
cor_path <- file.path(table_dir, "div30_first_urd_marker_correlations.tsv")
write.table(cor_df, cor_path, sep = "\t", row.names = FALSE, quote = FALSE)

cluster_summary <- aggregate(
  urd_pseudotime ~ paper_cluster_annotation,
  data = pt_df,
  FUN = function(x) median(x, na.rm = TRUE)
)
colnames(cluster_summary)[2] <- "median_urd_pseudotime"
cluster_summary$n_cells <- as.integer(table(pt_df$paper_cluster_annotation)[cluster_summary$paper_cluster_annotation])
cluster_path <- file.path(table_dir, "div30_first_urd_pseudotime_by_paper_cluster.tsv")
write.table(cluster_summary, cluster_path, sep = "\t", row.names = FALSE, quote = FALSE)

summary_df <- data.frame(
  key = c(
    "run_label",
    "input_dir",
    "root_label",
    "n_cells_urd",
    "n_root_cells_urd",
    "n_variable_genes",
    "knn",
    "sigma",
    "n_floods",
    "pseudotime_tsv",
    "correlation_tsv",
    "cluster_summary_tsv"
  ),
  value = c(
    run_label,
    input_dir,
    opt$`root-label`,
    nrow(pt_df),
    length(root_cells),
    length(urd@var.genes),
    knn,
    sigma_text,
    n_floods,
    pseudotime_path,
    cor_path,
    cluster_path
  )
)
summary_path <- file.path(table_dir, "div30_first_urd_summary.tsv")
write.table(summary_df, summary_path, sep = "\t", row.names = FALSE, quote = FALSE)

plot_base <- pt_df[is.finite(pt_df$UMAP_1) & is.finite(pt_df$UMAP_2), , drop = FALSE]
if (nrow(plot_base) > 0) {
  png(file.path(plot_dir, "div30_first_urd_umap_pseudotime.png"), width = 1800, height = 1500, res = 220)
  print(
    ggplot(plot_base, aes(x = UMAP_1, y = UMAP_2, color = urd_pseudotime)) +
      geom_point(size = 0.25, alpha = 0.8) +
      coord_equal() +
      scale_color_viridis_c(na.value = "grey85") +
      theme_void(base_size = 10) +
      labs(color = "URD pseudotime", title = "DIV30 first URD: paper Radial glia root")
  )
  dev.off()

  png(file.path(plot_dir, "div30_first_urd_pseudotime_by_paper_cluster.png"), width = 1800, height = 1200, res = 220)
  print(
    ggplot(pt_df, aes(x = paper_cluster_annotation, y = urd_pseudotime, fill = paper_cluster_annotation)) +
      geom_boxplot(outlier.size = 0.2) +
      theme_bw(base_size = 10) +
      theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "none") +
      labs(x = NULL, y = "URD pseudotime", title = "DIV30 first URD pseudotime by paper/manual annotation")
  )
  dev.off()
}

rds_path <- file.path(outdir, "div30_first_urd_object.rds")
saveRDS(urd, rds_path)
log_msg("Wrote pseudotime: ", pseudotime_path)
log_msg("Wrote summary: ", summary_path)
log_msg("Saved URD object: ", rds_path)
