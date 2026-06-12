#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
  library(Matrix)
  library(ggplot2)
  library(patchwork)
})

get_env <- function(name, default = "") {
  value <- Sys.getenv(name, unset = "")
  if (nzchar(value)) value else default
}

timestamp <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")

log_msg <- function(...) {
  cat(sprintf("[%s] %s\n", timestamp(), paste0(..., collapse = " ")))
  flush.console()
}

write_tsv <- function(df, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.table(df, path, sep = "\t", row.names = FALSE, quote = FALSE, na = "NA")
}

save_plot <- function(plot, stem, width = 8, height = 6, dpi = 300) {
  ggsave(paste0(stem, ".png"), plot, width = width, height = height, dpi = dpi)
  ggsave(paste0(stem, ".pdf"), plot, width = width, height = height)
}

project_root <- get_env("PROJECT_ROOT", "/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder")
repo_root <- get_env("REPO_ROOT", "/home/elcrespo/Desktop/githubprojects/mge_organoid_pipeline")
run_label <- get_env("XIANG_TSNE_RUN_LABEL", "xiang_2017_tsne_reproduction_v1")
out_dir <- get_env(
  "XIANG_TSNE_OUT_DIR",
  file.path(project_root, "results", "xiang_2017_tsne_reproduction", run_label)
)
seed <- as.integer(get_env("XIANG_TSNE_SEED", "1007"))
max_variable_genes_env <- get_env("XIANG_MAX_VARIABLE_GENES", "0")
max_variable_genes <- as.integer(max_variable_genes_env)
if (is.na(max_variable_genes)) max_variable_genes <- 0L
regression_model <- get_env("XIANG_REGRESSION_MODEL", "linear")

table_dir <- file.path(out_dir, "tables")
plot_dir <- file.path(out_dir, "plots")
rds_dir <- file.path(out_dir, "rds")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(rds_dir, recursive = TRUE, showWarnings = FALSE)

suppl_dir <- file.path(project_root, "data/raw/xiang_2018_geo_files/suppl")
matrix_path <- file.path(suppl_dir, "GSE98201_matrix.mtx.gz")
genes_path <- file.path(suppl_dir, "GSE98201_genes.tsv.gz")
barcodes_path <- file.path(suppl_dir, "GSE98201_barcodes.tsv.gz")
series_path <- file.path(
  project_root,
  "data/raw/xiang_2018_geo_files/matrix/GSE97882-GPL20301_series_matrix.txt.gz"
)

stopifnot(file.exists(matrix_path), file.exists(genes_path), file.exists(barcodes_path), file.exists(series_path))

log_msg("Output dir: ", out_dir)
log_msg("Reading Xiang GSE98201 matrix and metadata")

counts <- readMM(matrix_path)
genes <- fread(cmd = paste("zcat", shQuote(genes_path)), header = FALSE, col.names = c("ensembl_id", "gene_symbol"))
barcodes <- fread(cmd = paste("zcat", shQuote(barcodes_path)), header = FALSE, col.names = "barcode")

if (nrow(counts) != nrow(genes)) {
  stop("Gene table rows do not match matrix rows.")
}
if (ncol(counts) != nrow(barcodes)) {
  stop("Barcode rows do not match matrix columns.")
}

feature_names <- ifelse(
  !is.na(genes$gene_symbol) & nzchar(genes$gene_symbol),
  genes$gene_symbol,
  genes$ensembl_id
)
feature_names <- make.unique(feature_names)
rownames(counts) <- feature_names
colnames(counts) <- barcodes$barcode

suffix_map <- data.frame(
  barcode_suffix = as.character(1:8),
  sample = c(
    "hMGEO_d30_rep1",
    "hCO_d30_rep1",
    "hMGEO_d72_rep1",
    "hCO_d72_rep1",
    "hMGEO_d30_rep2",
    "hCO_d30_rep2",
    "hMGEO_d79_rep2",
    "hCO_d79_rep2"
  ),
  condition = c("hMGEO", "hCO", "hMGEO", "hCO", "hMGEO", "hCO", "hMGEO", "hCO"),
  timepoint = c("d30", "d30", "d72", "d72", "d30", "d30", "d79", "d79"),
  replicate = c("rep1", "rep1", "rep1", "rep1", "rep2", "rep2", "rep2", "rep2"),
  geo_accession = c(
    "GSM2589129",
    "GSM2589130",
    "GSM2589131",
    "GSM2589132",
    "GSM2684867",
    "GSM2684868",
    "GSM2684869",
    "GSM2684870"
  ),
  stringsAsFactors = FALSE
)
write_tsv(suffix_map, file.path(table_dir, "Xiang_suffix_sample_map.tsv"))

barcode_suffix <- sub("^.*-", "", barcodes$barcode)
metadata <- data.frame(
  barcode = barcodes$barcode,
  barcode_suffix = barcode_suffix,
  stringsAsFactors = FALSE
)
metadata <- merge(metadata, suffix_map, by = "barcode_suffix", all.x = TRUE, sort = FALSE)
metadata <- metadata[match(barcodes$barcode, metadata$barcode), , drop = FALSE]
rownames(metadata) <- metadata$barcode

if (any(is.na(metadata$sample))) {
  missing_suffixes <- paste(sort(unique(metadata$barcode_suffix[is.na(metadata$sample)])), collapse = ",")
  stop("Some barcode suffixes are missing from suffix_map: ", missing_suffixes)
}

metadata$sample <- factor(metadata$sample, levels = suffix_map$sample)
metadata$condition <- factor(metadata$condition, levels = c("hMGEO", "hCO"))
metadata$timepoint <- factor(metadata$timepoint, levels = c("d30", "d72", "d79"))
metadata$replicate <- factor(metadata$replicate, levels = c("rep1", "rep2"))
metadata$geo_accession <- factor(metadata$geo_accession, levels = suffix_map$geo_accession)
metadata$barcode_suffix <- factor(metadata$barcode_suffix, levels = as.character(1:8))

metadata_table <- metadata[, c("barcode", "sample", "condition", "timepoint", "replicate", "geo_accession", "barcode_suffix")]
write_tsv(metadata_table, file.path(table_dir, "Xiang_metadata_reconstructed.tsv"))

sample_validation <- as.data.frame(table(metadata$barcode_suffix, metadata$sample), stringsAsFactors = FALSE)
colnames(sample_validation) <- c("barcode_suffix", "sample", "n_cells")
sample_validation <- sample_validation[sample_validation$n_cells > 0, , drop = FALSE]
sample_validation <- merge(sample_validation, suffix_map, by = c("barcode_suffix", "sample"), all.x = TRUE, sort = FALSE)
sample_validation <- sample_validation[
  order(as.integer(as.character(sample_validation$barcode_suffix))),
  c("barcode_suffix", "sample", "condition", "timepoint", "replicate", "geo_accession", "n_cells")
]
write_tsv(sample_validation, file.path(table_dir, "Xiang_barcode_suffix_sample_validation.tsv"))

log_msg("Barcode suffix validation:")
print(sample_validation)

log_msg("Creating Seurat object with paper-style filters")
obj <- CreateSeuratObject(
  counts = counts,
  project = "Xiang2017_GSE98201",
  meta.data = metadata[, c("sample", "condition", "timepoint", "replicate", "geo_accession", "barcode_suffix")],
  min.cells = 4,
  min.features = 201
)
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")

qc_by_sample <- do.call(
  rbind,
  lapply(levels(obj$sample), function(one_sample) {
    cells <- colnames(obj)[obj$sample == one_sample]
    sample_counts <- GetAssayData(obj, assay = "RNA", layer = "counts")[, cells, drop = FALSE]
    data.frame(
      sample = one_sample,
      n_cells = length(cells),
      mean_reads_per_cell = NA_real_,
      median_genes_per_cell = median(obj$nFeature_RNA[cells]),
      total_reads = NA_real_,
      total_genes_detected = sum(Matrix::rowSums(sample_counts > 0) > 0),
      median_umi_counts_per_cell = median(obj$nCount_RNA[cells]),
      stringsAsFactors = FALSE
    )
  })
)
write_tsv(qc_by_sample, file.path(table_dir, "Xiang_figure_s3a_qc_by_sample.tsv"))

log_msg("Normalizing per cell and log2 transforming")
counts_filtered <- GetAssayData(obj, assay = "RNA", layer = "counts")
n_umi <- Matrix::colSums(counts_filtered)
log2_norm <- t(t(counts_filtered) / n_umi * 10000)
log2_norm@x <- log2(log2_norm@x + 1)
obj <- SetAssayData(obj, assay = "RNA", layer = "data", new.data = log2_norm)

gene_means <- Matrix::rowMeans(log2_norm)
gene_sq_means <- Matrix::rowMeans(log2_norm ^ 2)
gene_vars <- pmax(gene_sq_means - gene_means ^ 2, 0)
gene_dispersion <- gene_vars / pmax(gene_means, .Machine$double.eps)
variable_df <- data.frame(
  feature = rownames(log2_norm),
  mean = gene_means,
  variance = gene_vars,
  dispersion = gene_dispersion,
  stringsAsFactors = FALSE
)
variable_df <- variable_df[variable_df$mean > 0 & variable_df$dispersion > 0, , drop = FALSE]
variable_df <- variable_df[order(variable_df$dispersion, decreasing = TRUE), , drop = FALSE]
if (max_variable_genes > 0L && nrow(variable_df) > max_variable_genes) {
  variable_df$selected_for_pca <- seq_len(nrow(variable_df)) <= max_variable_genes
} else {
  variable_df$selected_for_pca <- TRUE
}
write_tsv(variable_df, file.path(table_dir, "Xiang_variable_genes_dispersion.tsv"))
variable_features <- variable_df$feature[variable_df$selected_for_pca]
VariableFeatures(obj) <- variable_features
log_msg("Variable genes with dispersion > 0: ", nrow(variable_df))
log_msg("Variable genes selected for regression/PCA: ", length(variable_features))

log_msg("Regressing sample and nUMI, model=", regression_model)
obj <- ScaleData(
  obj,
  features = variable_features,
  vars.to.regress = c("sample", "nCount_RNA"),
  model.use = regression_model,
  verbose = TRUE
)

log_msg("Running PCA, npcs=20")
set.seed(seed)
obj <- RunPCA(obj, features = variable_features, npcs = 20, approx = TRUE, verbose = FALSE)
pca_stdev <- obj[["pca"]]@stdev
pca_table <- data.frame(
  pc = paste0("PC", seq_along(pca_stdev)),
  stdev = pca_stdev,
  variance = pca_stdev ^ 2,
  stringsAsFactors = FALSE
)
write_tsv(pca_table, file.path(table_dir, "Xiang_pca_singular_values_modern.tsv"))

log_msg("Running tSNE on PC1-PC5 only")
set.seed(seed)
obj <- RunTSNE(
  obj,
  reduction = "pca",
  dims = 1:5,
  reduction.name = "xiang_tsne_pc1_5",
  reduction.key = "XiangTSNE_",
  seed.use = seed,
  check_duplicates = FALSE
)
tsne <- Embeddings(obj, "xiang_tsne_pc1_5")
tsne_table <- data.frame(
  barcode = rownames(tsne),
  tSNE_1 = tsne[, 1],
  tSNE_2 = tsne[, 2],
  obj@meta.data[rownames(tsne), c("sample", "condition", "timepoint", "replicate", "geo_accession", "barcode_suffix")],
  stringsAsFactors = FALSE
)
write_tsv(tsne_table, file.path(table_dir, "Xiang_tSNE_coordinates.tsv"))

for (group_col in c("sample", "condition", "timepoint", "replicate", "barcode_suffix")) {
  p <- DimPlot(obj, reduction = "xiang_tsne_pc1_5", group.by = group_col, pt.size = 0.2, raster = TRUE) +
    ggtitle(paste("Xiang 2017 tSNE PC1-PC5 by", group_col)) +
    theme_classic()
  save_plot(p, file.path(plot_dir, paste0("Xiang_tSNE_by_", group_col)), width = 8, height = 6)
}

marker_genes <- c(
  "NKX2-1", "DLX1", "DLX2", "GAD1", "GAD2", "TAC1",
  "PAX6", "NEUROG2", "TBR1", "NEUROD2", "BCL11B",
  "VIM", "NES", "HES1", "STMN2", "GAP43", "DCX", "GFAP"
)
present_markers <- marker_genes[marker_genes %in% rownames(obj)]
missing_markers <- setdiff(marker_genes, present_markers)
write_tsv(
  data.frame(marker = marker_genes, present = marker_genes %in% rownames(obj), stringsAsFactors = FALSE),
  file.path(table_dir, "Xiang_marker_presence.tsv")
)

if (length(present_markers) > 0) {
  marker_plots <- FeaturePlot(
    obj,
    features = present_markers,
    reduction = "xiang_tsne_pc1_5",
    cols = c("grey90", "firebrick"),
    pt.size = 0.15,
    raster = TRUE,
    combine = FALSE
  )
  marker_plot <- wrap_plots(marker_plots, ncol = 3)
  save_plot(marker_plot, file.path(plot_dir, "Xiang_tSNE_marker_featureplots"), width = 12, height = 18)
}

if (all(c("GFAP", "TBR1") %in% rownames(obj))) {
  gfap <- as.numeric(GetAssayData(obj, assay = "RNA", layer = "data")["GFAP", colnames(obj)])
  tbr1 <- as.numeric(GetAssayData(obj, assay = "RNA", layer = "data")["TBR1", colnames(obj)])
  doublet_df <- data.frame(
    barcode = colnames(obj),
    GFAP = gfap,
    TBR1 = tbr1,
    category = ifelse(gfap > 0 & tbr1 > 0, "GFAP+TBR1+",
      ifelse(gfap > 0, "GFAP_only", ifelse(tbr1 > 0, "TBR1_only", "neither"))
    ),
    stringsAsFactors = FALSE
  )
  write_tsv(as.data.frame(table(doublet_df$category), stringsAsFactors = FALSE), file.path(table_dir, "Xiang_s3c_GFAP_TBR1_category_counts.tsv"))
  p_doublet <- ggplot(doublet_df, aes(x = GFAP, y = TBR1)) +
    geom_point(alpha = 0.15, size = 0.2) +
    theme_classic() +
    labs(title = "Xiang S3C approximation: GFAP vs TBR1", x = "GFAP log2 normalized", y = "TBR1 log2 normalized")
  save_plot(p_doublet, file.path(plot_dir, "Xiang_s3c_GFAP_vs_TBR1"), width = 6, height = 5)
}

set.seed(seed)
pc_embed <- Embeddings(obj, "pca")[, 1:5, drop = FALSE]
pair_records <- list()
for (pair_type in c("within_replicate", "between_replicate")) {
  collected <- 0L
  rows <- vector("list", 1000)
  attempts <- 0L
  while (collected < 1000L && attempts < 200000L) {
    attempts <- attempts + 1L
    pair <- sample(seq_len(nrow(pc_embed)), 2)
    same_rep <- obj$replicate[pair[1]] == obj$replicate[pair[2]]
    if ((pair_type == "within_replicate" && same_rep) || (pair_type == "between_replicate" && !same_rep)) {
      collected <- collected + 1L
      rows[[collected]] <- data.frame(
        pair_type = pair_type,
        cell_1 = rownames(pc_embed)[pair[1]],
        cell_2 = rownames(pc_embed)[pair[2]],
        euclidean_pc1_5 = sqrt(sum((pc_embed[pair[1], ] - pc_embed[pair[2], ]) ^ 2)),
        stringsAsFactors = FALSE
      )
    }
  }
  pair_records[[pair_type]] <- do.call(rbind, rows[seq_len(collected)])
}
pair_df <- do.call(rbind, pair_records)
write_tsv(pair_df, file.path(table_dir, "Xiang_s3b_random_pair_distances_pc1_5.tsv"))

saveRDS(obj, file.path(rds_dir, "Xiang_2017_tSNE_reproduction_modern_approximation.rds"))

report <- c(
  "# Xiang 2017 tSNE Reproduction Report",
  "",
  "## Status",
  "",
  "This run is a modern approximation of the Xiang et al. Seurat v1.4 workflow.",
  "Old Seurat v1.4.0.14 was not used.",
  "",
  "## Inputs",
  "",
  paste0("- Matrix: `", matrix_path, "`"),
  paste0("- Genes: `", genes_path, "`"),
  paste0("- Barcodes: `", barcodes_path, "`"),
  paste0("- GEO series matrix: `", series_path, "`"),
  "",
  "## Sample Metadata",
  "",
  "Barcode suffixes `-1` through `-8` were present in the GSE98201 barcodes.",
  "The suffix-to-sample mapping follows the GEO sample order for the eight scRNA-seq samples in GSE98201.",
  "No separate Cell Ranger aggregation CSV was found in the local download, so this mapping should be treated as GEO-order verified rather than independently aggregation-file verified.",
  "",
  "## Processing",
  "",
  "- Filtering: genes detected in more than 3 cells; cells with more than 200 genes.",
  "- Normalization: per-cell counts per 10,000 followed by log2(x + 1).",
  paste0("- Variable genes: dispersion > 0; selected for PCA/regression: ", length(variable_features), "."),
  paste0("- Regression: `sample` and `nCount_RNA`, Seurat ScaleData model `", regression_model, "`."),
  "- PCA: 20 PCs.",
  "- tSNE: PC1-PC5 only.",
  paste0("- Seed: ", seed, "."),
  "",
  "## Unknown Paper Parameters",
  "",
  "- tSNE perplexity",
  "- tSNE theta",
  "- tSNE number of iterations",
  "- tSNE seed",
  "- exact Seurat v1 variable-gene details beyond dispersion > 0",
  "- exact Cell Ranger aggregation CSV / suffix mapping file",
  "",
  "## Missing Markers",
  "",
  if (length(missing_markers) == 0) "- None" else paste0("- ", paste(missing_markers, collapse = ", ")),
  "",
  "## Outputs",
  "",
  paste0("- Tables: `", table_dir, "`"),
  paste0("- Plots: `", plot_dir, "`"),
  paste0("- RDS: `", rds_dir, "`")
)
writeLines(report, file.path(out_dir, "Xiang_tSNE_reproduction_report.md"))

log_msg("Done. Outputs written to: ", out_dir)
