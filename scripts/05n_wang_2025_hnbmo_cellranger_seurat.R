#!/usr/bin/env Rscript

# Wang et al. 2025 hnbMO exploratory Seurat processing from rerun
# Cell Ranger filtered_feature_bc_matrix outputs.
#
# This branch is intentionally separate from the raw-H5 reconstruction and is
# not registered as a cross-study input. It compares merged vs Seurat-integrated
# UMAPs and records cells removed by the additional Seurat QC layer.

suppressPackageStartupMessages({
  library(Matrix)
  library(Seurat)
  library(ggplot2)
  library(patchwork)
})

options(future.globals.maxSize = 8 * 1024^3)

parse_args <- function(args) {
  out <- list(
    `project-root` = Sys.getenv("PROJECT_ROOT", unset = ""),
    `cellranger-dir` = "",
    outdir = "",
    `min-features` = "200",
    `max-features` = "Inf",
    `max-percent-mt` = "20",
    `gene-min-cells` = "3",
    nfeatures = "2000",
    dims = "10",
    resolution = "0.5",
    seed = "11",
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
    if (a == "-p") a <- "--project-root"
    if (a == "-o") a <- "--outdir"
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
  cat(
    paste(
      "Usage:",
      "  Rscript scripts/05n_wang_2025_hnbmo_cellranger_seurat.R --project-root <PROJECT_ROOT> [options]",
      "",
      "Options:",
      "  --cellranger-dir <path>  Directory containing per-sample Cell Ranger count outputs",
      "  --outdir <path>          Output directory",
      "  --min-features <int>     Additional Seurat QC lower nFeature_RNA filter (default: 200)",
      "  --max-features <num>     Additional Seurat QC upper nFeature_RNA filter (default: Inf)",
      "  --max-percent-mt <num>   Additional Seurat QC mitochondrial filter (default: 20)",
      "  --gene-min-cells <int>   Keep genes expressed in at least this many cells (default: 3)",
      "  --nfeatures <int>        HVGs for FindVariableFeatures (default: 2000)",
      "  --dims <int>             PCs used for neighbors/UMAP (default: 10)",
      "  --resolution <num>       Cluster resolution (default: 0.5)",
      "  --seed <int>             Random seed (default: 11)",
      sep = "\n"
    )
  )
}

write_tsv <- function(df, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.table(df, path, sep = "\t", row.names = FALSE, quote = FALSE, na = "NA")
}

as_numeric_inf <- function(x) {
  if (tolower(as.character(x)) %in% c("inf", "infinity", "none", "na")) return(Inf)
  as.numeric(x)
}

read_cellranger_metrics <- function(path, sample_id) {
  if (!file.exists(path)) {
    return(data.frame(sample_id = sample_id, metric = character(), value = character()))
  }
  x <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  data.frame(
    sample_id = sample_id,
    metric = names(x),
    value = as.character(x[1, ]),
    stringsAsFactors = FALSE
  )
}

clean_metric_number <- function(x) {
  as.numeric(gsub("[,%\"]", "", x))
}

make_sample_umap <- function(obj, reduction = "umap", title = NULL) {
  DimPlot(
    obj,
    reduction = reduction,
    group.by = "figure1c_sample",
    pt.size = 0.25
  ) +
    ggtitle(title) +
    theme_classic(base_size = 13) +
    theme(plot.title = element_text(face = "bold"))
}

make_cluster_umap <- function(obj, reduction = "umap", title = NULL) {
  DimPlot(
    obj,
    reduction = reduction,
    group.by = "seurat_clusters",
    label = TRUE,
    repel = TRUE,
    pt.size = 0.25
  ) +
    ggtitle(title) +
    theme_classic(base_size = 13) +
    theme(plot.title = element_text(face = "bold"))
}

make_metric_umap <- function(obj, metric, title = NULL) {
  emb <- as.data.frame(Embeddings(obj, reduction = "umap"))
  emb[[metric]] <- obj@meta.data[[metric]]
  ggplot(emb, aes(x = UMAP_1, y = UMAP_2, color = .data[[metric]])) +
    geom_point(size = 0.25, alpha = 0.85) +
    scale_color_viridis_c(option = "magma") +
    labs(x = "UMAP_1", y = "UMAP_2", color = metric, title = title) +
    theme_classic(base_size = 13) +
    theme(plot.title = element_text(face = "bold"))
}

opt <- parse_args(commandArgs(trailingOnly = TRUE))
if (isTRUE(opt$help)) {
  print_usage()
  quit(save = "no", status = 0)
}

if (!nzchar(opt$`project-root`)) stop("PROJECT_ROOT or --project-root is required")
project_root <- sub("/+$", "", opt$`project-root`)
cellranger_dir <- if (nzchar(opt$`cellranger-dir`)) {
  opt$`cellranger-dir`
} else {
  file.path(project_root, "results/liu_2025_hnbmo_cellranger_counts")
}
outdir <- if (nzchar(opt$outdir)) {
  opt$outdir
} else {
  file.path(project_root, "results/wang_2025_hnbmo_cellranger_seurat_exploratory")
}
plot_dir <- file.path(outdir, "plots")
table_dir <- file.path(outdir, "tables")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

min_features <- as.integer(opt$`min-features`)
max_features <- as_numeric_inf(opt$`max-features`)
max_percent_mt <- as.numeric(opt$`max-percent-mt`)
gene_min_cells <- as.integer(opt$`gene-min-cells`)
nfeatures <- as.integer(opt$nfeatures)
dims_n <- as.integer(opt$dims)
resolution <- as.numeric(opt$resolution)
seed <- as.integer(opt$seed)
set.seed(seed)

sample_info <- data.frame(
  sample_id = c("BF_H9_D36", "BF_H9_D63", "BFCO_IMR_D63"),
  figure1c_sample = c("H9 D36", "H9 D63", "IMR90-4 D63"),
  sample_geo_accession = c("GSM8721440", "GSM8721441", "GSM8721442"),
  disease_status = "healthy",
  day = c("D36", "D63", "D63"),
  day_numeric = c(36L, 63L, 63L),
  cell_line = c("H9", "H9", "IMR90-4"),
  stringsAsFactors = FALSE
)
sample_info$filtered_matrix_dir <- file.path(cellranger_dir, sample_info$sample_id, "outs", "filtered_feature_bc_matrix")
sample_info$metrics_summary_csv <- file.path(cellranger_dir, sample_info$sample_id, "outs", "metrics_summary.csv")

missing_dirs <- sample_info$filtered_matrix_dir[!dir.exists(sample_info$filtered_matrix_dir)]
if (length(missing_dirs) > 0L) stop("Missing Cell Ranger filtered matrix dir(s): ", paste(missing_dirs, collapse = ", "))

write_tsv(sample_info, file.path(table_dir, "sample_manifest_used.tsv"))

write_tsv(
  data.frame(
    key = c(
      "status",
      "input",
      "scope",
      "cellranger_rerun",
      "seurat_qc",
      "merged_analysis",
      "integrated_analysis",
      "not_registered"
    ),
    value = c(
      "Exploratory Wang et al. 2025 hnbMO branch; do not include in cross-study panel yet.",
      "Cell Ranger 6.1.2 rerun filtered_feature_bc_matrix outputs from healthy SRA FASTQs.",
      "Healthy H9 D36, H9 D63, and IMR90-4 D63 only.",
      "Uses Great Lakes cellranger/6.1.2, SC3Pv3, GRCh38-2020-A, no BAM.",
      paste0("Additional Seurat QC after Cell Ranger filtering: nFeature_RNA >= ", min_features,
             "; nFeature_RNA < ", ifelse(is.infinite(max_features), "Inf", as.character(max_features)),
             "; percent.mt < ", max_percent_mt,
             "; gene min.cells=", gene_min_cells, "."),
      paste0("Merged object: NormalizeData, FindVariableFeatures nfeatures=", nfeatures,
             ", ScaleData, RunPCA, FindNeighbors/FindClusters/RunUMAP dims=1:", dims_n,
             ", resolution=", resolution, "."),
      "Seurat anchor integration across samples, followed by PCA/UMAP/clustering on integrated assay.",
      "Outputs are exploratory and not wired into study-level/cross-study manifests."
    ),
    stringsAsFactors = FALSE
  ),
  file.path(table_dir, "analysis_assumptions.tsv")
)

message("Reading Cell Ranger filtered matrices from: ", cellranger_dir)

objects <- list()
qc_rows <- list()
metrics_rows <- list()
qc_cell_rows <- list()

for (i in seq_len(nrow(sample_info))) {
  info <- sample_info[i, , drop = FALSE]
  message("Reading ", info$sample_id, ": ", info$filtered_matrix_dir)
  mat <- Read10X(data.dir = info$filtered_matrix_dir, gene.column = 2, unique.features = TRUE)
  if (is.list(mat)) mat <- mat[["Gene Expression"]]

  n_features <- Matrix::colSums(mat > 0)
  n_counts <- Matrix::colSums(mat)
  mito_features <- grepl("^MT-", rownames(mat))
  mt_counts <- if (any(mito_features)) Matrix::colSums(mat[mito_features, , drop = FALSE]) else rep(0, ncol(mat))
  percent_mt <- ifelse(n_counts > 0, 100 * mt_counts / n_counts, 0)

  pass_min_features <- n_features >= min_features
  pass_max_features <- if (is.infinite(max_features)) rep(TRUE, length(n_features)) else n_features < max_features
  pass_percent_mt <- percent_mt < max_percent_mt
  keep <- pass_min_features & pass_max_features & pass_percent_mt

  qc_cell_rows[[info$sample_id]] <- data.frame(
    sample_id = info$sample_id,
    figure1c_sample = info$figure1c_sample,
    barcode = colnames(mat),
    nCount_RNA = as.numeric(n_counts),
    nFeature_RNA = as.numeric(n_features),
    percent.mt = as.numeric(percent_mt),
    pass_min_features = pass_min_features,
    pass_max_features = pass_max_features,
    pass_percent_mt = pass_percent_mt,
    seurat_qc_keep = keep,
    stringsAsFactors = FALSE
  )

  qc_rows[[info$sample_id]] <- data.frame(
    sample_id = info$sample_id,
    figure1c_sample = info$figure1c_sample,
    cellranger_filtered_cells = ncol(mat),
    removed_by_seurat_qc = sum(!keep),
    kept_after_seurat_qc = sum(keep),
    fail_min_features = sum(!pass_min_features),
    fail_max_features = sum(!pass_max_features),
    fail_percent_mt = sum(!pass_percent_mt),
    median_nFeature_before_qc = median(n_features),
    median_nCount_before_qc = median(n_counts),
    median_percent_mt_before_qc = median(percent_mt),
    median_nFeature_after_qc = median(n_features[keep]),
    median_nCount_after_qc = median(n_counts[keep]),
    median_percent_mt_after_qc = median(percent_mt[keep]),
    stringsAsFactors = FALSE
  )

  metrics_rows[[info$sample_id]] <- read_cellranger_metrics(info$metrics_summary_csv, info$sample_id)

  obj <- CreateSeuratObject(
    counts = mat[, keep, drop = FALSE],
    project = info$sample_id,
    min.cells = gene_min_cells,
    min.features = 0
  )
  obj <- RenameCells(obj, add.cell.id = info$sample_id)
  obj$sample_id <- info$sample_id
  obj$figure1c_sample <- factor(info$figure1c_sample, levels = sample_info$figure1c_sample)
  obj$sample_geo_accession <- info$sample_geo_accession
  obj$disease_status <- info$disease_status
  obj$day <- info$day
  obj$day_numeric <- info$day_numeric
  obj$cell_line <- info$cell_line
  obj$cellranger_version <- "6.1.2"
  obj$source_branch <- "wang_2025_hnbmo_cellranger_seurat_exploratory"

  objects[[info$sample_id]] <- obj
}

qc_summary <- do.call(rbind, qc_rows)
write_tsv(qc_summary, file.path(table_dir, "seurat_qc_filtering_by_sample.tsv"))
write_tsv(do.call(rbind, qc_cell_rows), file.path(table_dir, "seurat_qc_per_cell.tsv"))
cellranger_metrics <- do.call(rbind, metrics_rows)
write_tsv(cellranger_metrics, file.path(table_dir, "cellranger_metrics_summary_long.tsv"))

cell_count_total <- data.frame(
  stage = c("cellranger_filtered", "after_additional_seurat_qc"),
  n_cells = c(sum(qc_summary$cellranger_filtered_cells), sum(qc_summary$kept_after_seurat_qc)),
  stringsAsFactors = FALSE
)
write_tsv(cell_count_total, file.path(table_dir, "cell_count_totals.tsv"))

qc_long <- do.call(
  rbind,
  lapply(qc_cell_rows, function(df) {
    data.frame(
      figure1c_sample = df$figure1c_sample,
      seurat_qc_status = ifelse(df$seurat_qc_keep, "kept", "filtered"),
      nCount_RNA = df$nCount_RNA,
      nFeature_RNA = df$nFeature_RNA,
      percent.mt = df$percent.mt,
      stringsAsFactors = FALSE
    )
  })
)
qc_long$figure1c_sample <- factor(qc_long$figure1c_sample, levels = sample_info$figure1c_sample)

qc_count_long <- rbind(
  data.frame(figure1c_sample = qc_summary$figure1c_sample, status = "kept", n_cells = qc_summary$kept_after_seurat_qc),
  data.frame(figure1c_sample = qc_summary$figure1c_sample, status = "filtered", n_cells = qc_summary$removed_by_seurat_qc)
)
qc_count_long$figure1c_sample <- factor(qc_count_long$figure1c_sample, levels = sample_info$figure1c_sample)
qc_count_long$status <- factor(qc_count_long$status, levels = c("filtered", "kept"))

p_qc_counts <- ggplot(qc_count_long, aes(x = figure1c_sample, y = n_cells, fill = status)) +
  geom_col(width = 0.7) +
  scale_fill_manual(values = c(kept = "#4C78A8", filtered = "#E45756")) +
  labs(x = NULL, y = "Cells", fill = NULL, title = "Cells retained after additional Seurat QC") +
  theme_classic(base_size = 13) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))
ggsave(file.path(plot_dir, "qc_cells_kept_filtered_by_sample.png"), p_qc_counts, width = 7, height = 5, dpi = 300)
ggsave(file.path(plot_dir, "qc_cells_kept_filtered_by_sample.pdf"), p_qc_counts, width = 7, height = 5)

qc_metric_long <- rbind(
  data.frame(figure1c_sample = qc_long$figure1c_sample, seurat_qc_status = qc_long$seurat_qc_status, metric = "nCount_RNA", value = qc_long$nCount_RNA),
  data.frame(figure1c_sample = qc_long$figure1c_sample, seurat_qc_status = qc_long$seurat_qc_status, metric = "nFeature_RNA", value = qc_long$nFeature_RNA),
  data.frame(figure1c_sample = qc_long$figure1c_sample, seurat_qc_status = qc_long$seurat_qc_status, metric = "percent.mt", value = qc_long$percent.mt)
)
p_qc_violin <- ggplot(qc_metric_long, aes(x = figure1c_sample, y = value, fill = seurat_qc_status)) +
  geom_violin(scale = "width", trim = TRUE, color = "gray25", linewidth = 0.2) +
  facet_wrap(~metric, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = c(kept = "#4C78A8", filtered = "#E45756")) +
  labs(x = NULL, y = NULL, fill = NULL, title = "QC metrics before/after additional Seurat filtering") +
  theme_classic(base_size = 13) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))
ggsave(file.path(plot_dir, "qc_metric_violins_kept_vs_filtered.png"), p_qc_violin, width = 11, height = 5, dpi = 300)
ggsave(file.path(plot_dir, "qc_metric_violins_kept_vs_filtered.pdf"), p_qc_violin, width = 11, height = 5)

message("Cells after Cell Ranger filtering: ", sum(qc_summary$cellranger_filtered_cells))
message("Cells after additional Seurat QC: ", sum(qc_summary$kept_after_seurat_qc))

message("Running merged non-integrated Seurat analysis")
merged <- Reduce(function(x, y) merge(x, y), objects)
merged$figure1c_sample <- factor(merged$figure1c_sample, levels = sample_info$figure1c_sample)
DefaultAssay(merged) <- "RNA"
merged <- NormalizeData(merged, verbose = FALSE)
merged <- FindVariableFeatures(merged, selection.method = "vst", nfeatures = nfeatures, verbose = FALSE)
merged <- ScaleData(merged, verbose = FALSE)
merged <- RunPCA(merged, npcs = max(30, dims_n), verbose = FALSE)
merged <- FindNeighbors(merged, dims = seq_len(dims_n), verbose = FALSE)
merged <- FindClusters(merged, resolution = resolution, verbose = FALSE)
merged <- RunUMAP(merged, dims = seq_len(dims_n), reduction.name = "umap", reduction.key = "UMAP_", verbose = FALSE)
saveRDS(merged, file.path(outdir, "wang_2025_hnbmo_cellranger_merged_seurat.rds"))

p_merged_sample <- make_sample_umap(merged, title = "Merged, no integration")
p_merged_cluster <- make_cluster_umap(merged, title = "Merged clusters")
p_merged_nfeature <- make_metric_umap(merged, "nFeature_RNA", title = "Merged nFeature_RNA after QC")
p_merged_ncount <- make_metric_umap(merged, "nCount_RNA", title = "Merged nCount_RNA after QC")
ggsave(file.path(plot_dir, "merged_umap_by_sample.png"), p_merged_sample, width = 7, height = 5.5, dpi = 300)
ggsave(file.path(plot_dir, "merged_umap_by_sample.pdf"), p_merged_sample, width = 7, height = 5.5)
ggsave(file.path(plot_dir, "merged_umap_by_cluster.png"), p_merged_cluster, width = 7, height = 5.5, dpi = 300)
ggsave(file.path(plot_dir, "merged_umap_by_cluster.pdf"), p_merged_cluster, width = 7, height = 5.5)
ggsave(file.path(plot_dir, "merged_umap_nFeature_RNA_after_qc.png"), p_merged_nfeature, width = 7, height = 5.5, dpi = 300)
ggsave(file.path(plot_dir, "merged_umap_nFeature_RNA_after_qc.pdf"), p_merged_nfeature, width = 7, height = 5.5)
ggsave(file.path(plot_dir, "merged_umap_nCount_RNA_after_qc.png"), p_merged_ncount, width = 7, height = 5.5, dpi = 300)
ggsave(file.path(plot_dir, "merged_umap_nCount_RNA_after_qc.pdf"), p_merged_ncount, width = 7, height = 5.5)

message("Running Seurat anchor integration")
split_objects <- objects
split_objects <- lapply(split_objects, function(x) {
  DefaultAssay(x) <- "RNA"
  x <- NormalizeData(x, verbose = FALSE)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = nfeatures, verbose = FALSE)
  x
})
anchors <- FindIntegrationAnchors(object.list = split_objects, dims = seq_len(dims_n), verbose = FALSE)
integrated <- IntegrateData(anchorset = anchors, dims = seq_len(dims_n), verbose = FALSE)
DefaultAssay(integrated) <- "integrated"
integrated <- ScaleData(integrated, verbose = FALSE)
integrated <- RunPCA(integrated, npcs = max(30, dims_n), verbose = FALSE)
integrated <- FindNeighbors(integrated, dims = seq_len(dims_n), verbose = FALSE)
integrated <- FindClusters(integrated, resolution = resolution, verbose = FALSE)
integrated <- RunUMAP(integrated, dims = seq_len(dims_n), reduction.name = "umap", reduction.key = "UMAP_", verbose = FALSE)
integrated$figure1c_sample <- factor(integrated$figure1c_sample, levels = sample_info$figure1c_sample)
saveRDS(integrated, file.path(outdir, "wang_2025_hnbmo_cellranger_integrated_seurat.rds"))

p_integrated_sample <- make_sample_umap(integrated, title = "Seurat integrated")
p_integrated_cluster <- make_cluster_umap(integrated, title = "Integrated clusters")
p_integrated_nfeature <- make_metric_umap(integrated, "nFeature_RNA", title = "Integrated nFeature_RNA after QC")
p_integrated_ncount <- make_metric_umap(integrated, "nCount_RNA", title = "Integrated nCount_RNA after QC")
ggsave(file.path(plot_dir, "integrated_umap_by_sample.png"), p_integrated_sample, width = 7, height = 5.5, dpi = 300)
ggsave(file.path(plot_dir, "integrated_umap_by_sample.pdf"), p_integrated_sample, width = 7, height = 5.5)
ggsave(file.path(plot_dir, "integrated_umap_by_cluster.png"), p_integrated_cluster, width = 7, height = 5.5, dpi = 300)
ggsave(file.path(plot_dir, "integrated_umap_by_cluster.pdf"), p_integrated_cluster, width = 7, height = 5.5)
ggsave(file.path(plot_dir, "integrated_umap_nFeature_RNA_after_qc.png"), p_integrated_nfeature, width = 7, height = 5.5, dpi = 300)
ggsave(file.path(plot_dir, "integrated_umap_nFeature_RNA_after_qc.pdf"), p_integrated_nfeature, width = 7, height = 5.5)
ggsave(file.path(plot_dir, "integrated_umap_nCount_RNA_after_qc.png"), p_integrated_ncount, width = 7, height = 5.5, dpi = 300)
ggsave(file.path(plot_dir, "integrated_umap_nCount_RNA_after_qc.pdf"), p_integrated_ncount, width = 7, height = 5.5)

p_compare <- p_merged_sample + p_integrated_sample + plot_layout(ncol = 2)
ggsave(file.path(plot_dir, "merged_vs_integrated_umap_by_sample.png"), p_compare, width = 13, height = 5.5, dpi = 300)
ggsave(file.path(plot_dir, "merged_vs_integrated_umap_by_sample.pdf"), p_compare, width = 13, height = 5.5)

p_metric_compare <- (p_merged_nfeature + p_integrated_nfeature) / (p_merged_ncount + p_integrated_ncount)
ggsave(file.path(plot_dir, "merged_vs_integrated_umap_qc_metrics_after_qc.png"), p_metric_compare, width = 13, height = 11, dpi = 300)
ggsave(file.path(plot_dir, "merged_vs_integrated_umap_qc_metrics_after_qc.pdf"), p_metric_compare, width = 13, height = 11)

cluster_counts_merged <- as.data.frame.matrix(table(merged$seurat_clusters, merged$figure1c_sample))
cluster_counts_merged$seurat_cluster <- rownames(cluster_counts_merged)
cluster_counts_merged <- cluster_counts_merged[, c("seurat_cluster", sample_info$figure1c_sample)]
write_tsv(cluster_counts_merged, file.path(table_dir, "merged_cluster_counts_by_sample.tsv"))

cluster_counts_integrated <- as.data.frame.matrix(table(integrated$seurat_clusters, integrated$figure1c_sample))
cluster_counts_integrated$seurat_cluster <- rownames(cluster_counts_integrated)
cluster_counts_integrated <- cluster_counts_integrated[, c("seurat_cluster", sample_info$figure1c_sample)]
write_tsv(cluster_counts_integrated, file.path(table_dir, "integrated_cluster_counts_by_sample.tsv"))

message("Wrote exploratory outputs to: ", outdir)
