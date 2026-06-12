#!/usr/bin/env Rscript

# Reproduce Samarasinghe_2021 integration/UMAP as described in the paper.
# Deviations (explicit): starting from GEO filtered counts (no FASTQs, so no Cell Ranger),
# using Seurat >=4.x + SeuratWrappers instead of Seurat v3.2.0, and the provided
# GRCh38 reference is not re-run here.

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(Seurat)
  library(ggplot2)
})

required_pkgs <- c("SeuratWrappers", "rliger")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop(
    "Missing required module-provided R package(s): ",
    paste(missing_pkgs, collapse = ", "),
    "\nLoad a module stack or R_LIBS_USER that provides these packages before running this script.",
    "\nCurrent .libPaths():\n  ",
    paste(.libPaths(), collapse = "\n  ")
  )
}

suppressPackageStartupMessages({
  library(SeuratWrappers)
  library(rliger)
})

# Minimal arg parse (no optparse)
args <- commandArgs(trailingOnly = TRUE)
arg_list <- list()
i <- 1
while (i <= length(args)) {
  arg <- args[[i]]
  if (grepl("^--", arg)) {
    keyval <- sub("^--", "", arg)
    if (grepl("=", keyval, fixed = TRUE)) {
      kv <- strsplit(keyval, "=", fixed = TRUE)[[1]]
      arg_list[[kv[[1]]]] <- if (length(kv) > 1) kv[[2]] else ""
    } else {
      arg_list[[keyval]] <- if (i < length(args) && !grepl("^--", args[[i + 1]])) args[[i + 1]] else ""
      if (i < length(args) && !grepl("^--", args[[i + 1]])) i <- i + 1
    }
  }
  i <- i + 1
}
arg_value <- function(name, default = "") {
  value <- arg_list[[name]]
  if (!is.null(value) && nzchar(value)) value else default
}
arg_flag <- function(name, default = FALSE) {
  value <- arg_list[[name]]
  if (is.null(value) || !nzchar(value)) return(default)
  tolower(value) %in% c("1", "true", "t", "yes", "y")
}

project_root <- arg_value("project-root", Sys.getenv("PROJECT_ROOT"))
if (project_root == "") stop("PROJECT_ROOT or --project-root is required")
project_root <- sub("/+$", "", project_root)
control_only <- arg_flag("control-only", FALSE)

counts_path <- arg_value("counts", file.path(
  project_root,
  "data/raw/samarasinghe_2021_geo_files/suppl/GSE165577_Filtered_counts_all_samples.csv.gz"
))
if (!file.exists(counts_path)) stop("Counts file not found: ", counts_path)

outdir_default <- if (control_only) {
  file.path(project_root, "results/samarasinghe_2021_liger_control_only")
} else {
  file.path(project_root, "results/samarasinghe_2021_liger")
}
outdir <- arg_value("outdir", outdir_default)
plot_dir <- file.path(outdir, "plots")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

message("Reading counts: ", counts_path)
dt <- data.table::fread(cmd = paste("gzip -dc", shQuote(counts_path)))
genes <- dt[[1]]
dt[[1]] <- NULL
mat <- as.matrix(dt)
rownames(mat) <- genes
spmat <- Matrix(mat, sparse = TRUE)
rm(mat, dt); gc()

# Derive sample ID (orig.ident) from column names: prefix before last underscore.
cells <- colnames(spmat)
sample_id <- sub("_[^_]+$", "", cells)
condition <- ifelse(grepl("_rett", sample_id), "MECP2_mutant_Rett", "control_WT_MECP2")

if (control_only) {
  control_samples <- c("d56_ctrl", "d70_ctrl", "d100_ctrl_docked_2")
  keep <- sample_id %in% control_samples
  message("Control-only mode: retaining ", sum(keep), " of ", length(keep), " cells")
  if (!any(keep)) stop("No control cells found in counts matrix")
  spmat <- spmat[, keep, drop = FALSE]
  sample_id <- sample_id[keep]
  condition <- condition[keep]
  cells <- cells[keep]
}

seu <- CreateSeuratObject(counts = spmat, project = "Samarasinghe2021",
                          min.cells = 1, min.features = 1)
seu$orig.ident <- sample_id
seu$sample_id <- sample_id
seu$culture_day <- sub("^d([0-9]+).*$", "d\\1", sample_id)
seu$condition <- condition

# QC metrics
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")

# Filtering: nFeature_RNA > 500, nFeature_RNA < mean+3*sd, percent.mt < 10%
upper_thresh <- mean(seu$nFeature_RNA) + 3 * sd(seu$nFeature_RNA)
qc_keep <- seu$nFeature_RNA > 500 & seu$nFeature_RNA < upper_thresh & seu$percent.mt < 10
message("QC filter retained ", sum(qc_keep), " of ", length(qc_keep), " cells")
seu <- seu[, qc_keep]

# Normalize and HVG as in paper (default Seurat)
seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu)

# Scale per batch without centering (split.by = orig.ident, do.center = FALSE)
seu <- ScaleData(seu, split.by = "orig.ident", do.center = FALSE)

# LIGER integration. SeuratWrappers::RunOptimizeALS is not compatible with
# rliger 2.2.x, so call the current rliger Seurat methods directly while keeping
# the paper parameters and Seurat-scaled input.
seu <- rliger::runINMF(
  seu,
  k = 20,
  lambda = 5,
  datasetVar = "orig.ident",
  layer = "scale.data",
  reduction = "inmf"
)
seu <- rliger::quantileNorm(
  seu,
  reduction = "inmf",
  quantiles = 50,
  minCells = 20,
  nNeighbors = 20,
  center = FALSE
)
seu[["iNMF"]] <- CreateDimReducObject(
  embeddings = Embeddings(seu, "inmfNorm"),
  loadings = Loadings(seu, "inmfNorm"),
  assay = DefaultAssay(seu),
  key = "iNMF_"
)

# Clustering on iNMF
seu <- FindNeighbors(seu, reduction = "iNMF", dims = 1:20)
seu <- FindClusters(seu, resolution = 0.3)

# UMAP on iNMF dims 1:ncol(iNMF)
iNMF_dims <- 1:ncol(Embeddings(seu, "iNMF"))
seu <- RunUMAP(seu, dims = iNMF_dims, reduction = "iNMF")

# Save objects and plots
saveRDS(seu, file.path(outdir, "samarasinghe_2021_liger.rds"))

sample_counts <- aggregate(
  cell_count ~ sample_id + culture_day + condition,
  data = data.frame(
    sample_id = seu$sample_id,
    culture_day = seu$culture_day,
    condition = seu$condition,
    cell_count = 1L
  ),
  FUN = sum
)
write.csv(sample_counts, file.path(outdir, "sample_cell_counts_after_qc.csv"), row.names = FALSE)

p1 <- DimPlot(seu, reduction = "umap", group.by = "seurat_clusters", label = TRUE) +
  ggtitle("Samarasinghe 2021 UMAP (LIGER iNMF, res=0.3)")
ggsave(file.path(plot_dir, "umap_by_cluster_liger.png"), p1, width = 8, height = 6, dpi = 300)
ggsave(file.path(plot_dir, "umap_by_cluster_liger.pdf"), p1, width = 8, height = 6)

p_sample <- DimPlot(seu, reduction = "umap", group.by = "sample_id") +
  ggtitle("Samarasinghe 2021 UMAP by sample (LIGER iNMF)")
ggsave(file.path(plot_dir, "umap_by_sample_liger.png"), p_sample, width = 8, height = 6, dpi = 300)
ggsave(file.path(plot_dir, "umap_by_sample_liger.pdf"), p_sample, width = 8, height = 6)

p_condition <- DimPlot(seu, reduction = "umap", group.by = "condition") +
  ggtitle("Samarasinghe 2021 UMAP by condition (LIGER iNMF)")
ggsave(file.path(plot_dir, "umap_by_condition_liger.png"), p_condition, width = 8, height = 6, dpi = 300)
ggsave(file.path(plot_dir, "umap_by_condition_liger.pdf"), p_condition, width = 8, height = 6)

p_day <- DimPlot(seu, reduction = "umap", group.by = "culture_day") +
  ggtitle("Samarasinghe 2021 UMAP by culture day (LIGER iNMF)")
ggsave(file.path(plot_dir, "umap_by_culture_day_liger.png"), p_day, width = 8, height = 6, dpi = 300)
ggsave(file.path(plot_dir, "umap_by_culture_day_liger.pdf"), p_day, width = 8, height = 6)

control_cells <- Cells(seu)[seu$condition == "control_WT_MECP2"]
if (!control_only && length(control_cells) > 0) {
  p_control <- DimPlot(
    seu,
    reduction = "umap",
    cells = control_cells,
    group.by = "seurat_clusters",
    label = TRUE
  ) + ggtitle("Samarasinghe 2021 WT-control cells on all-sample LIGER UMAP")
  ggsave(file.path(plot_dir, "umap_control_only_by_cluster_liger.png"), p_control, width = 8, height = 6, dpi = 300)
  ggsave(file.path(plot_dir, "umap_control_only_by_cluster_liger.pdf"), p_control, width = 8, height = 6)
}

message("Skipping iNMF elbow plot; this reduction does not store PCA-style standard deviations.")

message("Done. Outputs in: ", outdir)
