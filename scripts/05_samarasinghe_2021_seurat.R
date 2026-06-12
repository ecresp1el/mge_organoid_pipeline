#!/usr/bin/env Rscript

# Build a Seurat object and UMAP for Samarasinghe_2021 (GSE165577) using the
# filtered counts matrix from GEO. This mirrors the Bershteyn/Walsh workflow
# structure: outputs go under results/samarasinghe_2021/.

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(Seurat)
  library(ggplot2)
})

parse_args <- function(args) {
  out <- list()
  i <- 1
  while (i <= length(args)) {
    key <- args[[i]]
    if (startsWith(key, "--")) {
      key <- sub("^--", "", key)
      if (grepl("=", key, fixed = TRUE)) {
        parts <- strsplit(key, "=", fixed = TRUE)[[1]]
        out[[parts[[1]]]] <- if (length(parts) > 1) parts[[2]] else ""
      } else {
        out[[key]] <- if (i < length(args)) args[[i + 1]] else ""
        i <- i + 1
      }
    }
    i <- i + 1
  }
  out
}

arg_list <- parse_args(commandArgs(trailingOnly = TRUE))
get_arg <- function(name, default = "") {
  if (!is.null(arg_list[[name]]) && nzchar(arg_list[[name]])) {
    return(arg_list[[name]])
  }
  env_name <- toupper(gsub("-", "_", name))
  value <- Sys.getenv(env_name, unset = "")
  if (nzchar(value)) value else default
}

project_root <- get_arg("project-root", Sys.getenv("PROJECT_ROOT", unset = ""))
if (project_root == "") {
  stop("PROJECT_ROOT or --project-root is required")
}
project_root <- sub("/+$", "", project_root)
counts_path <- get_arg("counts", "")
if (!nzchar(counts_path)) {
  counts_path <- file.path(project_root, "data/raw/samarasinghe_2021_geo_files/suppl/GSE165577_Filtered_counts_all_samples.csv.gz")
}
if (!file.exists(counts_path)) stop("Counts file not found: ", counts_path)

outdir <- get_arg("outdir", "")
if (!nzchar(outdir)) {
  outdir <- file.path(project_root, "results/samarasinghe_2021")
}
plot_dir <- file.path(outdir, "plots")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

message("Reading counts: ", counts_path)
# Stream gzip explicitly so this does not depend on the optional R.utils
# package in data.table::fread().
dt <- fread(cmd = paste("gzip -dc", shQuote(counts_path)))
gene_col <- dt[[1]]
dt[[1]] <- NULL

mat <- as.matrix(dt)
rownames(mat) <- gene_col

# Convert to sparse to reduce memory footprint
spmat <- Matrix(mat, sparse = TRUE)
rm(mat, dt); gc()

message("Applying pre-Seurat QC filter...")
n_features_by_cell <- Matrix::colSums(spmat > 0)
n_counts_by_cell <- Matrix::colSums(spmat)
mito_features <- grepl("^MT-", rownames(spmat))
if (any(mito_features)) {
  mt_counts_by_cell <- Matrix::colSums(spmat[mito_features, , drop = FALSE])
} else {
  mt_counts_by_cell <- rep(0, ncol(spmat))
}
percent_mt_by_cell <- ifelse(n_counts_by_cell > 0, 100 * mt_counts_by_cell / n_counts_by_cell, 0)
keep_cells <- n_features_by_cell > 500 & percent_mt_by_cell < 20
message("Cells before QC: ", ncol(spmat), "; after QC: ", sum(keep_cells))
spmat <- spmat[, keep_cells, drop = FALSE]
n_features_by_cell <- n_features_by_cell[keep_cells]
n_counts_by_cell <- n_counts_by_cell[keep_cells]
percent_mt_by_cell <- percent_mt_by_cell[keep_cells]

message("Creating Seurat object...")
seu <- CreateSeuratObject(counts = spmat, project = "Samarasinghe2021",
                          min.cells = 3, min.features = 0)

# Sample prefixes are embedded in cell names, e.g. d56_ctrl_AAAC... and
# d100_ctrl_docked_2_AAAC...
sample_id <- sub("_[^_]+$", "", Cells(seu))
seu$sample_id <- sample_id
seu$orig.ident <- sample_id
seu$culture_day <- sub("^d([0-9]+).*$", "d\\1", sample_id)
seu$condition <- sub("^d[0-9]+_([^_]+).*$", "\\1", sample_id)
seu$sample_detail <- sub("^d[0-9]+_[^_]+_?", "", sample_id)
seu$sample_detail[seu$sample_detail == seu$sample_id] <- ""
seu$nFeature_RNA_prefilter <- as.numeric(n_features_by_cell[Cells(seu)])
seu$nCount_RNA_prefilter <- as.numeric(n_counts_by_cell[Cells(seu)])
seu$percent.mt.prefilter <- as.numeric(percent_mt_by_cell[Cells(seu)])

# Calculate percent.mt
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")

message("Normalization / HVG / scaling / PCA / UMAP / clustering...")
seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 3000)
seu <- ScaleData(seu, vars.to.regress = c("percent.mt"))
seu <- RunPCA(seu, npcs = 50)
seu <- FindNeighbors(seu, dims = 1:30)
seu <- FindClusters(seu, resolution = 0.8)
seu <- RunUMAP(seu, dims = 1:30)

message("Saving Seurat object and plots...")
saveRDS(seu, file.path(outdir, "samarasinghe_2021_seurat.rds"))

p1 <- DimPlot(seu, reduction = "umap", group.by = "seurat_clusters", label = TRUE) +
  ggtitle("Samarasinghe 2021 UMAP by cluster (res=0.8)")
ggsave(file.path(plot_dir, "umap_by_cluster.png"), p1, width = 8, height = 6, dpi = 300)
ggsave(file.path(plot_dir, "umap_by_cluster.pdf"), p1, width = 8, height = 6)

p_sample <- DimPlot(seu, reduction = "umap", group.by = "sample_id") +
  ggtitle("Samarasinghe 2021 UMAP by sample")
ggsave(file.path(plot_dir, "umap_by_sample.png"), p_sample, width = 8, height = 6, dpi = 300)
ggsave(file.path(plot_dir, "umap_by_sample.pdf"), p_sample, width = 8, height = 6)

p2 <- ElbowPlot(seu, ndims = 50) + ggtitle("PCA elbow")
ggsave(file.path(plot_dir, "pca_elbow.png"), p2, width = 6, height = 4, dpi = 300)

message("Done. Outputs in: ", outdir)
