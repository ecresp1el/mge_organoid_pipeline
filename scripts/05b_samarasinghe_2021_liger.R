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

# Ensure required wrappers/liger are available; install into a temp library if missing.
lib_extra <- file.path(tempdir(), "r_libs")
dir.create(lib_extra, showWarnings = FALSE, recursive = TRUE)
.libPaths(c(lib_extra, .libPaths()))

need_pkgs <- c("SeuratWrappers", "rliger")
dep_pkgs <- c("S4Vectors","DelayedArray","HDF5Array","MatrixGenerics","rhdf5lib","rhdf5",
              "RcppPlanc","leidenAlg","sccore","pbmcapply","pROC","checkmate",
              "HighFive","hdf5r","hdf5r.Extra","bit","bit64","backports",
              "fastmatch","easy.utils")
all_pkgs <- unique(c(need_pkgs, dep_pkgs))

to_install <- all_pkgs[!all_pkgs %in% rownames(installed.packages())]
cran_repo <- "https://cloud.r-project.org"
if (length(to_install) > 0) {
  install.packages(to_install, repos = cran_repo, lib = lib_extra, quiet = TRUE)
}

# Install any remaining via Bioconductor if needed
remaining <- all_pkgs[!all_pkgs %in% rownames(installed.packages())]
if (length(remaining) > 0) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = cran_repo, lib = lib_extra, quiet = TRUE)
  }
  BiocManager::install(remaining, lib = lib_extra, ask = FALSE, update = FALSE)
}

suppressPackageStartupMessages({
  library(SeuratWrappers)
  library(rliger)
})

# Minimal arg parse (no optparse)
args <- commandArgs(trailingOnly = TRUE)
arg_list <- list()
for (a in args) {
  if (grepl("^--", a)) {
    kv <- strsplit(sub("^--", "", a), "=", fixed = TRUE)[[1]]
    key <- kv[1]; val <- if (length(kv) > 1) kv[2] else NA
    arg_list[[key]] <- val
  }
}
project_root <- if (!is.null(arg_list[["project-root"]])) arg_list[["project-root"]] else Sys.getenv("PROJECT_ROOT")
if (project_root == "") stop("PROJECT_ROOT or --project-root is required")
project_root <- sub("/+$", "", project_root)

counts_path <- if (!is.null(arg_list[["counts"]])) {
  arg_list[["counts"]]
} else {
  file.path(project_root, "data/raw/samarasinghe_2021_geo_files/suppl/GSE165577_Filtered_counts_all_samples.csv.gz")
}
if (!file.exists(counts_path)) stop("Counts file not found: ", counts_path)

outdir <- if (!is.null(arg_list[["outdir"]])) arg_list[["outdir"]] else file.path(project_root, "results/samarasinghe_2021_liger")
plot_dir <- file.path(outdir, "plots")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

message("Reading counts: ", counts_path)
dt <- data.table::fread(counts_path)
genes <- dt[[1]]
dt[[1]] <- NULL
mat <- as.matrix(dt)
rownames(mat) <- genes
spmat <- Matrix(mat, sparse = TRUE)
rm(mat, dt); gc()

# Derive sample ID (orig.ident) from column names: prefix before last underscore.
cells <- colnames(spmat)
sample_id <- sub("_[^_]+$", "", cells)

seu <- CreateSeuratObject(counts = spmat, project = "Samarasinghe2021",
                          min.cells = 1, min.features = 1)
seu$orig.ident <- sample_id

# QC metrics
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")

# Filtering: nFeature_RNA > 500, nFeature_RNA < mean+3*sd, percent.mt < 10%
upper_thresh <- mean(seu$nFeature_RNA) + 3 * sd(seu$nFeature_RNA)
seu <- subset(seu, subset = nFeature_RNA > 500 & nFeature_RNA < upper_thresh & percent.mt < 10)

# Normalize and HVG as in paper (default Seurat)
seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu)

# Scale per batch without centering (split.by = orig.ident, do.center = FALSE)
seu <- ScaleData(seu, split.by = "orig.ident", do.center = FALSE)

# LIGER integration via SeuratWrappers
seu <- RunOptimizeALS(seu, k = 20, lambda = 5, split.by = "orig.ident")
seu <- RunQuantileNorm(seu, split.by = "orig.ident")

# Clustering on iNMF
seu <- FindNeighbors(seu, reduction = "iNMF", dims = 1:20)
seu <- FindClusters(seu, resolution = 0.3)

# UMAP on iNMF dims 1:ncol(iNMF)
iNMF_dims <- 1:ncol(seu[["iNMF"]])
seu <- RunUMAP(seu, dims = iNMF_dims, reduction = "iNMF")

# Save objects and plots
saveRDS(seu, file.path(outdir, "samarasinghe_2021_liger.rds"))

p1 <- DimPlot(seu, reduction = "umap", group.by = "seurat_clusters", label = TRUE) +
  ggtitle("Samarasinghe 2021 UMAP (LIGER iNMF, res=0.3)")
ggsave(file.path(plot_dir, "umap_by_cluster_liger.png"), p1, width = 8, height = 6, dpi = 300)
ggsave(file.path(plot_dir, "umap_by_cluster_liger.pdf"), p1, width = 8, height = 6)

p2 <- ElbowPlot(seu, ndims = min(50, length(iNMF_dims))) + ggtitle("iNMF elbow (proxy)")
ggsave(file.path(plot_dir, "inmf_elbow.png"), p2, width = 6, height = 4, dpi = 300)

message("Done. Outputs in: ", outdir)
