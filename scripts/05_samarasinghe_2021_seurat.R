#!/usr/bin/env Rscript

# Build a Seurat object and UMAP for Samarasinghe_2021 (GSE165577) using the
# filtered counts matrix from GEO. This mirrors the Bershteyn/Walsh workflow
# structure: outputs go under results/samarasinghe_2021/.

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(Matrix)
  library(Seurat)
  library(ggplot2)
})

option_list <- list(
  make_option(c("-p", "--project-root"), type = "character", default = Sys.getenv("PROJECT_ROOT"),
              help = "Runtime workspace (PROJECT_ROOT). Required."),
  make_option(c("-c", "--counts"), type = "character", default = NULL,
              help = "Path to filtered counts CSV.gz (default: GEO filtered counts)."),
  make_option(c("-o", "--outdir"), type = "character", default = NULL,
              help = "Output directory (default: PROJECT_ROOT/results/samarasinghe_2021)")
)
opt <- parse_args(OptionParser(option_list = option_list))

if (opt$project_root == "") {
  stop("PROJECT_ROOT or --project-root is required")
}
project_root <- rtrim <- function(x) sub("/+$", "", x); project_root <- rtrim(opt$project_root)
counts_path <- opt$counts
if (is.null(counts_path)) {
  counts_path <- file.path(project_root, "data/raw/samarasinghe_2021_geo_files/suppl/GSE165577_Filtered_counts_all_samples.csv.gz")
}
if (!file.exists(counts_path)) stop("Counts file not found: ", counts_path)

outdir <- if (is.null(opt$outdir)) file.path(project_root, "results/samarasinghe_2021") else opt$outdir
plot_dir <- file.path(outdir, "plots")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

message("Reading counts: ", counts_path)
# data.table::fread handles gz; first column is gene symbols, header first cell is blank.
dt <- fread(counts_path)
gene_col <- dt[[1]]
dt[[1]] <- NULL

mat <- as.matrix(dt)
rownames(mat) <- gene_col

# Convert to sparse to reduce memory footprint
spmat <- Matrix(mat, sparse = TRUE)
rm(mat, dt); gc()

message("Creating Seurat object...")
seu <- CreateSeuratObject(counts = spmat, project = "Samarasinghe2021",
                          min.cells = 3, min.features = 200)

# Calculate percent.mt
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")

# Basic filtering (light) — adjust if needed
seu <- subset(seu, subset = nFeature_RNA > 500 & percent.mt < 20)

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

p2 <- ElbowPlot(seu, ndims = 50) + ggtitle("PCA elbow")
ggsave(file.path(plot_dir, "pca_elbow.png"), p2, width = 6, height = 4, dpi = 300)

message("Done. Outputs in: ", outdir)
