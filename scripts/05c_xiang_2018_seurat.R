#!/usr/bin/env Rscript

# Build a Seurat object and UMAP for Xiang_2018 scRNA (GSE98201 10x matrix).

suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
  library(Matrix)
  library(ggplot2)
})

# Minimal arg parse
args <- commandArgs(trailingOnly = TRUE)
arg_list <- list()
for (i in seq_along(args)) {
  a <- args[i]
  if (grepl("^--", a)) {
    if (grepl("=", a)) {
      kv <- strsplit(sub("^--", "", a), "=", fixed = TRUE)[[1]]
      arg_list[[kv[1]]] <- if (length(kv) > 1) kv[2] else NA
    } else {
      key <- sub("^--", "", a)
      val <- if (i < length(args)) args[i + 1] else NA
      arg_list[[key]] <- val
    }
  }
}
project_root <- Sys.getenv("PROJECT_ROOT")
if (!is.null(arg_list[["project-root"]])) project_root <- arg_list[["project-root"]]
if (project_root == "" || is.na(project_root)) stop("PROJECT_ROOT or --project-root is required")
project_root <- sub("/+$", "", project_root)

suppl_dir <- file.path(project_root, "data/raw/xiang_2018_geo_files/suppl")
mtx <- file.path(suppl_dir, "GSE98201_matrix.mtx.gz")
genes <- file.path(suppl_dir, "GSE98201_genes.tsv.gz")
barcodes <- file.path(suppl_dir, "GSE98201_barcodes.tsv.gz")
stopifnot(file.exists(mtx), file.exists(genes), file.exists(barcodes))

outdir <- file.path(project_root, "results/xiang_2018")
plot_dir <- file.path(outdir, "plots")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

message("Reading 10x matrix...")
# Seurat 4.1 expects default filenames; create a temp view with symlinks.
tmpdir <- file.path(outdir, "tmp_read10x")
dir.create(tmpdir, showWarnings = FALSE, recursive = TRUE)
file.symlink(mtx, file.path(tmpdir, "matrix.mtx.gz"))
file.symlink(genes, file.path(tmpdir, "genes.tsv.gz"))
file.symlink(genes, file.path(tmpdir, "features.tsv.gz"))
file.symlink(barcodes, file.path(tmpdir, "barcodes.tsv.gz"))

counts <- Read10X(data.dir = tmpdir, gene.column = 1, unique.features = TRUE)

unlink(tmpdir, recursive = TRUE, force = TRUE)

message("Creating Seurat object...")
seu <- CreateSeuratObject(counts, project = "Xiang2018", min.cells = 3, min.features = 200)

seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
seu <- subset(seu, subset = nFeature_RNA > 500 & percent.mt < 15)

seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 3000)
seu <- ScaleData(seu, vars.to.regress = "percent.mt")
seu <- RunPCA(seu, npcs = 50)
seu <- FindNeighbors(seu, dims = 1:30)
seu <- FindClusters(seu, resolution = 0.5)
seu <- RunUMAP(seu, dims = 1:30)

saveRDS(seu, file.path(outdir, "xiang_2018_seurat.rds"))

p1 <- DimPlot(seu, reduction = "umap", group.by = "seurat_clusters", label = TRUE) +
  ggtitle("Xiang 2018 UMAP (GSE98201)")
ggsave(file.path(plot_dir, "umap_by_cluster.png"), p1, width = 8, height = 6, dpi = 300)
ggsave(file.path(plot_dir, "umap_by_cluster.pdf"), p1, width = 8, height = 6)

p2 <- ElbowPlot(seu, ndims = 50) + ggtitle("PCA elbow")
ggsave(file.path(plot_dir, "pca_elbow.png"), p2, width = 6, height = 4, dpi = 300)

message("Done. Outputs in ", outdir)
