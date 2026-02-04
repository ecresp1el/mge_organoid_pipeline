#!/usr/bin/env Rscript

# Load Bershteyn_2023 provided Seurat object (GSE208672_Seurat_allsamples.rds.gz)
# and save UMAP plots and a copy of the object in results/bershteyn_2023.

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})

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

rds_path <- file.path(project_root, "data/raw/bershteyn_2023_geo_files/suppl/GSE208672_Seurat_allsamples.rds.gz")
stopifnot(file.exists(rds_path))

outdir <- file.path(project_root, "results/bershteyn_2023")
plot_dir <- file.path(outdir, "plots")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

message("Reading Seurat object: ", rds_path)
# File appears double-gzipped. Try nested gz connections first.
read_seurat_double_gz <- function(path) {
  # attempt gzcon(gzfile())
  con <- gzcon(gzfile(path, open = "rb"))
  on.exit(close(con), add = TRUE)
  tryCatch(readRDS(con), error = function(e) {
    # fallback: stream double gunzip to stdout into R
    cmd <- sprintf("gunzip -c %s | gunzip -c", shQuote(path))
    con2 <- pipe(cmd, open = "rb")
    on.exit(close(con2), add = TRUE)
    readRDS(con2)
  })
}
seu <- read_seurat_double_gz(rds_path)

# If UMAP missing, build quick neighbors/UMAP on existing PCA (dims 1:30)
if (!"umap" %in% Reductions(seu)) {
  if (!"pca" %in% Reductions(seu)) {
    seu <- RunPCA(seu, npcs = 50)
  }
  seu <- FindNeighbors(seu, dims = 1:30)
  seu <- FindClusters(seu, resolution = 0.5)
  seu <- RunUMAP(seu, dims = 1:30)
}

saveRDS(seu, file.path(outdir, "bershteyn_2023_seurat.rds"))

p1 <- DimPlot(seu, reduction = "umap", group.by = "seurat_clusters", label = TRUE) +
  ggtitle("Bershteyn 2023 UMAP (GSE208672)")
ggsave(file.path(plot_dir, "umap_by_cluster.png"), p1, width = 8, height = 6, dpi = 300)
ggsave(file.path(plot_dir, "umap_by_cluster.pdf"), p1, width = 8, height = 6)

message("Done. Outputs in ", outdir)
