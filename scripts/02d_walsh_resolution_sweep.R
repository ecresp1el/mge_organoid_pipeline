#!/usr/bin/env Rscript
# Walsh day-75: resolution sweep (fixed dims=1:20, k=20) starting from post-stress PCA checkpoint.

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
})

log_msg <- function(...) { cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste0(..., collapse = " "))); flush.console() }

args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args) >= 1) args[1] else Sys.getenv("PROJECT_ROOT", "")
if (project_root == "") stop("Provide PROJECT_ROOT as arg1 or env PROJECT_ROOT")
project_root <- normalizePath(project_root, mustWork = TRUE)

ckpt_dir <- file.path(project_root, "results", "walsh_day75", "checkpoints")
plots_dir <- file.path(project_root, "results", "walsh_day75", "resolution_sweep_plots")
out_dir <- file.path(project_root, "results", "walsh_day75")
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

pca_ckpt <- file.path(ckpt_dir, "walsh_pca_poststress.rds")
if (!file.exists(pca_ckpt)) stop("Missing PCA checkpoint: ", pca_ckpt)
obj <- readRDS(pca_ckpt)
log_msg("Loaded post-stress PCA checkpoint:", pca_ckpt)

res_vals <- c(1.8,1.6,1.4,1.2,1.0,0.8,0.6)
seeds <- list(neighbors=1004L, clusters=1005L, umap=1006L)
summary <- list()

for (res in res_vals) {
  log_msg("Running resolution=", res)
  set.seed(seeds$neighbors)
  obj <- FindNeighbors(obj, dims = 1:20, k.param = 20, algorithm = 1, verbose = FALSE)
  set.seed(seeds$clusters)
  obj <- FindClusters(obj, resolution = res, algorithm = 1, verbose = FALSE)
  set.seed(seeds$umap)
  red_name <- paste0("umap_res", gsub("\\.", "_", res))
  obj <- suppressWarnings(RunUMAP(obj, dims = 1:20, reduction.name = red_name, reduction.key = paste0("UMAPres", gsub("\\.", "_", res), "_"),
                 seed.use = seeds$umap, n.neighbors = 30, min.dist = 0.3, spread = 1,
                 metric = "cosine", umap.method = "uwot", return.model = FALSE, verbose = FALSE))

  p_cluster <- DimPlot(obj, reduction = red_name, group.by = "seurat_clusters", label = TRUE, repel = TRUE, pt.size = 0.4, shuffle = FALSE, order = TRUE, raster = TRUE) +
    ggtitle(paste0("UMAP (res=", res, ") by cluster"))
  p_domain <- DimPlot(obj, reduction = red_name, group.by = "domain", pt.size = 0.4, shuffle = FALSE, order = TRUE, raster = TRUE) +
    ggtitle(paste0("UMAP (res=", res, ") by domain"))

  f_clust <- file.path(plots_dir, paste0("umap_by_cluster_res", res, ".png"))
  f_domain <- file.path(plots_dir, paste0("umap_by_domain_res", res, ".png"))
  ggsave(f_clust, p_cluster, width = 8, height = 6, dpi = 300)
  ggsave(f_domain, p_domain, width = 8, height = 6, dpi = 300)

  n_clust <- length(unique(Idents(obj)))
  summary[[length(summary)+1]] <- data.frame(resolution=res, n_clusters=n_clust, cluster_plot=f_clust, domain_plot=f_domain)
  log_msg("res=", res, " produced ", n_clust, " clusters")
}

summary_df <- bind_rows(summary)
write.table(summary_df, file.path(out_dir, "walsh_day75_resolution_sweep_summary.tsv"), sep="\t", row.names = FALSE, quote = FALSE)
log_msg("Resolution sweep complete -> ", file.path(out_dir, "walsh_day75_resolution_sweep_summary.tsv"))
