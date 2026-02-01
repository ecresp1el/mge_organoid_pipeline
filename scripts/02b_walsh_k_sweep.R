#!/usr/bin/env Rscript
# Walsh day-75: neighbor k sweep post-stress using saved PCA checkpoint.

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
})

log_msg <- function(...) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste0(..., collapse = " ")))
  flush.console()
}

args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args) >= 1) args[1] else Sys.getenv("PROJECT_ROOT", "")
if (project_root == "") stop("Provide PROJECT_ROOT as arg1 or env PROJECT_ROOT")
project_root <- normalizePath(project_root, mustWork = TRUE)

ckpt_dir <- file.path(project_root, "results", "walsh_day75", "checkpoints")
plots_dir <- file.path(project_root, "results", "walsh_day75", "k_sweep_plots")
out_dir <- file.path(project_root, "results", "walsh_day75")
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

pca_ckpt <- file.path(ckpt_dir, "walsh_pca_poststress.rds")
if (!file.exists(pca_ckpt)) stop("Missing PCA checkpoint: ", pca_ckpt)

obj <- readRDS(pca_ckpt)
log_msg("Loaded post-stress PCA checkpoint:", pca_ckpt)

k_vals <- c(20,30,40,50)
summary <- list()
seeds <- list(neighbors=1004L, clusters=1005L, umap=1006L)

for (k in k_vals) {
  log_msg("Running k.param=", k)
  set.seed(seeds$neighbors)
  obj <- FindNeighbors(obj, dims = 1:20, k.param = k, algorithm = 1, verbose = FALSE)
  set.seed(seeds$clusters)
  obj <- FindClusters(obj, resolution = 2.0, algorithm = 1, verbose = FALSE)
  set.seed(seeds$umap)
  red_name <- paste0("umap_k", k)
  obj <- suppressWarnings(RunUMAP(obj, dims = 1:20, reduction.name = red_name, reduction.key = paste0("UMAP", k, "_"),
                 seed.use = seeds$umap, n.neighbors = 30, min.dist = 0.3, spread = 1,
                 metric = "cosine", umap.method = "uwot", return.model = FALSE, verbose = FALSE))

  # Plots
  p_cluster <- DimPlot(obj, reduction = red_name, group.by = "seurat_clusters", label = TRUE, repel = TRUE, pt.size = 0.4, shuffle = FALSE, order = TRUE, raster = TRUE) +
    ggtitle(paste0("UMAP (", red_name, ") by cluster"))
  p_domain <- DimPlot(obj, reduction = red_name, group.by = "domain", pt.size = 0.4, shuffle = FALSE, order = TRUE, raster = TRUE) +
    ggtitle(paste0("UMAP (", red_name, ") by domain"))

  f_clust <- file.path(plots_dir, paste0("umap_by_cluster_k", k, ".png"))
  f_domain <- file.path(plots_dir, paste0("umap_by_domain_k", k, ".png"))
  ggsave(f_clust, p_cluster, width = 8, height = 6, dpi = 300)
  ggsave(f_domain, p_domain, width = 8, height = 6, dpi = 300)

  n_clust <- length(unique(Idents(obj)))
  summary[[length(summary)+1]] <- data.frame(k_param=k, n_clusters=n_clust, cluster_plot=f_clust, domain_plot=f_domain)
  log_msg("k=", k, " produced ", n_clust, " clusters")
}

summary_df <- bind_rows(summary)
write.table(summary_df, file.path(out_dir, "walsh_day75_k_sweep_summary.tsv"), sep="\t", row.names = FALSE, quote = FALSE)
log_msg("Sweep complete. Summary -> ", file.path(out_dir, "walsh_day75_k_sweep_summary.tsv"))
