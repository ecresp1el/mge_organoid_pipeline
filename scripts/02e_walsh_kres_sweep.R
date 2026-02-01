#!/usr/bin/env Rscript
# Walsh day-75: 2D sweep over k.param and resolution (0.8, 0.6) starting from post-stress PCA checkpoint.

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
plots_dir <- file.path(project_root, "results", "walsh_day75", "kres_sweep_plots")
out_dir <- file.path(project_root, "results", "walsh_day75")
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

pca_ckpt <- file.path(ckpt_dir, "walsh_pca_poststress.rds")
if (!file.exists(pca_ckpt)) stop("Missing PCA checkpoint: ", pca_ckpt)
obj_base <- readRDS(pca_ckpt)
log_msg("Loaded post-stress PCA checkpoint:", pca_ckpt)

res_vals <- c(0.8, 0.6)
k_vals <- c(20, 30, 40, 50)
seeds <- list(neighbors=1004L, clusters=1005L, umap=1006L)
summary <- list()
found14 <- FALSE
found_pair <- NULL

for (res in res_vals) {
  for (k in k_vals) {
    log_msg("Running res=", res, " k=", k)
    obj <- obj_base
    set.seed(seeds$neighbors)
    obj <- FindNeighbors(obj, dims = 1:20, k.param = k, algorithm = 1, verbose = FALSE)
    set.seed(seeds$clusters)
    obj <- FindClusters(obj, resolution = res, algorithm = 1, verbose = FALSE)
    set.seed(seeds$umap)
    red_name <- paste0("umap_res", gsub("\\.", "_", res), "_k", k)
    obj <- suppressWarnings(RunUMAP(obj, dims = 1:20, reduction.name = red_name, reduction.key = paste0("UMAP", gsub("\\.", "_", res), "k", k, "_"),
                   seed.use = seeds$umap, n.neighbors = 30, min.dist = 0.3, spread = 1,
                   metric = "cosine", umap.method = "uwot", return.model = FALSE, verbose = FALSE))

    p_cluster <- DimPlot(obj, reduction = red_name, group.by = "seurat_clusters", label = TRUE, repel = TRUE, pt.size = 1.2, shuffle = FALSE, order = TRUE, raster = TRUE) +
      ggtitle(paste0("UMAP res=", res, " k=", k, " by cluster"))
    p_domain <- DimPlot(obj, reduction = red_name, group.by = "domain", pt.size = 1.2, shuffle = FALSE, order = TRUE, raster = TRUE) +
      ggtitle(paste0("UMAP res=", res, " k=", k, " by domain"))

    f_clust <- file.path(plots_dir, paste0("umap_by_cluster_res", res, "_k", k, ".png"))
    f_domain <- file.path(plots_dir, paste0("umap_by_domain_res", res, "_k", k, ".png"))
    ggsave(f_clust, p_cluster, width = 8, height = 6, dpi = 300)
    ggsave(f_domain, p_domain, width = 8, height = 6, dpi = 300)

    n_clust <- length(unique(Idents(obj)))
    summary[[length(summary)+1]] <- data.frame(resolution=res, k_param=k, n_clusters=n_clust, cluster_plot=f_clust, domain_plot=f_domain)
    log_msg("res=", res, " k=", k, " produced ", n_clust, " clusters")

    if (n_clust == 14) {
      found14 <- TRUE
      found_pair <- c(resolution=res, k_param=k)
      break
    }
  }
  if (found14) break
}

summary_df <- bind_rows(summary)
write.table(summary_df, file.path(out_dir, "walsh_day75_kres_sweep_summary.tsv"), sep="\t", row.names = FALSE, quote = FALSE)
if (found14) {
  log_msg("Found ~14 clusters at res=", found_pair["resolution"], " k=", found_pair["k_param"])
} else {
  log_msg("No combination reached 14 clusters; see summary for closest counts")
}
