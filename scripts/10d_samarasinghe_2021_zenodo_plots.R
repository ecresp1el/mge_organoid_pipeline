#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})

get_arg <- function(args, key, default = "") {
  hit <- which(args == key)
  if (length(hit) > 0 && hit[1] < length(args)) return(args[hit[1] + 1])
  prefix <- paste0(key, "=")
  hit <- grep(paste0("^", prefix), args)
  if (length(hit) > 0) return(sub(prefix, "", args[hit[1]]))
  default
}

write_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
}

save_plot <- function(plot, stem, width = 8, height = 6, dpi = 300) {
  dir.create(dirname(stem), recursive = TRUE, showWarnings = FALSE)
  ggsave(paste0(stem, ".png"), plot, width = width, height = height, dpi = dpi)
  ggsave(paste0(stem, ".pdf"), plot, width = width, height = height)
}

args <- commandArgs(trailingOnly = TRUE)
project_root <- get_arg(args, "--project-root", Sys.getenv("PROJECT_ROOT", unset = ""))
if (!nzchar(project_root)) {
  project_root <- "/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder"
}
project_root <- sub("/+$", "", project_root)

rdata <- get_arg(
  args,
  "--rdata",
  file.path(project_root, "data/raw/samarasinghe_2021_zenodo/Samarasinghe_2021_seurat_object.RData")
)
out_dir <- get_arg(
  args,
  "--out-dir",
  file.path(project_root, "results/samarasinghe_2021_zenodo_processed_object")
)
plot_dir <- file.path(out_dir, "plots")
table_dir <- file.path(out_dir, "tables")
stopifnot(file.exists(rdata))

message("Loading Zenodo RData: ", rdata)
loaded <- load(rdata)
seurat_names <- loaded[vapply(loaded, function(nm) inherits(get(nm), "Seurat"), logical(1))]
if (!length(seurat_names)) stop("No Seurat object found in RData")
obj <- get(seurat_names[[1]])

saveRDS(obj, file.path(out_dir, "samarasinghe_2021_zenodo_seurat.rds"))

meta <- obj[[]]
qc_upper <- mean(meta$nFeature_RNA) + 3 * sd(meta$nFeature_RNA)
qc_pass <- meta$nFeature_RNA > 500 & meta$nFeature_RNA < qc_upper & meta$percent.mt < 10
qc_summary <- data.frame(
  metric = c("n_cells", "nFeature_mean", "nFeature_sd", "nFeature_upper_mean_plus_3sd", "qc_pass_cells", "qc_fail_cells"),
  value = c(ncol(obj), mean(meta$nFeature_RNA), sd(meta$nFeature_RNA), qc_upper, sum(qc_pass), sum(!qc_pass)),
  stringsAsFactors = FALSE
)
write_tsv(qc_summary, file.path(table_dir, "zenodo_qc_summary_from_metadata.tsv"))

sample_counts <- as.data.frame(table(meta$orig.ident, meta$Time, meta$Genotype), stringsAsFactors = FALSE)
colnames(sample_counts) <- c("orig.ident", "Time", "Genotype", "n_cells")
sample_counts <- sample_counts[sample_counts$n_cells > 0, , drop = FALSE]
write_tsv(sample_counts, file.path(table_dir, "zenodo_sample_counts.tsv"))

cluster_counts <- as.data.frame(table(meta$seurat_clusters, meta$new.cluster.ids), stringsAsFactors = FALSE)
colnames(cluster_counts) <- c("seurat_clusters", "new.cluster.ids", "n_cells")
cluster_counts <- cluster_counts[cluster_counts$n_cells > 0, , drop = FALSE]
write_tsv(cluster_counts, file.path(table_dir, "zenodo_cluster_celltype_counts.tsv"))

save_plot(
  DimPlot(obj, reduction = "umap", group.by = "new.cluster.ids", label = TRUE, repel = TRUE) +
    ggtitle("Samarasinghe 2021 Zenodo UMAP by cell type"),
  file.path(plot_dir, "zenodo_umap_by_celltype")
)
save_plot(
  DimPlot(obj, reduction = "umap", group.by = "seurat_clusters", label = TRUE, repel = TRUE) +
    ggtitle("Samarasinghe 2021 Zenodo UMAP by Seurat cluster"),
  file.path(plot_dir, "zenodo_umap_by_cluster")
)
save_plot(
  DimPlot(obj, reduction = "umap", group.by = "Genotype") +
    ggtitle("Samarasinghe 2021 Zenodo UMAP by genotype"),
  file.path(plot_dir, "zenodo_umap_by_genotype")
)
save_plot(
  DimPlot(obj, reduction = "umap", group.by = "Time") +
    ggtitle("Samarasinghe 2021 Zenodo UMAP by time"),
  file.path(plot_dir, "zenodo_umap_by_time")
)
save_plot(
  DimPlot(obj, reduction = "umap", group.by = "orig.ident") +
    ggtitle("Samarasinghe 2021 Zenodo UMAP by original sample"),
  file.path(plot_dir, "zenodo_umap_by_orig_ident")
)

ctrl_cells <- rownames(meta)[meta$Genotype == "Ctrl"]
if (length(ctrl_cells) > 0) {
  save_plot(
    DimPlot(
      obj,
      reduction = "umap",
      cells = ctrl_cells,
      group.by = "new.cluster.ids",
      label = TRUE,
      repel = TRUE
    ) + ggtitle("Samarasinghe 2021 Zenodo Ctrl cells on full UMAP"),
    file.path(plot_dir, "zenodo_umap_ctrl_only_by_celltype")
  )
}

message("Done. Outputs in: ", out_dir)
