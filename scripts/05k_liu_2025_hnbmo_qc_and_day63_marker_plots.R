#!/usr/bin/env Rscript

# Supplemental-style plots for Liu 2025 hnbMO scRNA-seq:
#   - day-63 UMAP expression for CHAT and RBFOX3
#   - violin plots for nCount_RNA and nFeature_RNA by Figure 1C sample

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})

get_env <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (nzchar(value)) value else default
}

project_root <- get_env("PROJECT_ROOT", "/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder")
seurat_rds <- get_env(
  "LIU2025_HNBMO_SEURAT_RDS",
  file.path(project_root, "results/liu_2025_hnbmo_integrated/liu_2025_hnbmo_healthy_integrated_seurat.rds")
)
out_dir <- get_env(
  "LIU2025_HNBMO_SUPP_PLOT_DIR",
  file.path(project_root, "results/liu_2025_hnbmo_integrated/plots")
)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

message("Reading Seurat object: ", seurat_rds)
obj <- readRDS(seurat_rds)

if (!("figure1c_sample" %in% colnames(obj@meta.data))) {
  obj$figure1c_sample <- factor(
    paste(obj$cell_line, obj$day, sep = " "),
    levels = c("H9 D36", "H9 D63", "IMR90-4 D63")
  )
}

DefaultAssay(obj) <- "RNA"

if (!("umap" %in% names(obj@reductions))) {
  stop("UMAP reduction not found in object")
}

marker_genes <- c("CHAT", "RBFOX3")
present_markers <- intersect(marker_genes, rownames(obj))
missing_markers <- setdiff(marker_genes, present_markers)
if (length(missing_markers) > 0L) {
  warning("Missing marker(s): ", paste(missing_markers, collapse = ", "))
}
if (length(present_markers) == 0L) {
  stop("None of the requested markers were found: ", paste(marker_genes, collapse = ", "))
}

day63 <- subset(obj, subset = day == "D63")
message("Day 63 cells: ", ncol(day63))

p_day63_markers <- FeaturePlot(
  day63,
  features = present_markers,
  reduction = "umap",
  ncol = length(present_markers),
  order = TRUE
)

ggsave(
  file.path(out_dir, "supp_like_day63_umap_CHAT_RBFOX3.png"),
  p_day63_markers,
  width = 10,
  height = 5,
  dpi = 300
)
ggsave(
  file.path(out_dir, "supp_like_day63_umap_CHAT_RBFOX3.pdf"),
  p_day63_markers,
  width = 10,
  height = 5
)

qc_df <- obj@meta.data
qc_df$figure1c_sample <- factor(
  qc_df$figure1c_sample,
  levels = c("H9 D36", "H9 D63", "IMR90-4 D63")
)

metric_labels <- c(nCount_RNA = "nCount_RNA", nFeature_RNA = "nFeature_RNA")
qc_long <- do.call(
  rbind,
  lapply(names(metric_labels), function(metric) {
    data.frame(
      figure1c_sample = qc_df$figure1c_sample,
      metric = metric_labels[[metric]],
      value = as.numeric(qc_df[[metric]]),
      stringsAsFactors = FALSE
    )
  })
)

p_violin <- ggplot(qc_long, aes(x = figure1c_sample, y = value, fill = figure1c_sample)) +
  geom_violin(scale = "width", trim = TRUE, color = "gray25", linewidth = 0.25) +
  geom_boxplot(width = 0.12, outlier.shape = NA, color = "gray20", fill = "white", alpha = 0.75) +
  facet_wrap(~metric, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = c("H9 D36" = "#F8766D", "H9 D63" = "#00BA38", "IMR90-4 D63" = "#619CFF")) +
  labs(x = NULL, y = NULL, title = "scRNA-seq metrics by hnbMO sample") +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 25, hjust = 1),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold")
  )

ggsave(
  file.path(out_dir, "supp_like_violin_nCount_nFeature_by_sample.png"),
  p_violin,
  width = 9,
  height = 5,
  dpi = 300
)
ggsave(
  file.path(out_dir, "supp_like_violin_nCount_nFeature_by_sample.pdf"),
  p_violin,
  width = 9,
  height = 5
)

message("Wrote plots to: ", out_dir)
