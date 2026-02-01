#!/usr/bin/env Rscript
# Walsh day-75: annotate clusters with paper marker sets without re-running preprocessing.
# Inputs: results/walsh_day75/walsh_day75_final.rds (already post-stress, dims/UMAP done)
# Outputs: annotation tables and UMAP colored by walsh_group.

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
})

log_msg <- function(...) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste0(..., collapse = " ")))
  flush.console()
}

args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args) >= 1) args[1] else Sys.getenv("PROJECT_ROOT", "")
if (project_root == "") stop("Provide PROJECT_ROOT as arg1 or env PROJECT_ROOT")
project_root <- normalizePath(project_root, mustWork = TRUE)

final_rds <- file.path(project_root, "results", "walsh_day75", "walsh_day75_final.rds")
if (!file.exists(final_rds)) stop("Missing final RDS: ", final_rds)
obj <- readRDS(final_rds)
log_msg("Loaded final object:", final_rds)

# Define marker sets from paper
marker_sets <- list(
  `MGE_progenitors (Walsh 6/12)` = c("NKX2-1", "FOXG1"),
  `immature_interneurons (Walsh 0/1/2/4)` = c("LHX6", "GAD1", "ERBB4", "DCX"),
  `SST_interneurons (Walsh 5)` = c("SST", "LHX6", "GAD1"),
  `striatal_interneurons (Walsh 10)` = c("LHX8", "ACHE"),
  `CGE_interneurons (Walsh 8)` = c("VIP", "SCGN", "CCK"),
  `striatal_GP_neurons (Walsh 3/9?)` = c("RXRG", "FOXP2"),
  `astrocytes (paper legend)` = c("AQP4", "S100B"),
  `OPC_like (paper legend)` = c("PDGFRB", "S100B")
)

# Calculate module scores per cell
for (nm in names(marker_sets)) {
  feats <- marker_sets[[nm]]
  present <- feats[feats %in% rownames(obj)]
  if (length(present) == 0) {
    log_msg("No markers present for set", nm)
    obj[[paste0("score_", nm)]] <- 0
    next
  }
  obj <- AddModuleScore(obj, features = list(present), name = paste0("score_", nm), verbose = FALSE)
}

# Collect score columns (AddModuleScore adds columns starting with provided name)
score_cols <- grep("^score_", colnames(obj@meta.data), value = TRUE)
if (length(score_cols) == 0) stop("No module score columns found")

# Per-cell best group
score_mat <- obj@meta.data[, score_cols, drop = FALSE]
colnames(score_mat) <- sub("1$", "", colnames(score_mat)) # strip trailing 1 if present
max_idx <- apply(score_mat, 1, function(x) which.max(x))
best_group <- colnames(score_mat)[max_idx]
names(best_group) <- colnames(obj)
obj$walsh_group <- best_group

# Cluster-level majority vote
cluster_df <- obj@meta.data %>%
  mutate(cluster = as.character(seurat_clusters)) %>%
  group_by(cluster, walsh_group) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(cluster) %>%
  slice_max(n, with_ties = FALSE) %>%
  ungroup()

# Domain composition per walsh_group
domain_comp <- obj@meta.data %>%
  group_by(walsh_group, domain) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(walsh_group) %>%
  mutate(frac = n / sum(n)) %>%
  ungroup()

# Output crosswalk
crosswalk <- cluster_df %>%
  left_join(obj@meta.data %>% mutate(cluster = as.character(seurat_clusters)) %>%
              group_by(cluster) %>% summarise(cells = n(), .groups = "drop"), by = "cluster") %>%
  rename(walsh_group = walsh_group, our_cluster = cluster, cluster_cells = cells)

crosswalk_file <- file.path(project_root, "results", "walsh_day75", "walsh_cluster_annotation.tsv")
write.table(crosswalk, crosswalk_file, sep = "\t", quote = FALSE, row.names = FALSE)
log_msg("Cluster annotation saved ->", crosswalk_file)

# Domain composition table
domain_file <- file.path(project_root, "results", "walsh_day75", "walsh_domain_composition_by_group.tsv")
write.table(domain_comp, domain_file, sep = "\t", quote = FALSE, row.names = FALSE)
log_msg("Domain composition saved ->", domain_file)

# UMAP plot colored by walsh_group
plots_dir <- file.path(project_root, "results", "walsh_day75", "annotation_plots")
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
reduction_use <- if ("umap_sel" %in% names(obj@reductions)) "umap_sel" else "umap"
p_group <- DimPlot(obj, reduction = reduction_use, group.by = "walsh_group", label = FALSE, pt.size = 1.6, shuffle = FALSE, order = TRUE, raster = TRUE) +
  ggtitle("UMAP by Walsh-style group")
ggsave(file.path(plots_dir, "umap_by_walsh_group.pdf"), p_group, width = 8, height = 6)
ggsave(file.path(plots_dir, "umap_by_walsh_group.png"), p_group, width = 8, height = 6, dpi = 300)

# Save updated object with walsh_group metadata
saveRDS(obj, file.path(project_root, "results", "walsh_day75", "walsh_day75_final_annotated.rds"))
log_msg("Annotated object saved ->", file.path(project_root, "results", "walsh_day75", "walsh_day75_final_annotated.rds"))
