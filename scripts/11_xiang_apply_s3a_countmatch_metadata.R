#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
  library(ggplot2)
})

project_root <- Sys.getenv(
  "PROJECT_ROOT",
  "/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder"
)
run_label <- Sys.getenv("XIANG_TSNE_RUN_LABEL", "xiang_2017_tsne_preview_3000genes_v1")
out_dir <- file.path(project_root, "results", "xiang_2017_tsne_reproduction", run_label)
rds_in <- file.path(out_dir, "rds", "Xiang_2017_tSNE_reproduction_modern_approximation.rds")
rds_out <- file.path(out_dir, "rds", "Xiang_2017_tSNE_reproduction_modern_approximation_s3a_countmatch.rds")
plot_dir <- file.path(out_dir, "plots")
table_dir <- file.path(out_dir, "tables")

dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

message("[", Sys.time(), "] Loading object: ", rds_in)
obj <- readRDS(rds_in)

s3a_map <- data.table(
  barcode_suffix = as.character(1:8),
  sample_s3a_countmatch = c(
    "early_hCO_rep1", "early_hMGEO_rep1", "early_hMGEO_rep2", "early_hCO_rep2",
    "late_hMGEO_rep1", "late_hCO_rep1", "late_hMGEO_rep2", "late_hCO_rep2"
  ),
  condition_s3a_countmatch = c("hCO", "hMGEO", "hMGEO", "hCO", "hMGEO", "hCO", "hMGEO", "hCO"),
  timepoint_s3a_countmatch = c("d30", "d30", "d30", "d30", "d72", "d72", "d79", "d79"),
  replicate_s3a_countmatch = c("rep1", "rep1", "rep2", "rep2", "rep1", "rep1", "rep2", "rep2")
)

obj$barcode_suffix <- as.character(obj$barcode_suffix)
idx <- match(obj$barcode_suffix, s3a_map$barcode_suffix)
if (any(is.na(idx))) {
  stop("Unmapped barcode_suffix values: ", paste(sort(unique(obj$barcode_suffix[is.na(idx)])), collapse = ", "))
}

obj$sample_s3a_countmatch <- s3a_map$sample_s3a_countmatch[idx]
obj$condition_s3a_countmatch <- s3a_map$condition_s3a_countmatch[idx]
obj$timepoint_s3a_countmatch <- s3a_map$timepoint_s3a_countmatch[idx]
obj$replicate_s3a_countmatch <- s3a_map$replicate_s3a_countmatch[idx]

message("[", Sys.time(), "] Saving corrected object copy: ", rds_out)
saveRDS(obj, rds_out)

meta_cols <- c(
  "sample", "condition", "timepoint", "replicate", "barcode_suffix",
  "sample_s3a_countmatch", "condition_s3a_countmatch",
  "timepoint_s3a_countmatch", "replicate_s3a_countmatch"
)
meta_out <- data.table(barcode = colnames(obj), obj@meta.data[colnames(obj), meta_cols, drop = FALSE])
rm(obj)
gc()

fwrite(meta_out, file.path(table_dir, "Xiang_metadata_with_s3a_countmatch.tsv"), sep = "\t")
fwrite(s3a_map, file.path(table_dir, "Xiang_s3a_countmatch_suffix_map.tsv"), sep = "\t")

suffix_totals <- meta_out[
  ,
  .N,
  by = .(
    barcode_suffix,
    sample_s3a_countmatch,
    condition_s3a_countmatch,
    timepoint_s3a_countmatch,
    replicate_s3a_countmatch
  )
][order(as.integer(barcode_suffix))]
setnames(suffix_totals, "N", "n_cells")
fwrite(suffix_totals, file.path(table_dir, "Xiang_s3a_countmatch_suffix_totals.tsv"), sep = "\t")

condition_totals <- meta_out[, .N, by = condition_s3a_countmatch][order(condition_s3a_countmatch)]
setnames(condition_totals, "N", "n_cells")
fwrite(condition_totals, file.path(table_dir, "Xiang_s3a_countmatch_condition_totals.tsv"), sep = "\t")

sample_totals <- meta_out[
  ,
  .N,
  by = .(
    sample_s3a_countmatch,
    condition_s3a_countmatch,
    timepoint_s3a_countmatch,
    replicate_s3a_countmatch
  )
][order(condition_s3a_countmatch, timepoint_s3a_countmatch, replicate_s3a_countmatch, sample_s3a_countmatch)]
setnames(sample_totals, "N", "n_cells")
fwrite(sample_totals, file.path(table_dir, "Xiang_s3a_countmatch_sample_totals.tsv"), sep = "\t")

expected <- data.table(
  condition_s3a_countmatch = c("hMGEO", "hCO"),
  figure_5a_s3a_expected_cells = c(26949L, 32286L)
)
comparison <- merge(condition_totals, expected, by = "condition_s3a_countmatch", all = TRUE)
comparison[, difference_object_minus_expected := n_cells - figure_5a_s3a_expected_cells]
fwrite(comparison, file.path(table_dir, "Xiang_s3a_countmatch_condition_totals_vs_figure5a.tsv"), sep = "\t")

message("[", Sys.time(), "] Plotting from existing tSNE coordinate table")
tsne <- fread(file.path(table_dir, "Xiang_tSNE_coordinates.tsv"))
tsne[, barcode_suffix := as.character(barcode_suffix)]
plot_df <- merge(
  tsne[, .(barcode, tSNE_1, tSNE_2, barcode_suffix)],
  s3a_map,
  by = "barcode_suffix",
  all.x = TRUE
)
plot_df[, barcode_suffix := factor(barcode_suffix, levels = as.character(1:8))]
plot_df[, sample_s3a_countmatch := factor(sample_s3a_countmatch, levels = s3a_map$sample_s3a_countmatch)]
plot_df[, condition_s3a_countmatch := factor(condition_s3a_countmatch, levels = c("hMGEO", "hCO"))]
plot_df[, timepoint_s3a_countmatch := factor(timepoint_s3a_countmatch, levels = c("d30", "d72", "d79"))]
plot_df[, replicate_s3a_countmatch := factor(replicate_s3a_countmatch, levels = c("rep1", "rep2"))]

pal_suffix <- c(
  `1` = "#0072B2", `2` = "#D55E00", `3` = "#009E73", `4` = "#CC79A7",
  `5` = "#7A3DB8", `6` = "#E69F00", `7` = "#56B4E9", `8` = "#C2185B"
)
pal_sample <- c(
  early_hCO_rep1 = "#0072B2", early_hMGEO_rep1 = "#D55E00",
  early_hMGEO_rep2 = "#009E73", early_hCO_rep2 = "#CC79A7",
  late_hMGEO_rep1 = "#7A3DB8", late_hCO_rep1 = "#E69F00",
  late_hMGEO_rep2 = "#56B4E9", late_hCO_rep2 = "#C2185B"
)
pal_condition <- c(hMGEO = "#009E73", hCO = "#0072B2")
pal_timepoint <- c(d30 = "#D55E00", d72 = "#0072B2", d79 = "#7A3DB8")
pal_replicate <- c(rep1 = "#0072B2", rep2 = "#D55E00")

base_theme <- theme_classic(base_size = 16) +
  theme(
    plot.title = element_text(face = "bold", size = 20),
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 12)
  )

save_tsne_plot <- function(group_col, palette, label) {
  p <- ggplot(plot_df, aes(x = tSNE_1, y = tSNE_2, color = .data[[group_col]])) +
    geom_point(size = 0.42, alpha = 0.9, shape = 16, stroke = 0) +
    scale_color_manual(values = palette, drop = FALSE) +
    coord_equal() +
    base_theme +
    labs(title = paste("Xiang 2017 tSNE PC1-PC5 by", label, "S3A count-match"), x = "XiangTSNE_1", y = "XiangTSNE_2")
  stem <- file.path(plot_dir, paste0("Xiang_tSNE_by_", label, "_s3a_countmatch"))
  ggsave(paste0(stem, ".png"), p, width = 8, height = 6, dpi = 300)
  ggsave(paste0(stem, ".pdf"), p, width = 8, height = 6)
}

save_tsne_plot("barcode_suffix", pal_suffix, "suffix")
save_tsne_plot("sample_s3a_countmatch", pal_sample, "sample")
save_tsne_plot("condition_s3a_countmatch", pal_condition, "condition")
save_tsne_plot("timepoint_s3a_countmatch", pal_timepoint, "timepoint")
save_tsne_plot("replicate_s3a_countmatch", pal_replicate, "replicate")

message("[", Sys.time(), "] Corrected RDS: ", rds_out)
message("Condition totals from new metadata:")
print(condition_totals)
message("Sample totals from new metadata:")
print(sample_totals)
message("Condition totals vs Xiang Figure 5A/S3A:")
print(comparison)
message("[", Sys.time(), "] Done")
