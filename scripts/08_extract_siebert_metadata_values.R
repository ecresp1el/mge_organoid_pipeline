#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
})

get_env <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (nzchar(value)) value else default
}

write_tsv <- function(df, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.table(df, path, sep = "\t", row.names = FALSE, quote = FALSE, na = "NA")
}

count_values <- function(values) {
  values <- as.character(values)
  values[is.na(values) | !nzchar(values)] <- "NA"
  tab <- sort(table(values), decreasing = TRUE)
  data.frame(value = names(tab), n_cells = as.integer(tab), stringsAsFactors = FALSE)
}

project_root <- get_env("PROJECT_ROOT", "/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder")
seurat_rds <- get_env(
  "SIEBERT_SEURAT_RDS",
  file.path(project_root, "results/siebert_2026/siebert_2026_seurat.rds")
)
out_dir <- get_env(
  "SIEBERT_METADATA_OUT_DIR",
  file.path(project_root, "results/siebert_2026/metadata_value_inventory")
)

message("Reading Siebert Seurat object: ", seurat_rds)
obj <- readRDS(seurat_rds)
md <- obj@meta.data

target_columns <- c(
  "orig.ident",
  "sample",
  "cellLine",
  "donor_batch",
  "predictions",
  "predictions_class"
)
present_columns <- intersect(target_columns, colnames(md))

column_summary <- data.frame(
  column_name = present_columns,
  n_non_na = vapply(md[present_columns], function(x) sum(!is.na(x)), integer(1)),
  n_unique_non_na = vapply(md[present_columns], function(x) length(unique(as.character(x[!is.na(x)]))), integer(1)),
  stringsAsFactors = FALSE
)
write_tsv(column_summary, file.path(out_dir, "siebert_metadata_column_summary.tsv"))

value_counts <- do.call(
  rbind,
  lapply(present_columns, function(column_name) {
    counts <- count_values(md[[column_name]])
    counts$column_name <- column_name
    counts[, c("column_name", "value", "n_cells")]
  })
)
write_tsv(value_counts, file.path(out_dir, "siebert_metadata_value_counts.tsv"))

for (column_name in present_columns) {
  counts <- count_values(md[[column_name]])
  write_tsv(counts, file.path(out_dir, paste0("siebert_", column_name, "_counts.tsv")))
}

if (all(c("sample", "cellLine") %in% colnames(md))) {
  sample_cellline <- as.data.frame(table(md$sample, md$cellLine), stringsAsFactors = FALSE)
  colnames(sample_cellline) <- c("sample", "cellLine", "n_cells")
  sample_cellline <- sample_cellline[sample_cellline$n_cells > 0, , drop = FALSE]
  sample_cellline <- sample_cellline[order(sample_cellline$sample, sample_cellline$cellLine), , drop = FALSE]
  write_tsv(sample_cellline, file.path(out_dir, "siebert_sample_by_cellLine_counts.tsv"))
}

if (all(c("sample", "donor_batch") %in% colnames(md))) {
  sample_donor_batch <- as.data.frame(table(md$sample, md$donor_batch), stringsAsFactors = FALSE)
  colnames(sample_donor_batch) <- c("sample", "donor_batch", "n_cells")
  sample_donor_batch <- sample_donor_batch[sample_donor_batch$n_cells > 0, , drop = FALSE]
  sample_donor_batch <- sample_donor_batch[order(sample_donor_batch$sample, sample_donor_batch$donor_batch), , drop = FALSE]
  write_tsv(sample_donor_batch, file.path(out_dir, "siebert_sample_by_donor_batch_counts.tsv"))
}

message("Wrote Siebert metadata value inventory to: ", out_dir)
