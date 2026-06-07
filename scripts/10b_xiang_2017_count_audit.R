#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Matrix)
  library(data.table)
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

args <- commandArgs(trailingOnly = TRUE)
project_root <- get_arg(args, "--project-root", Sys.getenv("PROJECT_ROOT", unset = ""))
if (!nzchar(project_root)) {
  project_root <- "/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder"
}
project_root <- sub("/+$", "", project_root)

out_dir <- get_arg(
  args,
  "--out-dir",
  file.path(project_root, "results/xiang_2017_tsne_reproduction/audit_no_tsne")
)

suppl_dir <- file.path(project_root, "data/raw/xiang_2018_geo_files/suppl")
matrix_path <- file.path(suppl_dir, "GSE98201_matrix.mtx.gz")
genes_path <- file.path(suppl_dir, "GSE98201_genes.tsv.gz")
barcodes_path <- file.path(suppl_dir, "GSE98201_barcodes.tsv.gz")
stopifnot(file.exists(matrix_path), file.exists(genes_path), file.exists(barcodes_path))

message("Reading Xiang matrix for count audit only: ", matrix_path)
counts <- readMM(matrix_path)
genes <- fread(cmd = paste("zcat", shQuote(genes_path)), header = FALSE)
barcodes <- fread(cmd = paste("zcat", shQuote(barcodes_path)), header = FALSE, col.names = "barcode")
if (ncol(genes) == 1) {
  setnames(genes, "ensembl_id")
} else {
  setnames(genes, 1:2, c("ensembl_id", "gene_symbol"))
}
stopifnot(nrow(counts) == nrow(genes), ncol(counts) == nrow(barcodes))

barcode_suffix <- sub("^.*-", "", barcodes$barcode)
suffix_levels <- as.character(1:8)

# GEO sample order from GSE97882/GSE98201 records used by the Xiang tSNE
# reproduction script. This maps suffixes only; cell counts are audited below.
suffix_map <- data.frame(
  barcode_suffix = suffix_levels,
  sample = c(
    "hMGEO_d30_rep1",
    "hCO_d30_rep1",
    "hMGEO_d72_rep1",
    "hCO_d72_rep1",
    "hMGEO_d30_rep2",
    "hCO_d30_rep2",
    "hMGEO_d79_rep2",
    "hCO_d79_rep2"
  ),
  condition = c("hMGEO", "hCO", "hMGEO", "hCO", "hMGEO", "hCO", "hMGEO", "hCO"),
  timepoint = c("d30", "d30", "d72", "d72", "d30", "d30", "d79", "d79"),
  replicate = c("rep1", "rep1", "rep1", "rep1", "rep2", "rep2", "rep2", "rep2"),
  geo_accession = c(
    "GSM2589129",
    "GSM2589130",
    "GSM2589131",
    "GSM2589132",
    "GSM2684867",
    "GSM2684868",
    "GSM2684869",
    "GSM2684870"
  ),
  stringsAsFactors = FALSE
)

metadata <- data.frame(
  barcode = barcodes$barcode,
  barcode_suffix = barcode_suffix,
  stringsAsFactors = FALSE
)
metadata <- merge(metadata, suffix_map, by = "barcode_suffix", all.x = TRUE, sort = FALSE)
metadata <- metadata[match(barcodes$barcode, metadata$barcode), , drop = FALSE]

feature_names <- if ("gene_symbol" %in% names(genes)) {
  ifelse(!is.na(genes$gene_symbol) & nzchar(genes$gene_symbol), genes$gene_symbol, genes$ensembl_id)
} else {
  genes$ensembl_id
}
feature_names_unique <- make.unique(feature_names)

stage_rows <- list()
suffix_rows <- list()

suffix_count_df <- function(stage, keep_cells) {
  tab <- as.data.frame(table(factor(barcode_suffix[keep_cells], levels = suffix_levels)), stringsAsFactors = FALSE)
  colnames(tab) <- c("barcode_suffix", "n_cells")
  tab$stage <- stage
  tab <- merge(tab, suffix_map, by = "barcode_suffix", all.x = TRUE, sort = FALSE)
  tab <- tab[order(as.integer(tab$barcode_suffix)), ]
  tab[, c("stage", "barcode_suffix", "sample", "condition", "timepoint", "replicate", "geo_accession", "n_cells")]
}

add_stage <- function(stage_order, stage, keep_cells, n_genes, notes = "") {
  stage_rows[[length(stage_rows) + 1L]] <<- data.frame(
    stage_order = stage_order,
    stage = stage,
    n_cells = sum(keep_cells),
    n_genes = n_genes,
    notes = notes,
    stringsAsFactors = FALSE
  )
  suffix_rows[[length(suffix_rows) + 1L]] <<- suffix_count_df(stage, keep_cells)
}

all_cells <- rep(TRUE, ncol(counts))
add_stage(1, "raw_matrix_loaded", all_cells, nrow(counts), "MatrixMarket raw dimensions")
add_stage(2, "after_barcode_suffix_parsing", all_cells, nrow(counts), paste0("unmapped_suffixes=", sum(!barcode_suffix %in% suffix_levels)))
add_stage(3, "after_sample_metadata_assignment", all_cells, nrow(counts), paste0("missing_sample_metadata=", sum(is.na(metadata$sample))))

# Gene-symbol conversion is chronological before filtering in the script. It
# does not drop cells or genes because make.unique preserves one row per input.
add_stage(
  7,
  "after_gene_symbol_conversion",
  all_cells,
  length(feature_names_unique),
  paste0("duplicated_raw_symbols=", sum(duplicated(feature_names)), "; rows_preserved=true")
)

# Paper-style CreateSeuratObject filters decomposed into gene and cell effects:
# min.cells = 4, min.features = 201.
gene_keep <- Matrix::rowSums(counts > 0) >= 4
genes_after_gene_filter <- sum(gene_keep)
add_stage(
  4,
  "after_gene_filter_min_cells_ge_4",
  all_cells,
  genes_after_gene_filter,
  paste0("genes_removed=", sum(!gene_keep))
)

nfeature_after_gene_filter <- Matrix::colSums(counts[gene_keep, , drop = FALSE] > 0)
cell_keep <- nfeature_after_gene_filter >= 201
add_stage(
  5,
  "after_cell_filter_min_features_ge_201",
  cell_keep,
  genes_after_gene_filter,
  paste0("cells_removed=", sum(!cell_keep))
)
add_stage(
  6,
  "after_mito_count_doublet_filtering",
  cell_keep,
  genes_after_gene_filter,
  "none in paper-style Xiang reproduction after min.features filter"
)

variable_gene_file <- file.path(
  project_root,
  "results/xiang_2017_tsne_reproduction/xiang_2017_tsne_reproduction_v1/tables/Xiang_variable_genes_dispersion.tsv"
)
variable_genes <- if (file.exists(variable_gene_file)) {
  max(0L, length(readLines(variable_gene_file)) - 1L)
} else {
  NA_integer_
}
add_stage(
  8,
  "after_integration_or_feature_intersection",
  cell_keep,
  if (is.na(variable_genes)) genes_after_gene_filter else variable_genes,
  "no integration or cross-dataset feature intersection; n_genes is dispersion-positive variable-gene table size when available"
)

stage_counts <- do.call(rbind, stage_rows)
suffix_counts <- do.call(rbind, suffix_rows)
stage_counts <- stage_counts[order(stage_counts$stage_order), ]

user_expected <- data.frame(
  barcode_suffix = suffix_levels,
  user_expected_sample = c(
    "hMGEO_d30_rep1",
    "hCO_d30_rep1",
    "hMGEO_d72_rep1",
    "hCO_d72_rep1",
    "hMGEO_d30_rep2",
    "hCO_d30_rep2",
    "hMGEO_d79_rep2",
    "hCO_d79_rep2"
  ),
  user_expected_cells = c(2969, 6480, 5438, 6355, 9722, 10258, 8820, 9193),
  stringsAsFactors = FALSE
)
raw_suffix <- suffix_count_df("raw_matrix_loaded", all_cells)
comparison <- merge(
  user_expected,
  raw_suffix[, c("barcode_suffix", "sample", "geo_accession", "n_cells")],
  by = "barcode_suffix",
  all.x = TRUE,
  sort = FALSE
)
comparison <- comparison[order(as.integer(comparison$barcode_suffix)), ]
comparison$observed_minus_user_expected <- comparison$n_cells - comparison$user_expected_cells

write_tsv(stage_counts, file.path(out_dir, "xiang_audit_stage_counts.tsv"))
write_tsv(suffix_counts, file.path(out_dir, "xiang_audit_suffix_counts_by_stage.tsv"))
write_tsv(comparison, file.path(out_dir, "xiang_audit_user_expected_vs_observed_suffix_counts.tsv"))

message("Audit complete: ", out_dir)
print(stage_counts, row.names = FALSE)
