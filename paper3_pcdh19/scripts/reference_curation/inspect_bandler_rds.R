#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

parse_args <- function(args) {
  values <- list()
  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, "--") || i == length(args)) {
      stop("Arguments must be supplied as --name value pairs")
    }
    values[[substring(key, 3L)]] <- args[[i + 1L]]
    i <- i + 2L
  }
  required <- c("input", "study-dir", "summary")
  missing <- required[!required %in% names(values)]
  if (length(missing) > 0L) stop("Missing arguments: ", paste(missing, collapse = ", "))
  values
}

write_tsv <- function(data, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.table(data, file = path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
input_path <- normalizePath(args$input, mustWork = TRUE)
study_dir <- normalizePath(args[["study-dir"]], mustWork = FALSE)
summary_path <- args$summary
audit_dir <- file.path(study_dir, "audit")
metadata_dir <- file.path(study_dir, "metadata")
dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(metadata_dir, recursive = TRUE, showWarnings = FALSE)

connection <- if (grepl("\\.gz$", input_path, ignore.case = TRUE)) gzfile(input_path, open = "rb") else input_path
on.exit(if (inherits(connection, "connection")) close(connection), add = TRUE)
obj <- readRDS(connection)

object_class <- paste(class(obj), collapse = "|")
object_dim <- dim(obj)
if (is.null(object_dim)) object_dim <- integer(0)
dimension_text <- if (length(object_dim) > 0L) paste(object_dim, collapse = " x ") else "no_dim_attribute"

structure_lines <- capture.output({
  cat("path=", input_path, "\n", sep = "")
  cat("class=", object_class, "\n", sep = "")
  cat("dimensions=", dimension_text, "\n", sep = "")
  str(obj, max.level = 2L, give.attr = TRUE)
})
writeLines(structure_lines, file.path(audit_dir, "rds_structure.txt"))
writeLines(capture.output(sessionInfo()), file.path(audit_dir, "r_session_info.txt"))

metadata_columns <- character(0)
reduction_names <- character(0)
actual_type <- object_class

if (inherits(obj, "Seurat")) {
  if (!requireNamespace("SeuratObject", quietly = TRUE)) {
    stop("Object is Seurat but SeuratObject is unavailable")
  }
  metadata <- obj[[]]
  metadata_columns <- colnames(metadata)
  reduction_names <- names(obj@reductions)
  actual_type <- paste0("Seurat:", object_class)
} else if (inherits(obj, "SingleCellExperiment")) {
  if (!requireNamespace("SummarizedExperiment", quietly = TRUE) || !requireNamespace("SingleCellExperiment", quietly = TRUE)) {
    stop("Object is SingleCellExperiment but required packages are unavailable")
  }
  metadata_columns <- colnames(as.data.frame(SummarizedExperiment::colData(obj)))
  reduction_names <- SingleCellExperiment::reducedDimNames(obj)
  actual_type <- paste0("SingleCellExperiment:", object_class)
} else if (is.matrix(obj) || inherits(obj, "Matrix")) {
  actual_type <- paste0("expression_matrix:", object_class)
} else if (is.data.frame(obj)) {
  metadata_columns <- colnames(obj)
  actual_type <- paste0("data.frame:", object_class)
}

write_tsv(
  data.frame(column = metadata_columns, stringsAsFactors = FALSE),
  file.path(metadata_dir, "metadata_columns.tsv")
)
write_tsv(
  data.frame(reduction = reduction_names, stringsAsFactors = FALSE),
  file.path(audit_dir, "reduction_inventory.tsv")
)

annotation_pattern <- "cluster|class|subclass|cell.?type|type|subtype|identity|annotation|taxonomy|lineage|state"
annotation_columns <- metadata_columns[grepl(annotation_pattern, metadata_columns, ignore.case = TRUE)]
embedding_hits <- reduction_names[grepl("umap|tsne", reduction_names, ignore.case = TRUE)]

summary <- data.frame(
  paper = "Bandler2022",
  P0_file = basename(input_path),
  actual_object_type = actual_type,
  dimensions = dimension_text,
  author_embedding_present = if (length(embedding_hits) > 0L) paste0("yes:", paste(embedding_hits, collapse = "|")) else "no_saved_umap_or_tsne_detected",
  annotation_columns_present = if (length(annotation_columns) > 0L) paste0("yes:", paste(annotation_columns, collapse = "|")) else "no_cell_level_annotation_columns_detected",
  MGE_selectable = "yes_at_sample_level_CA301_is_WT_MGE",
  age_selectable = "yes_at_sample_level_CA301_is_E15",
  cell_level_labels_immediately_usable = if (length(annotation_columns) > 0L) "yes_pending_barcode_validation" else "no",
  next_minimal_action = if (length(annotation_columns) > 0L) {
    "Validate author label semantics and CA301 barcode identity before plotting saved embeddings."
  } else {
    "Find the smallest published barcode-to-embryonic-annotation mapping for CA301; do not reconstruct the full study yet."
  },
  stringsAsFactors = FALSE
)
write_tsv(summary, summary_path)
