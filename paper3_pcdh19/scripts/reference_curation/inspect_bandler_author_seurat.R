#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

parse_args <- function(args) {
  values <- list()
  index <- 1L
  while (index <= length(args)) {
    key <- args[[index]]
    if (!startsWith(key, "--") || index == length(args)) stop("Arguments must be --name value pairs")
    values[[substring(key, 3L)]] <- args[[index + 1L]]
    index <- index + 2L
  }
  required <- c("object", "ca301", "output-dir")
  missing <- required[!required %in% names(values)]
  if (length(missing)) stop("Missing arguments: ", paste(missing, collapse = ", "))
  values
}

write_tsv <- function(data, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- paste0(path, ".", Sys.getpid(), ".tmp")
  write.table(data, temporary, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  if (!file.rename(temporary, path)) stop("Could not atomically publish ", path)
}

read_nested_gzip_rds <- function(path) {
  payload <- readBin(path, what = "raw", n = file.info(path)$size)
  while (length(payload) >= 2L && identical(as.integer(payload[1:2]), c(31L, 139L))) {
    payload <- memDecompress(payload, type = "gzip")
  }
  connection <- rawConnection(payload, open = "rb")
  on.exit(close(connection))
  readRDS(connection)
}

canonical_barcode <- function(values) {
  found <- grepl("[ACGTN]+-[0-9]+$", values, perl = TRUE)
  result <- rep(NA_character_, length(values))
  result[found] <- sub(".*?([ACGTN]+-[0-9]+)$", "\\1", values[found], perl = TRUE)
  result
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
object_path <- normalizePath(args$object, mustWork = TRUE)
ca301_path <- normalizePath(args$ca301, mustWork = TRUE)
output_dir <- normalizePath(args[["output-dir"]], mustWork = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (!requireNamespace("SeuratObject", quietly = TRUE)) stop("SeuratObject is required")
object <- readRDS(object_path)
if (!inherits(object, "Seurat")) stop("Author artifact is not a Seurat object: ", paste(class(object), collapse = "|"))

metadata <- object[[]]
cell_ids <- colnames(object)
reduction_names <- names(object@reductions)
assay_names <- names(object@assays)

structure <- data.frame(
  field = c("path", "class", "features", "cells", "assays", "reductions", "metadata_columns"),
  value = c(
    object_path, paste(class(object), collapse = "|"), nrow(object), ncol(object),
    paste(assay_names, collapse = "|"), paste(reduction_names, collapse = "|"), ncol(metadata)
  )
)
write_tsv(structure, file.path(output_dir, "author_seurat_structure.tsv"))

column_rows <- lapply(colnames(metadata), function(column) {
  values <- metadata[[column]]
  observed <- as.character(values[!is.na(values)])
  data.frame(
    column = column,
    r_class = paste(class(values), collapse = "|"),
    n_unique_including_na = length(unique(as.character(values))),
    missing = sum(is.na(values)),
    examples = paste(utils::head(unique(observed), 12L), collapse = "|"),
    stringsAsFactors = FALSE
  )
})
write_tsv(do.call(rbind, column_rows), file.path(output_dir, "author_seurat_metadata_columns.tsv"))

count_rows <- list()
for (column in colnames(metadata)) {
  values <- as.character(metadata[[column]])
  values[is.na(values)] <- "<NA>"
  if (length(unique(values)) <= 500L) {
    counts <- sort(table(values), decreasing = TRUE)
    count_rows[[length(count_rows) + 1L]] <- data.frame(
      column = column, value = names(counts), cells = as.integer(counts), stringsAsFactors = FALSE
    )
  }
}
if (length(count_rows)) {
  write_tsv(do.call(rbind, count_rows), file.path(output_dir, "author_seurat_metadata_value_counts.tsv"))
}

reduction_rows <- lapply(reduction_names, function(name) {
  coordinates <- SeuratObject::Embeddings(object[[name]])
  data.frame(reduction = name, cells = nrow(coordinates), dimensions = ncol(coordinates), stringsAsFactors = FALSE)
})
write_tsv(
  if (length(reduction_rows)) do.call(rbind, reduction_rows) else data.frame(reduction=character(), cells=integer(), dimensions=integer()),
  file.path(output_dir, "author_seurat_reductions.tsv")
)

ca301 <- read_nested_gzip_rds(ca301_path)
ca301_ids <- colnames(ca301)
author_canonical <- canonical_barcode(cell_ids)
ca301_canonical <- canonical_barcode(ca301_ids)
overlap <- data.frame(
  comparison = c("exact_cell_id", "canonical_10x_barcode"),
  author_cells = c(length(cell_ids), sum(!is.na(author_canonical))),
  ca301_cells = c(length(ca301_ids), sum(!is.na(ca301_canonical))),
  overlap = c(length(intersect(cell_ids, ca301_ids)), length(intersect(author_canonical, ca301_canonical))),
  stringsAsFactors = FALSE
)
write_tsv(overlap, file.path(output_dir, "author_seurat_ca301_barcode_overlap.tsv"))

likely_columns <- grep(
  "orig.ident|sample|dataset|stage|age|collect|region|tissue|cluster|class|type|ident|state|assignment|lineage",
  colnames(metadata), ignore.case = TRUE, value = TRUE
)
target_pattern <- paste(c(
  "CA301", "MGE", "CGE", "LGE", "progen", "mitotic", "radial", "neuroblast", "precursor",
  "i_Six3", "Gucy1a3", "i_Ebf1", "i_Phlda1", "i_Nr2f2", "i_Nxph1",
  "astro", "oligo", "OPC", "microgl", "macroph", "vascular", "endothelial", "pericyte", "ependymal"
), collapse = "|")
hit_rows <- list()
for (column in likely_columns) {
  values <- as.character(metadata[[column]])
  values[is.na(values)] <- "<NA>"
  counts <- sort(table(values), decreasing = TRUE)
  selected <- grepl(target_pattern, names(counts), ignore.case = TRUE)
  if (any(selected)) {
    hit_rows[[length(hit_rows) + 1L]] <- data.frame(
      column = column,
      value = names(counts)[selected],
      cells = as.integer(counts[selected]),
      stringsAsFactors = FALSE
    )
  }
}
write_tsv(
  if (length(hit_rows)) do.call(rbind, hit_rows) else data.frame(column=character(), value=character(), cells=integer()),
  file.path(output_dir, "author_seurat_target_label_hits.tsv")
)

identifier_examples <- data.frame(
  source = c("author_seurat", "CA301_counts"),
  cells = c(length(cell_ids), length(ca301_ids)),
  examples = c(paste(head(cell_ids, 20L), collapse = "|"), paste(head(ca301_ids, 20L), collapse = "|")),
  stringsAsFactors = FALSE
)
write_tsv(identifier_examples, file.path(output_dir, "author_seurat_identifier_examples.tsv"))
writeLines(capture.output(sessionInfo()), file.path(output_dir, "author_seurat_R_session_info.txt"))
