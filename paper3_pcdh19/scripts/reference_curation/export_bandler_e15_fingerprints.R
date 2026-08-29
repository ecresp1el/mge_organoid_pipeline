#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

if (!requireNamespace("Matrix", quietly = TRUE)) {
  stop("Matrix is required")
}

parse_args <- function(args) {
  values <- list()
  i <- 1L
  while (i <= length(args)) {
    if (!startsWith(args[[i]], "--") || i == length(args)) {
      stop("Arguments must be supplied as --name value pairs")
    }
    values[[substring(args[[i]], 3L)]] <- args[[i + 1L]]
    i <- i + 2L
  }
  required <- c("source-root", "genes", "output", "inventory")
  missing <- required[!required %in% names(values)]
  if (length(missing) > 0L) stop("Missing arguments: ", paste(missing, collapse = ", "))
  values
}

read_nested_rds <- function(path) {
  payload <- readBin(path, what = "raw", n = file.info(path)$size)
  gzip_layers <- 0L
  while (length(payload) >= 2L && identical(as.integer(payload[1:2]), c(31L, 139L))) {
    payload <- memDecompress(payload, type = "gzip")
    gzip_layers <- gzip_layers + 1L
  }
  connection <- rawConnection(payload, open = "rb")
  on.exit(close(connection), add = TRUE)
  list(object = readRDS(connection), gzip_layers = gzip_layers)
}

write_tsv <- function(data, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- tempfile(pattern = paste0(".", basename(path), "."), tmpdir = dirname(path))
  on.exit(unlink(temporary), add = TRUE)
  write.table(data, temporary, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  if (!file.rename(temporary, path)) stop("Could not atomically publish ", path)
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
source_root <- normalizePath(args[["source-root"]], mustWork = TRUE)
genes <- strsplit(args$genes, ",", fixed = TRUE)[[1L]]
samples <- data.frame(
  sample = c("CA301", "CA302", "CA303"),
  gsm = c("GSM5684876", "GSM5684875", "GSM5684874"),
  region = c("MGE", "CGE", "LGE"),
  filename = c(
    "GSM5684876_CA301_filtered_RNA_counts.RDS.gz",
    "GSM5684875_CA302_filtered_RNA_counts.RDS.gz",
    "GSM5684874_CA303_filtered_RNA_counts.RDS.gz"
  ),
  stringsAsFactors = FALSE
)

fingerprints <- list()
inventory <- list()
for (i in seq_len(nrow(samples))) {
  sample_row <- samples[i, ]
  path <- file.path(source_root, sample_row$filename)
  if (!file.exists(path)) stop("Missing deposited matrix: ", path)
  loaded <- read_nested_rds(path)
  matrix <- loaded$object
  if (!(is.matrix(matrix) || inherits(matrix, "Matrix"))) {
    stop(sample_row$sample, " did not resolve to an expression matrix: ", paste(class(matrix), collapse = "|"))
  }
  if (is.null(rownames(matrix)) || is.null(colnames(matrix))) stop("Matrix identifiers are required for ", sample_row$sample)
  missing_genes <- setdiff(genes, rownames(matrix))
  if (length(missing_genes) > 0L) stop("Fingerprint genes absent from ", sample_row$sample, ": ", paste(missing_genes, collapse = ","))
  library_size <- Matrix::colSums(matrix)
  if (any(library_size <= 0)) stop("Zero-library cells found in ", sample_row$sample)
  selected <- as.matrix(matrix[genes, , drop = FALSE])
  normalized <- log1p(t(t(selected) / library_size * 10000))
  frame <- data.frame(
    sample = sample_row$sample,
    gsm = sample_row$gsm,
    region = sample_row$region,
    original_sample_order = seq_len(ncol(matrix)),
    cell_id = colnames(matrix),
    library_size = as.numeric(library_size),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  for (gene in genes) frame[[gene]] <- as.numeric(normalized[gene, ])
  fingerprints[[i]] <- frame
  inventory[[i]] <- data.frame(
    sample = sample_row$sample,
    gsm = sample_row$gsm,
    age = "E15.5 (GEO shorthand E15)",
    region = sample_row$region,
    genotype = "WT",
    deposited_cells = ncol(matrix),
    genes = nrow(matrix),
    object_class = paste(class(matrix), collapse = "|"),
    gzip_layers_removed = loaded$gzip_layers,
    first_cell_id = colnames(matrix)[1L],
    last_cell_id = colnames(matrix)[ncol(matrix)],
    stringsAsFactors = FALSE
  )
  rm(matrix, selected, normalized, loaded)
  invisible(gc())
}

write_tsv(do.call(rbind, fingerprints), args$output)
write_tsv(do.call(rbind, inventory), args$inventory)

