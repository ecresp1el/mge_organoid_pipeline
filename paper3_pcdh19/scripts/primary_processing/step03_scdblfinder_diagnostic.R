#!/usr/bin/env Rscript
# Instrumented, non-publishing Step 03 scDblFinder diagnostic.
#
# This script runs exactly one serial scDblFinder scores-only invocation. It
# never subsets, rewrites, or publishes the approved Step 02 object.

run_started <- Sys.time()

process_memory <- function() {
  path <- "/proc/self/status"
  if (!file.exists(path)) return("memory=unavailable")
  status <- readLines(path, warn = FALSE)
  values <- grep("^(VmRSS|VmHWM|VmSize|VmPeak):", status, value = TRUE)
  if (!length(values)) return("memory=unavailable")
  paste(trimws(values), collapse = "; ")
}

log_msg <- function(...) {
  elapsed <- as.numeric(difftime(Sys.time(), run_started, units = "mins"))
  cat(
    sprintf(
      "[%s] elapsed=%.2f min pid=%d %s | ",
      format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
      elapsed,
      Sys.getpid(),
      process_memory()
    ),
    ...,
    "\n",
    sep = ""
  )
  flush.console()
}

parse_options <- function(arguments) {
  if (length(arguments) %% 2L != 0L || any(!startsWith(arguments[seq(1L, length(arguments), 2L)], "--"))) {
    stop("Options must be supplied as --name value pairs")
  }
  keys <- sub("^--", "", arguments[seq(1L, length(arguments), 2L)])
  values <- arguments[seq(2L, length(arguments), 2L)]
  stats::setNames(as.list(values), keys)
}

require_option <- function(options, name) {
  value <- options[[name]]
  if (is.null(value) || !nzchar(value)) stop("Missing required option: --", name)
  value
}

log_msg("CHECKPOINT script start")
options <- parse_options(commandArgs(trailingOnly = TRUE))

bridge_h5 <- require_option(options, "bridge-h5")
metadata_tsv <- require_option(options, "metadata-tsv")
prepare_checks_tsv <- require_option(options, "prepare-checks-tsv")
input_h5ad <- require_option(options, "input-h5ad")
expected_input_sha256 <- require_option(options, "expected-input-sha256")
output_dir <- require_option(options, "output-dir")
expected_cells <- as.integer(require_option(options, "expected-cells"))
expected_genes <- as.integer(require_option(options, "expected-genes"))
capture_id <- require_option(options, "capture-id")
primary_seed <- as.integer(require_option(options, "primary-seed"))

for (path in c(bridge_h5, metadata_tsv, prepare_checks_tsv, input_h5ad)) {
  if (!file.exists(path)) stop("Missing diagnostic input: ", path)
}
if (dir.exists(output_dir) || file.exists(output_dir)) stop("Refusing existing diagnostic output directory: ", output_dir)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

log_msg("CHECKPOINT package loading BEGIN")
suppressPackageStartupMessages({
  library(BiocParallel)
  library(data.table)
  library(DropletUtils)
  library(scDblFinder)
  library(SingleCellExperiment)
})
log_msg("CHECKPOINT package loading END")

log_msg("CHECKPOINT sessionInfo/package versions BEGIN")
packages <- c("scDblFinder", "SingleCellExperiment", "DropletUtils", "data.table", "R6", "rhdf5", "BiocParallel")
for (package in packages) log_msg("PACKAGE ", package, "=", as.character(packageVersion(package)))
for (line in capture.output(sessionInfo())) log_msg("SESSION ", line)
log_msg("CHECKPOINT sessionInfo/package versions END")

log_msg("CHECKPOINT raw-count fingerprint verification BEGIN")
sha_output <- system2("sha256sum", input_h5ad, stdout = TRUE, stderr = TRUE)
sha_status <- attr(sha_output, "status")
if (!is.null(sha_status) && sha_status != 0L) stop("sha256sum failed: ", paste(sha_output, collapse = " "))
observed_input_sha256 <- strsplit(sha_output[[1L]], "[[:space:]]+")[[1L]][[1L]]
if (!identical(observed_input_sha256, expected_input_sha256)) {
  stop("Input SHA-256 mismatch: observed ", observed_input_sha256, "; expected ", expected_input_sha256)
}
prepare_checks <- data.table::fread(prepare_checks_tsv, sep = "\t", data.table = FALSE)
required_bridge_checks <- c("input_sha256", "bridge_data_identity", "bridge_indices_identity", "bridge_indptr_identity")
bridge_rows <- prepare_checks[match(required_bridge_checks, prepare_checks$check_id), , drop = FALSE]
if (any(is.na(bridge_rows$check_id)) || any(bridge_rows$status != "PASS")) {
  stop("Frozen bridge identity checks are incomplete or failed")
}
log_msg("CHECKPOINT raw-count fingerprint verification END sha256=", observed_input_sha256, "; bridge arrays exact=TRUE")

log_msg("CHECKPOINT input object loading BEGIN bridge=", bridge_h5)
sce <- DropletUtils::read10xCounts(bridge_h5, type = "HDF5", col.names = TRUE)
metadata <- data.table::fread(metadata_tsv, sep = "\t", data.table = FALSE)
log_msg("CHECKPOINT input object loading END")

log_msg("CHECKPOINT input dimensions BEGIN")
if (!identical(dim(sce), c(expected_genes, expected_cells))) {
  stop("Unexpected SCE dimensions: ", paste(dim(sce), collapse = "x"))
}
if (nrow(metadata) != expected_cells || anyDuplicated(metadata$cell_id)) stop("Diagnostic metadata row/ID validation failed")
if (!identical(as.character(colnames(sce)), as.character(metadata$cell_id))) stop("Diagnostic metadata order differs from SCE")
log_msg("CHECKPOINT input dimensions END genes=", nrow(sce), "; cells=", ncol(sce), "; metadata_rows=", nrow(metadata))

log_msg("CHECKPOINT capture_id validation BEGIN")
if (!identical(unique(as.character(metadata$capture_id)), capture_id) || !all(metadata$capture_id == capture_id)) {
  stop("Expected one constant capture_id=", capture_id)
}
for (field in colnames(metadata)) {
  if (field != "cell_id") colData(sce)[[field]] <- metadata[[field]]
}
log_msg("CHECKPOINT capture_id validation END capture_id=", capture_id, "; unique_captures=", length(unique(metadata$capture_id)))

bp <- BiocParallel::SerialParam(progressbar = TRUE)
log_msg("CHECKPOINT BPPARAM configured class=", class(bp)[[1L]], "; workers=", BiocParallel::bpnworkers(bp), "; progressbar=", BiocParallel::bpprogressbar(bp))
log_msg("CHECKPOINT beginning scDblFinder seed=", primary_seed, "; samples=capture_id; clusters=TRUE; dbr.sd=1; dbr=not supplied; verbose=TRUE; returnType=scores")
log_msg("ENTER scDblFinder")

t0 <- Sys.time()
pt0 <- proc.time()

set.seed(primary_seed)
out <- scDblFinder::scDblFinder(
  sce,
  samples = "capture_id",
  clusters = TRUE,
  dbr.sd = 1,
  verbose = TRUE,
  BPPARAM = bp,
  returnType = "scores"
)

log_msg(
  "EXIT scDblFinder; wall=",
  round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2),
  " min; CPU=",
  paste(round(proc.time() - pt0, 2), collapse = ",")
)
log_msg("CHECKPOINT scDblFinder returned")

log_msg("CHECKPOINT score extraction BEGIN")
scores <- as.data.frame(out)
scores$cell_id <- rownames(scores)
required_fields <- c("cell_id", "score", "class")
if (!all(required_fields %in% colnames(scores))) stop("Missing required scDblFinder score fields")
if (!identical(as.character(scores$cell_id), as.character(colnames(sce)))) stop("Returned score order differs from input cells")
result <- data.frame(
  cell_id = scores$cell_id,
  technical_sample_id = as.character(metadata$technical_sample_id),
  capture_id = as.character(metadata$capture_id),
  score = as.numeric(scores$score),
  class = as.character(scores$class),
  stringsAsFactors = FALSE
)
log_msg("CHECKPOINT score extraction END rows=", nrow(result), "; score_min=", min(result$score), "; score_max=", max(result$score))

log_msg("CHECKPOINT class summary BEGIN")
class_summary <- as.data.frame(table(result$class), stringsAsFactors = FALSE)
colnames(class_summary) <- c("class", "cells")
class_summary$pct <- 100 * class_summary$cells / nrow(result)
for (index in seq_len(nrow(class_summary))) {
  log_msg("CLASS ", class_summary$class[[index]], " cells=", class_summary$cells[[index]], " pct=", round(class_summary$pct[[index]], 6))
}
log_msg("CHECKPOINT class summary END")

log_msg("CHECKPOINT per-sample summary BEGIN")
result_dt <- data.table::as.data.table(result)
sample_summary <- result_dt[, .(
  cells = .N,
  called_doublets = sum(class == "doublet"),
  called_doublet_pct = 100 * mean(class == "doublet"),
  mean_score = mean(score),
  median_score = median(score)
), by = technical_sample_id]
for (index in seq_len(nrow(sample_summary))) {
  log_msg(
    "SAMPLE ", sample_summary$technical_sample_id[[index]],
    " cells=", sample_summary$cells[[index]],
    " doublets=", sample_summary$called_doublets[[index]],
    " pct=", round(sample_summary$called_doublet_pct[[index]], 6)
  )
}
log_msg("CHECKPOINT per-sample summary END samples=", nrow(sample_summary))

log_msg("CHECKPOINT output serialization BEGIN output_dir=", output_dir)
data.table::fwrite(result, file.path(output_dir, "diagnostic_scdblfinder_scores.tsv.gz"), sep = "\t", compress = "gzip")
data.table::fwrite(class_summary, file.path(output_dir, "diagnostic_class_summary.tsv"), sep = "\t")
data.table::fwrite(sample_summary, file.path(output_dir, "diagnostic_per_sample_summary.tsv"), sep = "\t")
data.table::fwrite(
  data.frame(
    setting = c("capture_id", "clusters", "dbr.sd", "dbr", "verbose", "BPPARAM", "returnType", "seed", "scDblFinder_invocations", "second_seed", "cells_removed", "published_as_step03"),
    value = c(capture_id, "TRUE", "1", "not supplied", "TRUE", "SerialParam(progressbar=TRUE)", "scores", primary_seed, "1", "none", "0", "FALSE")
  ),
  file.path(output_dir, "diagnostic_method_contract.tsv"),
  sep = "\t"
)
log_msg("CHECKPOINT output serialization END")

log_msg("CHECKPOINT output validation BEGIN")
reopened <- data.table::fread(file.path(output_dir, "diagnostic_scdblfinder_scores.tsv.gz"), sep = "\t", data.table = FALSE)
if (nrow(reopened) != expected_cells || anyDuplicated(reopened$cell_id)) stop("Serialized diagnostic score rows/IDs failed validation")
if (!identical(as.character(reopened$cell_id), as.character(metadata$cell_id))) stop("Serialized diagnostic cell order differs from input")
if (!all(is.finite(reopened$score)) || any(reopened$score < 0) || any(reopened$score > 1)) stop("Serialized scores are outside finite [0,1]")
if (!setequal(unique(reopened$class), c("singlet", "doublet"))) stop("Unexpected serialized class values")
log_msg("CHECKPOINT output validation END rows=", nrow(reopened), "; cells_removed=0; published_as_step03=FALSE")

writeLines(
  c(
    "DIAGNOSTIC_SUCCESS",
    paste0("cells=", expected_cells),
    paste0("genes=", expected_genes),
    "cells_removed=0",
    "published_as_step03=FALSE",
    paste0("completed=", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"))
  ),
  file.path(output_dir, "DIAGNOSTIC_SUCCESS.txt")
)
log_msg("CHECKPOINT script completion; diagnostic success; cells_removed=0; no Step 03 publication")
