#!/usr/bin/env Rscript
#' Run observable, non-filtering Step 03 scDblFinder detection in native R.

run_started <- Sys.time()

#' Return process resident memory in MiB when Linux exposes it.
process_memory <- function() {
  status <- "/proc/self/status"
  if (!file.exists(status)) return("unavailable")
  line <- grep("^VmRSS:", readLines(status, warn = FALSE), value = TRUE)
  if (!length(line)) return("unavailable")
  kib <- as.numeric(sub("^VmRSS:\\s+([0-9]+).*", "\\1", line[[1L]]))
  sprintf("%.1f MiB RSS", kib / 1024)
}

#' Emit one timestamped, immediately flushed checkpoint.
log_msg <- function(...) {
  elapsed <- round(as.numeric(difftime(Sys.time(), run_started, units = "mins")), 2)
  cat(
    sprintf("[%s] [elapsed=%.2f min] [memory=%s] ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), elapsed, process_memory()),
    ...,
    "\n",
    sep = ""
  )
  flush.console()
}

#' Parse named command-line options supplied as --name value pairs.
parse_options <- function(arguments) {
  if (length(arguments) %% 2L != 0L || any(!startsWith(arguments[seq(1L, length(arguments), 2L)], "--"))) {
    stop("Options must be supplied as --name value pairs")
  }
  keys <- sub("^--", "", arguments[seq(1L, length(arguments), 2L)])
  values <- arguments[seq(2L, length(arguments), 2L)]
  stats::setNames(as.list(values), keys)
}

#' Require a named command-line option.
require_option <- function(options, name) {
  value <- options[[name]]
  if (is.null(value) || !nzchar(value)) stop("Missing required option: --", name)
  value
}

#' Verify one preparation check by identifier.
require_prepare_check <- function(checks, check_id) {
  row <- checks[checks$check_id == check_id, , drop = FALSE]
  if (nrow(row) != 1L || row$status[[1L]] != "PASS") stop("Preparation check not PASS: ", check_id)
  invisible(row)
}

#' Write exact runtime versions and session details.
write_versions <- function(output_dir) {
  packages <- c("scDblFinder", "SingleCellExperiment", "DropletUtils", "data.table", "rhdf5", "BiocParallel")
  versions <- data.frame(
    component = c("R", packages),
    version = c(R.version.string, vapply(packages, function(x) as.character(utils::packageVersion(x)), character(1L)))
  )
  data.table::fwrite(versions, file.path(output_dir, "r_software_versions.tsv"), sep = "\t")
  capture.output(sessionInfo(), file = file.path(output_dir, "r_session_info.txt"))
}

#' Run the single approved observable scDblFinder invocation.
run_scdblfinder <- function(sce, seed) {
  bp <- BiocParallel::SerialParam(progressbar = TRUE)
  set.seed(as.integer(seed))
  log_msg("ENTER scDblFinder")
  t0 <- Sys.time()
  pt0 <- proc.time()
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
  out
}

#' Serialize and validate scores and calls without removing any cells.
write_results <- function(out, cell_ids, technical_sample_ids, output_dir) {
  results <- as.data.frame(out)
  if (!all(c("score", "class") %in% colnames(results))) stop("scDblFinder score return lacks score/class")
  if (nrow(results) != length(cell_ids)) stop("scDblFinder returned the wrong number of cells")
  if (!is.null(rownames(results)) && !identical(as.character(rownames(results)), as.character(cell_ids))) {
    stop("scDblFinder result order differs from input cell order")
  }
  per_cell <- data.frame(
    cell_id = as.character(cell_ids),
    primary_score = as.numeric(results$score),
    primary_class = as.character(results$class),
    stringsAsFactors = FALSE
  )
  data.table::fwrite(per_cell, file.path(output_dir, "scdblfinder_per_cell_results.tsv.gz"), sep = "\t", compress = "gzip")
  class_summary <- as.data.frame(table(class = per_cell$primary_class), stringsAsFactors = FALSE)
  class_summary$percent <- 100 * class_summary$Freq / nrow(per_cell)
  data.table::fwrite(class_summary, file.path(output_dir, "scdblfinder_class_summary.tsv"), sep = "\t")
  sample_frame <- data.table::data.table(
    technical_sample_id = as.character(technical_sample_ids),
    primary_score = per_cell$primary_score,
    primary_class = per_cell$primary_class
  )
  sample_summary <- sample_frame[, .(
    cells = .N,
    called_doublets = sum(primary_class == "doublet"),
    called_doublet_pct = 100 * mean(primary_class == "doublet"),
    mean_score = mean(primary_score),
    median_score = median(primary_score)
  ), by = technical_sample_id]
  data.table::fwrite(sample_summary, file.path(output_dir, "scdblfinder_per_sample_summary.tsv"), sep = "\t")
  contract <- data.frame(
    setting = c("samples", "clusters", "dbr.sd", "dbr", "other_model_parameters", "BPPARAM", "verbose", "returnType", "invocations", "cell_removal"),
    value = c("capture_id (constant GEX_1)", "TRUE", "1", "not supplied", "package defaults", "SerialParam(progressbar=TRUE)", "TRUE", "scores", "one", "none"),
    role = c(rep("approved scientific", 5), rep("observability/execution", 4), "review boundary")
  )
  data.table::fwrite(contract, file.path(output_dir, "scdblfinder_method_contract.tsv"), sep = "\t")
}

log_msg("CHECKPOINT script start")
options <- parse_options(commandArgs(trailingOnly = TRUE))
output_dir <- require_option(options, "output-dir")
if (dir.exists(output_dir) || file.exists(output_dir)) stop("Refusing existing R output directory: ", output_dir)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

log_msg("CHECKPOINT package loading begin")
suppressPackageStartupMessages({
  library(BiocParallel)
  library(data.table)
  library(DropletUtils)
  library(scDblFinder)
  library(SingleCellExperiment)
})
log_msg("CHECKPOINT package loading complete")

log_msg("CHECKPOINT sessionInfo/package versions begin")
write_versions(output_dir)
log_msg("R=", R.version.string, "; scDblFinder=", as.character(utils::packageVersion("scDblFinder")), "; BiocParallel=", as.character(utils::packageVersion("BiocParallel")))
log_msg("CHECKPOINT sessionInfo/package versions complete")

log_msg("CHECKPOINT input object loading begin")
sce <- DropletUtils::read10xCounts(require_option(options, "bridge-h5"), type = "HDF5", col.names = TRUE)
metadata <- data.table::fread(require_option(options, "metadata-tsv"), sep = "\t", data.table = FALSE)
log_msg("CHECKPOINT input object loading complete")

expected_cells <- as.integer(require_option(options, "expected-cells"))
expected_genes <- as.integer(require_option(options, "expected-genes"))
capture_id <- require_option(options, "capture-id")
log_msg("CHECKPOINT input dimensions begin")
stopifnot(identical(dim(sce), c(expected_genes, expected_cells)))
stopifnot(nrow(metadata) == expected_cells, identical(as.character(colnames(sce)), as.character(metadata$cell_id)))
log_msg("CHECKPOINT input dimensions complete; genes=", nrow(sce), "; cells=", ncol(sce))

log_msg("CHECKPOINT raw-count fingerprint verification begin")
prepare_checks <- data.table::fread(require_option(options, "prepare-checks-tsv"), sep = "\t", data.table = FALSE)
sha_output <- system2("sha256sum", require_option(options, "input-h5ad"), stdout = TRUE, stderr = TRUE)
sha_status <- attr(sha_output, "status")
if (!is.null(sha_status) && sha_status != 0L) stop("sha256sum failed: ", paste(sha_output, collapse = " "))
observed_sha256 <- strsplit(sha_output[[1L]], "[[:space:]]+")[[1L]][[1L]]
expected_sha256 <- require_option(options, "expected-input-sha256")
if (!identical(observed_sha256, expected_sha256)) stop("Input SHA-256 mismatch")
for (check_id in c("input_sha256", "bridge_data_identity", "bridge_indices_identity", "bridge_indptr_identity")) {
  require_prepare_check(prepare_checks, check_id)
}
log_msg("CHECKPOINT raw-count fingerprint verification complete; sha256=", observed_sha256, "; bridge arrays exact=TRUE")

log_msg("CHECKPOINT capture_id validation begin")
stopifnot(!anyDuplicated(metadata$cell_id), identical(unique(as.character(metadata$capture_id)), capture_id), all(metadata$capture_id == capture_id))
for (field in colnames(metadata)) if (field != "cell_id") colData(sce)[[field]] <- metadata[[field]]
log_msg("CHECKPOINT capture_id validation complete; capture_id=", capture_id, "; unique captures=", length(unique(metadata$capture_id)))

log_msg("CHECKPOINT beginning scDblFinder")
out <- run_scdblfinder(sce, require_option(options, "primary-seed"))
log_msg("CHECKPOINT scDblFinder returned")

log_msg("CHECKPOINT score extraction begin")
classes <- as.character(out$class)
scores <- as.numeric(out$score)
stopifnot(length(scores) == expected_cells, length(classes) == expected_cells)
log_msg("CHECKPOINT score extraction complete; rows=", length(scores))

log_msg("CHECKPOINT class summary begin")
for (line in capture.output(print(table(classes, useNA = "ifany")))) log_msg(line)
log_msg("CHECKPOINT class summary complete")

log_msg("CHECKPOINT per-sample summary begin")
sample_log <- data.table::data.table(technical_sample_id = metadata$technical_sample_id, score = scores, class = classes)[, .(
  cells = .N,
  called_doublets = sum(class == "doublet"),
  called_doublet_pct = 100 * mean(class == "doublet")
), by = technical_sample_id]
for (line in capture.output(print(sample_log))) log_msg(line)
log_msg("CHECKPOINT per-sample summary complete")

log_msg("CHECKPOINT output serialization begin")
write_results(out, colnames(sce), metadata$technical_sample_id, output_dir)
log_msg("CHECKPOINT output serialization complete")

log_msg("CHECKPOINT output validation begin")
saved <- data.table::fread(file.path(output_dir, "scdblfinder_per_cell_results.tsv.gz"), sep = "\t", data.table = FALSE)
stopifnot(
  nrow(saved) == expected_cells,
  identical(as.character(saved$cell_id), as.character(colnames(sce))),
  all(is.finite(saved$primary_score)),
  all(saved$primary_score >= 0 & saved$primary_score <= 1),
  all(saved$primary_class %in% c("singlet", "doublet"))
)
log_msg("CHECKPOINT output validation complete; cells retained=", nrow(saved), "; cells removed=0")
log_msg("CHECKPOINT script completion; no cells removed")
