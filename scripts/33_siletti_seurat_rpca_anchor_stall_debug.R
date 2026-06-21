#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
})

timestamp <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")

parse_args <- function(args) {
  out <- list(
    bridge_dir = NULL,
    transfer_dir = NULL,
    outdir = NULL,
    test_id = NULL,
    nfeatures = 500L,
    dims = 10L,
    k_filter = "100",
    reduction = "rpca",
    max_reference_cells = 2000L,
    max_query_cells = 2000L,
    subtype_label_column = "candidate_jia_group",
    query_class_col = "div90_broad_class",
    exclude_label = "Excluded / not assigned to Jia-style 9 groups",
    seed = 0L,
    k_anchor = 5L,
    k_score = 30L,
    max_features = 200L,
    nn_method = "annoy",
    n_trees = 50L,
    seurat_verbose = "true"
  )
  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, "--")) stop("Unknown argument: ", key, call. = FALSE)
    name <- gsub("-", "_", substring(key, 3L))
    if (!(name %in% names(out))) stop("Unknown argument: ", key, call. = FALSE)
    if (i == length(args)) stop("Missing value for argument: ", key, call. = FALSE)
    out[[name]] <- args[[i + 1L]]
    i <- i + 2L
  }
  int_names <- c("nfeatures", "dims", "max_reference_cells", "max_query_cells", "seed", "k_anchor", "k_score", "max_features", "n_trees")
  for (name in int_names) out[[name]] <- as.integer(out[[name]])
  out$seurat_verbose <- tolower(as.character(out$seurat_verbose)) %in% c("true", "t", "1", "yes", "y")
  if (toupper(as.character(out$k_filter)) == "NA") out$k_filter <- NA_integer_ else out$k_filter <- as.integer(out$k_filter)
  out$reduction <- tolower(as.character(out$reduction))
  if (!(out$reduction %in% c("rpca", "cca", "pcaproject", "lsiproject"))) {
    stop("Unsupported reduction: ", out$reduction, call. = FALSE)
  }
  out
}

read_tsv <- function(path) {
  utils::read.delim(path, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
}

write_tsv <- function(x, path) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  utils::write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
}

rss_mb <- function() {
  status_path <- "/proc/self/status"
  if (!file.exists(status_path)) return(NA_real_)
  line <- grep("^VmRSS:", readLines(status_path, warn = FALSE), value = TRUE)
  if (length(line) == 0) return(NA_real_)
  as.numeric(gsub("[^0-9]", "", line[[1]])) / 1024
}

progress_path <- NULL
progress_rows <- list()

mark_progress <- function(step, status, detail = "") {
  row <- data.frame(
    timestamp = timestamp(),
    step = step,
    status = status,
    detail = as.character(detail),
    rss_mb = round(rss_mb(), 3),
    stringsAsFactors = FALSE
  )
  progress_rows[[length(progress_rows) + 1L]] <<- row
  if (!is.null(progress_path)) {
    con <- file(progress_path, open = if (file.exists(progress_path)) "at" else "wt")
    on.exit(close(con), add = TRUE)
    utils::write.table(row, con, sep = "\t", quote = FALSE, row.names = FALSE, col.names = !file.exists(progress_path))
  }
  message("[R ", row$timestamp, "] ", step, " [", status, "]: ", detail, " rss_mb=", row$rss_mb)
  flush.console()
  invisible(NULL)
}

timed_step <- function(step, expr, detail = "") {
  mark_progress(step, "start", detail)
  elapsed_start <- proc.time()[["elapsed"]]
  value <- force(expr)
  elapsed <- proc.time()[["elapsed"]] - elapsed_start
  mark_progress(step, "end", paste0("elapsed_sec=", round(elapsed, 3)))
  value
}

load_bridge_object <- function(prefix, bridge_dir) {
  mat <- Matrix::readMM(file.path(bridge_dir, paste0(prefix, "_counts.mtx")))
  genes <- read_tsv(file.path(bridge_dir, paste0(prefix, "_genes.tsv")))$gene
  barcodes <- read_tsv(file.path(bridge_dir, paste0(prefix, "_barcodes.tsv")))$cell_id
  meta <- read_tsv(file.path(bridge_dir, paste0(prefix, "_metadata.tsv.gz")))
  if (nrow(mat) != length(genes)) stop(prefix, ": matrix rows do not match genes", call. = FALSE)
  if (ncol(mat) != length(barcodes)) stop(prefix, ": matrix columns do not match barcodes", call. = FALSE)
  rownames(mat) <- make.unique(as.character(genes))
  colnames(mat) <- as.character(barcodes)
  rownames(meta) <- as.character(meta$seurat_cell_id)
  CreateSeuratObject(counts = mat, meta.data = meta)
}

read_selected_features <- function(transfer_dir, nfeatures) {
  path <- file.path(transfer_dir, "fast_knn", "selected_transfer_features.tsv")
  if (!file.exists(path)) return(NULL)
  features <- read_tsv(path)$gene
  unique(as.character(features))[seq_len(min(length(unique(features)), nfeatures))]
}

stratified_cells <- function(meta, group_col, max_cells, seed) {
  if (max_cells <= 0L || nrow(meta) <= max_cells) return(rownames(meta))
  set.seed(seed)
  groups <- split(rownames(meta), as.character(meta[[group_col]]), drop = TRUE)
  groups <- groups[order(names(groups))]
  per_group <- ceiling(max_cells / length(groups))
  picked <- unlist(lapply(groups, function(x) if (length(x) <= per_group) x else sample(x, per_group)), use.names = FALSE)
  if (length(picked) > max_cells) picked <- sample(picked, max_cells)
  sort(picked)
}

main <- function() {
  opt <- parse_args(commandArgs(trailingOnly = TRUE))
  if (is.null(opt$bridge_dir) || is.null(opt$transfer_dir) || is.null(opt$outdir) || is.null(opt$test_id)) {
    stop("--bridge-dir, --transfer-dir, --outdir, and --test-id are required", call. = FALSE)
  }
  test_dir <- file.path(opt$outdir, paste0("test_", opt$test_id))
  dir.create(test_dir, recursive = TRUE, showWarnings = FALSE)
  progress_path <<- file.path(test_dir, "progress.tsv")
  if (file.exists(progress_path)) file.remove(progress_path)

  warnings_seen <- character()
  error_message <- NA_character_
  completed <- FALSE
  anchor_count <- NA_integer_
  anchor_elapsed <- NA_real_

  mark_progress("run", "start", paste0("test=", opt$test_id))
  mark_progress(
    "parameters",
    "info",
    paste0(
      "features=", opt$nfeatures,
      "; dims=1:", opt$dims,
      "; k.filter=", ifelse(is.na(opt$k_filter), "NA", opt$k_filter),
      "; reduction=", opt$reduction,
      "; nn.method=", opt$nn_method,
      "; max_reference_cells=", opt$max_reference_cells,
      "; max_query_cells=", opt$max_query_cells
    )
  )

  reference <- timed_step("load_reference", load_bridge_object("reference", opt$bridge_dir), opt$bridge_dir)
  query <- timed_step("load_query", load_bridge_object("query", opt$bridge_dir), opt$bridge_dir)
  initial_reference_cells <- ncol(reference)
  initial_query_cells <- ncol(query)
  shared_genes <- intersect(rownames(reference), rownames(query))
  mark_progress("object_sizes_initial", "info", paste0("reference=", initial_reference_cells, "; query=", initial_query_cells, "; shared_genes=", length(shared_genes)))

  if (!(opt$subtype_label_column %in% colnames(reference[[]]))) stop("Missing adult subtype column: ", opt$subtype_label_column, call. = FALSE)
  if (!(opt$query_class_col %in% colnames(query[[]]))) stop("Missing query class column: ", opt$query_class_col, call. = FALSE)
  reference[[opt$subtype_label_column]][, 1] <- as.character(reference[[opt$subtype_label_column]][, 1])
  query[[opt$query_class_col]][, 1] <- as.character(query[[opt$query_class_col]][, 1])

  if (!is.na(opt$exclude_label) && nzchar(opt$exclude_label) && toupper(opt$exclude_label) != "NONE") {
    keep <- rownames(reference[[]])[as.character(reference[[opt$subtype_label_column]][, 1]) != opt$exclude_label]
    reference <- subset(reference, cells = keep)
    mark_progress("filter_reference_label", "end", paste0("before=", initial_reference_cells, "; after=", ncol(reference), "; excluded=", opt$exclude_label))
  }

  ref_cells <- stratified_cells(reference[[]], opt$subtype_label_column, opt$max_reference_cells, opt$seed)
  query_cells <- stratified_cells(query[[]], opt$query_class_col, opt$max_query_cells, opt$seed + 1L)
  reference <- subset(reference, cells = ref_cells)
  query <- subset(query, cells = query_cells)
  mark_progress("stratified_downsample", "end", paste0("reference=", ncol(reference), "; query=", ncol(query)))

  ref_counts <- as.data.frame(table(reference[[opt$subtype_label_column]][, 1]), stringsAsFactors = FALSE)
  colnames(ref_counts) <- c("adult_subtype", "n_cells")
  query_counts <- as.data.frame(table(query[[opt$query_class_col]][, 1]), stringsAsFactors = FALSE)
  colnames(query_counts) <- c("div90_class", "n_cells")
  write_tsv(ref_counts, file.path(test_dir, "adult_subtype_downsample_counts.tsv"))
  write_tsv(query_counts, file.path(test_dir, "div90_class_downsample_counts.tsv"))

  selected_features <- read_selected_features(opt$transfer_dir, opt$nfeatures)
  if (is.null(selected_features)) {
    selected_features <- shared_genes[seq_len(min(length(shared_genes), opt$nfeatures))]
    feature_source <- "first_shared_genes_fallback"
  } else {
    feature_source <- "existing_fast_knn_selected_features"
  }
  features <- intersect(selected_features, intersect(rownames(reference), rownames(query)))
  features <- features[seq_len(min(length(features), opt$nfeatures))]
  if (length(features) < 50L) stop("Too few shared selected features: ", length(features), call. = FALSE)
  write_tsv(data.frame(feature = features), file.path(test_dir, "selected_features.tsv"))
  mark_progress("selected_features", "info", paste0("source=", feature_source, "; n=", length(features)))

  npcs <- min(50L, length(features) - 1L, ncol(reference) - 1L, ncol(query) - 1L)
  dims_use <- seq_len(min(opt$dims, npcs))
  VariableFeatures(reference) <- features
  VariableFeatures(query) <- features

  reference <- timed_step("normalize_reference", NormalizeData(reference, normalization.method = "LogNormalize", verbose = opt$seurat_verbose))
  query <- timed_step("normalize_query", NormalizeData(query, normalization.method = "LogNormalize", verbose = opt$seurat_verbose))
  reference <- timed_step("scale_reference", ScaleData(reference, features = features, verbose = opt$seurat_verbose))
  query <- timed_step("scale_query", ScaleData(query, features = features, verbose = opt$seurat_verbose))
  reference <- timed_step("pca_reference", RunPCA(reference, features = features, npcs = npcs, verbose = opt$seurat_verbose))
  query <- timed_step("pca_query", RunPCA(query, features = features, npcs = npcs, verbose = opt$seurat_verbose))

  object_sizes <- data.frame(
    object = c("reference", "query", "reference_scale_data", "query_scale_data"),
    size_bytes = c(
      as.numeric(object.size(reference)),
      as.numeric(object.size(query)),
      as.numeric(object.size(reference[["RNA"]]$scale.data)),
      as.numeric(object.size(query[["RNA"]]$scale.data))
    ),
    stringsAsFactors = FALSE
  )
  write_tsv(object_sizes, file.path(test_dir, "object_sizes.tsv"))

  anchor_args <- list(
    reference = reference,
    query = query,
    normalization.method = "LogNormalize",
    reference.assay = "RNA",
    query.assay = "RNA",
    reduction = opt$reduction,
    features = features,
    dims = dims_use,
    npcs = npcs,
    k.anchor = opt$k_anchor,
    k.filter = opt$k_filter,
    k.score = opt$k_score,
    max.features = opt$max_features,
    nn.method = opt$nn_method,
    n.trees = opt$n_trees,
    verbose = opt$seurat_verbose
  )
  if (opt$reduction %in% c("rpca", "pcaproject", "lsiproject")) {
    anchor_args$reference.reduction <- "pca"
  }
  mark_progress(
    "find_transfer_anchors",
    "start",
    paste0(
      "features=", length(features),
      "; dims=1:", max(dims_use),
      "; k.filter=", ifelse(is.na(opt$k_filter), "NA", opt$k_filter),
      "; reduction=", opt$reduction,
      "; nn.method=", opt$nn_method
    )
  )
  anchor_start <- proc.time()[["elapsed"]]
  anchors <- NULL
  result <- tryCatch(
    withCallingHandlers(
      {
        anchors <- do.call(FindTransferAnchors, anchor_args)
        anchors
      },
      warning = function(w) {
        warnings_seen <<- c(warnings_seen, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) {
      error_message <<- conditionMessage(e)
      NULL
    }
  )
  anchor_elapsed <- proc.time()[["elapsed"]] - anchor_start
  if (!is.null(result)) {
    completed <- TRUE
    anchor_count <- nrow(as.data.frame(result@anchors))
    mark_progress("find_transfer_anchors", "end", paste0("elapsed_sec=", round(anchor_elapsed, 3), "; anchors=", anchor_count))
  } else {
    mark_progress("find_transfer_anchors", "error", paste0("elapsed_sec=", round(anchor_elapsed, 3), "; error=", error_message))
  }

  if (length(warnings_seen) > 0L) writeLines(unique(warnings_seen), file.path(test_dir, "warnings.txt"))
  utils::capture.output(sessionInfo(), file = file.path(test_dir, "session_info.txt"))

  summary <- data.frame(
    test_id = opt$test_id,
    completed = completed,
    n_anchors = anchor_count,
    find_transfer_anchors_elapsed_sec = round(anchor_elapsed, 3),
    error = error_message,
    n_warnings = length(unique(warnings_seen)),
    reference_cells_initial = initial_reference_cells,
    query_cells_initial = initial_query_cells,
    reference_cells_used = ncol(reference),
    query_cells_used = ncol(query),
    n_features = length(features),
    npcs = npcs,
    dims = paste0("1:", max(dims_use)),
    k_filter = ifelse(is.na(opt$k_filter), "NA", as.character(opt$k_filter)),
    reduction = opt$reduction,
    nn_method = opt$nn_method,
    rss_mb_end = round(rss_mb(), 3),
    reference_object_bytes = as.numeric(object.size(reference)),
    query_object_bytes = as.numeric(object.size(query)),
    stringsAsFactors = FALSE
  )
  write_tsv(summary, file.path(test_dir, "test_summary.tsv"))
  mark_progress("run", ifelse(completed, "end", "failed"), paste0("completed=", completed))
  if (!completed) quit(status = 1L, save = "no")
}

main()
