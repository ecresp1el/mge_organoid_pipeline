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
    outdir = NULL,
    mode = NULL,
    max_reference_cells = 1000L,
    max_query_cells = 1000L,
    nfeatures = 2000L,
    npcs = 20L,
    dims = 20L,
    seed = 0L,
    subtype_label_column = "candidate_jia_group",
    query_class_col = "div90_broad_class",
    exclude_label = "Excluded / not assigned to Jia-style 9 groups",
    k_anchor = 5L,
    k_filter = "NA",
    k_score = 30L,
    max_features = 200L,
    nn_method = "annoy",
    n_trees = 10L,
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
  int_names <- c("max_reference_cells", "max_query_cells", "nfeatures", "npcs", "dims", "seed", "k_anchor", "k_score", "max_features", "n_trees")
  for (name in int_names) out[[name]] <- as.integer(out[[name]])
  if (toupper(as.character(out$k_filter)) == "NA") out$k_filter <- NA_integer_ else out$k_filter <- as.integer(out$k_filter)
  out$seurat_verbose <- tolower(as.character(out$seurat_verbose)) %in% c("true", "t", "1", "yes", "y")
  allowed <- c("integration_rpca", "transfer_rpca", "transfer_rpca_l2_false", "transfer_rpca_approx_false")
  if (is.null(out$mode) || !(out$mode %in% allowed)) stop("--mode must be one of: ", paste(allowed, collapse = ", "), call. = FALSE)
  out
}

read_tsv <- function(path) {
  utils::read.delim(path, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
}

write_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
}

rss_mb <- function() {
  status_path <- "/proc/self/status"
  if (!file.exists(status_path)) return(NA_real_)
  line <- grep("^VmRSS:", readLines(status_path, warn = FALSE), value = TRUE)
  if (length(line) == 0L) return(NA_real_)
  as.numeric(gsub("[^0-9]", "", line[[1]])) / 1024
}

progress_path <- NULL

mark_progress <- function(step, status, detail = "") {
  row <- data.frame(
    timestamp = timestamp(),
    step = step,
    status = status,
    detail = as.character(detail),
    rss_mb = round(rss_mb(), 3),
    stringsAsFactors = FALSE
  )
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
  CreateSeuratObject(counts = mat, meta.data = meta, min.cells = 0, min.features = 0)
}

layer_data <- function(object, assay, layer) {
  tryCatch(
    LayerData(object = object[[assay]], layer = layer),
    error = function(e) NULL
  )
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

make_stripped_from_data_layer <- function(object, assay, features) {
  data <- layer_data(object, assay, "data")
  data <- data[features, colnames(object), drop = FALSE]
  meta <- object[[]]
  stripped <- CreateSeuratObject(counts = data, assay = assay, meta.data = meta, min.cells = 0, min.features = 0)
  stripped <- SetAssayData(stripped, assay = assay, layer = "data", new.data = data)
  DefaultAssay(stripped) <- assay
  stripped
}

prepare_base_objects <- function(opt) {
  reference <- timed_step("load_reference", load_bridge_object("reference", opt$bridge_dir), opt$bridge_dir)
  query <- timed_step("load_query", load_bridge_object("query", opt$bridge_dir), opt$bridge_dir)
  if (!(opt$subtype_label_column %in% colnames(reference[[]]))) stop("Missing adult subtype column: ", opt$subtype_label_column, call. = FALSE)
  if (!(opt$query_class_col %in% colnames(query[[]]))) stop("Missing query class column: ", opt$query_class_col, call. = FALSE)
  reference$adult_subtype_label <- as.character(reference[[opt$subtype_label_column]][, 1])
  query$div90_class <- as.character(query[[opt$query_class_col]][, 1])

  if (!is.na(opt$exclude_label) && nzchar(opt$exclude_label) && toupper(opt$exclude_label) != "NONE") {
    before <- ncol(reference)
    reference <- subset(reference, cells = colnames(reference)[reference$adult_subtype_label != opt$exclude_label])
    mark_progress("filter_reference_label", "end", paste0("before=", before, "; after=", ncol(reference)))
  }
  reference <- subset(reference, cells = stratified_cells(reference[[]], "adult_subtype_label", opt$max_reference_cells, opt$seed))
  query <- subset(query, cells = stratified_cells(query[[]], "div90_class", opt$max_query_cells, opt$seed + 1L))
  mark_progress("stratified_downsample", "end", paste0("reference=", ncol(reference), "; query=", ncol(query)))

  shared <- intersect(rownames(reference), rownames(query))
  reference <- subset(reference, features = shared)
  query <- subset(query, features = shared)
  mark_progress("shared_genes_only", "end", paste0("shared=", length(shared)))

  reference <- timed_step("normalize_reference", NormalizeData(reference, normalization.method = "LogNormalize", verbose = opt$seurat_verbose))
  query <- timed_step("normalize_query", NormalizeData(query, normalization.method = "LogNormalize", verbose = opt$seurat_verbose))
  reference <- timed_step("find_variable_features_reference", FindVariableFeatures(reference, selection.method = "vst", nfeatures = opt$nfeatures, verbose = opt$seurat_verbose))
  query <- timed_step("find_variable_features_query", FindVariableFeatures(query, selection.method = "vst", nfeatures = opt$nfeatures, verbose = opt$seurat_verbose))

  stripped_reference <- make_stripped_from_data_layer(reference, "RNA", shared)
  stripped_query <- make_stripped_from_data_layer(query, "RNA", shared)
  VariableFeatures(stripped_reference) <- VariableFeatures(reference)
  VariableFeatures(stripped_query) <- VariableFeatures(query)
  list(reference = stripped_reference, query = stripped_query, shared = shared)
}

run_integration_rpca <- function(reference, query, opt, mode_dir) {
  object_list <- list(reference = reference, query = query)
  features <- timed_step("select_integration_features", SelectIntegrationFeatures(object.list = object_list, nfeatures = opt$nfeatures), paste0("nfeatures=", opt$nfeatures))
  write_tsv(data.frame(feature = features), file.path(mode_dir, "selected_integration_features.tsv"))
  object_list <- lapply(seq_along(object_list), function(i) {
    name <- names(object_list)[[i]]
    obj <- object_list[[i]]
    obj <- timed_step(paste0("scale_", name), ScaleData(obj, features = features, verbose = opt$seurat_verbose))
    obj <- timed_step(paste0("pca_", name), RunPCA(obj, features = features, npcs = opt$npcs, verbose = opt$seurat_verbose))
    obj
  })
  names(object_list) <- c("reference", "query")
  dims_use <- seq_len(min(opt$dims, opt$npcs))
  mark_progress("find_integration_anchors_rpca", "start", paste0("features=", length(features), "; dims=1:", max(dims_use)))
  elapsed_start <- proc.time()[["elapsed"]]
  anchors <- FindIntegrationAnchors(
    object.list = object_list,
    anchor.features = features,
    reduction = "rpca",
    dims = dims_use,
    k.anchor = opt$k_anchor,
    k.filter = opt$k_filter,
    k.score = opt$k_score,
    nn.method = opt$nn_method,
    n.trees = opt$n_trees,
    verbose = opt$seurat_verbose
  )
  elapsed <- proc.time()[["elapsed"]] - elapsed_start
  anchor_count <- nrow(as.data.frame(anchors@anchors))
  mark_progress("find_integration_anchors_rpca", "end", paste0("elapsed_sec=", round(elapsed, 3), "; anchors=", anchor_count))
  data.frame(mode = opt$mode, completed = TRUE, n_anchors = anchor_count, elapsed_sec = round(elapsed, 3), stringsAsFactors = FALSE)
}

run_transfer_rpca <- function(reference, query, opt, mode_dir) {
  object_list <- list(reference = reference, query = query)
  features <- timed_step("select_integration_features", SelectIntegrationFeatures(object.list = object_list, nfeatures = opt$nfeatures), paste0("nfeatures=", opt$nfeatures))
  write_tsv(data.frame(feature = features), file.path(mode_dir, "selected_integration_features.tsv"))
  reference <- timed_step("scale_reference", ScaleData(reference, features = features, verbose = opt$seurat_verbose))
  query <- timed_step("scale_query", ScaleData(query, features = features, verbose = opt$seurat_verbose))
  reference <- timed_step("pca_reference", RunPCA(reference, features = features, npcs = opt$npcs, verbose = opt$seurat_verbose))
  query <- timed_step("pca_query", RunPCA(query, features = features, npcs = opt$npcs, verbose = opt$seurat_verbose))
  dims_use <- seq_len(min(opt$dims, opt$npcs))
  anchor_args <- list(
    reference = reference,
    query = query,
    normalization.method = "LogNormalize",
    reference.assay = "RNA",
    query.assay = "RNA",
    reduction = "rpca",
    reference.reduction = "pca",
    features = features,
    dims = dims_use,
    npcs = opt$npcs,
    k.anchor = opt$k_anchor,
    k.filter = opt$k_filter,
    k.score = opt$k_score,
    max.features = opt$max_features,
    nn.method = opt$nn_method,
    n.trees = opt$n_trees,
    verbose = opt$seurat_verbose
  )
  if (opt$mode == "transfer_rpca_l2_false") {
    anchor_args$l2.norm <- FALSE
  }
  if (opt$mode == "transfer_rpca_approx_false") {
    anchor_args$approx.pca <- FALSE
  }
  mark_progress(
    "find_transfer_anchors_rpca",
    "start",
    paste0(
      "features=", length(features),
      "; dims=1:", max(dims_use),
      "; l2.norm=", ifelse(is.null(anchor_args$l2.norm), "default", anchor_args$l2.norm),
      "; approx.pca=", ifelse(is.null(anchor_args$approx.pca), "default", anchor_args$approx.pca)
    )
  )
  elapsed_start <- proc.time()[["elapsed"]]
  anchors <- do.call(FindTransferAnchors, anchor_args)
  elapsed <- proc.time()[["elapsed"]] - elapsed_start
  anchor_count <- nrow(as.data.frame(anchors@anchors))
  mark_progress("find_transfer_anchors_rpca", "end", paste0("elapsed_sec=", round(elapsed, 3), "; anchors=", anchor_count))
  data.frame(mode = opt$mode, completed = TRUE, n_anchors = anchor_count, elapsed_sec = round(elapsed, 3), stringsAsFactors = FALSE)
}

main <- function() {
  opt <- parse_args(commandArgs(trailingOnly = TRUE))
  if (is.null(opt$bridge_dir) || is.null(opt$outdir) || is.null(opt$mode)) {
    stop("--bridge-dir, --outdir, and --mode are required", call. = FALSE)
  }
  mode_dir <- file.path(opt$outdir, opt$mode)
  dir.create(mode_dir, recursive = TRUE, showWarnings = FALSE)
  progress_path <<- file.path(mode_dir, "progress.tsv")
  if (file.exists(progress_path)) file.remove(progress_path)
  utils::capture.output(sessionInfo(), file = file.path(mode_dir, "session_info.txt"))
  mark_progress("run", "start", opt$mode)
  mark_progress(
    "parameters",
    "info",
    paste0(
      "cells=", opt$max_reference_cells, "x", opt$max_query_cells,
      "; nfeatures=", opt$nfeatures,
      "; npcs=", opt$npcs,
      "; dims=1:", opt$dims,
      "; k.filter=", ifelse(is.na(opt$k_filter), "NA", opt$k_filter),
      "; nn.method=", opt$nn_method,
      "; n.trees=", opt$n_trees
    )
  )
  base <- prepare_base_objects(opt)
  object_sizes <- data.frame(
    object = c("reference", "query"),
    cells = c(ncol(base$reference), ncol(base$query)),
    genes = c(nrow(base$reference), nrow(base$query)),
    size_bytes = c(as.numeric(object.size(base$reference)), as.numeric(object.size(base$query))),
    stringsAsFactors = FALSE
  )
  write_tsv(object_sizes, file.path(mode_dir, "object_sizes.tsv"))
  result <- if (opt$mode == "integration_rpca") {
    run_integration_rpca(base$reference, base$query, opt, mode_dir)
  } else {
    run_transfer_rpca(base$reference, base$query, opt, mode_dir)
  }
  result$reference_cells <- ncol(base$reference)
  result$query_cells <- ncol(base$query)
  result$nfeatures <- opt$nfeatures
  result$npcs <- opt$npcs
  result$dims <- paste0("1:", opt$dims)
  result$k_filter <- ifelse(is.na(opt$k_filter), "NA", as.character(opt$k_filter))
  result$nn_method <- opt$nn_method
  write_tsv(result, file.path(mode_dir, "summary.tsv"))
  mark_progress("run", "end", "completed")
}

main()
