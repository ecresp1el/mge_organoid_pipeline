#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
})

timestamp <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")
log_msg <- function(...) message("[R ", timestamp(), "] ", paste0(..., collapse = ""))

parse_args <- function(args) {
  out <- list(
    study_id = NULL,
    label = NULL,
    seurat = NULL,
    h5ad = NULL,
    outdir = NULL
  )
  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, "--")) stop("Unknown argument: ", key, call. = FALSE)
    name <- substring(key, 3L)
    if (!(name %in% names(out))) stop("Unknown argument: ", key, call. = FALSE)
    if (i == length(args)) stop("Missing value for argument: ", key, call. = FALSE)
    out[[name]] <- args[[i + 1L]]
    i <- i + 2L
  }
  out
}

safe_class <- function(x) paste(class(x), collapse = ",")

safe_one_text <- function(x) {
  if (is.null(x) || length(x) == 0) return("")
  if (is.function(x)) return("<function>")
  if (is.environment(x)) return("<environment>")
  if (is.language(x)) return(gsub("[\r\n\t]+", " ", paste(deparse(x), collapse = " ")))
  if (isS4(x)) return(paste0("<S4:", safe_class(x), ">"))
  if (is.list(x) && !is.data.frame(x)) return(paste0("<list:", length(x), ">"))
  text <- tryCatch(as.character(x), error = function(e) paste0("<", safe_class(x), ">"))
  if (length(text) == 0) return("")
  text[is.na(text)] <- "<NA>"
  gsub("[\r\n\t]+", " ", text[[1]])
}

safe_text <- function(x, collapse = "; ") {
  if (is.null(x) || length(x) == 0) return("")
  if (is.list(x) && !is.data.frame(x)) {
    values <- vapply(x, safe_one_text, character(1))
  } else {
    values <- tryCatch(
      as.character(x),
      error = function(e) safe_one_text(x)
    )
  }
  if (length(values) == 0) return("")
  values[is.na(values)] <- "<NA>"
  values <- gsub("[\r\n\t]+", " ", values)
  paste(values, collapse = collapse)
}

safe_dim <- function(x) {
  d <- tryCatch(dim(x), error = function(e) NULL)
  if (is.null(d)) return(c(NA_integer_, NA_integer_))
  if (length(d) < 2) return(c(as.integer(d[[1]]), NA_integer_))
  c(as.integer(d[[1]]), as.integer(d[[2]]))
}

write_tsv <- function(df, path) {
  write.table(
    df,
    path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    na = ""
  )
}

slot_value <- function(obj, slot_name, default = NA_character_) {
  if (!(slot_name %in% slotNames(obj))) return(default)
  value <- tryCatch(slot(obj, slot_name), error = function(e) default)
  if (length(value) == 0) return(default)
  value
}

list_assay_layers <- function(assay_obj) {
  layers <- character(0)
  if (exists("Layers", where = asNamespace("SeuratObject"), inherits = FALSE)) {
    layers <- tryCatch(SeuratObject::Layers(assay_obj), error = function(e) character(0))
  }
  if (length(layers) == 0) {
    candidate_slots <- intersect(c("counts", "data", "scale.data"), slotNames(assay_obj))
    layers <- candidate_slots[vapply(candidate_slots, function(slot_name) {
      value <- tryCatch(slot(assay_obj, slot_name), error = function(e) NULL)
      !is.null(value) && length(value) > 0
    }, logical(1))]
  }
  unique(as.character(layers))
}

get_assay_layer_dim <- function(obj, assay_name, layer_name) {
  mat <- tryCatch(
    SeuratObject::GetAssayData(obj, assay = assay_name, layer = layer_name),
    error = function(e) NULL
  )
  if (is.null(mat)) return(c(NA_integer_, NA_integer_))
  safe_dim(mat)
}

summarize_top_values <- function(values, n = 10) {
  values <- as.character(values)
  values[is.na(values)] <- "<NA>"
  tab <- sort(table(values, useNA = "ifany"), decreasing = TRUE)
  if (length(tab) == 0) return("")
  top <- head(tab, n)
  safe_text(paste0(names(top), ": ", as.integer(top)))
}

opt <- parse_args(commandArgs(trailingOnly = TRUE))
required <- c("study_id", "seurat", "outdir")
missing <- required[vapply(required, function(name) is.null(opt[[name]]) || !nzchar(opt[[name]]), logical(1))]
if (length(missing) > 0) stop("Missing required argument(s): ", paste(missing, collapse = ", "), call. = FALSE)
if (!file.exists(opt$seurat)) stop("Seurat source does not exist: ", opt$seurat, call. = FALSE)

dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

log_msg("R executable: ", R.home("bin"))
log_msg("R version: ", paste(R.version$major, R.version$minor, sep = "."))
log_msg("Seurat version: ", as.character(utils::packageVersion("Seurat")))
log_msg("SeuratObject version: ", as.character(utils::packageVersion("SeuratObject")))
log_msg("Reading Seurat RDS: ", opt$seurat)

obj <- readRDS(opt$seurat)
if (!inherits(obj, "Seurat")) stop("Object is not a Seurat object: ", opt$seurat, call. = FALSE)
log_msg("Loaded Seurat object with ", ncol(obj), " cells and ", nrow(obj), " features")

assay_names <- Seurat::Assays(obj)
reduction_names <- Seurat::Reductions(obj)
graph_names <- if ("graphs" %in% slotNames(obj)) names(slot(obj, "graphs")) else character(0)
image_names <- if ("images" %in% slotNames(obj)) names(slot(obj, "images")) else character(0)
neighbor_names <- if ("neighbors" %in% slotNames(obj)) names(slot(obj, "neighbors")) else character(0)
misc_names <- if ("misc" %in% slotNames(obj)) names(slot(obj, "misc")) else character(0)
tool_names <- if ("tools" %in% slotNames(obj)) names(slot(obj, "tools")) else character(0)
command_names <- if ("commands" %in% slotNames(obj)) names(slot(obj, "commands")) else character(0)

object_summary <- data.frame(
  key = c(
    "study_id",
    "label",
    "seurat_path",
    "h5ad_path",
    "file_size_gb",
    "seurat_class",
    "n_cells",
    "n_features_default_assay",
    "default_assay",
    "project_name",
    "assays",
    "reductions",
    "graphs",
    "images",
    "neighbors",
    "misc_keys",
    "tool_keys",
    "commands",
    "meta_data_columns",
    "active_ident_levels"
  ),
  value = c(
    opt$study_id,
    ifelse(is.null(opt$label), "", opt$label),
    normalizePath(opt$seurat, mustWork = FALSE),
    ifelse(is.null(opt$h5ad), "", opt$h5ad),
    sprintf("%.3f", file.info(opt$seurat)$size / 1024^3),
    safe_class(obj),
    as.character(ncol(obj)),
    as.character(nrow(obj)),
    Seurat::DefaultAssay(obj),
    safe_text(slot_value(obj, "project.name", "")),
    safe_text(assay_names),
    safe_text(reduction_names),
    safe_text(graph_names),
    safe_text(image_names),
    safe_text(neighbor_names),
    safe_text(misc_names),
    safe_text(tool_names),
    safe_text(command_names),
    safe_text(colnames(obj@meta.data)),
    safe_text(levels(Seurat::Idents(obj)))
  ),
  stringsAsFactors = FALSE
)
write_tsv(object_summary, file.path(opt$outdir, "seurat_object_summary.tsv"))

assay_rows <- list()
assay_layer_rows <- list()
for (assay_name in assay_names) {
  assay_obj <- obj[[assay_name]]
  layers <- list_assay_layers(assay_obj)
  variable_features <- tryCatch(SeuratObject::VariableFeatures(obj, assay = assay_name), error = function(e) character(0))
  assay_rows[[length(assay_rows) + 1L]] <- data.frame(
    assay = assay_name,
    class = safe_class(assay_obj),
    n_features = nrow(assay_obj),
    n_cells = ncol(assay_obj),
    layers = safe_text(layers),
    n_variable_features = length(variable_features),
    stringsAsFactors = FALSE
  )
  if (length(variable_features) > 0) {
    write_tsv(
      data.frame(assay = assay_name, variable_feature = variable_features, stringsAsFactors = FALSE),
      file.path(opt$outdir, paste0("seurat_variable_features_", assay_name, ".tsv"))
    )
  }
  for (layer_name in layers) {
    layer_dim <- get_assay_layer_dim(obj, assay_name, layer_name)
    assay_layer_rows[[length(assay_layer_rows) + 1L]] <- data.frame(
      assay = assay_name,
      layer = layer_name,
      n_features = layer_dim[[1]],
      n_cells = layer_dim[[2]],
      stringsAsFactors = FALSE
    )
  }
}
write_tsv(do.call(rbind, assay_rows), file.path(opt$outdir, "seurat_assays.tsv"))
if (length(assay_layer_rows) > 0) {
  write_tsv(do.call(rbind, assay_layer_rows), file.path(opt$outdir, "seurat_assay_layers.tsv"))
}

reduction_rows <- list()
for (reduction_name in reduction_names) {
  reduction <- obj@reductions[[reduction_name]]
  emb <- tryCatch(SeuratObject::Embeddings(obj, reduction = reduction_name), error = function(e) NULL)
  loadings <- tryCatch(SeuratObject::Loadings(obj, reduction = reduction_name), error = function(e) NULL)
  emb_dim <- if (is.null(emb)) c(NA_integer_, NA_integer_) else safe_dim(emb)
  loadings_dim <- if (is.null(loadings)) c(NA_integer_, NA_integer_) else safe_dim(loadings)
  stdev <- slot_value(reduction, "stdev", numeric(0))
  reduction_rows[[length(reduction_rows) + 1L]] <- data.frame(
    reduction = reduction_name,
    class = safe_class(reduction),
    key = safe_text(slot_value(reduction, "key", "")),
    assay_used = safe_text(slot_value(reduction, "assay.used", "")),
    n_cells = emb_dim[[1]],
    n_components = emb_dim[[2]],
    n_loading_features = loadings_dim[[1]],
    n_loading_components = loadings_dim[[2]],
    n_stdev = length(stdev),
    stringsAsFactors = FALSE
  )
}
if (length(reduction_rows) > 0) {
  write_tsv(do.call(rbind, reduction_rows), file.path(opt$outdir, "seurat_reductions.tsv"))
}

if (length(graph_names) > 0) {
  graph_rows <- lapply(graph_names, function(graph_name) {
    graph <- obj@graphs[[graph_name]]
    graph_dim <- safe_dim(graph)
    data.frame(
      graph = graph_name,
      class = safe_class(graph),
      n_rows = graph_dim[[1]],
      n_cols = graph_dim[[2]],
      stringsAsFactors = FALSE
    )
  })
  write_tsv(do.call(rbind, graph_rows), file.path(opt$outdir, "seurat_graphs.tsv"))
}

meta <- obj@meta.data
metadata_rows <- list()
metadata_top_rows <- list()
for (column in colnames(meta)) {
  values <- meta[[column]]
  n_missing <- sum(is.na(values))
  n_unique <- length(unique(values[!is.na(values)]))
  numeric_values <- suppressWarnings(as.numeric(values))
  numeric_nonmissing <- sum(!is.na(numeric_values))
  is_numeric_like <- is.numeric(values) || (numeric_nonmissing > 0 && numeric_nonmissing == sum(!is.na(values)))
  metadata_rows[[length(metadata_rows) + 1L]] <- data.frame(
    column = column,
    class = safe_class(values),
    n_missing = as.integer(n_missing),
    n_unique = as.integer(n_unique),
    is_numeric_like = is_numeric_like,
    min = if (is_numeric_like && numeric_nonmissing > 0) min(numeric_values, na.rm = TRUE) else NA_real_,
    median = if (is_numeric_like && numeric_nonmissing > 0) median(numeric_values, na.rm = TRUE) else NA_real_,
    mean = if (is_numeric_like && numeric_nonmissing > 0) mean(numeric_values, na.rm = TRUE) else NA_real_,
    max = if (is_numeric_like && numeric_nonmissing > 0) max(numeric_values, na.rm = TRUE) else NA_real_,
    top_values = summarize_top_values(values, n = 10),
    stringsAsFactors = FALSE
  )
  top <- sort(table(as.character(values), useNA = "ifany"), decreasing = TRUE)
  if (length(top) > 0) {
    top <- head(top, 25)
    metadata_top_rows[[length(metadata_top_rows) + 1L]] <- data.frame(
      column = column,
      value = names(top),
      count = as.integer(top),
      stringsAsFactors = FALSE
    )
  }
}
write_tsv(do.call(rbind, metadata_rows), file.path(opt$outdir, "seurat_metadata_columns.tsv"))
if (length(metadata_top_rows) > 0) {
  write_tsv(do.call(rbind, metadata_top_rows), file.path(opt$outdir, "seurat_metadata_top_values.tsv"))
}

idents <- Seurat::Idents(obj)
ident_counts <- sort(table(as.character(idents), useNA = "ifany"), decreasing = TRUE)
write_tsv(
  data.frame(identity = names(ident_counts), n_cells = as.integer(ident_counts), stringsAsFactors = FALSE),
  file.path(opt$outdir, "seurat_identity_counts.tsv")
)

if (length(command_names) > 0) {
  command_rows <- lapply(command_names, function(command_name) {
    command <- obj@commands[[command_name]]
    params <- slot_value(command, "params", list())
    param_text <- ""
    if (length(params) > 0) {
      param_text <- safe_text(paste0(names(params), "=", vapply(params, function(x) safe_text(x, collapse = ","), character(1))))
    }
    data.frame(
      command = command_name,
      class = safe_class(command),
      name = safe_text(slot_value(command, "name", "")),
      time_stamp = safe_text(slot_value(command, "time.stamp", "")),
      assay_used = safe_text(slot_value(command, "assay.used", "")),
      call_string = safe_text(slot_value(command, "call.string", "")),
      params = param_text,
      stringsAsFactors = FALSE
    )
  })
  write_tsv(do.call(rbind, command_rows), file.path(opt$outdir, "seurat_commands.tsv"))
}

if (length(misc_names) > 0) {
  misc_rows <- lapply(misc_names, function(name) {
    value <- obj@misc[[name]]
    d <- safe_dim(value)
    data.frame(key = name, class = safe_class(value), n_rows = d[[1]], n_cols = d[[2]], stringsAsFactors = FALSE)
  })
  write_tsv(do.call(rbind, misc_rows), file.path(opt$outdir, "seurat_misc.tsv"))
}

if (length(tool_names) > 0) {
  tool_rows <- lapply(tool_names, function(name) {
    value <- obj@tools[[name]]
    d <- safe_dim(value)
    data.frame(key = name, class = safe_class(value), n_rows = d[[1]], n_cols = d[[2]], stringsAsFactors = FALSE)
  })
  write_tsv(do.call(rbind, tool_rows), file.path(opt$outdir, "seurat_tools.tsv"))
}

write_tsv(
  data.frame(
    study_id = opt$study_id,
    completed_at = timestamp(),
    outdir = normalizePath(opt$outdir, mustWork = FALSE),
    stringsAsFactors = FALSE
  ),
  file.path(opt$outdir, "seurat_inventory_complete.tsv")
)

log_msg("Finished Seurat inventory: ", opt$outdir)
