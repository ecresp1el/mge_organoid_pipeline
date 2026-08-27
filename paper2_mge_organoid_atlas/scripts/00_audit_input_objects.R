#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, width = 160)

parse_args <- function(args) {
  out <- list()
  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, "--") || i == length(args)) {
      stop("Expected --key value arguments; got: ", key)
    }
    out[[substring(key, 3L)]] <- args[[i + 1L]]
    i <- i + 2L
  }
  out
}

write_tsv <- function(x, path) {
  write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
}

write_tsv_gz <- function(x, path) {
  connection <- gzfile(path, "wt")
  on.exit(close(connection), add = TRUE)
  write.table(x, connection, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
}

collapse_examples <- function(x, limit = 5L) {
  x <- unique(as.character(x[!is.na(x)]))
  if (!length(x)) return("")
  paste(utils::head(x, limit), collapse = " | ")
}

split_priority <- function(value) {
  values <- trimws(strsplit(as.character(value), ";", fixed = TRUE)[[1L]])
  values[nzchar(values)]
}

first_available_column <- function(priority, available) {
  found <- priority[priority %in% available]
  if (length(found)) found[[1L]] else ""
}

display_cluster_label <- function(cluster_id, cluster_name) {
  cluster_id <- as.character(cluster_id)
  cluster_name <- as.character(cluster_name)
  if (!nzchar(cluster_name) || cluster_name == cluster_id) return(cluster_id)
  id_pattern <- paste0("^", gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", cluster_id), "([[:space:]]*[-.:)]|[[:space:]]+)")
  if (grepl(id_pattern, cluster_name)) cluster_name else paste(cluster_id, cluster_name, sep = " - ")
}

sha256_file <- function(path) {
  result <- system2("sha256sum", args = shQuote(path), stdout = TRUE, stderr = TRUE)
  status <- attr(result, "status")
  if (!is.null(status) && status != 0L) {
    stop("sha256sum failed for ", path, ": ", paste(result, collapse = " "))
  }
  strsplit(result[[1L]], "[[:space:]]+")[[1L]][[1L]]
}

assay_layers <- function(assay_object) {
  if ("Layers" %in% getNamespaceExports("SeuratObject")) {
    layers <- tryCatch(SeuratObject::Layers(assay_object), error = function(e) character())
    if (length(layers)) return(layers)
  }
  candidates <- intersect(c("counts", "data", "scale.data"), slotNames(assay_object))
  candidates[vapply(candidates, function(slot_name) {
    value <- slot(assay_object, slot_name)
    length(dim(value)) == 2L && all(dim(value) > 0L)
  }, logical(1))]
}

layer_dimensions <- function(assay_object, layer_name) {
  value <- tryCatch({
    if (inherits(assay_object, "Assay5")) {
      SeuratObject::LayerData(assay_object, layer = layer_name)
    } else {
      slot(assay_object, layer_name)
    }
  }, error = function(e) NULL)
  if (is.null(value) || length(dim(value)) != 2L) return(c(NA_integer_, NA_integer_))
  as.integer(dim(value))
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
manifest_path <- args[["manifest"]]
cluster_overrides_path <- args[["cluster-overrides"]]
outdir <- args[["outdir"]]
if (is.null(manifest_path) || is.null(cluster_overrides_path) || is.null(outdir)) {
  stop("Usage: 00_audit_input_objects.R --manifest INPUT.tsv --cluster-overrides OVERRIDES.tsv --outdir RUN_DIR")
}

dir.create(file.path(outdir, "tables", "features"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(outdir, "tables", "umap", "per_study"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(outdir, "provenance"), recursive = TRUE, showWarnings = FALSE)

manifest <- read.delim(manifest_path, check.names = FALSE)
required_columns <- c(
  "study_id", "display_name", "object_path", "input_role", "selection_status",
  "umap_reduction", "cluster_id_column", "cluster_name_priority", "sample_column_priority", "selection_note"
)
missing_columns <- setdiff(required_columns, colnames(manifest))
if (length(missing_columns)) stop("Manifest is missing columns: ", paste(missing_columns, collapse = ", "))
if (nrow(manifest) != 6L) stop("Paper 2 input registry must contain exactly six studies; found ", nrow(manifest))
if (anyDuplicated(manifest$study_id)) stop("study_id values must be unique")

cluster_overrides <- read.delim(cluster_overrides_path, check.names = FALSE)
override_columns <- c("study_id", "raw_cluster_id", "display_cluster_id", "display_cluster_name", "mapping_source")
missing_override_columns <- setdiff(override_columns, colnames(cluster_overrides))
if (length(missing_override_columns)) {
  stop("Cluster override table is missing columns: ", paste(missing_override_columns, collapse = ", "))
}
cluster_overrides$raw_cluster_id <- as.character(cluster_overrides$raw_cluster_id)

manifest$file_exists <- file.exists(manifest$object_path)
if (!all(manifest$file_exists)) {
  stop("Missing registered objects: ", paste(manifest$object_path[!manifest$file_exists], collapse = ", "))
}

info <- file.info(manifest$object_path)
manifest$file_size_bytes <- as.numeric(info$size)
manifest$file_modified <- format(info$mtime, "%Y-%m-%dT%H:%M:%S%z")
manifest$sha256 <- vapply(manifest$object_path, sha256_file, character(1))
write_tsv(manifest, file.path(outdir, "tables", "input_objects_audited.tsv"))

object_rows <- list()
layer_rows <- list()
metadata_rows <- list()
reduction_rows <- list()
umap_source_rows <- list()
umap_cluster_rows <- list()
feature_sets <- list()

for (index in seq_len(nrow(manifest))) {
  study_id <- manifest$study_id[[index]]
  object_path <- manifest$object_path[[index]]
  message("Reading ", study_id, ": ", object_path)
  object <- readRDS(object_path)
  if (!inherits(object, "Seurat")) {
    stop("Registered object is not a Seurat object: ", study_id, " (", paste(class(object), collapse = "/"), ")")
  }

  default_assay <- SeuratObject::DefaultAssay(object)
  assays <- names(object@assays)
  metadata <- object[[]]
  reductions <- names(object@reductions)
  features <- rownames(object[[default_assay]])
  feature_sets[[study_id]] <- unique(features)

  feature_connection <- gzfile(file.path(outdir, "tables", "features", paste0(study_id, "_default_assay_features.tsv.gz")), "wt")
  writeLines(c("feature", features), feature_connection)
  close(feature_connection)

  object_rows[[length(object_rows) + 1L]] <- data.frame(
    study_id = study_id,
    display_name = manifest$display_name[[index]],
    object_class = paste(class(object), collapse = "/"),
    n_features_default_assay = nrow(object),
    n_cells = ncol(object),
    default_assay = default_assay,
    assays = paste(assays, collapse = " | "),
    n_metadata_columns = ncol(metadata),
    reductions = paste(reductions, collapse = " | "),
    stringsAsFactors = FALSE
  )

  for (assay_name in assays) {
    assay_object <- object[[assay_name]]
    layers <- assay_layers(assay_object)
    if (!length(layers)) {
      layer_rows[[length(layer_rows) + 1L]] <- data.frame(
        study_id = study_id, assay = assay_name, assay_class = paste(class(assay_object), collapse = "/"),
        layer = NA_character_, n_features = nrow(assay_object), n_cells = ncol(assay_object), stringsAsFactors = FALSE
      )
    } else {
      for (layer_name in layers) {
        dims <- layer_dimensions(assay_object, layer_name)
        layer_rows[[length(layer_rows) + 1L]] <- data.frame(
          study_id = study_id, assay = assay_name, assay_class = paste(class(assay_object), collapse = "/"),
          layer = layer_name, n_features = dims[[1L]], n_cells = dims[[2L]], stringsAsFactors = FALSE
        )
      }
    }
  }

  for (column_name in colnames(metadata)) {
    values <- metadata[[column_name]]
    metadata_rows[[length(metadata_rows) + 1L]] <- data.frame(
      study_id = study_id,
      column = column_name,
      class = paste(class(values), collapse = "/"),
      n_missing = sum(is.na(values)),
      fraction_missing = mean(is.na(values)),
      n_unique_nonmissing = length(unique(values[!is.na(values)])),
      example_values = collapse_examples(values),
      stringsAsFactors = FALSE
    )
  }

  if (!length(reductions)) {
    reduction_rows[[length(reduction_rows) + 1L]] <- data.frame(
      study_id = study_id, reduction = NA_character_, reduction_class = NA_character_,
      n_cells = NA_integer_, n_dimensions = NA_integer_, stringsAsFactors = FALSE
    )
  } else {
    for (reduction_name in reductions) {
      reduction_object <- object[[reduction_name]]
      embeddings <- SeuratObject::Embeddings(reduction_object)
      reduction_rows[[length(reduction_rows) + 1L]] <- data.frame(
        study_id = study_id,
        reduction = reduction_name,
        reduction_class = paste(class(reduction_object), collapse = "/"),
        n_cells = nrow(embeddings),
        n_dimensions = ncol(embeddings),
        stringsAsFactors = FALSE
      )
    }
  }

  requested_reduction <- as.character(manifest$umap_reduction[[index]])
  reduction_used <- requested_reduction
  reduction_resolution <- "requested_reduction"
  if (!requested_reduction %in% reductions) {
    fallback <- grep("umap", reductions, ignore.case = TRUE, value = TRUE)
    if (!length(fallback)) {
      stop(study_id, " has no requested or fallback UMAP reduction; requested=", requested_reduction)
    }
    reduction_used <- fallback[[1L]]
    reduction_resolution <- "fallback_first_umap_named_reduction"
  }
  embeddings <- SeuratObject::Embeddings(object[[reduction_used]])
  if (ncol(embeddings) < 2L) stop(study_id, " UMAP reduction has fewer than two dimensions")
  cell_ids <- rownames(embeddings)
  if (is.null(cell_ids) || !all(cell_ids %in% rownames(metadata))) {
    stop(study_id, " UMAP cell IDs do not align to Seurat metadata row names")
  }
  metadata_umap <- metadata[cell_ids, , drop = FALSE]

  requested_cluster_column <- as.character(manifest$cluster_id_column[[index]])
  if (requested_cluster_column %in% colnames(metadata_umap)) {
    raw_cluster_id <- as.character(metadata_umap[[requested_cluster_column]])
    cluster_id_source <- requested_cluster_column
  } else {
    raw_cluster_id <- as.character(SeuratObject::Idents(object)[cell_ids])
    cluster_id_source <- "active.ident"
  }
  if (anyNA(raw_cluster_id) || any(!nzchar(raw_cluster_id))) {
    stop(study_id, " has missing raw cluster IDs in ", cluster_id_source)
  }

  cluster_name_column <- first_available_column(
    split_priority(manifest$cluster_name_priority[[index]]),
    colnames(metadata_umap)
  )
  if (nzchar(cluster_name_column)) {
    display_cluster_name <- as.character(metadata_umap[[cluster_name_column]])
    cluster_label_source <- paste0("metadata:", cluster_name_column)
  } else {
    display_cluster_name <- raw_cluster_id
    cluster_label_source <- "raw_cluster_id_no_name_column_available"
  }
  display_cluster_id <- raw_cluster_id
  cluster_label_source_values <- rep(cluster_label_source, length(raw_cluster_id))

  study_overrides <- cluster_overrides[cluster_overrides$study_id == study_id, , drop = FALSE]
  if (nrow(study_overrides)) {
    override_index <- match(raw_cluster_id, study_overrides$raw_cluster_id)
    if (anyNA(override_index)) {
      missing_raw <- sort(unique(raw_cluster_id[is.na(override_index)]))
      stop(study_id, " cluster override table does not cover raw clusters: ", paste(missing_raw, collapse = ", "))
    }
    display_cluster_id <- as.character(study_overrides$display_cluster_id[override_index])
    display_cluster_name <- as.character(study_overrides$display_cluster_name[override_index])
    cluster_label_source_values <- as.character(study_overrides$mapping_source[override_index])
  }
  display_cluster_name[is.na(display_cluster_name) | !nzchar(display_cluster_name)] <- raw_cluster_id[is.na(display_cluster_name) | !nzchar(display_cluster_name)]
  display_cluster_labels <- mapply(display_cluster_label, display_cluster_id, display_cluster_name, USE.NAMES = FALSE)

  sample_column <- first_available_column(
    split_priority(manifest$sample_column_priority[[index]]),
    colnames(metadata_umap)
  )
  if (nzchar(sample_column)) {
    sample_id <- as.character(metadata_umap[[sample_column]])
    sample_source <- paste0("metadata:", sample_column)
  } else {
    sample_id <- rep(study_id, length(cell_ids))
    sample_source <- "study_id_fallback"
  }

  plot_transform <- if (study_id == "varela_div90") {
    "UMAP1_plot=UMAP1_original;UMAP2_plot=-UMAP2_original"
  } else {
    "UMAP1_plot=UMAP1_original;UMAP2_plot=UMAP2_original"
  }
  umap_table <- data.frame(
    study_id = study_id,
    display_name = manifest$display_name[[index]],
    cell_id_original = cell_ids,
    atlas_candidate_cell_id = paste(study_id, cell_ids, sep = "__"),
    sample_id = sample_id,
    sample_source = sample_source,
    requested_reduction = requested_reduction,
    reduction_used = reduction_used,
    UMAP1_original = as.numeric(embeddings[, 1L]),
    UMAP2_original = as.numeric(embeddings[, 2L]),
    plot_transform = plot_transform,
    raw_cluster_id = raw_cluster_id,
    cluster_id_source = cluster_id_source,
    display_cluster_id = display_cluster_id,
    display_cluster_name = display_cluster_name,
    display_cluster_label = display_cluster_labels,
    cluster_label_source = cluster_label_source_values,
    stringsAsFactors = FALSE
  )
  write_tsv_gz(
    umap_table,
    file.path(outdir, "tables", "umap", "per_study", paste0(study_id, "_umap_cluster_inventory.tsv.gz"))
  )

  count_columns <- c(
    "study_id", "display_name", "raw_cluster_id", "cluster_id_source",
    "display_cluster_id", "display_cluster_name", "display_cluster_label", "cluster_label_source"
  )
  umap_counts <- aggregate(
    rep(1L, nrow(umap_table)),
    by = umap_table[count_columns],
    FUN = sum
  )
  colnames(umap_counts)[ncol(umap_counts)] <- "n_cells"
  umap_cluster_rows[[length(umap_cluster_rows) + 1L]] <- umap_counts
  umap_source_rows[[length(umap_source_rows) + 1L]] <- data.frame(
    study_id = study_id,
    display_name = manifest$display_name[[index]],
    requested_reduction = requested_reduction,
    reduction_used = reduction_used,
    reduction_resolution = reduction_resolution,
    n_cells_with_umap = nrow(umap_table),
    raw_cluster_id_source = cluster_id_source,
    cluster_name_column_available = if (nzchar(cluster_name_column)) cluster_name_column else "",
    cluster_override_applied = nrow(study_overrides) > 0L,
    n_raw_clusters = length(unique(raw_cluster_id)),
    n_display_clusters = length(unique(display_cluster_labels)),
    display_cluster_labels = paste(sort(unique(display_cluster_labels)), collapse = " | "),
    sample_source = sample_source,
    plot_transform = plot_transform,
    all_cells_retained_for_audit_plot = TRUE,
    stringsAsFactors = FALSE
  )

  rm(object, metadata, metadata_umap, features, embeddings, umap_table)
  invisible(gc())
}

object_table <- do.call(rbind, object_rows)
layer_table <- do.call(rbind, layer_rows)
metadata_table <- do.call(rbind, metadata_rows)
reduction_table <- do.call(rbind, reduction_rows)
umap_source_table <- do.call(rbind, umap_source_rows)
umap_cluster_table <- do.call(rbind, umap_cluster_rows)

write_tsv(object_table, file.path(outdir, "tables", "object_summary.tsv"))
write_tsv(layer_table, file.path(outdir, "tables", "assay_layer_inventory.tsv"))
write_tsv(metadata_table, file.path(outdir, "tables", "metadata_column_inventory.tsv"))
write_tsv(reduction_table, file.path(outdir, "tables", "reduction_inventory.tsv"))
write_tsv(umap_source_table, file.path(outdir, "tables", "umap_source_inventory.tsv"))
write_tsv(umap_cluster_table, file.path(outdir, "tables", "umap_cluster_counts.tsv"))

pair_rows <- list()
study_ids <- names(feature_sets)
for (i in seq_along(study_ids)) {
  for (j in i:length(study_ids)) {
    left <- study_ids[[i]]
    right <- study_ids[[j]]
    overlap <- length(intersect(feature_sets[[left]], feature_sets[[right]]))
    pair_rows[[length(pair_rows) + 1L]] <- data.frame(
      study_id_1 = left,
      study_id_2 = right,
      n_features_1 = length(feature_sets[[left]]),
      n_features_2 = length(feature_sets[[right]]),
      n_shared = overlap,
      fraction_of_study_1 = overlap / length(feature_sets[[left]]),
      fraction_of_study_2 = overlap / length(feature_sets[[right]]),
      stringsAsFactors = FALSE
    )
  }
}
write_tsv(do.call(rbind, pair_rows), file.path(outdir, "tables", "feature_overlap_pairwise.tsv"))

shared_all <- Reduce(intersect, feature_sets)
union_all <- Reduce(union, feature_sets)
write_tsv(data.frame(
  n_studies = length(feature_sets),
  n_feature_union = length(union_all),
  n_features_shared_all = length(shared_all),
  stringsAsFactors = FALSE
), file.path(outdir, "tables", "feature_overlap_summary.tsv"))

shared_connection <- gzfile(file.path(outdir, "tables", "features", "features_shared_all_six_studies.tsv.gz"), "wt")
writeLines(c("feature", sort(shared_all)), shared_connection)
close(shared_connection)

writeLines(capture.output(sessionInfo()), file.path(outdir, "provenance", "sessionInfo.txt"))
message("Audit tables written to ", outdir)
