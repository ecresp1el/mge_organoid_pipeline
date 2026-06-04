#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
})

parse_args <- function(args) {
  out <- list()
  i <- 1
  while (i <= length(args)) {
    key <- args[[i]]
    if (startsWith(key, "--")) {
      key <- sub("^--", "", key)
      if (grepl("=", key, fixed = TRUE)) {
        parts <- strsplit(key, "=", fixed = TRUE)[[1]]
        out[[parts[[1]]]] <- if (length(parts) > 1) parts[[2]] else ""
      } else {
        out[[key]] <- if (i < length(args)) args[[i + 1]] else ""
        i <- i + 1
      }
    }
    i <- i + 1
  }
  out
}

arg_list <- parse_args(commandArgs(trailingOnly = TRUE))

get_arg <- function(name, default = "") {
  value <- Sys.getenv(name, unset = "")
  cli_name <- tolower(gsub("_", "-", name))
  cli_name <- sub("^cross-study-inventory-", "", cli_name)
  if (!is.null(arg_list[[cli_name]]) && nzchar(arg_list[[cli_name]])) {
    value <- arg_list[[cli_name]]
  }
  if (!nzchar(value)) default else value
}

timestamp <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")

log_msg <- function(...) {
  cat(sprintf("[%s] %s\n", timestamp(), paste0(..., collapse = " ")))
  flush.console()
}

collapse_csv <- function(x) {
  x <- x[!is.na(x) & nzchar(as.character(x))]
  if (length(x) == 0) "" else paste(as.character(x), collapse = ",")
}

safe_chr <- function(expr, default = "") {
  tryCatch(as.character(expr), error = function(e) default)
}

safe_int <- function(expr, default = NA_integer_) {
  tryCatch(as.integer(expr), error = function(e) default)
}

resolve_path <- function(path, project_root) {
  if (!nzchar(path)) return(path)
  if (startsWith(path, "/")) path else file.path(project_root, path)
}

write_tsv <- function(df, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.table(df, path, sep = "\t", row.names = FALSE, quote = FALSE, na = "NA")
}

class_string <- function(x) paste(class(x), collapse = ",")

unique_examples <- function(x, max_n = 5) {
  values <- unique(as.character(x[!is.na(x)]))
  collapse_csv(head(values, max_n))
}

count_values <- function(values) {
  values <- as.character(values)
  values[is.na(values) | !nzchar(values)] <- "NA"
  tab <- sort(table(values), decreasing = TRUE)
  data.frame(
    value = names(tab),
    n_cells = as.integer(tab),
    stringsAsFactors = FALSE
  )
}

project_root <- get_arg("PROJECT_ROOT", "/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder")
results_dirname <- get_arg("CROSS_STUDY_INVENTORY_RESULTS_DIRNAME", "cross_study_prediction_inventory")
run_label <- get_arg("CROSS_STUDY_INVENTORY_RUN_LABEL", "cross_study_prediction_inventory_v1")
shi_reference_rds <- get_arg(
  "SHI_REFERENCE_RDS",
  file.path(project_root, "results/shi_2019_paper_qc/shi_2019_seurat.rds")
)
div30_obs_path <- get_arg(
  "DIV30_SHI_TRANSFER_OBS",
  file.path(project_root, "results/shi_reference_div30_seurat_label_transfer/shi_reference_div30_seurat_label_transfer_v1/tables/div30_shi_seurat_label_transfer_obs.tsv.gz")
)
div90_obs_path <- get_arg(
  "DIV90_SHI_TRANSFER_OBS",
  file.path(project_root, "results/shi_reference_div90_seurat_label_transfer/shi_reference_div90_seurat_label_transfer_v1/tables/div90_shi_seurat_label_transfer_obs.tsv.gz")
)

run_dir <- file.path(project_root, "results", results_dirname, run_label)
table_dir <- file.path(run_dir, "tables")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

target_studies <- data.frame(
  study_id = c("varela_div30", "varela_div90", "walsh", "bershteyn_2025", "bershteyn_2023"),
  study_label = c("Varela DIV30", "Varela DIV90", "Walsh", "Bershteyn 2025", "Bershteyn 2023"),
  object_path = c(
    "results/varela_this_paper/varela_this_paper_seurat.rds",
    "/nfs/turbo/umms-parent/Manny_test/ventral_sosrs_output/umap_props_output/clustered_day90_with_cluster_names_2.rds",
    "results/walsh_day75/walsh_day75_final_annotated.rds",
    "results/bershteyn_2025/bershteyn_2025_seurat.rds",
    "results/bershteyn_2023/bershteyn_2023_seurat.rds"
  ),
  assay = c("RNA", "RNA", "RNA", "RNA", "RNA"),
  reduction = c("umap", "umap", "umap_sel", "umap", "umap"),
  stringsAsFactors = FALSE
)
target_studies$object_path_resolved <- vapply(target_studies$object_path, resolve_path, character(1), project_root = project_root)

log_msg("Project root: ", project_root)
log_msg("Run dir: ", run_dir)
log_msg("Loading Shi reference feature names: ", shi_reference_rds)
reference_features <- character()
reference_status <- "missing_reference"
if (file.exists(shi_reference_rds)) {
  reference_obj <- readRDS(shi_reference_rds)
  reference_assay <- if ("RNA" %in% Assays(reference_obj)) "RNA" else DefaultAssay(reference_obj)
  reference_features <- rownames(reference_obj[[reference_assay]])
  reference_status <- "ok"
  rm(reference_obj)
  gc()
}
log_msg("Reference status: ", reference_status, "; n_features=", length(reference_features))

object_rows <- list()
metadata_rows <- list()
value_count_rows <- list()
sample_count_rows <- list()
cluster_count_rows <- list()
readiness_rows <- list()

sample_priority <- c("sample_id", "orig.ident", "sample", "samples")
cluster_priority <- c("seurat_clusters", "cluster_id", "cluster_number_name", "predicted.GEcluster", "celltype", "process")
metadata_count_candidates <- c(
  "orig.ident", "sample_id", "sample", "samples", "domain", "type", "Phase", "old.ident",
  "seurat_clusters", "cluster_id", "cluster_number_name", "celltype", "process",
  "predicted.GEtype", "predicted.GEcluster", "predicted.GEgws", "predicted.macaclass", "predicted.musclass"
)

for (idx in seq_len(nrow(target_studies))) {
  study <- target_studies[idx, ]
  log_msg("Inspecting ", study$study_id, ": ", study$object_path_resolved)
  if (!file.exists(study$object_path_resolved)) {
    object_rows[[length(object_rows) + 1]] <- data.frame(
      study_id = study$study_id,
      study_label = study$study_label,
      status = "missing_object",
      object_path = study$object_path_resolved,
      object_class = "",
      default_assay = "",
      assays_available = "",
      assay = study$assay,
      assay_class = "",
      assay_layers_available = "",
      reductions_available = "",
      requested_reduction = study$reduction,
      n_cells = NA_integer_,
      n_features = NA_integer_,
      n_metadata_columns = NA_integer_,
      n_reference_features = length(reference_features),
      n_shared_reference_features = NA_integer_,
      shared_reference_fraction = NA_real_,
      stringsAsFactors = FALSE
    )
    next
  }

  obj <- readRDS(study$object_path_resolved)
  assays_available <- safe_chr(collapse_csv(Assays(obj)))
  default_assay <- safe_chr(DefaultAssay(obj))
  assay_use <- if (study$assay %in% Assays(obj)) study$assay else default_assay
  assay_obj <- obj[[assay_use]]
  features <- safe_chr(rownames(assay_obj))
  cells <- safe_chr(Cells(obj))
  reductions_available <- safe_chr(collapse_csv(Reductions(obj)))
  layer_names <- tryCatch(Layers(assay_obj), error = function(e) character())
  metadata <- obj@meta.data
  metadata_cols <- colnames(metadata)
  requested_reduction_ok <- study$reduction %in% Reductions(obj)
  shared_features <- intersect(features, reference_features)
  sample_col <- sample_priority[sample_priority %in% metadata_cols]
  sample_col <- if (length(sample_col) > 0) sample_col[[1]] else ""
  cluster_col <- cluster_priority[cluster_priority %in% metadata_cols]
  cluster_col <- if (length(cluster_col) > 0) cluster_col[[1]] else ""
  sample_cols_present <- intersect(sample_priority, metadata_cols)
  cluster_cols_present <- unique(c(
    intersect(cluster_priority, metadata_cols),
    grep("(_snn_res\\.|seurat_clusters|cluster|celltype|process|predicted\\.GE)", metadata_cols, value = TRUE)
  ))
  prediction_cols_present <- grep("shi.*prediction|prediction_score|predicted_shi", metadata_cols, value = TRUE, ignore.case = TRUE)

  object_rows[[length(object_rows) + 1]] <- data.frame(
    study_id = study$study_id,
    study_label = study$study_label,
    status = "ok",
    object_path = study$object_path_resolved,
    object_class = class_string(obj),
    default_assay = default_assay,
    assays_available = assays_available,
    assay = assay_use,
    assay_class = class_string(assay_obj),
    assay_layers_available = collapse_csv(layer_names),
    reductions_available = reductions_available,
    requested_reduction = study$reduction,
    n_cells = length(cells),
    n_features = length(features),
    n_metadata_columns = length(metadata_cols),
    n_reference_features = length(reference_features),
    n_shared_reference_features = length(shared_features),
    shared_reference_fraction = if (length(reference_features) > 0) length(shared_features) / length(reference_features) else NA_real_,
    stringsAsFactors = FALSE
  )

  readiness_missing <- character()
  if (!nzchar(assay_use)) readiness_missing <- c(readiness_missing, "RNA assay")
  if (!requested_reduction_ok) readiness_missing <- c(readiness_missing, paste0("requested reduction ", study$reduction))
  if (!nzchar(sample_col)) readiness_missing <- c(readiness_missing, "sample metadata column")
  if (!nzchar(cluster_col)) readiness_missing <- c(readiness_missing, "cluster metadata column")
  if (length(shared_features) < 1000) readiness_missing <- c(readiness_missing, ">=1000 shared Shi reference features")

  readiness_rows[[length(readiness_rows) + 1]] <- data.frame(
    study_id = study$study_id,
    study_label = study$study_label,
    ready_for_seurat_label_transfer = length(shared_features) >= 1000 && nzchar(assay_use),
    ready_for_sample_level_score_plots = nzchar(sample_col),
    ready_for_cluster_summaries = nzchar(cluster_col),
    primary_sample_col = sample_col,
    sample_cols_present = collapse_csv(sample_cols_present),
    primary_cluster_col = cluster_col,
    cluster_cols_present = collapse_csv(cluster_cols_present),
    requested_reduction = study$reduction,
    requested_reduction_ok = requested_reduction_ok,
    existing_shi_prediction_metadata_cols = collapse_csv(prediction_cols_present),
    n_shared_reference_features = length(shared_features),
    missing_for_requested_workflow = collapse_csv(readiness_missing),
    stringsAsFactors = FALSE
  )

  for (col in metadata_cols) {
    x <- metadata[[col]]
    metadata_rows[[length(metadata_rows) + 1]] <- data.frame(
      study_id = study$study_id,
      study_label = study$study_label,
      column_name = col,
      column_class = class_string(x),
      n_non_na = sum(!is.na(x)),
      n_unique_non_na = length(unique(x[!is.na(x)])),
      example_values = unique_examples(x, max_n = 5),
      stringsAsFactors = FALSE
    )
  }

  count_cols <- unique(c(
    intersect(metadata_count_candidates, metadata_cols),
    grep("(_snn_res\\.|seurat_clusters|cluster|celltype|process|predicted\\.GE|walsh_group)", metadata_cols, value = TRUE)
  ))
  for (col in count_cols) {
    x <- metadata[[col]]
    n_unique <- length(unique(x[!is.na(x)]))
    if (n_unique <= 200) {
      counts <- count_values(x)
      counts$study_id <- study$study_id
      counts$study_label <- study$study_label
      counts$column_name <- col
      value_count_rows[[length(value_count_rows) + 1]] <- counts[, c("study_id", "study_label", "column_name", "value", "n_cells")]
    }
  }

  if (nzchar(sample_col)) {
    counts <- count_values(metadata[[sample_col]])
    counts$study_id <- study$study_id
    counts$study_label <- study$study_label
    counts$sample_col <- sample_col
    sample_count_rows[[length(sample_count_rows) + 1]] <- counts[, c("study_id", "study_label", "sample_col", "value", "n_cells")]
  }
  if (nzchar(cluster_col)) {
    counts <- count_values(metadata[[cluster_col]])
    counts$study_id <- study$study_id
    counts$study_label <- study$study_label
    counts$cluster_col <- cluster_col
    cluster_count_rows[[length(cluster_count_rows) + 1]] <- counts[, c("study_id", "study_label", "cluster_col", "value", "n_cells")]
  }

  rm(obj, assay_obj, metadata)
  gc()
}

bind_or_empty <- function(rows) {
  if (length(rows) == 0) data.frame() else do.call(rbind, rows)
}

object_summary <- bind_or_empty(object_rows)
metadata_columns <- bind_or_empty(metadata_rows)
metadata_value_counts <- bind_or_empty(value_count_rows)
sample_counts <- bind_or_empty(sample_count_rows)
cluster_counts <- bind_or_empty(cluster_count_rows)
prediction_readiness <- bind_or_empty(readiness_rows)

write_tsv(object_summary, file.path(table_dir, "cross_study_prediction_object_summary.tsv"))
write_tsv(metadata_columns, file.path(table_dir, "cross_study_prediction_metadata_columns.tsv"))
write_tsv(metadata_value_counts, file.path(table_dir, "cross_study_prediction_metadata_value_counts.tsv"))
write_tsv(sample_counts, file.path(table_dir, "cross_study_prediction_primary_sample_counts.tsv"))
write_tsv(cluster_counts, file.path(table_dir, "cross_study_prediction_primary_cluster_counts.tsv"))
write_tsv(prediction_readiness, file.path(table_dir, "cross_study_prediction_readiness.tsv"))

prediction_table_summary_rows <- list()
prediction_sample_rows <- list()
prediction_cluster_rows <- list()
prediction_label_rows <- list()

prediction_tables <- data.frame(
  timepoint = c("DIV30", "DIV90"),
  path = c(div30_obs_path, div90_obs_path),
  stringsAsFactors = FALSE
)

for (idx in seq_len(nrow(prediction_tables))) {
  rec <- prediction_tables[idx, ]
  log_msg("Inspecting existing Shi-transfer obs table for ", rec$timepoint, ": ", rec$path)
  if (!file.exists(rec$path)) {
    prediction_table_summary_rows[[length(prediction_table_summary_rows) + 1]] <- data.frame(
      timepoint = rec$timepoint,
      status = "missing",
      path = rec$path,
      n_cells = NA_integer_,
      sample_col = "",
      n_samples = NA_integer_,
      cluster_col = "",
      n_clusters = NA_integer_,
      n_prediction_score_cols = NA_integer_,
      n_week_score_cols = NA_integer_,
      label_col = "",
      week_label_col = "",
      stringsAsFactors = FALSE
    )
    next
  }
  obs <- read.delim(gzfile(rec$path), sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  sample_col <- if ("orig.ident" %in% colnames(obs)) "orig.ident" else ""
  cluster_col <- if ("seurat_clusters" %in% colnames(obs)) "seurat_clusters" else ""
  score_cols <- grep("^shi_seurat_full_prediction_score", colnames(obs), value = TRUE)
  week_score_cols <- grep("^shi_seurat_full_week_prediction_score", colnames(obs), value = TRUE)
  label_col <- if ("shi_seurat_full_predicted_shi_label" %in% colnames(obs)) "shi_seurat_full_predicted_shi_label" else ""
  week_label_col <- if ("shi_seurat_full_predicted_shi_week_label" %in% colnames(obs)) "shi_seurat_full_predicted_shi_week_label" else ""

  prediction_table_summary_rows[[length(prediction_table_summary_rows) + 1]] <- data.frame(
    timepoint = rec$timepoint,
    status = "ok",
    path = rec$path,
    n_cells = nrow(obs),
    sample_col = sample_col,
    n_samples = if (nzchar(sample_col)) length(unique(obs[[sample_col]])) else NA_integer_,
    cluster_col = cluster_col,
    n_clusters = if (nzchar(cluster_col)) length(unique(obs[[cluster_col]])) else NA_integer_,
    n_prediction_score_cols = length(score_cols),
    n_week_score_cols = length(week_score_cols),
    label_col = label_col,
    week_label_col = week_label_col,
    stringsAsFactors = FALSE
  )
  if (nzchar(sample_col)) {
    counts <- count_values(obs[[sample_col]])
    counts$timepoint <- rec$timepoint
    counts$sample_col <- sample_col
    prediction_sample_rows[[length(prediction_sample_rows) + 1]] <- counts[, c("timepoint", "sample_col", "value", "n_cells")]
  }
  if (nzchar(cluster_col)) {
    counts <- count_values(obs[[cluster_col]])
    counts$timepoint <- rec$timepoint
    counts$cluster_col <- cluster_col
    prediction_cluster_rows[[length(prediction_cluster_rows) + 1]] <- counts[, c("timepoint", "cluster_col", "value", "n_cells")]
  }
  if (nzchar(label_col)) {
    counts <- count_values(obs[[label_col]])
    counts$timepoint <- rec$timepoint
    counts$label_col <- label_col
    prediction_label_rows[[length(prediction_label_rows) + 1]] <- counts[, c("timepoint", "label_col", "value", "n_cells")]
  }
  rm(obs)
  gc()
}

write_tsv(bind_or_empty(prediction_table_summary_rows), file.path(table_dir, "div30_div90_existing_shi_prediction_table_summary.tsv"))
write_tsv(bind_or_empty(prediction_sample_rows), file.path(table_dir, "div30_div90_existing_shi_prediction_sample_counts.tsv"))
write_tsv(bind_or_empty(prediction_cluster_rows), file.path(table_dir, "div30_div90_existing_shi_prediction_cluster_counts.tsv"))
write_tsv(bind_or_empty(prediction_label_rows), file.path(table_dir, "div30_div90_existing_shi_prediction_label_counts.tsv"))

completion <- data.frame(
  run_label = run_label,
  status = "complete",
  completed_at = timestamp(),
  n_target_studies = nrow(target_studies),
  n_object_summary_rows = nrow(object_summary),
  n_metadata_column_rows = nrow(metadata_columns),
  n_metadata_value_count_rows = nrow(metadata_value_counts),
  n_prediction_readiness_rows = nrow(prediction_readiness),
  stringsAsFactors = FALSE
)
write_tsv(completion, file.path(table_dir, "cross_study_prediction_inventory_complete.tsv"))

log_msg("Complete: ", run_dir)
