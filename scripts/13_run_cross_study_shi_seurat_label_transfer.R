#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
})

RUN_LABEL_DEFAULT <- "cross_study_shi_seurat_label_transfer_v1"
PROJECT_ROOT_DEFAULT <- "/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder"
REFERENCE_DEFAULT <- file.path(PROJECT_ROOT_DEFAULT, "results/shi_2019_paper_qc/shi_2019_seurat.rds")
MAJOR_LABELS <- c(
  "MGE", "LGE", "CGE", "progenitor", "Excitatory IPC", "Excitatory neuron",
  "Thalamic neurons", "Microglia", "OPC", "Endothelial"
)

timestamp <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")
log_msg <- function(...) message("[", timestamp(), "] ", paste0(..., collapse = ""))

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  out <- list(
    project_root = Sys.getenv("PROJECT_ROOT", PROJECT_ROOT_DEFAULT),
    reference = Sys.getenv("SHI_REFERENCE_RDS", REFERENCE_DEFAULT),
    outdir = "",
    run_label = Sys.getenv("CROSS_STUDY_SHI_RUN_LABEL", RUN_LABEL_DEFAULT),
    study_id = character(),
    reference_label_col = Sys.getenv("SHI_REFERENCE_LABEL_COL", ""),
    reference_week_col = Sys.getenv("SHI_REFERENCE_WEEK_COL", ""),
    reference_labels_tsv = Sys.getenv(
      "SHI_REFERENCE_LABELS_TSV",
      file.path(
        Sys.getenv("PROJECT_ROOT", PROJECT_ROOT_DEFAULT),
        "results/shi_reference_div30_seurat_label_transfer/shi_reference_div30_seurat_label_transfer_v1/seurat/shi_reference_labels_for_seurat.tsv"
      )
    ),
    dims = "50",
    min_shared_features = "500",
    normalization_method = "LogNormalize",
    force_rerun = Sys.getenv("CROSS_STUDY_SHI_FORCE_RERUN", "false"),
    reuse_existing = Sys.getenv("CROSS_STUDY_SHI_REUSE_EXISTING", "true"),
    seed = "0"
  )
  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    value <- if (i < length(args)) args[[i + 1L]] else ""
    if (key == "--project-root") {
      out$project_root <- value; i <- i + 2L
    } else if (key == "--reference") {
      out$reference <- value; i <- i + 2L
    } else if (key == "--outdir") {
      out$outdir <- value; i <- i + 2L
    } else if (key == "--run-label") {
      out$run_label <- value; i <- i + 2L
    } else if (key == "--study-id") {
      out$study_id <- c(out$study_id, unlist(strsplit(value, "[,;[:space:]]+"))); i <- i + 2L
    } else if (key == "--reference-label-col") {
      out$reference_label_col <- value; i <- i + 2L
    } else if (key == "--reference-week-col") {
      out$reference_week_col <- value; i <- i + 2L
    } else if (key == "--reference-labels-tsv") {
      out$reference_labels_tsv <- value; i <- i + 2L
    } else if (key == "--dims") {
      out$dims <- value; i <- i + 2L
    } else if (key == "--min-shared-features") {
      out$min_shared_features <- value; i <- i + 2L
    } else if (key == "--normalization-method") {
      out$normalization_method <- value; i <- i + 2L
    } else if (key == "--force-rerun") {
      out$force_rerun <- value; i <- i + 2L
    } else if (key == "--reuse-existing") {
      out$reuse_existing <- value; i <- i + 2L
    } else if (key == "--seed") {
      out$seed <- value; i <- i + 2L
    } else {
      stop("Unknown argument: ", key, call. = FALSE)
    }
  }
  out$study_id <- out$study_id[nzchar(out$study_id)]
  if (!nzchar(out$outdir)) {
    out$outdir <- file.path(
      out$project_root,
      "results/cross_study_shi_seurat_label_transfer",
      out$run_label
    )
  }
  out
}

to_bool <- function(x, default = FALSE) {
  if (is.null(x) || !nzchar(x)) return(default)
  tolower(trimws(as.character(x))) %in% c("1", "true", "t", "yes", "y")
}

resolve_path <- function(path, project_root) {
  if (!nzchar(path)) return("")
  if (grepl("^/", path)) return(path)
  file.path(project_root, path)
}

write_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
}

write_tsv_gz <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  con <- gzfile(path, open = "wt")
  on.exit(close(con), add = TRUE)
  utils::write.table(x, con, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
}

read_tsv_maybe_gz <- function(path) {
  con <- if (grepl("\\.gz$", path, ignore.case = TRUE)) gzfile(path, "rt") else file(path, "rt")
  on.exit(close(con), add = TRUE)
  utils::read.delim(con, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
}

sanitize_token <- function(x) {
  y <- gsub("[^A-Za-z0-9]+", "_", trimws(as.character(x)))
  y <- gsub("_+", "_", y)
  y <- gsub("^_|_$", "", y)
  ifelse(nzchar(y), y, "value")
}

parse_gw_numeric <- function(x) {
  text <- as.character(x)
  out <- rep(NA_real_, length(text))
  for (i in seq_along(text)) {
    m <- gregexpr("[0-9]+", text[[i]], perl = TRUE)[[1]]
    if (length(m) > 0 && m[[1]] > 0) {
      nums <- as.numeric(regmatches(text[[i]], list(m))[[1]])
      out[[i]] <- mean(nums, na.rm = TRUE)
    }
  }
  out
}

default_studies <- function(project_root) {
  data.frame(
    study_id = c(
      "varela_div30", "varela_div90", "siebert_2026", "walsh",
      "bershteyn_2025", "bershteyn_2023", "samarasinghe_2021"
    ),
    study_label = c(
      "This Study, DIV 30", "This Study, DIV 90", "Siebert et al. 2026",
      "Walsh et al. 2025", "Bershteyn et al. 2025",
      "Bershteyn et al. 2023", "Samarasinghe et al. 2021"
    ),
    object_path = c(
      "results/varela_this_paper/varela_this_paper_seurat.rds",
      "/nfs/turbo/umms-parent/Manny_test/ventral_sosrs_output/umap_props_output/clustered_day90_with_cluster_names_2.rds",
      "results/siebert_2026/siebert_2026_seurat.rds",
      "results/walsh_day75/walsh_day75_final_annotated.rds",
      "results/bershteyn_2025/bershteyn_2025_seurat.rds",
      "results/bershteyn_2023/bershteyn_2023_seurat.rds",
      "results/samarasinghe_2021_zenodo_processed_object/samarasinghe_2021_zenodo_seurat.rds"
    ),
    reduction = c("umap", "umap", "umap", "umap_sel", "umap", "umap", "umap"),
    sample_col = c("orig.ident", "orig.ident", "orig.ident|sample", "sample_id", "sample", "orig.ident", "orig.ident"),
    cluster_col = rep("seurat_clusters", 7),
    existing_obs_path = c(
      file.path(project_root, "results/shi_reference_div30_seurat_label_transfer/shi_reference_div30_seurat_label_transfer_v1/tables/div30_shi_seurat_label_transfer_obs.tsv.gz"),
      file.path(project_root, "results/shi_reference_div90_seurat_label_transfer/shi_reference_div90_seurat_label_transfer_v1/tables/div90_shi_seurat_label_transfer_obs.tsv.gz"),
      "", "", "", "", ""
    ),
    existing_seurat_dir = c(
      file.path(project_root, "results/shi_reference_div30_seurat_label_transfer/shi_reference_div30_seurat_label_transfer_v1/seurat"),
      file.path(project_root, "results/shi_reference_div90_seurat_label_transfer/shi_reference_div90_seurat_label_transfer_v1/seurat"),
      "", "", "", "", ""
    ),
    existing_prefix = c("div30", "div90", "", "", "", "", ""),
    stringsAsFactors = FALSE
  )
}

sample_label_values <- function(study_id, sample) {
  sample <- as.character(sample)
  if (identical(study_id, "walsh")) {
    out <- sample
    out[out == "GSM7979671"] <- "MEL1_dFB_d75"
    out[out == "GSM7979672"] <- "MEL1_vFB_d75"
    return(out)
  }
  if (identical(study_id, "siebert_2026")) {
    out <- sample
    out[grepl("^Young_", out)] <- paste0(out[grepl("^Young_", out)], " (younger group, likely DIV51)")
    out[grepl("^Old_", out)] <- paste0(out[grepl("^Old_", out)], " (older group, likely DIV164)")
    return(out)
  }
  sample
}

join_layers_if_needed <- function(obj, assay) {
  if (
    inherits(obj[[assay]], "Assay5") &&
      exists("JoinLayers", where = asNamespace("SeuratObject"), inherits = FALSE)
  ) {
    log_msg("Joining Assay5 layers for assay ", assay)
    obj <- SeuratObject::JoinLayers(obj, assay = assay)
  }
  obj
}

get_data_layer <- function(obj, assay) {
  last_error <- NULL
  for (layer in c("data", "logcounts", "counts")) {
    mat <- tryCatch(
      SeuratObject::GetAssayData(obj, assay = assay, layer = layer),
      error = function(e) {
        last_error <<- e
        NULL
      }
    )
    if (!is.null(mat)) return(layer)
  }
  stop("Could not extract RNA matrix for assay ", assay, ": ", conditionMessage(last_error), call. = FALSE)
}

ensure_log_normalized <- function(obj, assay) {
  Seurat::DefaultAssay(obj) <- assay
  layer <- get_data_layer(obj, assay)
  if (identical(layer, "counts")) {
    log_msg("No data/logcounts layer found for assay ", assay, "; running NormalizeData")
    obj <- Seurat::NormalizeData(obj, assay = assay, normalization.method = "LogNormalize", verbose = FALSE)
  } else {
    log_msg("Using existing ", layer, " layer for assay ", assay)
  }
  obj
}

choose_existing_column <- function(meta, explicit, candidates, kind) {
  if (nzchar(explicit)) {
    if (!(explicit %in% colnames(meta))) stop("Requested ", kind, " column not found: ", explicit, call. = FALSE)
    return(explicit)
  }
  for (candidate in candidates) {
    if (candidate %in% colnames(meta)) return(candidate)
  }
  hits <- grep(paste(c("shi", "week", "gw", "label", "major", "type"), collapse = "|"), colnames(meta), value = TRUE, ignore.case = TRUE)
  stop(
    "Could not detect Shi reference ", kind, " metadata column. Pass --reference-", kind,
    "-col. Candidate metadata columns: ", paste(hits, collapse = ", "),
    call. = FALSE
  )
}

find_existing_column <- function(meta, explicit, candidates, kind) {
  if (nzchar(explicit)) {
    if (!(explicit %in% colnames(meta))) stop("Requested ", kind, " column not found: ", explicit, call. = FALSE)
    return(explicit)
  }
  for (candidate in candidates) {
    if (candidate %in% colnames(meta)) return(candidate)
  }
  ""
}

unique_key_map <- function(df, key_col, value_col) {
  if (!(key_col %in% colnames(df)) || !(value_col %in% colnames(df))) {
    return(stats::setNames(character(0), character(0)))
  }
  keys <- as.character(df[[key_col]])
  values <- as.character(df[[value_col]])
  keep <- !is.na(keys) & nzchar(keys) & !is.na(values) & nzchar(values)
  key_counts <- table(keys[keep])
  unique_keys <- names(key_counts)[key_counts == 1L]
  keep <- keep & keys %in% unique_keys
  stats::setNames(values[keep], keys[keep])
}

attach_reference_labels_from_tsv <- function(reference, labels_tsv) {
  if (!file.exists(labels_tsv)) {
    stop(
      "Shi reference has no major label metadata column and labels TSV is missing: ",
      labels_tsv,
      call. = FALSE
    )
  }
  log_msg("Attaching Shi major labels from TSV: ", labels_tsv)
  labels <- read_tsv_maybe_gz(labels_tsv)
  if (!("shi_label" %in% colnames(labels))) {
    stop("Reference labels TSV lacks required shi_label column: ", labels_tsv, call. = FALSE)
  }
  maps <- list(
    reference_obs_name = unique_key_map(labels, "reference_obs_name", "shi_label"),
    reference_cell_id = unique_key_map(labels, "reference_cell_id", "shi_label"),
    reference_raw_cell_id = unique_key_map(labels, "reference_raw_cell_id", "shi_label")
  )
  week_maps <- list(
    reference_obs_name = unique_key_map(labels, "reference_obs_name", "shi_week_label"),
    reference_cell_id = unique_key_map(labels, "reference_cell_id", "shi_week_label"),
    reference_raw_cell_id = unique_key_map(labels, "reference_raw_cell_id", "shi_week_label")
  )
  cells <- colnames(reference)
  meta <- reference@meta.data
  result <- rep(NA_character_, length(cells))
  week_result <- rep(NA_character_, length(cells))
  source <- rep("", length(cells))
  for (i in seq_along(cells)) {
    candidates <- list(reference_obs_name = cells[[i]])
    if ("raw_cell_id" %in% colnames(meta)) {
      candidates$reference_raw_cell_id <- as.character(meta[i, "raw_cell_id"])
    }
    for (candidate_name in names(candidates)) {
      key <- candidates[[candidate_name]]
      if (!is.na(key) && nzchar(key) && key %in% names(maps[[candidate_name]])) {
        result[[i]] <- unname(maps[[candidate_name]][[key]])
        source[[i]] <- candidate_name
        if (key %in% names(week_maps[[candidate_name]])) {
          week_result[[i]] <- unname(week_maps[[candidate_name]][[key]])
        }
        break
      }
    }
  }
  reference$shi_label <- result
  reference$shi_label_source <- source
  if (!("shi_week_label" %in% colnames(reference@meta.data)) && any(!is.na(week_result) & nzchar(week_result))) {
    reference$shi_week_label <- week_result
  }
  n_labelled <- sum(!is.na(reference$shi_label) & nzchar(as.character(reference$shi_label)))
  if (n_labelled < 50) {
    stop("Only ", n_labelled, " Shi reference cells received major labels from TSV", call. = FALSE)
  }
  reference
}

select_first_metadata_col <- function(meta, raw_col_spec) {
  candidates <- unlist(strsplit(raw_col_spec, "\\|", perl = TRUE))
  for (candidate in candidates) {
    if (nzchar(candidate) && candidate %in% colnames(meta)) return(candidate)
  }
  ""
}

prediction_score_cols <- function(df) {
  setdiff(grep("^prediction\\.score\\.", colnames(df), value = TRUE), "prediction.score.max")
}

rename_score_columns <- function(df, prefix, score_cols) {
  out <- df[, c("cell_id", score_cols), drop = FALSE]
  new_names <- colnames(out)
  for (col in score_cols) {
    label <- sub("^prediction\\.score\\.", "", col)
    new_names[new_names == col] <- paste0(prefix, sanitize_token(label))
  }
  colnames(out) <- new_names
  out
}

expected_gw_from_scores <- function(scores, prefix) {
  score_cols <- setdiff(colnames(scores), "cell_id")
  labels <- sub(paste0("^", prefix), "", score_cols)
  gw_numeric <- parse_gw_numeric(labels)
  if (any(!is.finite(gw_numeric))) {
    stop("Could not parse numeric GW labels from columns: ", paste(score_cols[!is.finite(gw_numeric)], collapse = ", "), call. = FALSE)
  }
  observed <- sort(unique(gw_numeric))
  even_map <- stats::setNames(seq_along(observed), as.character(observed))
  mat <- as.matrix(scores[, score_cols, drop = FALSE])
  storage.mode(mat) <- "double"
  denom <- rowSums(mat, na.rm = TRUE)
  expected_numeric <- as.numeric(mat %*% gw_numeric) / denom
  expected_even <- as.numeric(mat %*% as.numeric(even_map[as.character(gw_numeric)])) / denom
  expected_numeric[!is.finite(expected_numeric)] <- NA_real_
  expected_even[!is.finite(expected_even)] <- NA_real_
  data.frame(
    cell_id = scores$cell_id,
    shi_seurat_full_expected_shi_gw_numeric = expected_numeric,
    shi_seurat_full_expected_shi_gw_even = expected_even,
    stringsAsFactors = FALSE
  )
}

validate_score_range <- function(df, cols, study_id) {
  for (col in cols) {
    values <- suppressWarnings(as.numeric(df[[col]]))
    bad <- !is.na(values) & (values < -1e-8 | values > 1 + 1e-8)
    if (any(bad)) stop(study_id, " has prediction scores outside [0,1] in ", col, call. = FALSE)
  }
}

build_obs_table <- function(study, query, predictions, label_scores, week_predictions, week_scores) {
  study_id <- study$study_id[[1]]
  meta <- query@meta.data
  if (!(study$reduction[[1]] %in% names(query@reductions))) {
    stop(study_id, " missing UMAP reduction: ", study$reduction[[1]], call. = FALSE)
  }
  emb <- Seurat::Embeddings(query, reduction = study$reduction[[1]])
  if (ncol(emb) < 2) stop(study_id, " reduction has fewer than two dimensions: ", study$reduction[[1]], call. = FALSE)
  cells <- Reduce(intersect, list(colnames(query), rownames(meta), rownames(emb), predictions$cell_id, week_predictions$cell_id))
  if (length(cells) != ncol(query)) {
    stop(study_id, " cell ID join mismatch: matched ", length(cells), " of ", ncol(query), " query cells", call. = FALSE)
  }
  sample_col <- select_first_metadata_col(meta, study$sample_col[[1]])
  cluster_col <- select_first_metadata_col(meta, study$cluster_col[[1]])
  if (!nzchar(sample_col)) stop(study_id, " missing sample metadata column(s): ", study$sample_col[[1]], call. = FALSE)
  if (!nzchar(cluster_col)) stop(study_id, " missing cluster metadata column: ", study$cluster_col[[1]], call. = FALSE)

  predictions <- predictions[match(cells, predictions$cell_id), , drop = FALSE]
  week_predictions <- week_predictions[match(cells, week_predictions$cell_id), , drop = FALSE]
  label_scores <- label_scores[match(cells, label_scores$cell_id), , drop = FALSE]
  week_scores <- week_scores[match(cells, week_scores$cell_id), , drop = FALSE]

  sample <- as.character(meta[cells, sample_col])
  out <- data.frame(
    cell_id = cells,
    study_id = study_id,
    study_label = study$study_label[[1]],
    sample = sample,
    sample_label = sample_label_values(study_id, sample),
    cluster = as.character(meta[cells, cluster_col]),
    umap_1 = emb[cells, 1],
    umap_2 = emb[cells, 2],
    shi_seurat_full_predicted_shi_label = as.character(predictions$predicted.id),
    shi_seurat_full_prediction_score = as.numeric(predictions$prediction.score.max),
    shi_seurat_full_uncertainty_score = 1 - as.numeric(predictions$prediction.score.max),
    shi_seurat_full_predicted_shi_week_label = as.character(week_predictions$predicted.id),
    shi_seurat_full_week_prediction_score = as.numeric(week_predictions$prediction.score.max),
    shi_seurat_full_week_uncertainty_score = 1 - as.numeric(week_predictions$prediction.score.max),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  out <- merge(out, label_scores, by = "cell_id", all.x = TRUE, sort = FALSE)
  out <- merge(out, week_scores, by = "cell_id", all.x = TRUE, sort = FALSE)
  out <- merge(out, expected_gw_from_scores(week_scores, "shi_seurat_full_week_prediction_score_"), by = "cell_id", all.x = TRUE, sort = FALSE)
  out[match(cells, out$cell_id), , drop = FALSE]
}

copy_existing_raw_exports <- function(study, seurat_dir) {
  prefix <- study$existing_prefix[[1]]
  source_dir <- study$existing_seurat_dir[[1]]
  if (!nzchar(prefix) || !dir.exists(source_dir)) return(FALSE)
  suffixes <- c(
    "_shi_seurat_full_predictions.tsv.gz",
    "_shi_seurat_full_prediction_scores.tsv.gz",
    "_shi_seurat_full_week_predictions.tsv.gz",
    "_shi_seurat_full_week_prediction_scores.tsv.gz",
    "_shi_seurat_full_transfer_diagnostics.tsv"
  )
  ok <- TRUE
  for (suffix in suffixes) {
    src <- file.path(source_dir, paste0(prefix, suffix))
    dst <- file.path(seurat_dir, paste0(study$study_id[[1]], suffix))
    if (file.exists(src)) {
      file.copy(src, dst, overwrite = TRUE)
    } else {
      ok <- FALSE
    }
  }
  ok
}

augment_reused_obs_from_query <- function(study, obs, project_root) {
  needs_query <- !all(c("umap_1", "umap_2") %in% colnames(obs)) ||
    !("sample" %in% colnames(obs)) ||
    !("cluster" %in% colnames(obs))
  if (!needs_query) return(obs)
  study_id <- study$study_id[[1]]
  object_path <- resolve_path(study$object_path[[1]], project_root)
  if (!file.exists(object_path)) {
    stop(study_id, " reused prediction table needs target metadata/UMAP, but target object is missing: ", object_path, call. = FALSE)
  }
  log_msg("Loading query metadata/UMAP for reused ", study_id, ": ", object_path)
  query <- readRDS(object_path)
  if (!inherits(query, "Seurat")) stop(study_id, " target is not a Seurat object", call. = FALSE)
  if (!(study$reduction[[1]] %in% names(query@reductions))) {
    stop(study_id, " missing UMAP reduction: ", study$reduction[[1]], call. = FALSE)
  }
  emb <- Seurat::Embeddings(query, reduction = study$reduction[[1]])
  meta <- query@meta.data
  match_col <- ""
  for (candidate in c("cell_id", "cell_id_for_join")) {
    if (candidate %in% colnames(obs)) {
      values <- as.character(obs[[candidate]])
      if (all(values %in% rownames(emb)) && all(values %in% rownames(meta))) {
        match_col <- candidate
        break
      }
    }
  }
  if (!nzchar(match_col)) {
    stop(study_id, " reused prediction table cell IDs cannot be matched to target object UMAP metadata", call. = FALSE)
  }
  cells <- as.character(obs[[match_col]])
  if (!all(c("umap_1", "umap_2") %in% colnames(obs))) {
    obs$umap_1 <- emb[cells, 1]
    obs$umap_2 <- emb[cells, 2]
  }
  sample_col <- select_first_metadata_col(meta, study$sample_col[[1]])
  cluster_col <- select_first_metadata_col(meta, study$cluster_col[[1]])
  if (!nzchar(sample_col)) stop(study_id, " missing sample metadata column(s): ", study$sample_col[[1]], call. = FALSE)
  if (!nzchar(cluster_col)) stop(study_id, " missing cluster metadata column: ", study$cluster_col[[1]], call. = FALSE)
  if (!("sample" %in% colnames(obs))) obs$sample <- as.character(meta[cells, sample_col])
  if (!("cluster" %in% colnames(obs))) obs$cluster <- as.character(meta[cells, cluster_col])
  obs
}

reuse_existing_obs <- function(study, table_per_study_dir, seurat_dir, project_root) {
  path <- study$existing_obs_path[[1]]
  if (!nzchar(path) || !file.exists(path)) return(NULL)
  log_msg("Reusing existing Seurat prediction obs table for ", study$study_id[[1]], ": ", path)
  obs <- read_tsv_maybe_gz(path)
  obs <- augment_reused_obs_from_query(study, obs, project_root)
  required <- c(
    "cell_id", "umap_1", "umap_2", "shi_seurat_full_predicted_shi_label",
    "shi_seurat_full_prediction_score", "shi_seurat_full_predicted_shi_week_label",
    "shi_seurat_full_week_prediction_score"
  )
  missing <- setdiff(required, colnames(obs))
  if (length(missing) > 0) stop("Existing obs table for ", study$study_id[[1]], " is missing: ", paste(missing, collapse = ", "), call. = FALSE)
  if (!("sample" %in% colnames(obs)) && "orig.ident" %in% colnames(obs)) obs$sample <- obs[["orig.ident"]]
  if (!("cluster" %in% colnames(obs)) && "seurat_clusters" %in% colnames(obs)) obs$cluster <- obs[["seurat_clusters"]]
  if (!("sample" %in% colnames(obs))) obs$sample <- ""
  if (!("cluster" %in% colnames(obs))) obs$cluster <- ""
  obs$study_id <- study$study_id[[1]]
  obs$study_label <- study$study_label[[1]]
  obs$sample_label <- sample_label_values(study$study_id[[1]], obs$sample)
  if (!("shi_seurat_full_uncertainty_score" %in% colnames(obs))) {
    obs$shi_seurat_full_uncertainty_score <- 1 - as.numeric(obs$shi_seurat_full_prediction_score)
  }
  if (!("shi_seurat_full_week_uncertainty_score" %in% colnames(obs))) {
    obs$shi_seurat_full_week_uncertainty_score <- 1 - as.numeric(obs$shi_seurat_full_week_prediction_score)
  }
  if ("shi_seurat_full_expected_shi_week_numeric" %in% colnames(obs) && !("shi_seurat_full_expected_shi_gw_numeric" %in% colnames(obs))) {
    obs$shi_seurat_full_expected_shi_gw_numeric <- obs$shi_seurat_full_expected_shi_week_numeric
  }
  week_score_cols <- grep("^shi_seurat_full_week_prediction_score_", colnames(obs), value = TRUE)
  if (!("shi_seurat_full_expected_shi_gw_numeric" %in% colnames(obs)) && length(week_score_cols) > 0) {
    obs <- merge(
      obs,
      expected_gw_from_scores(obs[, c("cell_id", week_score_cols), drop = FALSE], "shi_seurat_full_week_prediction_score_"),
      by = "cell_id",
      all.x = TRUE,
      sort = FALSE
    )
  }
  label_score_cols <- grep("^shi_seurat_full_prediction_score_", colnames(obs), value = TRUE)
  validate_score_range(obs, c(label_score_cols, week_score_cols, "shi_seurat_full_prediction_score", "shi_seurat_full_week_prediction_score"), study$study_id[[1]])
  out_path <- file.path(table_per_study_dir, paste0(study$study_id[[1]], "_shi_seurat_label_transfer_obs.tsv.gz"))
  write_tsv_gz(obs, out_path)
  raw_copied <- copy_existing_raw_exports(study, seurat_dir)
  diagnostics <- data.frame(
    study_id = study$study_id[[1]],
    status = "ok",
    transfer_source = "reused_existing_varela_tables",
    n_query_cells = nrow(obs),
    n_reference_cells = NA_integer_,
    n_shared_features = NA_integer_,
    n_anchors = NA_integer_,
    label_score_columns_exported = paste(label_score_cols, collapse = ","),
    week_score_columns_exported = paste(week_score_cols, collapse = ","),
    n_missing_prediction = sum(is.na(obs$shi_seurat_full_predicted_shi_label) | !nzchar(as.character(obs$shi_seurat_full_predicted_shi_label))),
    umap_reduction = study$reduction[[1]],
    sample_column_used = study$sample_col[[1]],
    cluster_column_used = study$cluster_col[[1]],
    raw_seurat_exports_copied = raw_copied,
    stringsAsFactors = FALSE
  )
  write_tsv(diagnostics, file.path(seurat_dir, paste0(study$study_id[[1]], "_shi_seurat_full_transfer_diagnostics.tsv")))
  diagnostics
}

run_transfer_one <- function(study, reference, label_col, week_col, opt, table_per_study_dir, seurat_dir) {
  study_id <- study$study_id[[1]]
  object_path <- resolve_path(study$object_path[[1]], opt$project_root)
  if (!file.exists(object_path)) stop(study_id, " target object not found: ", object_path, call. = FALSE)
  log_msg("Loading query ", study_id, ": ", object_path)
  query <- readRDS(object_path)
  if (!inherits(query, "Seurat")) stop(study_id, " target is not a Seurat object", call. = FALSE)
  if (!("RNA" %in% Seurat::Assays(query))) stop(study_id, " target lacks RNA assay", call. = FALSE)
  query <- join_layers_if_needed(query, "RNA")
  if (!(study$reduction[[1]] %in% names(query@reductions))) {
    stop(study_id, " missing UMAP reduction: ", study$reduction[[1]], call. = FALSE)
  }

  reference <- join_layers_if_needed(reference, "RNA")
  reference <- ensure_log_normalized(reference, "RNA")
  query <- ensure_log_normalized(query, "RNA")
  ref_features <- rownames(reference[["RNA"]])
  query_features <- rownames(query[["RNA"]])
  shared_features <- intersect(ref_features, query_features)
  min_shared <- as.integer(opt$min_shared_features)
  dims_n <- as.integer(opt$dims)
  if (length(shared_features) < min_shared) {
    stop(study_id, " has too few shared RNA features with Shi reference: ", length(shared_features), " < ", min_shared, call. = FALSE)
  }
  if (length(shared_features) <= dims_n) {
    stop(study_id, " has only ", length(shared_features), " shared features; dims=", dims_n, " requires more features", call. = FALSE)
  }
  if (ncol(reference) <= dims_n || ncol(query) <= dims_n) {
    stop(study_id, " has too few cells for dims=", dims_n, " (reference=", ncol(reference), ", query=", ncol(query), ")", call. = FALSE)
  }
  dims_use <- seq_len(dims_n)

  log_msg("Scaling reference and running PCA for ", study_id, " with ", length(shared_features), " shared features")
  reference <- Seurat::ScaleData(reference, assay = "RNA", features = shared_features, verbose = FALSE)
  reference <- Seurat::RunPCA(reference, assay = "RNA", features = shared_features, npcs = dims_n, verbose = FALSE)

  log_msg("Finding transfer anchors for ", study_id)
  anchors <- Seurat::FindTransferAnchors(
    reference = reference,
    query = query,
    normalization.method = opt$normalization_method,
    reference.assay = "RNA",
    query.assay = "RNA",
    features = shared_features,
    reference.reduction = "pca",
    reduction = "pcaproject",
    dims = dims_use,
    verbose = TRUE
  )

  log_msg("TransferData major Shi labels for ", study_id)
  predictions <- as.data.frame(Seurat::TransferData(
    anchorset = anchors,
    refdata = reference@meta.data[[label_col]],
    dims = dims_use,
    verbose = TRUE
  ), stringsAsFactors = FALSE)
  predictions <- data.frame(cell_id = rownames(predictions), predictions, check.names = FALSE)
  label_score_cols_raw <- prediction_score_cols(predictions)
  if (length(label_score_cols_raw) == 0) stop(study_id, " TransferData returned no major-label score columns", call. = FALSE)
  label_scores <- rename_score_columns(predictions, "shi_seurat_full_prediction_score_", label_score_cols_raw)

  log_msg("TransferData Shi gestational-week labels for ", study_id)
  week_predictions <- as.data.frame(Seurat::TransferData(
    anchorset = anchors,
    refdata = reference@meta.data[[week_col]],
    dims = dims_use,
    verbose = TRUE
  ), stringsAsFactors = FALSE)
  week_predictions <- data.frame(cell_id = rownames(week_predictions), week_predictions, check.names = FALSE)
  week_score_cols_raw <- prediction_score_cols(week_predictions)
  if (length(week_score_cols_raw) == 0) stop(study_id, " TransferData returned no week score columns", call. = FALSE)
  week_scores <- rename_score_columns(week_predictions, "shi_seurat_full_week_prediction_score_", week_score_cols_raw)

  obs <- build_obs_table(study, query, predictions, label_scores, week_predictions, week_scores)
  label_score_cols <- grep("^shi_seurat_full_prediction_score_", colnames(obs), value = TRUE)
  week_score_cols <- grep("^shi_seurat_full_week_prediction_score_", colnames(obs), value = TRUE)
  validate_score_range(obs, c(label_score_cols, week_score_cols, "shi_seurat_full_prediction_score", "shi_seurat_full_week_prediction_score"), study_id)
  gw_range <- range(parse_gw_numeric(sub("^shi_seurat_full_week_prediction_score_", "", week_score_cols)), na.rm = TRUE)
  expected <- obs$shi_seurat_full_expected_shi_gw_numeric
  if (any(!is.na(expected) & (expected < gw_range[[1]] - 1e-8 | expected > gw_range[[2]] + 1e-8))) {
    stop(study_id, " expected GW numeric is outside observed Shi GW score range", call. = FALSE)
  }

  write_tsv_gz(predictions, file.path(seurat_dir, paste0(study_id, "_shi_seurat_full_predictions.tsv.gz")))
  write_tsv_gz(label_scores, file.path(seurat_dir, paste0(study_id, "_shi_seurat_full_prediction_scores.tsv.gz")))
  write_tsv_gz(week_predictions, file.path(seurat_dir, paste0(study_id, "_shi_seurat_full_week_predictions.tsv.gz")))
  write_tsv_gz(week_scores, file.path(seurat_dir, paste0(study_id, "_shi_seurat_full_week_prediction_scores.tsv.gz")))
  write_tsv_gz(obs, file.path(table_per_study_dir, paste0(study_id, "_shi_seurat_label_transfer_obs.tsv.gz")))

  zero_labels <- setdiff(MAJOR_LABELS, unique(as.character(obs$shi_seurat_full_predicted_shi_label)))
  if (length(zero_labels) > 0) warning(study_id, " has zero predicted cells for Shi labels: ", paste(zero_labels, collapse = ", "))
  diagnostics <- data.frame(
    study_id = study_id,
    status = "ok",
    transfer_source = "new_seurat_transfer",
    n_query_cells = ncol(query),
    n_reference_cells = ncol(reference),
    n_shared_features = length(shared_features),
    n_anchors = nrow(anchors@anchors),
    label_score_columns_exported = paste(label_score_cols, collapse = ","),
    week_score_columns_exported = paste(week_score_cols, collapse = ","),
    n_missing_prediction = sum(is.na(obs$shi_seurat_full_predicted_shi_label) | !nzchar(as.character(obs$shi_seurat_full_predicted_shi_label))),
    umap_reduction = study$reduction[[1]],
    sample_column_used = select_first_metadata_col(query@meta.data, study$sample_col[[1]]),
    cluster_column_used = select_first_metadata_col(query@meta.data, study$cluster_col[[1]]),
    raw_seurat_exports_copied = FALSE,
    stringsAsFactors = FALSE
  )
  write_tsv(diagnostics, file.path(seurat_dir, paste0(study_id, "_shi_seurat_full_transfer_diagnostics.tsv")))
  rm(query, anchors, predictions, week_predictions, label_scores, week_scores, obs)
  gc()
  diagnostics
}

main <- function() {
  opt <- parse_args()
  set.seed(as.integer(opt$seed))
  seurat_dir <- file.path(opt$outdir, "seurat", "per_study")
  table_per_study_dir <- file.path(opt$outdir, "tables", "per_study")
  for (path in c(seurat_dir, table_per_study_dir, file.path(opt$outdir, "plots/umap_grids"), file.path(opt$outdir, "plots/summary"), file.path(opt$outdir, "h5ad"), file.path(opt$outdir, "diagnostics"))) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }

  studies <- default_studies(opt$project_root)
  if (length(opt$study_id) > 0L) studies <- studies[studies$study_id %in% opt$study_id, , drop = FALSE]
  if (nrow(studies) == 0L) stop("No target studies selected", call. = FALSE)
  write_tsv(studies, file.path(opt$outdir, "tables/cross_study_shi_seurat_label_transfer_studies.tsv"))

  force_rerun <- to_bool(opt$force_rerun)
  reuse_existing <- to_bool(opt$reuse_existing, default = TRUE)
  diagnostics <- list()
  pending_studies <- studies
  if (reuse_existing && !force_rerun) {
    reusable <- nzchar(studies$existing_obs_path) & file.exists(studies$existing_obs_path)
    if (any(reusable)) {
      for (idx in which(reusable)) {
        diagnostics[[length(diagnostics) + 1L]] <- reuse_existing_obs(
          studies[idx, , drop = FALSE],
          table_per_study_dir,
          seurat_dir,
          opt$project_root
        )
      }
      pending_studies <- studies[!reusable, , drop = FALSE]
    }
  }
  if (nrow(pending_studies) == 0L) {
    summary <- do.call(rbind, diagnostics)
    write_tsv(summary, file.path(opt$outdir, "diagnostics/cross_study_shi_transfer_diagnostics_summary.tsv"))
    write_tsv(summary, file.path(opt$outdir, "tables/cross_study_shi_transfer_diagnostics_summary.tsv"))
    print(summary)
    log_msg("Finished cross-study Shi Seurat label-transfer exports using reused existing tables only: ", opt$outdir)
    return(invisible(summary))
  }

  if (!file.exists(opt$reference)) stop("Missing Shi reference Seurat object: ", opt$reference, call. = FALSE)
  log_msg("Loading Shi reference: ", opt$reference)
  reference <- readRDS(opt$reference)
  if (!inherits(reference, "Seurat")) stop("Reference is not a Seurat object", call. = FALSE)
  if (!("RNA" %in% Seurat::Assays(reference))) stop("Shi reference lacks RNA assay", call. = FALSE)
  label_col <- find_existing_column(
    reference@meta.data,
    opt$reference_label_col,
    c("shi_label", "shi_transfer_label", "Major types", "major_type", "major_cell_type", "cell_type", "celltype"),
    "label"
  )
  if (!nzchar(label_col)) {
    reference <- attach_reference_labels_from_tsv(reference, opt$reference_labels_tsv)
    label_col <- find_existing_column(
      reference@meta.data,
      opt$reference_label_col,
      c("shi_label", "shi_transfer_label", "Major types", "major_type", "major_cell_type", "cell_type", "celltype"),
      "label"
    )
  }
  if (!nzchar(label_col)) {
    choose_existing_column(
      reference@meta.data,
      opt$reference_label_col,
      c("shi_label", "shi_transfer_label", "Major types", "major_type", "major_cell_type", "cell_type", "celltype"),
      "label"
    )
  }
  week_col <- find_existing_column(
    reference@meta.data,
    opt$reference_week_col,
    c("week_label", "shi_week_label", "shi_transfer_week_label", "sample_week_label", "GW", "gw", "week"),
    "week"
  )
  if (!nzchar(week_col)) {
    choose_existing_column(
      reference@meta.data,
      opt$reference_week_col,
      c("week_label", "shi_week_label", "shi_transfer_week_label", "sample_week_label", "GW", "gw", "week"),
      "week"
    )
  }
  log_msg("Reference major label column: ", label_col)
  log_msg("Reference GW label column: ", week_col)
  keep <- !is.na(reference@meta.data[[label_col]]) & nzchar(as.character(reference@meta.data[[label_col]])) &
    !is.na(reference@meta.data[[week_col]]) & nzchar(as.character(reference@meta.data[[week_col]]))
  if (sum(keep) < 50) stop("Too few Shi reference cells have both major label and GW metadata: ", sum(keep), call. = FALSE)
  reference <- subset(reference, cells = colnames(reference)[keep])
  reference@meta.data[[label_col]] <- factor(as.character(reference@meta.data[[label_col]]), levels = unique(c(MAJOR_LABELS, sort(unique(as.character(reference@meta.data[[label_col]]))))))
  reference@meta.data[[week_col]] <- factor(as.character(reference@meta.data[[week_col]]))
  week_numeric <- parse_gw_numeric(levels(reference@meta.data[[week_col]]))
  if (any(!is.finite(week_numeric))) {
    stop("Could not parse numeric GW values from reference week labels: ", paste(levels(reference@meta.data[[week_col]])[!is.finite(week_numeric)], collapse = ", "), call. = FALSE)
  }
  write_tsv(as.data.frame(table(reference@meta.data[[label_col]])), file.path(opt$outdir, "diagnostics/shi_reference_labels_used_by_seurat.tsv"))
  write_tsv(as.data.frame(table(reference@meta.data[[week_col]])), file.path(opt$outdir, "diagnostics/shi_reference_weeks_used_by_seurat.tsv"))

  for (idx in seq_len(nrow(pending_studies))) {
    study <- pending_studies[idx, , drop = FALSE]
    diagnostics[[length(diagnostics) + 1L]] <- run_transfer_one(study, reference, label_col, week_col, opt, table_per_study_dir, seurat_dir)
  }
  summary <- do.call(rbind, diagnostics)
  write_tsv(summary, file.path(opt$outdir, "diagnostics/cross_study_shi_transfer_diagnostics_summary.tsv"))
  write_tsv(summary, file.path(opt$outdir, "tables/cross_study_shi_transfer_diagnostics_summary.tsv"))
  print(summary)
  log_msg("Finished cross-study Shi Seurat label-transfer exports: ", opt$outdir)
}

main()
