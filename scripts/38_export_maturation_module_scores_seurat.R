#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
})

timestamp <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")
log_msg <- function(...) message("[MaturationScoreExport ", timestamp(), "] ", paste0(..., collapse = ""))

parse_args <- function(args) {
  opt <- list(
    "project-root" = Sys.getenv("PROJECT_ROOT", unset = ""),
    "div30-rds" = "",
    "div90-rds" = "",
    "marker-csv" = "",
    "outdir" = "",
    "assay" = "RNA",
    "ctrl" = "50",
    "nbin" = "24",
    "seed" = "0"
  )
  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, "--")) stop("Unknown argument: ", key, call. = FALSE)
    name <- substring(key, 3L)
    if (!(name %in% names(opt))) stop("Unknown argument: ", key, call. = FALSE)
    if (i == length(args)) stop("Missing value for argument: ", key, call. = FALSE)
    opt[[name]] <- args[[i + 1L]]
    i <- i + 2L
  }
  opt
}

write_tsv <- function(df, path) {
  con <- if (grepl("\\.gz$", path)) gzfile(path, open = "wt") else file(path, open = "wt")
  on.exit(close(con), add = TRUE)
  write.table(df, con, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
}

trim_trailing_slash <- function(x) sub("/+$", "", x)

read_jia_markers <- function(path) {
  markers <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE, fileEncoding = "UTF-8-BOM")
  names(markers) <- sub("^\ufeff", "", names(markers))
  required <- c("gene", "cluster")
  missing <- setdiff(required, names(markers))
  if (length(missing) > 0) stop("Jia marker CSV is missing column(s): ", paste(missing, collapse = ", "), call. = FALSE)
  markers <- markers[!is.na(markers$gene) & nzchar(markers$gene) & !is.na(markers$cluster) & nzchar(markers$cluster), , drop = FALSE]
  markers$gene <- as.character(markers$gene)
  markers$cluster <- as.character(markers$cluster)
  markers
}

find_genes_case_insensitive <- function(requested, available) {
  requested <- unique(as.character(requested))
  available_upper <- toupper(as.character(available))
  idx <- match(toupper(requested), available_upper)
  found <- as.character(available[idx[!is.na(idx)]])
  names(found) <- requested[!is.na(idx)]
  missing <- requested[is.na(idx)]
  list(found = unname(unique(found)), missing = unique(missing))
}

join_assay_layers_if_available <- function(obj, assay) {
  if (exists("JoinLayers", where = asNamespace("SeuratObject"), inherits = FALSE)) {
    obj <- tryCatch(SeuratObject::JoinLayers(obj, assay = assay), error = function(e) obj)
  } else if (exists("JoinLayers", where = asNamespace("Seurat"), inherits = FALSE)) {
    obj <- tryCatch(Seurat::JoinLayers(obj, assay = assay), error = function(e) obj)
  }
  obj
}

choose_umap_reduction <- function(obj) {
  reductions <- Seurat::Reductions(obj)
  if ("umap" %in% reductions) return("umap")
  hit <- reductions[grepl("umap", reductions, ignore.case = TRUE)]
  if (length(hit) > 0) return(hit[[1]])
  stop("No UMAP reduction found. Reductions: ", paste(reductions, collapse = ", "), call. = FALSE)
}

safe_metadata_vector <- function(meta, col) {
  if (!(col %in% colnames(meta))) return(NULL)
  values <- meta[[col]]
  if (is.factor(values)) values <- as.character(values)
  if (is.list(values) && !is.data.frame(values)) values <- vapply(values, function(x) paste(as.character(x), collapse = ";"), character(1))
  values
}

extract_cluster_numbers <- function(values) {
  values <- as.character(values)
  values[is.na(values)] <- ""
  out <- sub("^\\s*([0-9]+).*$", "\\1", values)
  out[!grepl("^\\s*[0-9]+", values)] <- NA_character_
  out
}

choose_visual_cluster <- function(meta) {
  candidates <- c(
    "seurat_clusters",
    "cluster_number",
    "cluster",
    "cluster_id",
    "RNA_snn_res.0.5",
    "RNA_snn_res.0.2",
    "integrated_snn_res.0.5",
    "SCT_snn_res.0.5",
    "cluster_number_name"
  )
  for (col in candidates) {
    if (!(col %in% colnames(meta))) next
    nums <- extract_cluster_numbers(meta[[col]])
    if (sum(!is.na(nums)) > 0) return(list(col = col, values = nums))
  }
  list(col = "", values = rep(NA_character_, nrow(meta)))
}

score_module <- function(obj, requested_genes, score_name, assay, ctrl, nbin, seed) {
  matched <- find_genes_case_insensitive(requested_genes, rownames(obj))
  if (length(matched$found) == 0) {
    obj@meta.data[[score_name]] <- NA_real_
    return(list(obj = obj, found = character(0), missing = matched$missing))
  }

  prefix <- paste0("tmp_", score_name, "_")
  before <- colnames(obj@meta.data)
  set.seed(seed)
  obj <- Seurat::AddModuleScore(
    object = obj,
    features = list(matched$found),
    assay = assay,
    name = prefix,
    ctrl = ctrl,
    nbin = nbin,
    seed = seed,
    search = FALSE
  )
  after <- colnames(obj@meta.data)
  new_cols <- setdiff(after, before)
  if (length(new_cols) == 0) {
    expected <- paste0(prefix, "1")
    if (expected %in% after) new_cols <- expected
  }
  if (length(new_cols) == 0) stop("AddModuleScore did not create a metadata column for ", score_name, call. = FALSE)
  obj@meta.data[[score_name]] <- as.numeric(obj@meta.data[[new_cols[[1]]]])
  for (col in new_cols) {
    if (col != score_name && col %in% colnames(obj@meta.data)) obj@meta.data[[col]] <- NULL
  }
  list(obj = obj, found = matched$found, missing = matched$missing)
}

build_gene_report_row <- function(dataset, score, requested, found, missing, assay, ctrl, nbin, seed) {
  data.frame(
    dataset = dataset,
    score = score,
    scoring_method = "Seurat::AddModuleScore",
    assay = assay,
    ctrl = ctrl,
    nbin = nbin,
    seed = seed,
    n_genes_requested = length(unique(requested)),
    n_genes_found = length(unique(found)),
    genes_requested = paste(unique(requested), collapse = ", "),
    genes_found = paste(unique(found), collapse = ", "),
    genes_missing = paste(unique(missing), collapse = ", "),
    stringsAsFactors = FALSE
  )
}

score_dataset <- function(dataset, rds_path, marker_programs, predefined_sets, assay, ctrl, nbin, seed) {
  log_msg("Reading ", dataset, " Seurat object: ", rds_path)
  obj <- readRDS(rds_path)
  if (!inherits(obj, "Seurat")) stop("Input is not a Seurat object: ", rds_path, call. = FALSE)
  if (!(assay %in% Seurat::Assays(obj))) {
    assay <- Seurat::DefaultAssay(obj)
    log_msg(dataset, ": requested assay not present; using default assay ", assay)
  }
  DefaultAssay(obj) <- assay
  obj <- join_assay_layers_if_available(obj, assay)
  log_msg(dataset, ": ", ncol(obj), " cells x ", nrow(obj), " genes; assay=", assay)

  reports <- list()
  score_names <- character(0)
  for (program in c("RGC1", "RGC2", "IPC")) {
    score_name <- paste0("jia_score_", program)
    result <- score_module(obj, marker_programs[[program]], score_name, assay, ctrl, nbin, seed)
    obj <- result$obj
    reports[[length(reports) + 1L]] <- build_gene_report_row(dataset, score_name, marker_programs[[program]], result$found, result$missing, assay, ctrl, nbin, seed)
    score_names <- c(score_names, score_name)
    log_msg(dataset, " ", score_name, ": found ", length(result$found), "/", length(unique(marker_programs[[program]])), " genes")
  }
  obj@meta.data[["jia_score_RGC1_RGC2_mean"]] <- rowMeans(obj@meta.data[, c("jia_score_RGC1", "jia_score_RGC2"), drop = FALSE], na.rm = FALSE)

  for (score_name in names(predefined_sets)) {
    result <- score_module(obj, predefined_sets[[score_name]], score_name, assay, ctrl, nbin, seed)
    obj <- result$obj
    reports[[length(reports) + 1L]] <- build_gene_report_row(dataset, score_name, predefined_sets[[score_name]], result$found, result$missing, assay, ctrl, nbin, seed)
    score_names <- c(score_names, score_name)
    log_msg(dataset, " ", score_name, ": found ", length(result$found), "/", length(unique(predefined_sets[[score_name]])), " genes")
  }

  reduction <- choose_umap_reduction(obj)
  emb <- Seurat::Embeddings(obj, reduction = reduction)
  emb <- emb[colnames(obj), , drop = FALSE]
  meta <- obj@meta.data
  visual_cluster <- choose_visual_cluster(meta)
  cluster_values <- visual_cluster$values
  excluded <- dataset == "DIV90" & cluster_values %in% c("6", "7")

  out <- data.frame(
    dataset = dataset,
    cell_id = colnames(obj),
    umap_1 = as.numeric(emb[, 1]),
    umap_2 = as.numeric(emb[, 2]),
    div90_visualization_cluster_col = ifelse(dataset == "DIV90", visual_cluster$col, ""),
    div90_visualization_cluster = ifelse(is.na(cluster_values), "", cluster_values),
    div90_excluded_stressed_cluster = excluded,
    plot_include = !excluded,
    stringsAsFactors = FALSE
  )

  metadata_cols <- c("orig.ident", "sample", "seurat_clusters", "cluster_number_name", "paper_cluster_annotation", "RNA_snn_res.0.5", "RNA_snn_res.0.2")
  for (col in metadata_cols) {
    values <- safe_metadata_vector(meta, col)
    if (!is.null(values)) out[[col]] <- values
  }
  for (col in c("jia_score_RGC1", "jia_score_RGC2", "jia_score_IPC", "jia_score_RGC1_RGC2_mean", names(predefined_sets))) {
    out[[col]] <- as.numeric(meta[[col]])
  }

  list(obs = out, gene_report = do.call(rbind, reports))
}

main <- function() {
  opt <- parse_args(commandArgs(trailingOnly = TRUE))
  project_root <- trim_trailing_slash(opt[["project-root"]])
  if (!nzchar(project_root)) stop("--project-root or PROJECT_ROOT is required", call. = FALSE)
  div30_rds <- if (nzchar(opt[["div30-rds"]])) opt[["div30-rds"]] else file.path(project_root, "results/varela_this_paper/varela_this_paper_seurat.rds")
  div90_rds <- if (nzchar(opt[["div90-rds"]])) opt[["div90-rds"]] else "/nfs/turbo/umms-parent/Manny_test/ventral_sosrs_output/umap_props_output/clustered_day90_with_cluster_names_2.rds"
  marker_csv <- if (nzchar(opt[["marker-csv"]])) opt[["marker-csv"]] else file.path(project_root, "reference/Jia_et_al_2026_Science_3_progs.csv")
  outdir <- if (nzchar(opt$outdir)) opt$outdir else file.path(project_root, "results/maturation_scores/maturation_scores_v1/tables")
  assay <- opt$assay
  ctrl <- as.integer(opt$ctrl)
  nbin <- as.integer(opt$nbin)
  seed <- as.integer(opt$seed)

  for (path in c(div30_rds, div90_rds, marker_csv)) {
    if (!file.exists(path)) stop("Input file not found: ", path, call. = FALSE)
  }
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  markers <- read_jia_markers(marker_csv)
  marker_programs <- lapply(c("RGC1", "RGC2", "IPC"), function(program) unique(markers$gene[markers$cluster == program]))
  names(marker_programs) <- c("RGC1", "RGC2", "IPC")
  predefined_sets <- list(
    immature_module_score = c("DCX", "STMN2", "STMN4", "SOX11", "TUBB3", "TUBB2B", "ELAVL4", "GAP43", "CXCR4", "ACKR3"),
    mature_module_score = c("RBFOX3", "SNAP25", "SYT1", "SYN1", "SYN2", "DLG4", "VAMP2", "SLC12A5", "GAD1", "GAD2", "SLC6A1", "ERBB4")
  )

  write_tsv(markers, file.path(outdir, "jia_program_markers_selected.tsv"))
  module_spec <- do.call(rbind, c(
    lapply(names(marker_programs), function(name) data.frame(source = "Jia_et_al", score = paste0("jia_score_", name), gene = marker_programs[[name]], stringsAsFactors = FALSE)),
    lapply(names(predefined_sets), function(name) data.frame(source = "predefined", score = name, gene = predefined_sets[[name]], stringsAsFactors = FALSE))
  ))
  write_tsv(module_spec, file.path(outdir, "maturation_score_module_gene_sets_requested.tsv"))

  log_msg("Scoring with Seurat::AddModuleScore ctrl=", ctrl, ", nbin=", nbin, ", seed=", seed)
  div30 <- score_dataset("DIV30", div30_rds, marker_programs, predefined_sets, assay, ctrl, nbin, seed)
  gc()
  div90 <- score_dataset("DIV90", div90_rds, marker_programs, predefined_sets, assay, ctrl, nbin, seed)

  obs <- rbind(div30$obs, div90$obs)
  report <- rbind(div30$gene_report, div90$gene_report)
  write_tsv(obs, file.path(outdir, "div30_div90_maturation_scores_obs.tsv.gz"))
  write_tsv(report, file.path(outdir, "maturation_score_gene_report.tsv"))
  write_tsv(subset(report, dataset == "DIV30"), file.path(outdir, "DIV30_predefined_maturation_gene_set_report.tsv"))
  write_tsv(subset(report, dataset == "DIV90"), file.path(outdir, "DIV90_predefined_maturation_gene_set_report.tsv"))

  filter_summary <- do.call(rbind, lapply(split(obs, obs$dataset), function(frame) {
    data.frame(
      dataset = frame$dataset[[1]],
      n_plot = sum(frame$plot_include),
      n_excluded = sum(!frame$plot_include),
      n_total = nrow(frame),
      stringsAsFactors = FALSE
    )
  }))
  write_tsv(filter_summary, file.path(outdir, "maturation_scores_plot_filter_summary.tsv"))
  log_msg("Wrote score table: ", file.path(outdir, "div30_div90_maturation_scores_obs.tsv.gz"))
}

main()
