#!/usr/bin/env Rscript

# Publish two explicit DIV30 progenitor correction outputs:
#   1. sample-only Seurat integration, using the completed CCA integration
#   2. sample integration followed by cell-line regression in integrated space
#
# The script also redraws UMAPs with fully opaque points so the audit plots do
# not look washed out.

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(ggplot2)
  library(patchwork)
})

timestamp <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")
log_msg <- function(...) {
  cat(sprintf("[%s] %s\n", timestamp(), paste0(..., collapse = "")))
  flush.console()
}

parse_args <- function(args) {
  out <- list(
    `project-root` = Sys.getenv("PROJECT_ROOT"),
    `source-dir` = Sys.getenv("SOURCE_DIR", ""),
    `sample-outdir` = Sys.getenv("SAMPLE_OUTDIR", ""),
    `sample-cellline-outdir` = Sys.getenv("SAMPLE_CELLLINE_OUTDIR", ""),
    `sample-map` = Sys.getenv("SAMPLE_MAP", ""),
    `sample-col` = Sys.getenv("SAMPLE_COL", "orig.ident"),
    dims = Sys.getenv("DIMS", "50"),
    npcs = Sys.getenv("NPCS", "60"),
    resolution = Sys.getenv("RESOLUTION", "0.8"),
    seed = Sys.getenv("SEED", "7"),
    `point-size` = Sys.getenv("POINT_SIZE", "0.28"),
    help = FALSE
  )
  i <- 1L
  while (i <= length(args)) {
    a <- args[[i]]
    if (a %in% c("--help", "-h")) {
      out$help <- TRUE
      i <- i + 1L
      next
    }
    if (!startsWith(a, "--")) stop("Unknown argument: ", a, call. = FALSE)
    key <- substring(a, 3L)
    if (!(key %in% names(out))) stop("Unknown argument: ", a, call. = FALSE)
    if (i == length(args)) stop("Missing value for argument: ", a, call. = FALSE)
    out[[key]] <- args[[i + 1L]]
    i <- i + 2L
  }
  out
}

print_usage <- function() {
  cat(paste(
    "Usage:",
    "  Rscript scripts/32_div30_jia_progenitor_publish_sample_and_cellline_outputs.R --project-root <PROJECT_ROOT>",
    "",
    "Inputs default to:",
    "  PROJECT_ROOT/results/div30_jia_progenitor_reclustering/div30_jia_progenitor_strong_integration_v1",
    "",
    "Outputs default to:",
    "  PROJECT_ROOT/results/div30_jia_progenitor_reclustering/div30_jia_progenitor_sample_integration_v1",
    "  PROJECT_ROOT/results/div30_jia_progenitor_reclustering/div30_jia_progenitor_sample_plus_cellline_correction_v1",
    sep = "\n"
  ))
}

trim_trailing_slash <- function(x) sub("/+$", "", x)
as_int <- function(x, name) {
  value <- suppressWarnings(as.integer(x))
  if (is.na(value)) stop(name, " must be an integer; got ", x, call. = FALSE)
  value
}
as_num <- function(x, name) {
  value <- suppressWarnings(as.numeric(x))
  if (is.na(value)) stop(name, " must be numeric; got ", x, call. = FALSE)
  value
}
write_tsv <- function(x, path) {
  con <- if (grepl("\\.gz$", path)) gzfile(path, open = "wt") else file(path, open = "wt")
  on.exit(close(con), add = TRUE)
  write.table(x, con, sep = "\t", row.names = FALSE, quote = FALSE, na = "")
}
save_plot_pair <- function(plot, stem, width, height, dpi = 300) {
  ggplot2::ggsave(paste0(stem, ".png"), plot = plot, width = width, height = height, dpi = dpi, bg = "white")
  ggplot2::ggsave(paste0(stem, ".pdf"), plot = plot, width = width, height = height, bg = "white")
}
safe_name <- function(x) gsub("[^A-Za-z0-9_+-]+", "_", x)

cell_line_from_label <- function(x) {
  out <- sub("_rep[0-9]+$", "", x)
  out <- sub("_old$", "", out)
  out
}

sample_key_variants <- function(x) {
  unique(c(
    x,
    sub("^9583-", "9853-", x),
    sub("^9853-", "9583-", x)
  ))
}

sample_map_expanded <- function(path) {
  if (!file.exists(path)) stop("Sample map not found: ", path, call. = FALSE)
  sample_map <- read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
  needed <- c("DIV", "run_sample_id", "biological_label")
  if (!all(needed %in% colnames(sample_map))) {
    stop("Sample map must contain: ", paste(needed, collapse = ", "), call. = FALSE)
  }
  sample_map <- sample_map[sample_map$DIV == "DIV30", needed, drop = FALSE]
  sample_map$cell_line <- cell_line_from_label(sample_map$biological_label)
  expanded <- do.call(rbind, lapply(seq_len(nrow(sample_map)), function(i) {
    data.frame(
      sample_join_key = sample_key_variants(sample_map$run_sample_id[[i]]),
      run_sample_id_from_sample_map = sample_map$run_sample_id[[i]],
      biological_label = sample_map$biological_label[[i]],
      cell_line = sample_map$cell_line[[i]],
      stringsAsFactors = FALSE
    )
  }))
  expanded[!duplicated(expanded$sample_join_key), , drop = FALSE]
}

attach_sample_map <- function(obj, sample_col, sample_map_path) {
  if (!(sample_col %in% colnames(obj@meta.data))) stop("Sample column missing from object: ", sample_col, call. = FALSE)
  expanded <- sample_map_expanded(sample_map_path)
  idx <- match(as.character(obj@meta.data[[sample_col]]), expanded$sample_join_key)
  obj$run_sample_id_from_sample_map <- expanded$run_sample_id_from_sample_map[idx]
  obj$biological_label <- expanded$biological_label[idx]
  obj$cell_line_from_sample_map <- expanded$cell_line[idx]
  missing <- unique(as.character(obj@meta.data[[sample_col]])[is.na(obj$cell_line_from_sample_map)])
  if (length(missing) > 0L) {
    stop("Could not map sample(s) to cell line: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  obj$cell_line_from_sample_map <- factor(
    as.character(obj$cell_line_from_sample_map),
    levels = c("H9", "79B", "2E", sort(setdiff(unique(as.character(obj$cell_line_from_sample_map)), c("H9", "79B", "2E"))))
  )
  obj
}

embedding_df <- function(obj, reduction, sample_col, cluster_col = NULL) {
  coords <- as.data.frame(Embeddings(obj, reduction))
  coords$cell_id <- rownames(coords)
  colnames(coords)[seq_len(2)] <- c("UMAP_1", "UMAP_2")
  meta <- obj@meta.data
  meta$cell_id <- rownames(meta)
  keep <- unique(c(
    "cell_id", sample_col, "cell_line_from_sample_map", "biological_label",
    "div30_parent_cluster", "div30_parent_progenitor_label", cluster_col,
    "nFeature_RNA", "nCount_RNA", "percent.mt", "S.Score", "G2M.Score", "CC.Difference",
    "jia_score_RGC1", "jia_score_RGC2", "jia_score_IPC"
  ))
  keep <- keep[keep %in% colnames(meta)]
  merge(coords, meta[, keep, drop = FALSE], by = "cell_id", all.x = TRUE)
}

plain_umap <- function(df, color_col, title, point_size) {
  ggplot(df, aes(x = UMAP_1, y = UMAP_2, color = .data[[color_col]])) +
    geom_point(size = point_size, alpha = 1, stroke = 0) +
    coord_equal() +
    labs(title = title, x = "UMAP_1", y = "UMAP_2", color = NULL) +
    guides(color = guide_legend(override.aes = list(size = 2.4, alpha = 1))) +
    theme_classic(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold", size = 11),
      legend.key.height = unit(0.34, "cm")
    )
}

feature_umap <- function(df, color_col, title, point_size) {
  ggplot(df, aes(x = UMAP_1, y = UMAP_2, color = .data[[color_col]])) +
    geom_point(size = point_size, alpha = 1, stroke = 0) +
    coord_equal() +
    scale_color_viridis_c(option = "magma") +
    labs(title = title, x = "UMAP_1", y = "UMAP_2", color = color_col) +
    theme_classic(base_size = 10) +
    theme(plot.title = element_text(face = "bold", size = 11))
}

write_counts <- function(obj, sample_col, cluster_col, table_dir, prefix) {
  meta <- obj@meta.data
  write_tsv(
    as.data.frame(table(sample = meta[[sample_col]], cell_line = meta$cell_line_from_sample_map), stringsAsFactors = FALSE),
    file.path(table_dir, paste0(prefix, "_sample_by_cell_line_counts.tsv"))
  )
  if (!is.null(cluster_col) && cluster_col %in% colnames(meta)) {
    write_tsv(
      as.data.frame(table(cluster = meta[[cluster_col]], sample = meta[[sample_col]]), stringsAsFactors = FALSE),
      file.path(table_dir, paste0(prefix, "_cluster_by_sample_counts.tsv"))
    )
    write_tsv(
      as.data.frame(table(cluster = meta[[cluster_col]], cell_line = meta$cell_line_from_sample_map), stringsAsFactors = FALSE),
      file.path(table_dir, paste0(prefix, "_cluster_by_cell_line_counts.tsv"))
    )
  }
}

copy_if_exists <- function(from, to) {
  if (file.exists(from)) {
    dir.create(dirname(to), recursive = TRUE, showWarnings = FALSE)
    file.copy(from, to, overwrite = TRUE)
  }
}

add_cell_line_regression_dummies <- function(obj, cell_line_col) {
  lines <- factor(as.character(obj@meta.data[[cell_line_col]]))
  mm <- stats::model.matrix(~ lines)
  if (ncol(mm) <= 1L) stop("Cell-line regression needs at least two cell lines.", call. = FALSE)
  dummy_cols <- character()
  for (j in 2:ncol(mm)) {
    col <- paste0("regress_cell_line_", safe_name(colnames(mm)[[j]]))
    obj[[col]] <- as.numeric(mm[, j])
    dummy_cols <- c(dummy_cols, col)
  }
  list(object = obj, dummy_cols = dummy_cols)
}

rerun_integrated_space <- function(obj, dims, npcs, resolution, seed, vars_to_regress = NULL) {
  DefaultAssay(obj) <- "integrated"
  set.seed(seed)
  obj <- ScaleData(obj, assay = "integrated", vars.to.regress = vars_to_regress, verbose = FALSE)
  obj <- RunPCA(obj, assay = "integrated", npcs = npcs, reduction.name = "sample_cellline_pca", reduction.key = "SCLPCA_", seed.use = seed, verbose = FALSE)
  obj <- FindNeighbors(
    obj,
    reduction = "sample_cellline_pca",
    dims = seq_len(dims),
    graph.name = c("sample_cellline_nn", "sample_cellline_snn"),
    verbose = FALSE
  )
  obj <- FindClusters(
    obj,
    graph.name = "sample_cellline_snn",
    resolution = resolution,
    random.seed = seed,
    cluster.name = "div30_progenitor_sample_cellline_cluster",
    verbose = FALSE
  )
  obj <- RunUMAP(
    obj,
    reduction = "sample_cellline_pca",
    dims = seq_len(dims),
    reduction.name = "sample_cellline_umap",
    reduction.key = "SCLUMAP_",
    seed.use = seed,
    verbose = FALSE
  )
  obj
}

make_sample_only_outputs <- function(uncorrected, sample_integrated, cfg) {
  dir.create(cfg$sample_table_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(cfg$sample_plot_dir, recursive = TRUE, showWarnings = FALSE)

  copy_if_exists(file.path(cfg$source_dir, "div30_jia_progenitor_uncorrected_baseline_seurat.rds"), file.path(cfg$sample_outdir, "div30_jia_progenitor_uncorrected_baseline_seurat.rds"))
  copy_if_exists(file.path(cfg$source_dir, "div30_jia_progenitor_strong_integrated_seurat.rds"), file.path(cfg$sample_outdir, "div30_jia_progenitor_sample_integrated_seurat.rds"))

  unc_df <- embedding_df(uncorrected, "uncorrected_umap", cfg$sample_col, "div30_progenitor_uncorrected_cluster")
  int_df <- embedding_df(sample_integrated, "integrated_umap", cfg$sample_col, "div30_progenitor_integrated_cluster")

  save_plot_pair(
    plain_umap(unc_df, cfg$sample_col, "Before integration: sample", cfg$point_size) +
      plain_umap(int_df, cfg$sample_col, "Sample-integrated: sample", cfg$point_size) +
      plot_layout(ncol = 2),
    file.path(cfg$sample_plot_dir, "uncorrected_vs_sample_integrated_umap_by_sample"),
    14, 5.8
  )
  save_plot_pair(
    plain_umap(unc_df, "cell_line_from_sample_map", "Before integration: cell line", cfg$point_size) +
      plain_umap(int_df, "cell_line_from_sample_map", "Sample-integrated: cell line", cfg$point_size) +
      plot_layout(ncol = 2),
    file.path(cfg$sample_plot_dir, "uncorrected_vs_sample_integrated_umap_by_cell_line"),
    12, 5.8
  )
  save_plot_pair(
    plain_umap(unc_df, "div30_progenitor_uncorrected_cluster", "Before integration: reclusters", cfg$point_size) +
      plain_umap(int_df, "div30_progenitor_integrated_cluster", "Sample-integrated: reclusters", cfg$point_size) +
      plot_layout(ncol = 2),
    file.path(cfg$sample_plot_dir, "uncorrected_vs_sample_integrated_umap_by_recluster"),
    14, 5.8
  )

  write_counts(uncorrected, cfg$sample_col, "div30_progenitor_uncorrected_cluster", cfg$sample_table_dir, "uncorrected")
  write_counts(sample_integrated, cfg$sample_col, "div30_progenitor_integrated_cluster", cfg$sample_table_dir, "sample_integrated")
  write_tsv(
    data.frame(
      status = "complete",
      completed_at = timestamp(),
      correction = "sample_only_seurat_cca_integration",
      source_dir = cfg$source_dir,
      note = "UMAPs redrawn with alpha=1 opaque points.",
      stringsAsFactors = FALSE
    ),
    file.path(cfg$sample_table_dir, "sample_integration_complete.tsv")
  )
}

make_sample_cellline_outputs <- function(sample_integrated, cfg) {
  dir.create(cfg$sample_cellline_table_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(cfg$sample_cellline_plot_dir, recursive = TRUE, showWarnings = FALSE)

  before_df <- embedding_df(sample_integrated, "integrated_umap", cfg$sample_col, "div30_progenitor_integrated_cluster")
  dummy <- add_cell_line_regression_dummies(sample_integrated, "cell_line_from_sample_map")
  obj <- dummy$object
  log_msg("Running sample-integrated plus cell-line regression with vars: ", paste(dummy$dummy_cols, collapse = ", "))
  obj <- rerun_integrated_space(
    obj,
    dims = cfg$dims,
    npcs = cfg$npcs,
    resolution = cfg$resolution,
    seed = cfg$seed,
    vars_to_regress = dummy$dummy_cols
  )
  after_df <- embedding_df(obj, "sample_cellline_umap", cfg$sample_col, "div30_progenitor_sample_cellline_cluster")

  save_plot_pair(
    plain_umap(before_df, cfg$sample_col, "Sample-integrated: sample", cfg$point_size) +
      plain_umap(after_df, cfg$sample_col, "Sample + cell-line corrected: sample", cfg$point_size) +
      plot_layout(ncol = 2),
    file.path(cfg$sample_cellline_plot_dir, "sample_integrated_vs_sample_cellline_corrected_umap_by_sample"),
    14, 5.8
  )
  save_plot_pair(
    plain_umap(before_df, "cell_line_from_sample_map", "Sample-integrated: cell line", cfg$point_size) +
      plain_umap(after_df, "cell_line_from_sample_map", "Sample + cell-line corrected: cell line", cfg$point_size) +
      plot_layout(ncol = 2),
    file.path(cfg$sample_cellline_plot_dir, "sample_integrated_vs_sample_cellline_corrected_umap_by_cell_line"),
    12, 5.8
  )
  save_plot_pair(
    plain_umap(before_df, "div30_progenitor_integrated_cluster", "Sample-integrated: reclusters", cfg$point_size) +
      plain_umap(after_df, "div30_progenitor_sample_cellline_cluster", "Sample + cell-line corrected: reclusters", cfg$point_size) +
      plot_layout(ncol = 2),
    file.path(cfg$sample_cellline_plot_dir, "sample_integrated_vs_sample_cellline_corrected_umap_by_recluster"),
    14, 5.8
  )

  for (metric in c("nFeature_RNA", "nCount_RNA", "percent.mt", "S.Score", "G2M.Score", "CC.Difference")) {
    if (metric %in% colnames(after_df)) {
      save_plot_pair(
        feature_umap(before_df, metric, paste("Sample-integrated:", metric), cfg$point_size) +
          feature_umap(after_df, metric, paste("Sample + cell-line corrected:", metric), cfg$point_size) +
          plot_layout(ncol = 2),
        file.path(cfg$sample_cellline_plot_dir, paste0("sample_integrated_vs_sample_cellline_corrected_umap_", safe_name(metric))),
        14, 5.8
      )
    }
  }

  saveRDS(obj, file.path(cfg$sample_cellline_outdir, "div30_jia_progenitor_sample_plus_cellline_corrected_seurat.rds"))
  write_counts(sample_integrated, cfg$sample_col, "div30_progenitor_integrated_cluster", cfg$sample_cellline_table_dir, "sample_integrated")
  write_counts(obj, cfg$sample_col, "div30_progenitor_sample_cellline_cluster", cfg$sample_cellline_table_dir, "sample_cellline_corrected")
  write_tsv(
    data.frame(
      status = "complete",
      completed_at = timestamp(),
      correction = "sample_seurat_cca_integration_plus_cell_line_regression",
      source_sample_integrated_object = file.path(cfg$source_dir, "div30_jia_progenitor_strong_integrated_seurat.rds"),
      cell_line_column = "cell_line_from_sample_map",
      regression_variables = paste(dummy$dummy_cols, collapse = ","),
      note = "Cell-line correction is linear regression on the sample-integrated assay before PCA/neighbors/clustering/UMAP.",
      stringsAsFactors = FALSE
    ),
    file.path(cfg$sample_cellline_table_dir, "sample_plus_cellline_correction_complete.tsv")
  )
}

opt <- parse_args(commandArgs(trailingOnly = TRUE))
if (isTRUE(opt$help)) {
  print_usage()
  quit(save = "no", status = 0)
}
if (!nzchar(opt$`project-root`)) stop("PROJECT_ROOT or --project-root is required", call. = FALSE)
project_root <- trim_trailing_slash(opt$`project-root`)
cfg <- list(
  project_root = project_root,
  source_dir = if (nzchar(opt$`source-dir`)) opt$`source-dir` else file.path(project_root, "results/div30_jia_progenitor_reclustering/div30_jia_progenitor_strong_integration_v1"),
  sample_outdir = if (nzchar(opt$`sample-outdir`)) opt$`sample-outdir` else file.path(project_root, "results/div30_jia_progenitor_reclustering/div30_jia_progenitor_sample_integration_v1"),
  sample_cellline_outdir = if (nzchar(opt$`sample-cellline-outdir`)) opt$`sample-cellline-outdir` else file.path(project_root, "results/div30_jia_progenitor_reclustering/div30_jia_progenitor_sample_plus_cellline_correction_v1"),
  sample_map = if (nzchar(opt$`sample-map`)) opt$`sample-map` else file.path(getwd(), "metadata/div30_div90_sample_id_to_biolabel_map.tsv"),
  sample_col = opt$`sample-col`,
  dims = as_int(opt$dims, "dims"),
  npcs = max(as_int(opt$npcs, "npcs"), as_int(opt$dims, "dims")),
  resolution = as_num(opt$resolution, "resolution"),
  seed = as_int(opt$seed, "seed"),
  point_size = as_num(opt$`point-size`, "point-size")
)
cfg$sample_table_dir <- file.path(cfg$sample_outdir, "tables")
cfg$sample_plot_dir <- file.path(cfg$sample_outdir, "plots")
cfg$sample_cellline_table_dir <- file.path(cfg$sample_cellline_outdir, "tables")
cfg$sample_cellline_plot_dir <- file.path(cfg$sample_cellline_outdir, "plots")

options(future.globals.maxSize = max(getOption("future.globals.maxSize", 0), 120 * 1024^3))
if (requireNamespace("future", quietly = TRUE)) future::plan("sequential")

uncorrected_path <- file.path(cfg$source_dir, "div30_jia_progenitor_uncorrected_baseline_seurat.rds")
sample_integrated_path <- file.path(cfg$source_dir, "div30_jia_progenitor_strong_integrated_seurat.rds")
if (!file.exists(uncorrected_path)) stop("Missing uncorrected object: ", uncorrected_path, call. = FALSE)
if (!file.exists(sample_integrated_path)) stop("Missing sample-integrated object: ", sample_integrated_path, call. = FALSE)

log_msg("Loading uncorrected object: ", uncorrected_path)
uncorrected <- readRDS(uncorrected_path)
log_msg("Loading sample-integrated object: ", sample_integrated_path)
sample_integrated <- readRDS(sample_integrated_path)

uncorrected <- attach_sample_map(uncorrected, cfg$sample_col, cfg$sample_map)
sample_integrated <- attach_sample_map(sample_integrated, cfg$sample_col, cfg$sample_map)

make_sample_only_outputs(uncorrected, sample_integrated, cfg)
make_sample_cellline_outputs(sample_integrated, cfg)

log_msg("Complete.")
log_msg("Sample-only outputs: ", cfg$sample_outdir)
log_msg("Sample + cell-line outputs: ", cfg$sample_cellline_outdir)
