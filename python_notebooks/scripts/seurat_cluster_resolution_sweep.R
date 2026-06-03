#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(ggplot2)
})

timestamp <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")
log_msg <- function(...) message("[R ", timestamp(), "] ", paste0(..., collapse = ""))

parse_args <- function(args) {
  out <- list(
    study_id = NULL,
    label = NULL,
    seurat = NULL,
    outdir = NULL,
    graph = "RNA_snn",
    reduction = "umap",
    resolutions = "0.1,0.2,0.3,0.4,0.5,0.6,0.8,1.0,1.2",
    algorithm = "1",
    random_seed = "0",
    overwrite = "false"
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

write_tsv <- function(df, path) {
  write.table(df, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
}

safe_res_tag <- function(resolution) {
  gsub("[^0-9A-Za-z]+", "p", sprintf("%.3f", as.numeric(resolution)))
}

parse_resolutions <- function(raw) {
  values <- strsplit(gsub("[:; ]+", ",", raw), ",", fixed = FALSE)[[1]]
  values <- values[nzchar(values)]
  resolutions <- as.numeric(values)
  if (any(is.na(resolutions))) stop("Invalid resolution list: ", raw, call. = FALSE)
  unique(resolutions)
}

opt <- parse_args(commandArgs(trailingOnly = TRUE))
required <- c("study_id", "seurat", "outdir")
missing <- required[vapply(required, function(name) is.null(opt[[name]]) || !nzchar(opt[[name]]), logical(1))]
if (length(missing) > 0) stop("Missing required argument(s): ", paste(missing, collapse = ", "), call. = FALSE)
if (!file.exists(opt$seurat)) stop("Seurat source does not exist: ", opt$seurat, call. = FALSE)

resolutions <- parse_resolutions(opt$resolutions)
algorithm <- as.integer(opt$algorithm)
random_seed <- as.integer(opt$random_seed)
overwrite <- tolower(opt$overwrite) %in% c("1", "true", "yes", "y")

dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)
plot_dir <- file.path(opt$outdir, "plots")
table_dir <- file.path(opt$outdir, "tables")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

complete_marker <- file.path(opt$outdir, "resolution_sweep_complete.tsv")
if (file.exists(complete_marker) && !overwrite) {
  log_msg("Using existing resolution sweep marker: ", complete_marker)
  quit(status = 0)
}

log_msg("R executable: ", R.home("bin"))
log_msg("R version: ", paste(R.version$major, R.version$minor, sep = "."))
log_msg("Seurat version: ", as.character(utils::packageVersion("Seurat")))
log_msg("SeuratObject version: ", as.character(utils::packageVersion("SeuratObject")))
log_msg("Reading Seurat RDS: ", opt$seurat)
obj <- readRDS(opt$seurat)
if (!inherits(obj, "Seurat")) stop("Object is not a Seurat object: ", opt$seurat, call. = FALSE)
log_msg("Loaded Seurat object with ", ncol(obj), " cells and ", nrow(obj), " features")

graph_names <- if ("graphs" %in% slotNames(obj)) names(slot(obj, "graphs")) else character(0)
reduction_names <- Seurat::Reductions(obj)
if (!(opt$graph %in% graph_names)) {
  stop("Missing graph '", opt$graph, "'. Available graphs: ", paste(graph_names, collapse = ", "), call. = FALSE)
}
if (!(opt$reduction %in% reduction_names)) {
  stop("Missing reduction '", opt$reduction, "'. Available reductions: ", paste(reduction_names, collapse = ", "), call. = FALSE)
}

run_parameters <- data.frame(
  key = c(
    "study_id",
    "label",
    "seurat_path",
    "outdir",
    "graph",
    "reduction",
    "resolutions",
    "algorithm",
    "random_seed",
    "n_cells",
    "n_features"
  ),
  value = c(
    opt$study_id,
    ifelse(is.null(opt$label), "", opt$label),
    normalizePath(opt$seurat, mustWork = FALSE),
    normalizePath(opt$outdir, mustWork = FALSE),
    opt$graph,
    opt$reduction,
    paste(resolutions, collapse = ","),
    as.character(algorithm),
    as.character(random_seed),
    as.character(ncol(obj)),
    as.character(nrow(obj))
  ),
  stringsAsFactors = FALSE
)
write_tsv(run_parameters, file.path(table_dir, "resolution_sweep_run_parameters.tsv"))

cell_ids <- colnames(obj)
assignment_wide <- data.frame(cell_id = cell_ids, stringsAsFactors = FALSE)
count_rows <- list()
plot_rows <- list()

for (resolution in resolutions) {
  res_tag <- safe_res_tag(resolution)
  cluster_col <- paste0("resolution_sweep_res_", res_tag)
  log_msg("Running FindClusters resolution=", resolution, " cluster_col=", cluster_col)

  obj <- Seurat::FindClusters(
    object = obj,
    graph.name = opt$graph,
    resolution = resolution,
    algorithm = algorithm,
    random.seed = random_seed,
    cluster.name = cluster_col,
    verbose = FALSE
  )

  clusters <- as.character(obj@meta.data[[cluster_col]])
  assignment_wide[[cluster_col]] <- clusters
  counts <- sort(table(clusters), decreasing = TRUE)
  count_rows[[length(count_rows) + 1L]] <- data.frame(
    study_id = opt$study_id,
    resolution = resolution,
    cluster_column = cluster_col,
    cluster = names(counts),
    n_cells = as.integer(counts),
    stringsAsFactors = FALSE
  )

  plot_path <- file.path(plot_dir, paste0("umap_", cluster_col, ".png"))
  log_msg("Saving UMAP plot: ", plot_path)
  p <- Seurat::DimPlot(
    obj,
    reduction = opt$reduction,
    group.by = cluster_col,
    label = TRUE,
    repel = TRUE,
    raster = TRUE
  ) +
    ggplot2::ggtitle(paste0(opt$study_id, " ", cluster_col)) +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 11))
  ggplot2::ggsave(plot_path, plot = p, width = 7, height = 6, dpi = 180)

  plot_rows[[length(plot_rows) + 1L]] <- data.frame(
    study_id = opt$study_id,
    resolution = resolution,
    cluster_column = cluster_col,
    plot_path = normalizePath(plot_path, mustWork = FALSE),
    stringsAsFactors = FALSE
  )
}

assignment_long <- do.call(
  rbind,
  lapply(names(assignment_wide)[names(assignment_wide) != "cell_id"], function(cluster_col) {
    resolution <- sub("^resolution_sweep_res_", "", cluster_col)
    data.frame(
      study_id = opt$study_id,
      cell_id = assignment_wide$cell_id,
      cluster_column = cluster_col,
      resolution_tag = resolution,
      cluster = assignment_wide[[cluster_col]],
      stringsAsFactors = FALSE
    )
  })
)

cluster_counts <- do.call(rbind, count_rows)
plot_manifest <- do.call(rbind, plot_rows)

write_tsv(assignment_wide, file.path(table_dir, "resolution_sweep_cluster_assignments_wide.tsv"))
write_tsv(assignment_long, file.path(table_dir, "resolution_sweep_cluster_assignments_long.tsv"))
write_tsv(cluster_counts, file.path(table_dir, "resolution_sweep_cluster_counts.tsv"))
write_tsv(plot_manifest, file.path(table_dir, "resolution_sweep_plot_manifest.tsv"))

write_tsv(
  data.frame(
    study_id = opt$study_id,
    completed_at = timestamp(),
    outdir = normalizePath(opt$outdir, mustWork = FALSE),
    n_resolutions = length(resolutions),
    n_plots = nrow(plot_manifest),
    stringsAsFactors = FALSE
  ),
  complete_marker
)

log_msg("Finished resolution sweep: ", opt$outdir)
