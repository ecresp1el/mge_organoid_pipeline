#!/usr/bin/env Rscript

# First-pass URD pseudotime for DIV30 cells.
#
# This script is deliberately organized as a small, auditable pipeline:
#
#   1. Read a plain Matrix Market input bundle produced by
#      python_notebooks/scripts/export_div30_first_urd_inputs.py.
#   2. Create an URD object with explicit gene/cell/count filters.
#   3. Select highly variable genes from URD's normalized logupx matrix.
#   4. Run URD PCA, diffusion map, and flood pseudotime.
#   5. Export pseudotime, imposed parameters, marker correlations, summaries,
#      plots, and the final URD object.
#
# No Seurat object is opened here. The root definition and all sidecar labels
# are inherited from the input metadata table, so the exact cells and labels used
# for lineage reconstruction are recoverable from the input bundle.

log_msg <- function(...) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste0(..., collapse = " ")))
  flush.console()
}

parse_args <- function(args) {
  out <- list(
    `project-root` = Sys.getenv("PROJECT_ROOT"),
    `input-dir` = NULL,
    outdir = NULL,
    `run-label` = Sys.getenv("RUN_LABEL", "div30_first_urd_paper_radial_glia_v1"),
    `root-label` = Sys.getenv("ROOT_LABEL", "Radial glia"),
    `pseudotime-name` = Sys.getenv("URD_PSEUDOTIME_NAME", "paper_radial_glia_root"),
    seed = Sys.getenv("URD_SEED", Sys.getenv("SEED", "7")),
    knn = Sys.getenv("URD_KNN", "100"),
    sigma = Sys.getenv("URD_SIGMA", "local"),
    `n-floods` = Sys.getenv("URD_N_FLOODS", "20"),
    `num-variable-genes` = Sys.getenv("URD_NUM_VARIABLE_GENES", "3000"),
    `min-genes` = Sys.getenv("URD_MIN_GENES", "500"),
    `min-cells` = Sys.getenv("URD_MIN_CELLS", "3"),
    `min-counts` = Sys.getenv("URD_MIN_COUNTS", "10"),
    `pca-mp-factor` = Sys.getenv("URD_PCA_MP_FACTOR", "2"),
    `flood-minimum-cells` = Sys.getenv("URD_FLOOD_MINIMUM_CELLS", "2"),
    `flood-max-frac-na` = Sys.getenv("URD_FLOOD_MAX_FRAC_NA", "0.4"),
    `flood-stability-div` = Sys.getenv("URD_FLOOD_STABILITY_DIV", ""),
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
    "  Rscript scripts/14_div30_first_urd.R --project-root <PROJECT_ROOT> [--input-dir <dir>] [--outdir <dir>]",
    "",
    "Inputs in --input-dir:",
    "  div30_first_urd_counts.mtx",
    "  div30_first_urd_features.tsv",
    "  div30_first_urd_cell_metadata.tsv",
    "",
    "Core imposed parameters:",
    "  --root-label <label>              Metadata label used to select root cells",
    "  --pseudotime-name <name>          Name assigned to flood pseudotime",
    "  --min-genes <int>                 URD createURD cell filter",
    "  --min-cells <int>                 URD createURD gene filter",
    "  --min-counts <int>                URD createURD gene/count filter",
    "  --num-variable-genes <int>        Top variable genes stored in urd@var.genes",
    "  --pca-mp-factor <numeric>         calcPCA Marchenko-Pastur factor",
    "  --knn <int>                       calcDM k nearest neighbors",
    "  --sigma <local|NULL|numeric>      calcDM sigma",
    "  --n-floods <int>                  floodPseudotime replicate count",
    "  --flood-minimum-cells <int>       floodPseudotime minimum.cells.flooded",
    "  --flood-max-frac-na <numeric>     floodPseudotimeProcess max.frac.NA",
    "  --flood-stability-div <int>       floodPseudotimeProcess stability.div",
    "",
    "Outputs under --outdir:",
    "  tables/div30_first_urd_parameters.tsv",
    "  tables/div30_first_urd_pseudotime.tsv",
    "  tables/div30_first_urd_summary.tsv",
    "  tables/div30_first_urd_marker_correlations.tsv",
    "  tables/div30_first_urd_pseudotime_by_paper_cluster.tsv",
    "  plots/*.png",
    "  div30_first_urd_object.rds",
    sep = "\n"
  ))
}

required <- c("Matrix", "URD", "ggplot2")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) {
  stop("Missing required R package(s): ", paste(missing, collapse = ", "), call. = FALSE)
}

suppressPackageStartupMessages({
  library(Matrix)
  library(URD)
  library(ggplot2)
})

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

parse_sigma <- function(sigma_text) {
  sigma_num <- suppressWarnings(as.numeric(sigma_text))
  if (!is.na(sigma_num)) return(sigma_num)
  if (toupper(sigma_text) == "NULL") return(NULL)
  sigma_text
}

build_config <- function(opt) {
  if (!nzchar(opt$`project-root`)) stop("PROJECT_ROOT or --project-root is required", call. = FALSE)
  project_root <- trim_trailing_slash(opt$`project-root`)
  run_label <- opt$`run-label`
  input_dir <- if (is.null(opt$`input-dir`) || !nzchar(opt$`input-dir`)) {
    file.path(project_root, "results/div30_first_urd", run_label, "inputs")
  } else {
    opt$`input-dir`
  }
  outdir <- if (is.null(opt$outdir) || !nzchar(opt$outdir)) {
    file.path(project_root, "results/div30_first_urd", run_label)
  } else {
    opt$outdir
  }
  n_floods <- as_int(opt$`n-floods`, "n-floods")
  stability_text <- opt$`flood-stability-div`
  stability_div <- if (nzchar(stability_text)) as_int(stability_text, "flood-stability-div") else min(10L, n_floods)
  list(
    project_root = project_root,
    input_dir = input_dir,
    outdir = outdir,
    plot_dir = file.path(outdir, "plots"),
    table_dir = file.path(outdir, "tables"),
    run_label = run_label,
    root_label = opt$`root-label`,
    seed = as_int(opt$seed, "seed"),
    min_genes = as_int(opt$`min-genes`, "min-genes"),
    min_cells = as_int(opt$`min-cells`, "min-cells"),
    min_counts = as_int(opt$`min-counts`, "min-counts"),
    num_variable_genes = as_int(opt$`num-variable-genes`, "num-variable-genes"),
    pca_mp_factor = as_num(opt$`pca-mp-factor`, "pca-mp-factor"),
    knn = as_int(opt$knn, "knn"),
    sigma_text = opt$sigma,
    sigma = parse_sigma(opt$sigma),
    n_floods = n_floods,
    flood_minimum_cells = as_int(opt$`flood-minimum-cells`, "flood-minimum-cells"),
    flood_max_frac_na = as_num(opt$`flood-max-frac-na`, "flood-max-frac-na"),
    flood_stability_div = stability_div,
    pseudotime_name = opt$`pseudotime-name`
  )
}

parameter_table <- function(cfg) {
  data.frame(
    stage = c(
      "input", "input", "input", "input",
      "root", "root",
      "filter", "filter", "filter",
      "variable_genes",
      "pca",
      "diffusion_map", "diffusion_map",
      "flood", "flood", "flood_process", "flood_process",
      "output"
    ),
    parameter = c(
      "project_root", "input_dir", "outdir", "run_label",
      "root_label", "pseudotime_name",
      "min_genes", "min_cells", "min_counts",
      "num_variable_genes",
      "pca_mp_factor",
      "knn", "sigma",
      "n_floods", "minimum_cells_flooded", "max_frac_NA", "stability_div",
      "seed"
    ),
    value = as.character(c(
      cfg$project_root, cfg$input_dir, cfg$outdir, cfg$run_label,
      cfg$root_label, cfg$pseudotime_name,
      cfg$min_genes, cfg$min_cells, cfg$min_counts,
      cfg$num_variable_genes,
      cfg$pca_mp_factor,
      cfg$knn, cfg$sigma_text,
      cfg$n_floods, cfg$flood_minimum_cells, cfg$flood_max_frac_na, cfg$flood_stability_div,
      cfg$seed
    )),
    meaning = c(
      "NFS project directory used for default paths.",
      "Matrix Market input bundle consumed by this R script.",
      "Run directory receiving tables, plots, and RDS output.",
      "Run label used to namespace inputs and outputs.",
      "Cells with urd_root_candidate TRUE are expected to have this metadata label.",
      "Name assigned to the URD pseudotime vector.",
      "createURD keeps cells with at least this many detected genes.",
      "createURD keeps genes detected in at least this many cells.",
      "createURD keeps genes with at least this many total counts.",
      "Number of highest-variance log-normalized genes stored in urd@var.genes.",
      "calcPCA Marchenko-Pastur multiplier used to estimate significant PCs.",
      "calcDM nearest-neighbor count.",
      "calcDM sigma setting; local means local adaptive sigma.",
      "Number of flood pseudotime replicates.",
      "floodPseudotime minimum.cells.flooded.",
      "floodPseudotimeProcess max.frac.NA.",
      "floodPseudotimeProcess stability.div; defaults to min(10, n_floods).",
      "Random seed set before PCA, diffusion map, and flood pseudotime."
    ),
    stringsAsFactors = FALSE
  )
}

input_paths <- function(input_dir) {
  list(
    counts = file.path(input_dir, "div30_first_urd_counts.mtx"),
    features = file.path(input_dir, "div30_first_urd_features.tsv"),
    metadata = file.path(input_dir, "div30_first_urd_cell_metadata.tsv")
  )
}

read_urd_input_bundle <- function(cfg) {
  paths <- input_paths(cfg$input_dir)
  missing_files <- unlist(paths)[!file.exists(unlist(paths))]
  if (length(missing_files) > 0) {
    stop("Missing URD input file(s): ", paste(missing_files, collapse = ", "), call. = FALSE)
  }

  log_msg("Reading counts: ", paths$counts)
  counts <- Matrix::readMM(paths$counts)
  if (!inherits(counts, "dgCMatrix")) counts <- as(counts, "dgCMatrix")

  features <- read.delim(paths$features, stringsAsFactors = FALSE, check.names = FALSE)
  metadata <- read.delim(paths$metadata, stringsAsFactors = FALSE, check.names = FALSE)
  if (!("feature_id" %in% colnames(features))) stop("Feature table needs feature_id column", call. = FALSE)
  if (!("cell_id" %in% colnames(metadata))) stop("Metadata table needs cell_id column", call. = FALSE)
  if (!("urd_root_candidate" %in% colnames(metadata))) {
    stop("Metadata table needs urd_root_candidate column", call. = FALSE)
  }
  if (nrow(counts) != nrow(features)) stop("Counts rows do not match features", call. = FALSE)
  if (ncol(counts) != nrow(metadata)) stop("Counts columns do not match metadata rows", call. = FALSE)

  rownames(counts) <- make.unique(as.character(features$feature_id))
  colnames(counts) <- as.character(metadata$cell_id)
  rownames(metadata) <- metadata$cell_id
  metadata$urd_root_candidate <- as.logical(metadata$urd_root_candidate)
  if (anyNA(metadata$urd_root_candidate)) {
    stop("urd_root_candidate contains NA after logical conversion", call. = FALSE)
  }

  list(counts = counts, features = features, metadata = metadata, paths = paths)
}

create_filtered_urd <- function(bundle, cfg) {
  log_msg(
    "Creating URD object with min.genes=", cfg$min_genes,
    " min.cells=", cfg$min_cells,
    " min.counts=", cfg$min_counts
  )
  urd <- createURD(
    count.data = bundle$counts,
    meta = bundle$metadata,
    min.cells = cfg$min_cells,
    min.genes = cfg$min_genes,
    min.counts = cfg$min_counts,
    verbose = TRUE
  )

  cells <- rownames(urd@meta)
  root_cells <- cells[as.logical(urd@meta[cells, "urd_root_candidate"])]
  if (length(root_cells) == 0) {
    stop("No root candidate cells remained after URD filtering", call. = FALSE)
  }
  log_msg("Root cells retained after filtering: ", length(root_cells))
  list(urd = urd, root_cells = root_cells)
}

select_variable_genes <- function(urd, cfg) {
  log_msg("Selecting top variable genes from URD logupx.data")
  log_data <- urd@logupx.data
  gene_means <- Matrix::rowMeans(log_data)
  log_data_sq <- log_data
  log_data_sq@x <- log_data_sq@x^2
  gene_vars <- Matrix::rowMeans(log_data_sq) - gene_means^2
  gene_vars <- gene_vars[is.finite(gene_vars)]
  gene_vars <- sort(gene_vars, decreasing = TRUE)
  urd@var.genes <- names(gene_vars)[seq_len(min(cfg$num_variable_genes, length(gene_vars)))]
  log_msg("Variable genes stored in URD object: ", length(urd@var.genes))
  urd
}

run_urd_geometry <- function(urd, cfg) {
  set.seed(cfg$seed)
  log_msg("Running calcPCA with mp.factor=", cfg$pca_mp_factor)
  urd <- calcPCA(urd, mp.factor = cfg$pca_mp_factor)
  log_msg("Running calcDM with knn=", cfg$knn, " sigma=", cfg$sigma_text)
  calcDM(urd, knn = cfg$knn, sigma = cfg$sigma)
}

run_flood_pseudotime <- function(urd, root_cells, cfg) {
  set.seed(cfg$seed)
  log_msg(
    "Running floodPseudotime n=", cfg$n_floods,
    " minimum.cells.flooded=", cfg$flood_minimum_cells
  )
  floods <- floodPseudotime(
    urd,
    root.cells = root_cells,
    n = cfg$n_floods,
    minimum.cells.flooded = cfg$flood_minimum_cells,
    verbose = TRUE
  )
  log_msg(
    "Processing flood pseudotime max.frac.NA=", cfg$flood_max_frac_na,
    " stability.div=", cfg$flood_stability_div
  )
  floodPseudotimeProcess(
    urd,
    floods,
    floods.name = cfg$pseudotime_name,
    max.frac.NA = cfg$flood_max_frac_na,
    pseudotime.fun = mean,
    stability.div = cfg$flood_stability_div
  )
}

extract_pseudotime <- function(object, name) {
  pt <- object@pseudotime
  if (is.data.frame(pt) || is.matrix(pt)) return(as.numeric(pt[rownames(object@meta), name]))
  if (is.list(pt) && name %in% names(pt)) return(as.numeric(pt[[name]][rownames(object@meta)]))
  stop("Could not extract pseudotime named ", name, call. = FALSE)
}

find_gene <- function(gene, available) {
  hit <- which(toupper(available) == toupper(gene))
  if (length(hit) == 0) return(NA_integer_)
  hit[[1]]
}

build_pseudotime_table <- function(urd, cfg) {
  pt_values <- extract_pseudotime(urd, cfg$pseudotime_name)
  out_meta <- urd@meta[rownames(urd@meta), , drop = FALSE]
  pt_df <- data.frame(
    cell_id = rownames(out_meta),
    urd_pseudotime = pt_values,
    urd_pseudotime_name = cfg$pseudotime_name,
    urd_root_label = cfg$root_label,
    stringsAsFactors = FALSE
  )
  pt_df <- cbind(pt_df, out_meta)

  marker_genes <- c("DLX2", "ASCL1", "DCX", "STMN2", "TUBB3", "MAP2", "RBFOX3", "SYT1", "SNAP25")
  available_genes <- rownames(urd@logupx.data)
  for (gene in marker_genes) {
    idx <- find_gene(gene, available_genes)
    col <- paste0("logupx_", gene)
    if (is.na(idx)) {
      pt_df[[col]] <- NA_real_
    } else {
      pt_df[[col]] <- as.numeric(urd@logupx.data[idx, pt_df$cell_id])
    }
  }
  neuron_cols <- paste0("logupx_", c("DCX", "STMN2", "TUBB3", "MAP2", "RBFOX3", "SYT1", "SNAP25"))
  neuron_cols <- neuron_cols[neuron_cols %in% colnames(pt_df)]
  pt_df$neuronal_maturation_score <- rowMeans(pt_df[, neuron_cols, drop = FALSE], na.rm = TRUE)
  pt_df
}

safe_cor <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3) return(NA_real_)
  suppressWarnings(cor(x[ok], y[ok], method = "spearman"))
}

build_correlation_table <- function(pt_df) {
  score_cols <- intersect(
    c(
      "jia_score_RGC1",
      "jia_score_RGC2",
      "jia_score_IPC",
      "logupx_DLX2",
      "logupx_ASCL1",
      "logupx_DCX",
      "neuronal_maturation_score",
      "shi_seurat_full_prediction_score_MGE",
      "shi_seurat_full_prediction_score_progenitor",
      "shi_seurat_full_expected_shi_week_numeric"
    ),
    colnames(pt_df)
  )
  data.frame(
    feature = score_cols,
    spearman_with_urd_pseudotime = vapply(
      score_cols,
      function(col) safe_cor(pt_df$urd_pseudotime, as.numeric(pt_df[[col]])),
      numeric(1)
    ),
    stringsAsFactors = FALSE
  )
}

build_cluster_summary <- function(pt_df) {
  if (!("paper_cluster_annotation" %in% colnames(pt_df))) {
    return(data.frame())
  }
  cluster_summary <- aggregate(
    urd_pseudotime ~ paper_cluster_annotation,
    data = pt_df,
    FUN = function(x) median(x, na.rm = TRUE)
  )
  colnames(cluster_summary)[2] <- "median_urd_pseudotime"
  cluster_summary$n_cells <- as.integer(table(pt_df$paper_cluster_annotation)[cluster_summary$paper_cluster_annotation])
  cluster_summary
}

write_tables <- function(urd, pt_df, root_cells, cfg) {
  pseudotime_path <- file.path(cfg$table_dir, "div30_first_urd_pseudotime.tsv")
  parameter_path <- file.path(cfg$table_dir, "div30_first_urd_parameters.tsv")
  cor_path <- file.path(cfg$table_dir, "div30_first_urd_marker_correlations.tsv")
  cluster_path <- file.path(cfg$table_dir, "div30_first_urd_pseudotime_by_paper_cluster.tsv")
  summary_path <- file.path(cfg$table_dir, "div30_first_urd_summary.tsv")

  write.table(parameter_table(cfg), parameter_path, sep = "\t", row.names = FALSE, quote = FALSE)
  write.table(pt_df, pseudotime_path, sep = "\t", row.names = FALSE, quote = FALSE)
  write.table(build_correlation_table(pt_df), cor_path, sep = "\t", row.names = FALSE, quote = FALSE)
  write.table(build_cluster_summary(pt_df), cluster_path, sep = "\t", row.names = FALSE, quote = FALSE)

  summary_df <- data.frame(
    key = c(
      "run_label",
      "input_dir",
      "root_label",
      "n_cells_urd",
      "n_root_cells_urd",
      "n_variable_genes",
      "knn",
      "sigma",
      "n_floods",
      "parameter_tsv",
      "pseudotime_tsv",
      "correlation_tsv",
      "cluster_summary_tsv"
    ),
    value = c(
      cfg$run_label,
      cfg$input_dir,
      cfg$root_label,
      nrow(pt_df),
      length(root_cells),
      length(urd@var.genes),
      cfg$knn,
      cfg$sigma_text,
      cfg$n_floods,
      parameter_path,
      pseudotime_path,
      cor_path,
      cluster_path
    )
  )
  write.table(summary_df, summary_path, sep = "\t", row.names = FALSE, quote = FALSE)
  list(
    parameter = parameter_path,
    pseudotime = pseudotime_path,
    correlation = cor_path,
    cluster = cluster_path,
    summary = summary_path
  )
}

write_plots <- function(pt_df, cfg) {
  if (!all(c("UMAP_1", "UMAP_2") %in% colnames(pt_df))) return(invisible(NULL))
  plot_base <- pt_df[is.finite(pt_df$UMAP_1) & is.finite(pt_df$UMAP_2), , drop = FALSE]
  if (nrow(plot_base) == 0) return(invisible(NULL))

  png(file.path(cfg$plot_dir, "div30_first_urd_umap_pseudotime.png"), width = 1800, height = 1500, res = 220)
  print(
    ggplot(plot_base, aes(x = UMAP_1, y = UMAP_2, color = urd_pseudotime)) +
      geom_point(size = 0.25, alpha = 0.8) +
      coord_equal() +
      scale_color_viridis_c(na.value = "grey85") +
      theme_void(base_size = 10) +
      labs(color = "URD pseudotime", title = "DIV30 first URD: paper Radial glia root")
  )
  dev.off()

  if ("paper_cluster_annotation" %in% colnames(pt_df)) {
    png(file.path(cfg$plot_dir, "div30_first_urd_pseudotime_by_paper_cluster.png"), width = 1800, height = 1200, res = 220)
    print(
      ggplot(pt_df, aes(x = paper_cluster_annotation, y = urd_pseudotime, fill = paper_cluster_annotation)) +
        geom_boxplot(outlier.size = 0.2) +
        theme_bw(base_size = 10) +
        theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "none") +
        labs(x = NULL, y = "URD pseudotime", title = "DIV30 first URD pseudotime by paper/manual annotation")
    )
    dev.off()
  }
  invisible(NULL)
}

run_pipeline <- function(cfg) {
  dir.create(cfg$plot_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(cfg$table_dir, recursive = TRUE, showWarnings = FALSE)

  log_msg("Writing imposed parameter table")
  write.table(
    parameter_table(cfg),
    file.path(cfg$table_dir, "div30_first_urd_parameters.tsv"),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )

  bundle <- read_urd_input_bundle(cfg)
  created <- create_filtered_urd(bundle, cfg)
  rm(bundle)
  invisible(gc())

  urd <- select_variable_genes(created$urd, cfg)
  urd <- run_urd_geometry(urd, cfg)
  urd <- run_flood_pseudotime(urd, created$root_cells, cfg)

  pt_df <- build_pseudotime_table(urd, cfg)
  paths <- write_tables(urd, pt_df, created$root_cells, cfg)
  write_plots(pt_df, cfg)

  rds_path <- file.path(cfg$outdir, "div30_first_urd_object.rds")
  saveRDS(urd, rds_path)
  log_msg("Wrote parameters: ", paths$parameter)
  log_msg("Wrote pseudotime: ", paths$pseudotime)
  log_msg("Wrote summary: ", paths$summary)
  log_msg("Saved URD object: ", rds_path)
  invisible(list(urd = urd, pseudotime = pt_df, paths = paths, rds = rds_path))
}

opt <- parse_args(commandArgs(trailingOnly = TRUE))
if (isTRUE(opt$help)) {
  print_usage()
  quit(save = "no", status = 0)
}

cfg <- build_config(opt)
run_pipeline(cfg)
