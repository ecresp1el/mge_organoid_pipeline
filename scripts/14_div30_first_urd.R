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
#   5. Save resumable URD checkpoints after expensive geometry stages.
#   6. Export pseudotime, imposed parameters, marker correlations, summaries,
#      plots, and the final URD object.
#
# No Seurat object is opened here. The root definition and all sidecar labels
# are inherited from the input metadata table, so the exact cells and labels used
# for lineage reconstruction are recoverable from the input bundle.

log_msg <- function(...) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste0(..., collapse = " ")))
  flush.console()
}

# Stage timings are written after every completed/failed major step. This is
# separate from the final parameter table so partial long runs still leave a
# progress trail if they stop before the URD object is saved.
.stage_timings <- data.frame(
  stage = character(),
  status = character(),
  start_time = character(),
  end_time = character(),
  elapsed_seconds = numeric(),
  memory = character(),
  stringsAsFactors = FALSE
)

memory_summary <- function() {
  stats <- gc()
  paste(
    sprintf("Ncells_used=%s", format(stats["Ncells", "used"], scientific = FALSE)),
    sprintf("Vcells_used_mb=%.1f", stats["Vcells", "used"] * 8 / 1024^2),
    sprintf("Vcells_max_mb=%.1f", stats["Vcells", "max used"] * 8 / 1024^2),
    sep = " "
  )
}

write_stage_timings <- function(table_dir) {
  if (!dir.exists(table_dir)) dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
  write.table(
    .stage_timings,
    file.path(table_dir, "div30_first_urd_stage_timings.tsv"),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
}

run_stage <- function(stage, table_dir, expr) {
  started <- Sys.time()
  log_msg("STAGE_START ", stage, " memory=", memory_summary())
  result <- tryCatch(
    force(expr),
    error = function(e) {
      ended <- Sys.time()
      .stage_timings <<- rbind(
        .stage_timings,
        data.frame(
          stage = stage,
          status = "failed",
          start_time = format(started, "%Y-%m-%d %H:%M:%S"),
          end_time = format(ended, "%Y-%m-%d %H:%M:%S"),
          elapsed_seconds = as.numeric(difftime(ended, started, units = "secs")),
          memory = memory_summary(),
          stringsAsFactors = FALSE
        )
      )
      write_stage_timings(table_dir)
      log_msg("STAGE_FAIL ", stage, " elapsed_seconds=", round(as.numeric(difftime(ended, started, units = "secs")), 2), " error=", conditionMessage(e))
      stop(e)
    }
  )
  ended <- Sys.time()
  .stage_timings <<- rbind(
    .stage_timings,
    data.frame(
      stage = stage,
      status = "completed",
      start_time = format(started, "%Y-%m-%d %H:%M:%S"),
      end_time = format(ended, "%Y-%m-%d %H:%M:%S"),
      elapsed_seconds = as.numeric(difftime(ended, started, units = "secs")),
      memory = memory_summary(),
      stringsAsFactors = FALSE
    )
  )
  write_stage_timings(table_dir)
  log_msg("STAGE_END ", stage, " elapsed_seconds=", round(as.numeric(difftime(ended, started, units = "secs")), 2), " memory=", memory_summary())
  result
}

as_bool <- function(x, name) {
  value <- tolower(trimws(as.character(x)))
  if (value %in% c("true", "t", "1", "yes", "y")) return(TRUE)
  if (value %in% c("false", "f", "0", "no", "n")) return(FALSE)
  stop(name, " must be true/false; got ", x, call. = FALSE)
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
    resume = Sys.getenv("URD_RESUME", "true"),
    `force-recompute` = Sys.getenv("URD_FORCE_RECOMPUTE", "false"),
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
    "  --resume <true|false>             Reuse saved checkpoints when present",
    "  --force-recompute <true|false>    Ignore checkpoints and recompute stages",
    "",
    "Outputs under --outdir:",
    "  tables/div30_first_urd_parameters.tsv",
    "  tables/div30_first_urd_pseudotime.tsv",
    "  tables/div30_first_urd_summary.tsv",
    "  tables/div30_first_urd_marker_correlations.tsv",
    "  tables/div30_first_urd_pseudotime_by_paper_cluster.tsv",
    "  plots/*.png",
    "  tables/div30_first_urd_stage_timings.tsv",
    "  tables/div30_first_urd_checkpoint_manifest.tsv",
    "  checkpoints/urd_after_*.rds",
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
    checkpoint_dir = file.path(outdir, "checkpoints"),
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
    pseudotime_name = opt$`pseudotime-name`,
    resume = as_bool(opt$resume, "resume"),
    force_recompute = as_bool(opt$`force-recompute`, "force-recompute")
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
      "checkpoint", "checkpoint",
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
      "resume", "force_recompute",
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
      cfg$resume, cfg$force_recompute,
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
      "Whether existing stage checkpoints should be reused if present.",
      "Whether to ignore existing stage checkpoints and recompute from inputs.",
      "Random seed set before PCA, diffusion map, and flood pseudotime."
    ),
    stringsAsFactors = FALSE
  )
}

checkpoint_path <- function(cfg, stage) {
  file.path(cfg$checkpoint_dir, paste0("urd_after_", stage, ".rds"))
}

checkpoint_manifest_path <- function(cfg) {
  file.path(cfg$table_dir, "div30_first_urd_checkpoint_manifest.tsv")
}

write_checkpoint_manifest <- function(cfg) {
  stages <- c("filter", "variable_genes", "pca", "diffusion_map", "flood_pseudotime", "final")
  paths <- c(
    checkpoint_path(cfg, "filter"),
    checkpoint_path(cfg, "variable_genes"),
    checkpoint_path(cfg, "pca"),
    checkpoint_path(cfg, "diffusion_map"),
    checkpoint_path(cfg, "flood_pseudotime"),
    file.path(cfg$outdir, "div30_first_urd_object.rds")
  )
  info <- lapply(paths, function(path) {
    if (!file.exists(path)) {
      return(data.frame(path = path, exists = FALSE, size_bytes = NA_real_, modified_time = NA_character_))
    }
    file_info <- file.info(path)
    data.frame(
      path = path,
      exists = TRUE,
      size_bytes = as.numeric(file_info$size),
      modified_time = format(file_info$mtime, "%Y-%m-%d %H:%M:%S"),
      stringsAsFactors = FALSE
    )
  })
  manifest <- data.frame(stage = stages, do.call(rbind, info), row.names = NULL, check.names = FALSE)
  write.table(manifest, checkpoint_manifest_path(cfg), sep = "\t", row.names = FALSE, quote = FALSE)
  invisible(manifest)
}

save_checkpoint <- function(urd, cfg, stage) {
  dir.create(cfg$checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
  path <- checkpoint_path(cfg, stage)
  tmp_path <- paste0(path, ".tmp")
  log_msg("Saving checkpoint ", stage, ": ", path)
  saveRDS(urd, tmp_path)
  if (!file.rename(tmp_path, path)) {
    stop("Failed to move checkpoint into place: ", tmp_path, " -> ", path, call. = FALSE)
  }
  write_checkpoint_manifest(cfg)
  invisible(path)
}

load_checkpoint <- function(cfg, stage) {
  path <- checkpoint_path(cfg, stage)
  if (!cfg$resume || cfg$force_recompute || !file.exists(path)) return(NULL)
  log_msg("Loading checkpoint ", stage, ": ", path)
  readRDS(path)
}

root_cells_from_urd <- function(urd) {
  cells <- rownames(urd@meta)
  if (!("urd_root_candidate" %in% colnames(urd@meta))) {
    stop("URD metadata does not contain urd_root_candidate", call. = FALSE)
  }
  root_cells <- cells[as.logical(urd@meta[cells, "urd_root_candidate"])]
  if (length(root_cells) == 0) {
    stop("No root candidate cells are present in the URD object", call. = FALSE)
  }
  root_cells
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

run_urd_pca <- function(urd, cfg) {
  set.seed(cfg$seed)
  log_msg("Running calcPCA with mp.factor=", cfg$pca_mp_factor)
  calcPCA(urd, mp.factor = cfg$pca_mp_factor)
}

run_urd_diffusion_map <- function(urd, cfg) {
  set.seed(cfg$seed)
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
  dir.create(cfg$checkpoint_dir, recursive = TRUE, showWarnings = FALSE)

  log_msg("Writing imposed parameter table")
  write.table(
    parameter_table(cfg),
    file.path(cfg$table_dir, "div30_first_urd_parameters.tsv"),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
  write_checkpoint_manifest(cfg)

  final_rds_path <- file.path(cfg$outdir, "div30_first_urd_object.rds")
  if (cfg$resume && !cfg$force_recompute && file.exists(final_rds_path)) {
    log_msg("Loading completed URD object because resume=true: ", final_rds_path)
    urd <- readRDS(final_rds_path)
  } else {
    urd <- load_checkpoint(cfg, "flood_pseudotime")
    if (is.null(urd)) {
      urd <- load_checkpoint(cfg, "diffusion_map")
    }
    if (is.null(urd)) {
      urd <- load_checkpoint(cfg, "pca")
      if (is.null(urd)) {
        urd <- load_checkpoint(cfg, "variable_genes")
        if (is.null(urd)) {
          urd <- load_checkpoint(cfg, "filter")
          if (is.null(urd)) {
            bundle <- run_stage("read_input_bundle", cfg$table_dir, read_urd_input_bundle(cfg))
            created <- run_stage("create_filtered_urd", cfg$table_dir, create_filtered_urd(bundle, cfg))
            urd <- created$urd
            rm(bundle, created)
            invisible(gc())
            run_stage("checkpoint_after_filter", cfg$table_dir, save_checkpoint(urd, cfg, "filter"))
          }
          urd <- run_stage("select_variable_genes", cfg$table_dir, select_variable_genes(urd, cfg))
          run_stage("checkpoint_after_variable_genes", cfg$table_dir, save_checkpoint(urd, cfg, "variable_genes"))
        }
        urd <- run_stage("calc_pca", cfg$table_dir, run_urd_pca(urd, cfg))
        run_stage("checkpoint_after_pca", cfg$table_dir, save_checkpoint(urd, cfg, "pca"))
      }
      urd <- run_stage("calc_diffusion_map", cfg$table_dir, run_urd_diffusion_map(urd, cfg))
      run_stage("checkpoint_after_diffusion_map", cfg$table_dir, save_checkpoint(urd, cfg, "diffusion_map"))
    }
    root_cells <- root_cells_from_urd(urd)
    urd <- run_stage("flood_pseudotime", cfg$table_dir, run_flood_pseudotime(urd, root_cells, cfg))
    run_stage("checkpoint_after_flood_pseudotime", cfg$table_dir, save_checkpoint(urd, cfg, "flood_pseudotime"))
  }

  root_cells <- root_cells_from_urd(urd)
  pt_df <- run_stage("build_pseudotime_table", cfg$table_dir, build_pseudotime_table(urd, cfg))
  paths <- run_stage("write_tables", cfg$table_dir, write_tables(urd, pt_df, root_cells, cfg))
  run_stage("write_plots", cfg$table_dir, write_plots(pt_df, cfg))

  rds_path <- final_rds_path
  run_stage("save_urd_object", cfg$table_dir, saveRDS(urd, rds_path))
  write_checkpoint_manifest(cfg)
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
