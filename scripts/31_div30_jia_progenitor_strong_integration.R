#!/usr/bin/env Rscript

# DIV30 Jia-style progenitor reclustering with sample-aware integration.
#
# This is a follow-up to scripts/28_div30_jia_progenitor_reclustering.R.
# It keeps the same progenitor-only input cells, creates an uncorrected baseline
# UMAP, then runs a stronger Seurat CCA anchor integration across samples/batches
# so the before/after UMAPs can be reviewed side by side.

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
    `seurat-rds` = NULL,
    `jia-scores` = NULL,
    outdir = NULL,
    `run-label` = Sys.getenv("RUN_LABEL", "div30_jia_progenitor_strong_integration_v1"),
    assay = Sys.getenv("SEURAT_ASSAY", "RNA"),
    `source-cluster-col` = Sys.getenv("SOURCE_CLUSTER_COL", ""),
    `progenitor-clusters` = Sys.getenv("PROGENITOR_CLUSTERS", "0,3,6,7"),
    `batch-col` = Sys.getenv("BATCH_COL", "orig.ident"),
    `nfeatures` = Sys.getenv("NFEATURES", "3000"),
    `integration-nfeatures` = Sys.getenv("INTEGRATION_NFEATURES", "5000"),
    dims = Sys.getenv("DIMS", "50"),
    npcs = Sys.getenv("NPCS", "60"),
    resolution = Sys.getenv("RESOLUTION", "0.8"),
    seed = Sys.getenv("SEED", "7"),
    `integration-reduction` = Sys.getenv("INTEGRATION_REDUCTION", "cca"),
    `k-anchor` = Sys.getenv("K_ANCHOR", "10"),
    `k-weight` = Sys.getenv("K_WEIGHT", "200"),
    `future-globals-max-gb` = Sys.getenv("FUTURE_GLOBALS_MAX_GB", "120"),
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
    "  Rscript scripts/31_div30_jia_progenitor_strong_integration.R --project-root <PROJECT_ROOT> [options]",
    "",
    "Purpose:",
    "  Recluster the same DIV30 progenitor cells before and after strong sample-aware integration.",
    "",
    "Key defaults:",
    "  --batch-col orig.ident",
    "  --progenitor-clusters 0,3,6,7",
    "  --integration-nfeatures 5000",
    "  --integration-reduction cca",
    "  --dims 50",
    "  --k-anchor 10",
    "  --k-weight 200",
    "  --future-globals-max-gb 120",
    sep = "\n"
  ))
}

trim_trailing_slash <- function(x) sub("/+$", "", x)
split_csv <- function(x) {
  vals <- trimws(strsplit(x, ",", fixed = TRUE)[[1]])
  vals[nzchar(vals)]
}
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
read_tsv <- function(path) {
  con <- if (grepl("\\.gz$", path)) gzfile(path, open = "rt") else file(path, open = "rt")
  on.exit(close(con), add = TRUE)
  read.delim(con, stringsAsFactors = FALSE, check.names = FALSE)
}
save_plot_pair <- function(plot, stem, width, height, dpi = 300) {
  ggplot2::ggsave(paste0(stem, ".png"), plot = plot, width = width, height = height, dpi = dpi, bg = "white")
  ggplot2::ggsave(paste0(stem, ".pdf"), plot = plot, width = width, height = height, bg = "white")
}
safe_name <- function(x) gsub("[^A-Za-z0-9_+-]+", "_", x)

join_assay_layers_if_available <- function(obj, assay) {
  if (exists("JoinLayers", where = asNamespace("SeuratObject"), inherits = FALSE)) {
    log_msg("Joining Seurat v5 assay layers for assay ", assay)
    return(SeuratObject::JoinLayers(obj, assay = assay))
  }
  if (exists("JoinLayers", where = asNamespace("Seurat"), inherits = FALSE)) {
    log_msg("Joining Seurat v5 assay layers for assay ", assay)
    return(Seurat::JoinLayers(obj, assay = assay))
  }
  obj
}

find_cluster_col <- function(obj, requested) {
  meta_cols <- colnames(obj@meta.data)
  if (nzchar(requested)) {
    if (!(requested %in% meta_cols)) stop("Requested source cluster column not found: ", requested, call. = FALSE)
    return(requested)
  }
  for (candidate in c("seurat_clusters", "RNA_snn_res.0.2")) {
    if (candidate %in% meta_cols) return(candidate)
  }
  stop("No source cluster column found. Tried seurat_clusters and RNA_snn_res.0.2.", call. = FALSE)
}

find_batch_col <- function(obj, requested) {
  meta_cols <- colnames(obj@meta.data)
  if (nzchar(requested) && !toupper(requested) %in% c("AUTO", "AUTOMATIC")) {
    if (!(requested %in% meta_cols)) stop("Requested batch column not found: ", requested, call. = FALSE)
    return(requested)
  }
  for (candidate in c("orig.ident", "sample_id", "sample", "Sample", "batch", "cell_line", "line", "donor")) {
    if (candidate %in% meta_cols && length(unique(stats::na.omit(obj@meta.data[[candidate]]))) > 1L) {
      return(candidate)
    }
  }
  stop("No usable batch/sample column found.", call. = FALSE)
}

score_columns <- function() c("jia_score_RGC1", "jia_score_RGC2", "jia_score_IPC")

build_config <- function(opt) {
  if (!nzchar(opt$`project-root`)) stop("PROJECT_ROOT or --project-root is required", call. = FALSE)
  project_root <- trim_trailing_slash(opt$`project-root`)
  seurat_rds <- if (is.null(opt$`seurat-rds`) || !nzchar(opt$`seurat-rds`)) {
    file.path(project_root, "results/varela_this_paper/varela_this_paper_seurat.rds")
  } else {
    opt$`seurat-rds`
  }
  jia_scores <- if (is.null(opt$`jia-scores`) || !nzchar(opt$`jia-scores`)) {
    file.path(project_root, "results/jia_program_div30_scoring/jia_program_div30_scoring_v1/tables/div30_jia_program_scores_obs.tsv")
  } else {
    opt$`jia-scores`
  }
  outdir <- if (is.null(opt$outdir) || !nzchar(opt$outdir)) {
    file.path(project_root, "results/div30_jia_progenitor_reclustering", opt$`run-label`)
  } else {
    opt$outdir
  }
  dims <- as_int(opt$dims, "dims")
  npcs <- max(as_int(opt$npcs, "npcs"), dims)
  list(
    project_root = project_root,
    seurat_rds = seurat_rds,
    jia_scores = jia_scores,
    outdir = outdir,
    table_dir = file.path(outdir, "tables"),
    plot_dir = file.path(outdir, "plots"),
    run_label = opt$`run-label`,
    assay = opt$assay,
    source_cluster_col_requested = opt$`source-cluster-col`,
    progenitor_clusters = split_csv(opt$`progenitor-clusters`),
    batch_col_requested = opt$`batch-col`,
    nfeatures = as_int(opt$nfeatures, "nfeatures"),
    integration_nfeatures = as_int(opt$`integration-nfeatures`, "integration-nfeatures"),
    dims = dims,
    npcs = npcs,
    resolution = as_num(opt$resolution, "resolution"),
    seed = as_int(opt$seed, "seed"),
    integration_reduction = opt$`integration-reduction`,
    k_anchor = as_int(opt$`k-anchor`, "k-anchor"),
    k_weight = as_int(opt$`k-weight`, "k-weight"),
    future_globals_max_gb = as_num(opt$`future-globals-max-gb`, "future-globals-max-gb"),
    uncorrected_cluster_col = "div30_progenitor_uncorrected_cluster",
    integrated_cluster_col = "div30_progenitor_integrated_cluster",
    uncorrected_umap = "uncorrected_umap",
    integrated_umap = "integrated_umap"
  )
}

parameter_table <- function(cfg, source_cluster_col, batch_col, n_parent_cells, n_prog_cells, effective_k_weight) {
  data.frame(
    parameter = c(
      "project_root", "seurat_rds", "jia_scores", "outdir", "run_label",
      "assay", "source_cluster_col", "batch_col", "progenitor_clusters",
      "n_parent_cells", "n_progenitor_cells", "uncorrected_workflow",
      "integration_method", "integration_reduction", "integration_nfeatures",
      "npcs", "dims", "resolution", "k_anchor", "k_weight_requested",
      "k_weight_effective", "future_globals_max_gb", "seed", "jia_scores_use"
    ),
    value = as.character(c(
      cfg$project_root, cfg$seurat_rds, cfg$jia_scores, cfg$outdir, cfg$run_label,
      cfg$assay, source_cluster_col, batch_col, paste(cfg$progenitor_clusters, collapse = ","),
      n_parent_cells, n_prog_cells,
      paste0("NormalizeData; FindVariableFeatures nfeatures=", cfg$nfeatures,
             "; ScaleData no regression; PCA/UMAP dims=1:", cfg$dims),
      "Seurat FindIntegrationAnchors + IntegrateData across batch_col",
      cfg$integration_reduction, cfg$integration_nfeatures,
      cfg$npcs, paste0("1:", cfg$dims), cfg$resolution, cfg$k_anchor,
      cfg$k_weight, effective_k_weight, cfg$future_globals_max_gb, cfg$seed,
      "Post hoc metadata overlays only; not used for integration or clustering"
    )),
    stringsAsFactors = FALSE
  )
}

subset_progenitors <- function(obj, cfg, source_cluster_col) {
  source_clusters <- as.character(obj@meta.data[[source_cluster_col]])
  names(source_clusters) <- colnames(obj)
  keep_cells <- names(source_clusters)[source_clusters %in% cfg$progenitor_clusters]
  if (length(keep_cells) == 0L) {
    stop("No cells found for progenitor clusters: ", paste(cfg$progenitor_clusters, collapse = ","), call. = FALSE)
  }
  log_msg("Subsetting progenitor clusters ", paste(cfg$progenitor_clusters, collapse = ","), ": ", length(keep_cells), " cells")
  prog <- subset(obj, cells = keep_cells)
  prog <- join_assay_layers_if_available(prog, cfg$assay)
  prog$div30_parent_cluster <- as.character(prog@meta.data[[source_cluster_col]])
  prog$div30_parent_progenitor_label <- ifelse(
    prog$div30_parent_cluster == "6",
    "Inhibitory progenitors",
    "Radial glia"
  )
  DefaultAssay(prog) <- cfg$assay
  prog
}

attach_jia_scores <- function(prog, cfg) {
  if (!file.exists(cfg$jia_scores)) {
    log_msg("Jia score table not found; continuing without Jia overlays: ", cfg$jia_scores)
    return(list(object = prog, score_cols = character()))
  }
  scores <- read_tsv(cfg$jia_scores)
  if (!("cell_id" %in% colnames(scores))) stop("Jia score table lacks cell_id: ", cfg$jia_scores, call. = FALSE)
  score_cols <- score_columns()[score_columns() %in% colnames(scores)]
  idx <- match(colnames(prog), scores$cell_id)
  for (cc in score_cols) {
    prog[[cc]] <- as.numeric(scores[[cc]][idx])
  }
  list(object = prog, score_cols = score_cols)
}

run_uncorrected <- function(prog, cfg) {
  DefaultAssay(prog) <- cfg$assay
  set.seed(cfg$seed)
  log_msg("Running uncorrected baseline PCA/UMAP")
  prog <- NormalizeData(prog, assay = cfg$assay, verbose = FALSE)
  prog <- FindVariableFeatures(prog, assay = cfg$assay, selection.method = "vst", nfeatures = cfg$nfeatures, verbose = FALSE)
  prog <- ScaleData(prog, assay = cfg$assay, verbose = FALSE)
  prog <- RunPCA(prog, assay = cfg$assay, npcs = cfg$npcs, reduction.name = "uncorrected_pca", reduction.key = "UNCPCA_", seed.use = cfg$seed, verbose = FALSE)
  prog <- FindNeighbors(
    prog,
    reduction = "uncorrected_pca",
    dims = seq_len(cfg$dims),
    graph.name = c("uncorrected_nn", "uncorrected_snn"),
    verbose = FALSE
  )
  prog <- FindClusters(
    prog,
    graph.name = "uncorrected_snn",
    resolution = cfg$resolution,
    random.seed = cfg$seed,
    cluster.name = cfg$uncorrected_cluster_col,
    verbose = FALSE
  )
  prog <- RunUMAP(
    prog,
    reduction = "uncorrected_pca",
    dims = seq_len(cfg$dims),
    reduction.name = cfg$uncorrected_umap,
    reduction.key = "UNCUMAP_",
    seed.use = cfg$seed,
    verbose = FALSE
  )
  prog
}

run_integrated <- function(prog, cfg, batch_col) {
  batch_values <- as.character(prog@meta.data[[batch_col]])
  if (any(is.na(batch_values) | !nzchar(batch_values))) {
    stop("Batch column contains missing/empty values: ", batch_col, call. = FALSE)
  }
  n_by_batch <- sort(table(batch_values), decreasing = TRUE)
  if (length(n_by_batch) < 2L) stop("Need at least two batches for integration in ", batch_col, call. = FALSE)
  min_batch <- min(as.integer(n_by_batch))
  effective_k_weight <- min(cfg$k_weight, max(1L, min_batch - 1L))
  if (effective_k_weight != cfg$k_weight) {
    log_msg("Reducing k.weight from ", cfg$k_weight, " to ", effective_k_weight, " because the smallest batch has ", min_batch, " cells")
  }

  log_msg("Splitting progenitor object by ", batch_col, " for integration")
  split_objects <- SplitObject(prog, split.by = batch_col)
  split_objects <- split_objects[vapply(split_objects, ncol, numeric(1)) > 0]

  split_objects <- lapply(split_objects, function(x) {
    DefaultAssay(x) <- cfg$assay
    x <- NormalizeData(x, assay = cfg$assay, verbose = FALSE)
    x <- FindVariableFeatures(x, assay = cfg$assay, selection.method = "vst", nfeatures = cfg$integration_nfeatures, verbose = FALSE)
    if (tolower(cfg$integration_reduction) == "rpca") {
      x <- ScaleData(x, assay = cfg$assay, features = VariableFeatures(x), verbose = FALSE)
      x <- RunPCA(x, assay = cfg$assay, features = VariableFeatures(x), npcs = cfg$npcs, seed.use = cfg$seed, verbose = FALSE)
    }
    x
  })

  log_msg("Selecting ", cfg$integration_nfeatures, " integration features")
  features <- SelectIntegrationFeatures(object.list = split_objects, nfeatures = cfg$integration_nfeatures)
  log_msg("Finding integration anchors with reduction=", cfg$integration_reduction, ", dims=1:", cfg$dims, ", k.anchor=", cfg$k_anchor)
  anchors <- FindIntegrationAnchors(
    object.list = split_objects,
    anchor.features = features,
    reduction = cfg$integration_reduction,
    dims = seq_len(cfg$dims),
    k.anchor = cfg$k_anchor,
    verbose = FALSE
  )
  log_msg("Integrating data with k.weight=", effective_k_weight)
  integrated <- IntegrateData(
    anchorset = anchors,
    dims = seq_len(cfg$dims),
    k.weight = effective_k_weight,
    verbose = FALSE
  )

  DefaultAssay(integrated) <- "integrated"
  set.seed(cfg$seed)
  log_msg("Running integrated PCA/UMAP/clustering")
  integrated <- ScaleData(integrated, assay = "integrated", verbose = FALSE)
  integrated <- RunPCA(integrated, assay = "integrated", npcs = cfg$npcs, reduction.name = "integrated_pca", reduction.key = "INTPCA_", seed.use = cfg$seed, verbose = FALSE)
  integrated <- FindNeighbors(
    integrated,
    reduction = "integrated_pca",
    dims = seq_len(cfg$dims),
    graph.name = c("integrated_nn", "integrated_snn"),
    verbose = FALSE
  )
  integrated <- FindClusters(
    integrated,
    graph.name = "integrated_snn",
    resolution = cfg$resolution,
    random.seed = cfg$seed,
    cluster.name = cfg$integrated_cluster_col,
    verbose = FALSE
  )
  integrated <- RunUMAP(
    integrated,
    reduction = "integrated_pca",
    dims = seq_len(cfg$dims),
    reduction.name = cfg$integrated_umap,
    reduction.key = "INTUMAP_",
    seed.use = cfg$seed,
    verbose = FALSE
  )
  list(object = integrated, effective_k_weight = effective_k_weight, batch_counts = n_by_batch)
}

cluster_counts <- function(obj, cluster_col, batch_col) {
  tab <- as.data.frame.matrix(table(obj@meta.data[[cluster_col]], obj@meta.data[[batch_col]]))
  tab$cluster <- rownames(tab)
  tab <- tab[, c("cluster", setdiff(colnames(tab), "cluster")), drop = FALSE]
  rownames(tab) <- NULL
  tab
}

batch_entropy <- function(obj, cluster_col, batch_col, label) {
  tab <- table(obj@meta.data[[cluster_col]], obj@meta.data[[batch_col]])
  rows <- lapply(rownames(tab), function(cluster) {
    counts <- as.numeric(tab[cluster, ])
    props <- counts / sum(counts)
    props_nonzero <- props[props > 0]
    entropy <- -sum(props_nonzero * log(props_nonzero))
    max_entropy <- log(length(counts))
    data.frame(
      analysis = label,
      cluster = cluster,
      n_cells = sum(counts),
      n_batches_present = sum(counts > 0),
      normalized_batch_entropy = ifelse(max_entropy > 0, entropy / max_entropy, NA_real_),
      max_batch_fraction = max(props),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

jia_score_summary <- function(obj, cluster_col, score_cols, label) {
  if (length(score_cols) == 0L) {
    return(data.frame(analysis = label, status = "no_jia_scores_available", stringsAsFactors = FALSE))
  }
  clusters <- as.character(obj@meta.data[[cluster_col]])
  rows <- list()
  for (cluster in sort(unique(clusters))) {
    meta <- obj@meta.data[clusters == cluster, , drop = FALSE]
    for (score_col in score_cols) {
      x <- as.numeric(meta[[score_col]])
      rows[[length(rows) + 1L]] <- data.frame(
        analysis = label,
        cluster = cluster,
        score = score_col,
        n_cells = nrow(meta),
        mean = mean(x, na.rm = TRUE),
        median = median(x, na.rm = TRUE),
        q25 = as.numeric(stats::quantile(x, 0.25, na.rm = TRUE)),
        q75 = as.numeric(stats::quantile(x, 0.75, na.rm = TRUE)),
        pct_non_na = mean(!is.na(x)),
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}

cell_table <- function(uncorrected, integrated, cfg, batch_col) {
  unc <- as.data.frame(Embeddings(uncorrected, cfg$uncorrected_umap))
  int <- as.data.frame(Embeddings(integrated, cfg$integrated_umap))
  unc$cell_id <- rownames(unc)
  int$cell_id <- rownames(int)
  colnames(unc)[seq_len(2)] <- c("uncorrected_umap_1", "uncorrected_umap_2")
  colnames(int)[seq_len(2)] <- c("integrated_umap_1", "integrated_umap_2")
  meta <- integrated@meta.data
  meta$cell_id <- rownames(meta)
  keep <- unique(c(
    "cell_id", batch_col, "div30_parent_cluster", "div30_parent_progenitor_label",
    cfg$uncorrected_cluster_col, cfg$integrated_cluster_col, score_columns()
  ))
  keep <- keep[keep %in% colnames(meta)]
  out <- merge(meta[, keep, drop = FALSE], unc, by = "cell_id", all.x = TRUE)
  out <- merge(out, int, by = "cell_id", all.x = TRUE)
  out
}

dim_plot <- function(obj, reduction, group_by, title, label = FALSE) {
  DimPlot(
    obj,
    reduction = reduction,
    group.by = group_by,
    label = label,
    repel = TRUE,
    raster = TRUE,
    pt.size = 0.25
  ) +
    ggtitle(title) +
    theme_classic(base_size = 12) +
    theme(plot.title = element_text(face = "bold", size = 11))
}

feature_plot <- function(obj, reduction, features, title) {
  p <- FeaturePlot(
    obj,
    reduction = reduction,
    features = features,
    ncol = length(features),
    order = TRUE,
    raster = TRUE
  )
  p + plot_annotation(title = title) &
    theme(plot.title = element_text(face = "bold", size = 11))
}

plot_outputs <- function(uncorrected, integrated, cfg, batch_col, score_cols) {
  p_unc_batch <- dim_plot(uncorrected, cfg$uncorrected_umap, batch_col, "Before integration: sample/batch")
  p_int_batch <- dim_plot(integrated, cfg$integrated_umap, batch_col, "After strong integration: sample/batch")
  save_plot_pair(p_unc_batch, file.path(cfg$plot_dir, "uncorrected_umap_by_batch"), 7, 5.8)
  save_plot_pair(p_int_batch, file.path(cfg$plot_dir, "integrated_umap_by_batch"), 7, 5.8)
  save_plot_pair(p_unc_batch + p_int_batch + plot_layout(ncol = 2), file.path(cfg$plot_dir, "uncorrected_vs_integrated_umap_by_batch"), 14, 5.8)

  p_unc_parent <- dim_plot(uncorrected, cfg$uncorrected_umap, "div30_parent_cluster", "Before integration: source cluster", TRUE)
  p_int_parent <- dim_plot(integrated, cfg$integrated_umap, "div30_parent_cluster", "After integration: source cluster", TRUE)
  save_plot_pair(p_unc_parent + p_int_parent + plot_layout(ncol = 2), file.path(cfg$plot_dir, "uncorrected_vs_integrated_umap_by_source_cluster"), 14, 5.8)

  p_unc_cluster <- dim_plot(uncorrected, cfg$uncorrected_umap, cfg$uncorrected_cluster_col, "Before integration: reclusters", TRUE)
  p_int_cluster <- dim_plot(integrated, cfg$integrated_umap, cfg$integrated_cluster_col, "After integration: reclusters", TRUE)
  save_plot_pair(p_unc_cluster, file.path(cfg$plot_dir, "uncorrected_umap_by_recluster"), 7, 5.8)
  save_plot_pair(p_int_cluster, file.path(cfg$plot_dir, "integrated_umap_by_recluster"), 7, 5.8)
  save_plot_pair(p_unc_cluster + p_int_cluster + plot_layout(ncol = 2), file.path(cfg$plot_dir, "uncorrected_vs_integrated_umap_by_recluster"), 14, 5.8)

  for (metric in c("nFeature_RNA", "nCount_RNA", "percent.mt", "S.Score", "G2M.Score", "CC.Difference")) {
    if (metric %in% colnames(integrated@meta.data)) {
      p_unc_metric <- FeaturePlot(uncorrected, reduction = cfg$uncorrected_umap, features = metric, order = TRUE, raster = TRUE) +
        ggtitle(paste("Before integration:", metric)) +
        theme(plot.title = element_text(face = "bold", size = 11))
      p_int_metric <- FeaturePlot(integrated, reduction = cfg$integrated_umap, features = metric, order = TRUE, raster = TRUE) +
        ggtitle(paste("After integration:", metric)) +
        theme(plot.title = element_text(face = "bold", size = 11))
      save_plot_pair(p_unc_metric + p_int_metric + plot_layout(ncol = 2), file.path(cfg$plot_dir, paste0("uncorrected_vs_integrated_umap_", safe_name(metric))), 14, 5.8)
    }
  }

  if (length(score_cols) > 0L) {
    save_plot_pair(
      feature_plot(uncorrected, cfg$uncorrected_umap, score_cols, "Before integration: Jia scores"),
      file.path(cfg$plot_dir, "uncorrected_jia_score_umap_overlays"),
      11, 4.2
    )
    save_plot_pair(
      feature_plot(integrated, cfg$integrated_umap, score_cols, "After integration: Jia scores"),
      file.path(cfg$plot_dir, "integrated_jia_score_umap_overlays"),
      11, 4.2
    )
  }
}

opt <- parse_args(commandArgs(trailingOnly = TRUE))
if (isTRUE(opt$help)) {
  print_usage()
  quit(save = "no", status = 0)
}

cfg <- build_config(opt)
dir.create(cfg$table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(cfg$plot_dir, recursive = TRUE, showWarnings = FALSE)
options(future.globals.maxSize = max(getOption("future.globals.maxSize", 0), cfg$future_globals_max_gb * 1024^3))
if (requireNamespace("future", quietly = TRUE)) future::plan("sequential")

log_msg("Loading Seurat object: ", cfg$seurat_rds)
obj <- readRDS(cfg$seurat_rds)
if (!inherits(obj, "Seurat")) stop("Input is not a Seurat object: ", cfg$seurat_rds, call. = FALSE)
if (!(cfg$assay %in% Assays(obj))) stop("Assay not found: ", cfg$assay, call. = FALSE)

source_cluster_col <- find_cluster_col(obj, cfg$source_cluster_col_requested)
prog <- subset_progenitors(obj, cfg, source_cluster_col)
rm(obj)
invisible(gc())

batch_col <- find_batch_col(prog, cfg$batch_col_requested)
prog[[batch_col]] <- as.character(prog@meta.data[[batch_col]])
write_tsv(
  data.frame(batch = names(table(prog@meta.data[[batch_col]])), n_cells = as.integer(table(prog@meta.data[[batch_col]])), stringsAsFactors = FALSE),
  file.path(cfg$table_dir, "div30_progenitor_cells_by_batch.tsv")
)

jia <- attach_jia_scores(prog, cfg)
prog <- jia$object
score_cols <- jia$score_cols

uncorrected <- run_uncorrected(prog, cfg)
integration <- run_integrated(prog, cfg, batch_col)
integrated <- integration$object

for (cc in c(cfg$uncorrected_cluster_col, score_cols)) {
  if (cc %in% colnames(uncorrected@meta.data)) integrated[[cc]] <- uncorrected@meta.data[colnames(integrated), cc]
}

write_tsv(
  parameter_table(cfg, source_cluster_col, batch_col, n_parent_cells = NA_integer_, n_prog_cells = ncol(prog), effective_k_weight = integration$effective_k_weight),
  file.path(cfg$table_dir, "div30_jia_progenitor_strong_integration_run_parameters.tsv")
)
write_tsv(cluster_counts(uncorrected, cfg$uncorrected_cluster_col, batch_col), file.path(cfg$table_dir, "uncorrected_cluster_counts_by_batch.tsv"))
write_tsv(cluster_counts(integrated, cfg$integrated_cluster_col, batch_col), file.path(cfg$table_dir, "integrated_cluster_counts_by_batch.tsv"))
write_tsv(
  rbind(
    batch_entropy(uncorrected, cfg$uncorrected_cluster_col, batch_col, "uncorrected"),
    batch_entropy(integrated, cfg$integrated_cluster_col, batch_col, "integrated")
  ),
  file.path(cfg$table_dir, "batch_entropy_by_cluster_before_after.tsv")
)
write_tsv(
  rbind(
    jia_score_summary(uncorrected, cfg$uncorrected_cluster_col, score_cols, "uncorrected"),
    jia_score_summary(integrated, cfg$integrated_cluster_col, score_cols, "integrated")
  ),
  file.path(cfg$table_dir, "jia_score_summary_by_cluster_before_after.tsv")
)
write_tsv(cell_table(uncorrected, integrated, cfg, batch_col), file.path(cfg$table_dir, "div30_jia_progenitor_strong_integration_cells.tsv.gz"))

plot_outputs(uncorrected, integrated, cfg, batch_col, score_cols)

saveRDS(uncorrected, file.path(cfg$outdir, "div30_jia_progenitor_uncorrected_baseline_seurat.rds"))
saveRDS(integrated, file.path(cfg$outdir, "div30_jia_progenitor_strong_integrated_seurat.rds"))
write_tsv(
  data.frame(
    status = "complete",
    completed_at = timestamp(),
    batch_col = batch_col,
    n_cells = ncol(integrated),
    uncorrected_object = file.path(cfg$outdir, "div30_jia_progenitor_uncorrected_baseline_seurat.rds"),
    integrated_object = file.path(cfg$outdir, "div30_jia_progenitor_strong_integrated_seurat.rds"),
    before_after_sample_umap = file.path(cfg$plot_dir, "uncorrected_vs_integrated_umap_by_batch.png"),
    stringsAsFactors = FALSE
  ),
  file.path(cfg$table_dir, "div30_jia_progenitor_strong_integration_complete.tsv")
)

log_msg("Complete. Outputs: ", cfg$outdir)
