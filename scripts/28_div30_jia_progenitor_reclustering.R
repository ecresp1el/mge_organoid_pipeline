#!/usr/bin/env Rscript

# DIV30 Jia-style Phase 1 progenitor reclustering.
#
# Scope:
#   - Use only the current DIV30 progenitor Seurat clusters 0, 3, 6, and 7.
#   - Recluster the progenitor-only object with the existing DIV30 Seurat-style
#     workflow: NormalizeData, FindVariableFeatures, ScaleData, RunPCA,
#     FindNeighbors, FindClusters, RunUMAP.
#   - Run FindAllMarkers on the new unsupervised clusters.
#   - Tentatively annotate clusters as VZ-RGC/RGC1, SVZ-RGC/RGC2, or IPC from
#     marker genes.
#   - Attach Jia scores only after clustering/marker identification for
#     validation overlays. Jia scores are not used for cell selection,
#     reclustering, PCA, neighbors, UMAP, or marker calls.
#
# Stop point:
#   Phase 1 only. This script intentionally does not perform lineage-committed
#   progenitor reclustering.

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(ggplot2)
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
    `run-label` = Sys.getenv("RUN_LABEL", "div30_jia_progenitor_reclustering_phase1_v1"),
    assay = Sys.getenv("SEURAT_ASSAY", "RNA"),
    `source-cluster-col` = Sys.getenv("SOURCE_CLUSTER_COL", ""),
    `progenitor-clusters` = Sys.getenv("PROGENITOR_CLUSTERS", "0,3,6,7"),
    `nfeatures` = Sys.getenv("NFEATURES", "3000"),
    dims = Sys.getenv("DIMS", "30"),
    npcs = Sys.getenv("NPCS", "50"),
    resolution = Sys.getenv("RESOLUTION", "0.8"),
    seed = Sys.getenv("SEED", "7"),
    `cluster-col` = Sys.getenv("PROGENITOR_CLUSTER_COL", "div30_progenitor_cluster"),
    `umap-reduction` = Sys.getenv("PROGENITOR_UMAP_REDUCTION", "progenitor_umap"),
    `min-pct` = Sys.getenv("MARKER_MIN_PCT", "0.1"),
    `logfc-threshold` = Sys.getenv("MARKER_LOGFC_THRESHOLD", "0.25"),
    `top-n-markers` = Sys.getenv("TOP_N_MARKERS", "50"),
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
    "  Rscript scripts/28_div30_jia_progenitor_reclustering.R --project-root <PROJECT_ROOT> [options]",
    "",
    "Default inputs:",
    "  --seurat-rds PROJECT_ROOT/results/varela_this_paper/varela_this_paper_seurat.rds",
    "  --jia-scores PROJECT_ROOT/results/jia_program_div30_scoring/jia_program_div30_scoring_v1/tables/div30_jia_program_scores_obs.tsv",
    "",
    "Core options:",
    "  --progenitor-clusters 0,3,6,7",
    "  --source-cluster-col <col>       Default: seurat_clusters, falling back to RNA_snn_res.0.2",
    "  --nfeatures 3000",
    "  --dims 30",
    "  --resolution 0.8",
    "",
    "Outputs:",
    "  tables/div30_jia_progenitor_phase1_run_parameters.tsv",
    "  tables/div30_jia_progenitor_phase1_cluster_counts.tsv",
    "  tables/div30_jia_progenitor_phase1_all_markers.tsv.gz",
    "  tables/div30_jia_progenitor_phase1_top_markers.tsv",
    "  tables/div30_jia_progenitor_phase1_candidate_marker_expression_by_cluster.tsv",
    "  tables/div30_jia_progenitor_phase1_tentative_annotations.tsv",
    "  tables/div30_jia_progenitor_phase1_jia_score_summary_by_cluster.tsv",
    "  tables/div30_jia_progenitor_phase1_cells.tsv.gz",
    "  plots/*.png and *.pdf",
    "  div30_jia_progenitor_phase1_seurat.rds",
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

get_assay_matrix <- function(obj, assay, layer) {
  tryCatch(
    SeuratObject::GetAssayData(obj, assay = assay, layer = layer),
    error = function(e) SeuratObject::GetAssayData(obj, assay = assay, slot = layer)
  )
}

find_genes_case_insensitive <- function(requested, available) {
  idx <- match(toupper(requested), toupper(available))
  found <- available[idx[!is.na(idx)]]
  names(found) <- requested[!is.na(idx)]
  found
}

candidate_marker_sets <- function() {
  list(
    VZ_RGC_RGC1 = c("HES1", "VIM", "NES"),
    SVZ_RGC_RGC2 = c("FBLN7", "CACNA1E", "DACH1"),
    IPC = c("DLX1", "DLX2", "ASCL1")
  )
}

marker_spec_table <- function() {
  sets <- candidate_marker_sets()
  do.call(rbind, lapply(names(sets), function(set_name) {
    data.frame(marker_set = set_name, gene = sets[[set_name]], stringsAsFactors = FALSE)
  }))
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
  npcs <- as_int(opt$npcs, "npcs")
  if (dims < 2L) stop("dims must be >= 2", call. = FALSE)
  if (npcs < dims) npcs <- dims

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
    nfeatures = as_int(opt$nfeatures, "nfeatures"),
    dims = dims,
    npcs = npcs,
    resolution = as_num(opt$resolution, "resolution"),
    seed = as_int(opt$seed, "seed"),
    cluster_col = opt$`cluster-col`,
    umap_reduction = opt$`umap-reduction`,
    min_pct = as_num(opt$`min-pct`, "min-pct"),
    logfc_threshold = as_num(opt$`logfc-threshold`, "logfc-threshold"),
    top_n_markers = as_int(opt$`top-n-markers`, "top-n-markers")
  )
}

parameter_table <- function(cfg, source_cluster_col, n_parent_cells = NA_integer_, n_prog_cells = NA_integer_) {
  data.frame(
    parameter = c(
      "project_root",
      "seurat_rds",
      "jia_scores",
      "outdir",
      "run_label",
      "assay",
      "source_cluster_col",
      "progenitor_clusters",
      "n_parent_cells",
      "n_progenitor_cells",
      "normalization",
      "variable_feature_method",
      "nfeatures",
      "scale_data",
      "npcs",
      "dims",
      "find_neighbors_graph",
      "resolution",
      "cluster_col",
      "umap_reduction",
      "seed",
      "marker_test",
      "marker_min_pct",
      "marker_logfc_threshold",
      "top_n_markers",
      "phase_stop"
    ),
    value = as.character(c(
      cfg$project_root,
      cfg$seurat_rds,
      cfg$jia_scores,
      cfg$outdir,
      cfg$run_label,
      cfg$assay,
      source_cluster_col,
      paste(cfg$progenitor_clusters, collapse = ","),
      n_parent_cells,
      n_prog_cells,
      "NormalizeData default LogNormalize",
      "FindVariableFeatures selection.method=vst",
      cfg$nfeatures,
      "ScaleData default, no regression",
      cfg$npcs,
      paste0("1:", cfg$dims),
      "FindNeighbors on PCA",
      cfg$resolution,
      cfg$cluster_col,
      cfg$umap_reduction,
      cfg$seed,
      "FindAllMarkers after unsupervised reclustering",
      cfg$min_pct,
      cfg$logfc_threshold,
      cfg$top_n_markers,
      "Phase 1 only: no lineage-committed progenitor reclustering"
    )),
    stringsAsFactors = FALSE
  )
}

recluster_progenitors <- function(obj, cfg, source_cluster_col) {
  source_clusters <- as.character(obj@meta.data[[source_cluster_col]])
  names(source_clusters) <- colnames(obj)
  keep_cells <- names(source_clusters)[source_clusters %in% cfg$progenitor_clusters]
  if (length(keep_cells) == 0) {
    stop("No cells found for progenitor clusters: ", paste(cfg$progenitor_clusters, collapse = ","), call. = FALSE)
  }

  log_msg("Subsetting progenitor clusters ", paste(cfg$progenitor_clusters, collapse = ","), ": ", length(keep_cells), " cells")
  prog <- subset(obj, cells = keep_cells)
  prog$div30_parent_cluster <- as.character(prog@meta.data[[source_cluster_col]])
  prog$div30_parent_progenitor_label <- ifelse(
    prog$div30_parent_cluster == "6",
    "Inhibitory progenitors",
    "Radial glia"
  )
  DefaultAssay(prog) <- cfg$assay

  set.seed(cfg$seed)
  log_msg("Running NormalizeData")
  prog <- NormalizeData(prog, assay = cfg$assay, verbose = FALSE)
  log_msg("Running FindVariableFeatures")
  prog <- FindVariableFeatures(prog, assay = cfg$assay, selection.method = "vst", nfeatures = cfg$nfeatures, verbose = FALSE)
  log_msg("Running ScaleData")
  prog <- ScaleData(prog, assay = cfg$assay, verbose = FALSE)
  log_msg("Running RunPCA")
  prog <- RunPCA(prog, assay = cfg$assay, npcs = cfg$npcs, seed.use = cfg$seed, verbose = FALSE)
  log_msg("Running FindNeighbors")
  prog <- FindNeighbors(
    prog,
    reduction = "pca",
    dims = seq_len(cfg$dims),
    graph.name = c("progenitor_nn", "progenitor_snn"),
    verbose = FALSE
  )
  log_msg("Running FindClusters")
  prog <- FindClusters(
    prog,
    graph.name = "progenitor_snn",
    resolution = cfg$resolution,
    random.seed = cfg$seed,
    cluster.name = cfg$cluster_col,
    verbose = FALSE
  )
  Idents(prog) <- prog@meta.data[[cfg$cluster_col]]
  log_msg("Running RunUMAP")
  prog <- RunUMAP(
    prog,
    reduction = "pca",
    dims = seq_len(cfg$dims),
    reduction.name = cfg$umap_reduction,
    reduction.key = "PROGUMAP_",
    seed.use = cfg$seed,
    verbose = FALSE
  )
  prog
}

run_markers <- function(prog, cfg) {
  Idents(prog) <- prog@meta.data[[cfg$cluster_col]]
  markers <- FindAllMarkers(
    prog,
    assay = cfg$assay,
    only.pos = TRUE,
    min.pct = cfg$min_pct,
    logfc.threshold = cfg$logfc_threshold,
    verbose = FALSE
  )
  if (nrow(markers) == 0) return(markers)
  if ("avg_logFC" %in% colnames(markers) && !("avg_log2FC" %in% colnames(markers))) {
    markers$avg_log2FC <- markers$avg_logFC
  }
  markers <- markers[order(markers$cluster, -markers$avg_log2FC, markers$p_val_adj), , drop = FALSE]
  rownames(markers) <- NULL
  markers
}

top_markers <- function(markers, top_n) {
  if (nrow(markers) == 0) return(markers)
  pieces <- split(markers, markers$cluster)
  out <- do.call(rbind, lapply(pieces, function(df) head(df, top_n)))
  rownames(out) <- NULL
  out
}

candidate_marker_summary <- function(prog, cfg) {
  spec <- marker_spec_table()
  found <- find_genes_case_insensitive(spec$gene, rownames(prog))
  spec$gene_found <- unname(found[spec$gene])
  spec$gene_present <- !is.na(spec$gene_found)

  clusters <- as.character(prog@meta.data[[cfg$cluster_col]])
  names(clusters) <- colnames(prog)
  expr <- get_assay_matrix(prog, cfg$assay, "data")

  rows <- list()
  for (i in seq_len(nrow(spec))) {
    if (!isTRUE(spec$gene_present[[i]])) {
      next
    }
    gene <- spec$gene_found[[i]]
    values <- as.numeric(expr[gene, colnames(prog)])
    names(values) <- colnames(prog)
    for (cluster in sort(unique(clusters))) {
      cells <- names(clusters)[clusters == cluster]
      x <- values[cells]
      rows[[length(rows) + 1L]] <- data.frame(
        cluster = cluster,
        marker_set = spec$marker_set[[i]],
        gene = spec$gene[[i]],
        gene_found = gene,
        n_cells = length(cells),
        mean_expression = mean(x, na.rm = TRUE),
        median_expression = median(x, na.rm = TRUE),
        pct_expressed = mean(x > 0, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    }
  }

  summary <- if (length(rows) > 0) do.call(rbind, rows) else {
    data.frame(
      cluster = character(), marker_set = character(), gene = character(), gene_found = character(),
      n_cells = integer(), mean_expression = numeric(), median_expression = numeric(), pct_expressed = numeric(),
      stringsAsFactors = FALSE
    )
  }
  list(spec = spec, summary = summary)
}

set_scores_by_cluster <- function(marker_summary) {
  if (nrow(marker_summary) == 0) {
    return(data.frame(cluster = character(), marker_set = character(), marker_set_score = numeric(), stringsAsFactors = FALSE))
  }
  aggregate(
    marker_summary$mean_expression,
    by = list(cluster = marker_summary$cluster, marker_set = marker_summary$marker_set),
    FUN = mean,
    na.rm = TRUE
  ) |>
    setNames(c("cluster", "marker_set", "marker_set_score"))
}

tentative_annotations <- function(marker_summary, cluster_counts) {
  scores <- set_scores_by_cluster(marker_summary)
  if (nrow(scores) == 0) {
    return(data.frame(
      cluster = cluster_counts$cluster,
      n_cells = cluster_counts$n_cells,
      tentative_identity = "unassigned_no_candidate_markers_found",
      marker_support = "",
      stringsAsFactors = FALSE
    ))
  }
  labels <- c(
    VZ_RGC_RGC1 = "RGC1 / VZ-RGC candidate",
    SVZ_RGC_RGC2 = "RGC2 / SVZ-RGC candidate",
    IPC = "IPC candidate"
  )
  rows <- lapply(split(scores, scores$cluster), function(df) {
    df <- df[order(-df$marker_set_score), , drop = FALSE]
    top <- df[1, , drop = FALSE]
    data.frame(
      cluster = as.character(top$cluster),
      tentative_identity = labels[[top$marker_set]],
      top_marker_set = top$marker_set,
      top_marker_set_score = top$marker_set_score,
      marker_support = paste(sprintf("%s=%.4f", df$marker_set, df$marker_set_score), collapse = "; "),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out <- merge(cluster_counts, out, by = "cluster", all.x = TRUE)
  out <- out[order(as.numeric(as.character(out$cluster))), , drop = FALSE]
  rownames(out) <- NULL
  out
}

attach_jia_scores <- function(prog, cfg) {
  if (!file.exists(cfg$jia_scores)) {
    log_msg("Jia score table not found; skipping validation overlays: ", cfg$jia_scores)
    return(list(object = prog, score_cols = character(), scores = data.frame()))
  }
  log_msg("Reading Jia scores for post hoc validation: ", cfg$jia_scores)
  scores <- read_tsv(cfg$jia_scores)
  if (!("cell_id" %in% colnames(scores))) {
    stop("Jia score table must contain cell_id column: ", cfg$jia_scores, call. = FALSE)
  }
  cols <- score_columns()
  score_cols <- cols[cols %in% colnames(scores)]
  missing <- setdiff(cols, score_cols)
  if (length(missing) > 0) {
    log_msg("Missing Jia score columns: ", paste(missing, collapse = ", "))
  }
  if (length(score_cols) == 0) {
    log_msg("No requested Jia score columns found; skipping validation overlays.")
    return(list(object = prog, score_cols = character(), scores = scores))
  }
  idx <- match(colnames(prog), scores$cell_id)
  for (cc in score_cols) {
    prog[[cc]] <- as.numeric(scores[[cc]][idx])
  }
  list(object = prog, score_cols = score_cols, scores = scores)
}

jia_score_summary <- function(prog, cfg, score_cols) {
  if (length(score_cols) == 0) {
    return(data.frame(status = "no_jia_scores_available", stringsAsFactors = FALSE))
  }
  clusters <- as.character(prog@meta.data[[cfg$cluster_col]])
  rows <- list()
  for (cluster in sort(unique(clusters))) {
    meta <- prog@meta.data[clusters == cluster, , drop = FALSE]
    for (score_col in score_cols) {
      x <- as.numeric(meta[[score_col]])
      rows[[length(rows) + 1L]] <- data.frame(
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

cell_table <- function(prog, cfg) {
  coords <- as.data.frame(Embeddings(prog, cfg$umap_reduction))
  coords$cell_id <- rownames(coords)
  colnames(coords)[seq_len(2)] <- c("progenitor_umap_1", "progenitor_umap_2")
  meta <- prog@meta.data
  meta$cell_id <- rownames(meta)
  keep <- unique(c(
    "cell_id",
    "div30_parent_cluster",
    "div30_parent_progenitor_label",
    cfg$cluster_col,
    score_columns()
  ))
  keep <- keep[keep %in% colnames(meta)]
  out <- merge(meta[, keep, drop = FALSE], coords, by = "cell_id", all.x = TRUE)
  out
}

plot_outputs <- function(prog, cfg, marker_spec, score_cols) {
  cluster_plot <- DimPlot(
    prog,
    reduction = cfg$umap_reduction,
    group.by = cfg$cluster_col,
    label = TRUE,
    repel = TRUE,
    raster = TRUE
  ) +
    ggtitle("DIV30 progenitor reclusters") +
    theme(plot.title = element_text(size = 11))
  save_plot_pair(cluster_plot, file.path(cfg$plot_dir, "div30_progenitor_umap_by_recluster"), 7, 6)

  parent_plot <- DimPlot(
    prog,
    reduction = cfg$umap_reduction,
    group.by = "div30_parent_cluster",
    label = TRUE,
    repel = TRUE,
    raster = TRUE
  ) +
    ggtitle("DIV30 progenitor UMAP by source cluster") +
    theme(plot.title = element_text(size = 11))
  save_plot_pair(parent_plot, file.path(cfg$plot_dir, "div30_progenitor_umap_by_source_cluster"), 7, 6)

  genes <- unique(stats::na.omit(marker_spec$gene_found))
  if (length(genes) > 0) {
    marker_plot <- FeaturePlot(
      prog,
      reduction = cfg$umap_reduction,
      features = genes,
      ncol = 3,
      order = TRUE,
      raster = TRUE
    )
    save_plot_pair(marker_plot, file.path(cfg$plot_dir, "div30_progenitor_candidate_marker_feature_grid"), 10, 9)

    dot_plot <- DotPlot(
      prog,
      assay = cfg$assay,
      features = genes,
      group.by = cfg$cluster_col
    ) +
      RotatedAxis() +
      ggtitle("Candidate VZ-RGC, SVZ-RGC, and IPC markers by recluster") +
      theme(plot.title = element_text(size = 11))
    save_plot_pair(dot_plot, file.path(cfg$plot_dir, "div30_progenitor_candidate_marker_dotplot"), 9, 4.8)
  }

  if (length(score_cols) > 0) {
    score_plot <- FeaturePlot(
      prog,
      reduction = cfg$umap_reduction,
      features = score_cols,
      ncol = length(score_cols),
      order = TRUE,
      raster = TRUE
    )
    save_plot_pair(score_plot, file.path(cfg$plot_dir, "div30_progenitor_jia_score_umap_overlays"), 11, 4)

    violin_plot <- VlnPlot(
      prog,
      features = score_cols,
      group.by = cfg$cluster_col,
      pt.size = 0,
      ncol = length(score_cols)
    ) +
      ggtitle("Post hoc Jia score validation by recluster") +
      theme(plot.title = element_text(size = 11))
    save_plot_pair(violin_plot, file.path(cfg$plot_dir, "div30_progenitor_jia_score_violins_by_recluster"), 11, 4)
  }
}

run_pipeline <- function(cfg) {
  dir.create(cfg$table_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(cfg$plot_dir, recursive = TRUE, showWarnings = FALSE)

  if (!file.exists(cfg$seurat_rds)) stop("Seurat RDS not found: ", cfg$seurat_rds, call. = FALSE)
  log_msg("Reading DIV30 Seurat object: ", cfg$seurat_rds)
  obj <- readRDS(cfg$seurat_rds)
  if (!inherits(obj, "Seurat")) stop("Input is not a Seurat object: ", cfg$seurat_rds, call. = FALSE)
  if (!(cfg$assay %in% Assays(obj))) stop("Assay not found: ", cfg$assay, ". Available: ", paste(Assays(obj), collapse = ","), call. = FALSE)
  source_cluster_col <- find_cluster_col(obj, cfg$source_cluster_col_requested)

  write_tsv(
    parameter_table(cfg, source_cluster_col, n_parent_cells = ncol(obj), n_prog_cells = NA_integer_),
    file.path(cfg$table_dir, "div30_jia_progenitor_phase1_run_parameters.tsv")
  )

  prog <- recluster_progenitors(obj, cfg, source_cluster_col)
  rm(obj)
  gc()

  cluster_counts <- as.data.frame(sort(table(as.character(prog@meta.data[[cfg$cluster_col]])), decreasing = TRUE), stringsAsFactors = FALSE)
  colnames(cluster_counts) <- c("cluster", "n_cells")
  cluster_counts <- cluster_counts[order(as.numeric(as.character(cluster_counts$cluster))), , drop = FALSE]
  write_tsv(cluster_counts, file.path(cfg$table_dir, "div30_jia_progenitor_phase1_cluster_counts.tsv"))

  write_tsv(
    parameter_table(cfg, source_cluster_col, n_parent_cells = NA_integer_, n_prog_cells = ncol(prog)),
    file.path(cfg$table_dir, "div30_jia_progenitor_phase1_run_parameters.tsv")
  )

  log_msg("Running FindAllMarkers on new progenitor clusters")
  markers <- run_markers(prog, cfg)
  write_tsv(markers, file.path(cfg$table_dir, "div30_jia_progenitor_phase1_all_markers.tsv.gz"))
  marker_cols <- intersect(c("cluster", "gene", "avg_log2FC", "pct.1", "pct.2", "p_val_adj"), colnames(markers))
  write_tsv(top_markers(markers, cfg$top_n_markers)[, marker_cols, drop = FALSE], file.path(cfg$table_dir, "div30_jia_progenitor_phase1_top_markers.tsv"))

  marker_info <- candidate_marker_summary(prog, cfg)
  write_tsv(marker_info$spec, file.path(cfg$table_dir, "div30_jia_progenitor_phase1_candidate_marker_panel.tsv"))
  write_tsv(marker_info$summary, file.path(cfg$table_dir, "div30_jia_progenitor_phase1_candidate_marker_expression_by_cluster.tsv"))

  annotations <- tentative_annotations(marker_info$summary, cluster_counts)
  write_tsv(annotations, file.path(cfg$table_dir, "div30_jia_progenitor_phase1_tentative_annotations.tsv"))
  annotation_map <- annotations$tentative_identity
  names(annotation_map) <- annotations$cluster
  prog$div30_progenitor_phase1_tentative_identity <- unname(annotation_map[as.character(prog@meta.data[[cfg$cluster_col]])])

  jia <- attach_jia_scores(prog, cfg)
  prog <- jia$object
  score_summary <- jia_score_summary(prog, cfg, jia$score_cols)
  write_tsv(score_summary, file.path(cfg$table_dir, "div30_jia_progenitor_phase1_jia_score_summary_by_cluster.tsv"))

  write_tsv(cell_table(prog, cfg), file.path(cfg$table_dir, "div30_jia_progenitor_phase1_cells.tsv.gz"))
  plot_outputs(prog, cfg, marker_info$spec, jia$score_cols)

  saveRDS(prog, file.path(cfg$outdir, "div30_jia_progenitor_phase1_seurat.rds"))
  write_tsv(
    data.frame(
      run_label = cfg$run_label,
      completed_at = timestamp(),
      outdir = normalizePath(cfg$outdir, mustWork = FALSE),
      n_cells = ncol(prog),
      n_clusters = length(unique(as.character(prog@meta.data[[cfg$cluster_col]]))),
      phase_stop = "Phase 1 complete; do not proceed to lineage-committed progenitor reclustering until reviewed.",
      stringsAsFactors = FALSE
    ),
    file.path(cfg$table_dir, "div30_jia_progenitor_phase1_complete.tsv")
  )
  log_msg("DIV30 Jia-style progenitor Phase 1 complete: ", cfg$outdir)
}

opt <- parse_args(commandArgs(trailingOnly = TRUE))
if (isTRUE(opt$help)) {
  print_usage()
  quit(save = "no", status = 0)
}
cfg <- build_config(opt)
run_pipeline(cfg)
