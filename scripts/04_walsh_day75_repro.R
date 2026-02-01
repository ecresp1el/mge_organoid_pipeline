#!/usr/bin/env Rscript

# Walsh day-75 reproduction (GSM7979671 dFB, GSM7979672 vFB) with explicit, paper-matched parameters.
# Enforces R 4.1.x and Seurat 4.2.0; fails otherwise.
# Deterministic: explicit seeds for HVG/PCA/Neighbors/Clustering/UMAP.
# Outputs: checkpoints, stress tables, plots, and a one-page reproduction report under PROJECT_ROOT/results/walsh_day75_repro.

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(data.table)
  library(ggplot2)
})

require_version <- function(pkg, expected) {
  v <- packageVersion(pkg)
  if (v != expected) {
    stop(sprintf("Version check failed for %s: found %s, expected %s", pkg, v, expected))
  }
}

check_R_version <- function() {
  v <- getRversion()
  if (!(v$major == "4" && v$minor == "1.0")) {
    stop(sprintf("R version must be 4.1.x; found %s", v))
  }
}

log_msg <- function(...) {
  msg <- sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste0(..., collapse = " "))
  cat(msg, "\n")
  flush.console()
}

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  res <- list(
    project_root = Sys.getenv("PROJECT_ROOT", ""),
    hypoxia_genes = Sys.getenv("HYPOXIA_GENES", ""),
    glycolysis_genes = Sys.getenv("GLYCOLYSIS_GENES", "")
  )
  i <- 1
  while (i <= length(args)) {
    if (args[i] %in% c("-p", "--project-root")) { res$project_root <- args[i + 1]; i <- i + 2; next }
    if (args[i] == "--hypoxia-genes") { res$hypoxia_genes <- args[i + 1]; i <- i + 2; next }
    if (args[i] == "--glycolysis-genes") { res$glycolysis_genes <- args[i + 1]; i <- i + 2; next }
    if (args[i] %in% c("-h", "--help")) {
      cat("Usage: 04_walsh_day75_repro.R --project-root <path> --hypoxia-genes <file> --glycolysis-genes <file>\n")
      quit(status = 0)
    }
    i <- i + 1
  }
  if (res$project_root == "") stop("PROJECT_ROOT is required via env or --project-root")
  if (res$hypoxia_genes == "" || res$glycolysis_genes == "") stop("Both hypoxia and glycolysis gene files are required")
  res
}

read_gene_list <- function(path) {
  if (!file.exists(path)) stop(sprintf("Gene list file not found: %s", path))
  genes <- scan(path, what = character(), quiet = TRUE)
  # Normalize PKM variants
  genes <- gsub("PKM[`’']+", "PKM", genes, perl = TRUE)
  unique(genes[genes != ""])
}

read_counts_csv <- function(path) {
  if (grepl("\\.gz$", path, ignore.case = TRUE)) {
    dt <- fread(cmd = paste("gunzip -c", shQuote(path)))
  } else {
    dt <- fread(path)
  }
  gene_col <- dt[[1]]
  mat <- as.matrix(dt[, -1, with = FALSE])
  rownames(mat) <- gene_col
  m <- Matrix(mat, sparse = TRUE)
  # Validate integer counts
  if (!all(m@x == round(m@x))) stop("Counts are not integer-valued")
  m
}

load_day75_counts <- function(project_root) {
  tarfile <- file.path(project_root, "data", "raw", "walsh_2025_geo_files", "suppl", "GSE250482_RAW.tar")
  if (!file.exists(tarfile)) stop("Missing tar: ", tarfile)
  outdir <- file.path(project_root, "data", "processed", "walsh_day75_repro", "raw_extract")
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  files_to_extract <- c("GSM7979671_dFB_d75.csv.gz", "GSM7979672_vFB_d75.csv.gz")
  utils::untar(tarfile, files = files_to_extract, exdir = outdir)
  list(
    dFB = file.path(outdir, "GSM7979671_dFB_d75.csv.gz"),
    vFB = file.path(outdir, "GSM7979672_vFB_d75.csv.gz")
  )
}

build_seurat <- function(counts, sample_id, domain_label) {
  obj <- CreateSeuratObject(counts = counts, project = sample_id, min.cells = 0, min.features = 0)
  obj$sample_id <- sample_id
  obj$domain <- domain_label
  obj
}

qc_filter <- function(obj) {
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  subset(obj, subset = nFeature_RNA >= 1000 & nFeature_RNA <= 5000 & percent.mt < 15)
}

log_norm_settings <- function() {
  log_msg("Normalization: NormalizeData(LogNormalize, scale.factor=1e4), HVG: vst, nfeatures=5000")
}

run_pipeline <- function(obj, seeds, hypoxia_genes, glycolysis_genes, out_prefix = NULL) {
  set.seed(seeds$hvg)
  obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 1e4, verbose = FALSE)
  obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 5000, verbose = FALSE)

  set.seed(seeds$cc)
  obj <- CellCycleScoring(obj, s.features = cc.genes.updated.2019$s.genes, g2m.features = cc.genes.updated.2019$g2m.genes, set.ident = FALSE)
  log_msg("ScaleData regressors: S.Score, G2M.Score")
  obj <- ScaleData(obj, features = rownames(obj), vars.to.regress = c("S.Score", "G2M.Score"), verbose = FALSE)

  set.seed(seeds$pca)
  obj <- RunPCA(obj, features = VariableFeatures(obj), verbose = FALSE)

  set.seed(seeds$neighbors)
  obj <- FindNeighbors(obj, dims = 1:20, k.param = 20, algorithm = 1, verbose = FALSE)

  set.seed(seeds$clusters)
  obj <- FindClusters(obj, resolution = 2.0, algorithm = 1, verbose = FALSE)

  set.seed(seeds$umap)
  obj <- RunUMAP(
    obj, dims = 1:20, seed.use = seeds$umap,
    n.neighbors = 30, min.dist = 0.3, spread = 1,
    metric = "cosine", umap.method = "umap-learn", verbose = FALSE
  )

  # Stress scoring
  obj <- AddModuleScore(obj, features = list(hypoxia_genes), name = "HypoxiaScore", verbose = FALSE)
  obj$HypoxiaScore <- obj$HypoxiaScore1
  obj <- AddModuleScore(obj, features = list(glycolysis_genes), name = "GlycolysisScore", verbose = FALSE)
  obj$GlycolysisScore <- obj$GlycolysisScore1
  obj
}

stress_summary <- function(obj) {
  meta <- obj@meta.data
  df <- meta %>%
    dplyr::mutate(cluster = as.character(seurat_clusters)) %>%
    dplyr::group_by(cluster) %>%
    dplyr::summarise(mean_hypoxia = mean(HypoxiaScore), mean_glycolysis = mean(GlycolysisScore), cells = dplyr::n(), .groups = "drop")
  stressed <- df %>% dplyr::filter(mean_hypoxia > 0.5 | mean_glycolysis > 0.5)
  list(table = df, stressed_clusters = stressed$cluster)
}

save_plots <- function(obj, outdir, markers, seed) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  common_args <- list(reduction = "umap", pt.size = 0.4)

  set.seed(seed)
  p1 <- DimPlot(obj, group.by = "seurat_clusters", label = TRUE, repel = TRUE, shuffle = FALSE, order = TRUE, raster = TRUE, !!!common_args) +
    ggtitle("UMAP by Louvain cluster (res=2.0, dims 1:20)")
  ggsave(file.path(outdir, "umap_by_cluster.pdf"), p1, width = 8, height = 6)
  ggsave(file.path(outdir, "umap_by_cluster.png"), p1, width = 8, height = 6, dpi = 300)

  if ("domain" %in% colnames(obj@meta.data)) {
    p2 <- DimPlot(obj, group.by = "domain", shuffle = FALSE, order = TRUE, raster = TRUE, !!!common_args) +
      ggtitle("UMAP by domain (vFB vs dFB)")
    ggsave(file.path(outdir, "umap_by_domain.pdf"), p2, width = 8, height = 6)
    ggsave(file.path(outdir, "umap_by_domain.png"), p2, width = 8, height = 6, dpi = 300)
  }
  if ("sample_id" %in% colnames(obj@meta.data)) {
    p3 <- DimPlot(obj, group.by = "sample_id", shuffle = FALSE, order = TRUE, raster = TRUE, !!!common_args) +
      ggtitle("UMAP by sample")
    ggsave(file.path(outdir, "umap_by_sample.pdf"), p3, width = 8, height = 6)
    ggsave(file.path(outdir, "umap_by_sample.png"), p3, width = 8, height = 6, dpi = 300)
  }

  # Marker feature plots
  markers_present <- markers[markers %in% rownames(obj)]
  if (length(markers_present) > 0) {
    fp <- FeaturePlot(obj, features = markers_present, reduction = "umap", order = TRUE, raster = TRUE, pt.size = 0.3) +
      ggtitle("Marker genes")
    ggsave(file.path(outdir, "umap_marker_features.pdf"), fp, width = 10, height = 8)
    ggsave(file.path(outdir, "umap_marker_features.png"), fp, width = 10, height = 8, dpi = 300)
  }
}

write_report <- function(path, info) {
  con <- file(path, "w")
  on.exit(close(con))
  writeLines(info, con)
}

main <- function() {
  check_R_version()
  require_version("Seurat", "4.2.0")
  args <- parse_args()

  project_root <- normalizePath(args$project_root, mustWork = TRUE)
  results_root <- file.path(project_root, "results", "walsh_day75_repro")
  plots_dir <- file.path(results_root, "plots")
  dir.create(results_root, recursive = TRUE, showWarnings = FALSE)

  seeds <- list(
    hvg = 1001L,
    cc = 1002L,
    pca = 1003L,
    neighbors = 1004L,
    clusters = 1005L,
    umap = 1006L,
    plots = 1007L
  )

  hypoxia_genes <- read_gene_list(args$hypoxia_genes)
  glycolysis_genes <- read_gene_list(args$glycolysis_genes)

  log_msg("Project root:", project_root)
  log_norm_settings()
  log_msg("HVG method: vst, nfeatures=5000")
  log_msg("Cell cycle regression vars: S.Score, G2M.Score")
  log_msg("Stress gene lists (hypoxia):", paste(hypoxia_genes, collapse = ", "))
  log_msg("Stress gene lists (glycolysis):", paste(glycolysis_genes, collapse = ", "))

  paths <- load_day75_counts(project_root)
  counts_dfb <- read_counts_csv(paths$dFB)
  counts_vfb <- read_counts_csv(paths$vFB)

  # Validate gene naming and orientation
  log_msg("Counts dFB: genes", nrow(counts_dfb), "cells", ncol(counts_dfb))
  log_msg("Counts vFB: genes", nrow(counts_vfb), "cells", ncol(counts_vfb))
  if (length(intersect(rownames(counts_dfb), rownames(counts_vfb))) < 1000) {
    stop("Gene sets between samples have low overlap; check orientation/naming")
  }
  log_msg("Mito genes identified by pattern: ^MT-")

  obj_dfb <- build_seurat(counts_dfb, "GSM7979671", "dFB_domain")
  obj_vfb <- build_seurat(counts_vfb, "GSM7979672", "vFB_domain")
  combined <- merge(obj_dfb, y = obj_vfb, add.cell.ids = c("GSM7979671", "GSM7979672"), project = "Walsh_d75_repro")

  combined_qc <- qc_filter(combined)
  saveRDS(combined_qc, file.path(results_root, "walsh_day75_postQC.rds"))
  log_msg("Post-QC cells:", ncol(combined_qc), "features:", nrow(combined_qc))

  # Initial pipeline
  obj_initial <- run_pipeline(combined_qc, seeds, hypoxia_genes, glycolysis_genes)
  saveRDS(obj_initial, file.path(results_root, "walsh_day75_pre_stress_filter.rds"))

  stress_info <- stress_summary(obj_initial)
  write.csv(stress_info$table, file.path(results_root, "walsh_day75_stress_scores_by_cluster.csv"), row.names = FALSE)
  log_msg("Clusters flagged for stress removal (mean > 0.5):", paste(stress_info$stressed_clusters, collapse = ", "))

  if (length(stress_info$stressed_clusters) > 0) {
    obj_filtered <- subset(obj_initial, idents = setdiff(levels(obj_initial), stress_info$stressed_clusters))
  } else {
    obj_filtered <- obj_initial
  }
  saveRDS(obj_filtered, file.path(results_root, "walsh_day75_post_stress_prerun.rds"))
  log_msg("Cells after stress-cluster removal:", ncol(obj_filtered))

  # Rerun pipeline after stress removal
  obj_final <- run_pipeline(obj_filtered, seeds, hypoxia_genes, glycolysis_genes)
  saveRDS(obj_final, file.path(results_root, "walsh_day75_final.rds"))

  # Log gene presence for stress lists
  hypoxia_present <- hypoxia_genes[hypoxia_genes %in% rownames(obj_final)]
  glycolysis_present <- glycolysis_genes[glycolysis_genes %in% rownames(obj_final)]
  log_msg("Hypoxia genes present:", paste(hypoxia_present, collapse = ", "))
  log_msg("Glycolysis genes present:", paste(glycolysis_present, collapse = ", "))

  markers <- c("NKX2-1", "DLX2", "GAD1", "GAD2", "SST")
  save_plots(obj_final, plots_dir, markers, seeds$plots)

  # Report
  report_lines <- c(
    sprintf("R version: %s", getRversion()),
    sprintf("Seurat version: %s", packageVersion("Seurat")),
    sprintf("Project root: %s", project_root),
    sprintf("Counts source: %s", file.path(project_root, "data", "raw", "walsh_2025_geo_files", "suppl", "GSE250482_RAW.tar")),
    sprintf("Samples: GSM7979671 (dFB), GSM7979672 (vFB)"),
    sprintf("MT pattern: ^MT- ; QC filters: 1000-5000 genes, <15%% MT"),
    sprintf("Normalization: LogNormalize, scale.factor=1e4"),
    sprintf("HVG: vst, nfeatures=5000"),
    sprintf("Cell cycle regression: S.Score, G2M.Score"),
    sprintf("Neighbors: dims=1:20, k.param=20, algorithm=1 (exact)"),
    sprintf("Clustering: Louvain (algorithm=1), resolution=2.0"),
    sprintf("UMAP: dims=1:20, method=umap-learn, metric=cosine, n.neighbors=30, min.dist=0.3, spread=1, seed=%s", seeds$umap),
    sprintf("Seeds: hvg=%s, cc=%s, pca=%s, neighbors=%s, clusters=%s, umap=%s, plots=%s",
            seeds$hvg, seeds$cc, seeds$pca, seeds$neighbors, seeds$clusters, seeds$umap, seeds$plots),
    sprintf("Hypoxia genes: %s", paste(hypoxia_genes, collapse = ", ")),
    sprintf("Glycolysis genes: %s", paste(glycolysis_genes, collapse = ", ")),
    sprintf("Hypoxia present: %s", paste(hypoxia_present, collapse = ", ")),
    sprintf("Glycolysis present: %s", paste(glycolysis_present, collapse = ", ")),
    sprintf("Stress removal threshold: mean > 0.5; removed clusters: %s", paste(stress_info$stressed_clusters, collapse = ", ")),
    sprintf("Markers plotted: %s", paste(markers, collapse = ", ")),
    sprintf("Outputs: %s", paste(c(
      "walsh_day75_postQC.rds",
      "walsh_day75_pre_stress_filter.rds",
      "walsh_day75_post_stress_prerun.rds",
      "walsh_day75_final.rds",
      "walsh_day75_stress_scores_by_cluster.csv",
      "plots/*.pdf/png"
    ), collapse = "; "))
  )
  write_report(file.path(results_root, "walsh_day75_reproduction_report.txt"), report_lines)

  sink(file.path(results_root, "sessionInfo_walsh_day75_repro.txt"))
  print(sessionInfo())
  sink()
  log_msg("Done. Results in", results_root)
}

main()
