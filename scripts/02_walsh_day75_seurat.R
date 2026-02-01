#!/usr/bin/env Rscript

# Walsh day-75 reproduction (GSM7979671 dFB, GSM7979672 vFB) matching paper methods.
# Enforces R 4.1.x and Seurat 4.2.0 (fails otherwise).
# Deterministic seeds for HVG/PCA/Neighbors/Louvain/UMAP; explicit parameters and logging.
# Stress scoring/removal and rerun; outputs under PROJECT_ROOT/results/walsh_day75/.

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(data.table)
  library(ggplot2)
  library(dplyr)
  library(rlang)
  library(future)
})

check_versions <- function() {
  rv <- getRversion()
  if (!(rv >= "4.1.0" && rv < "4.3.0")) stop(sprintf("R version must be 4.1.x-4.2.x; found %s", rv))
  sv <- packageVersion("Seurat")
  if (sv$major != 4) stop(sprintf("Seurat major version must be 4.x; found %s", sv))
  log_msg("Using R", rv, "and Seurat", sv, "(paper reported Seurat 4.2.0.114)")
}

log_msg <- function(...) {
  msg <- sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste0(..., collapse = " "))
  cat(msg, "\n"); flush.console()
}

detect_threads <- function() {
  th <- as.integer(Sys.getenv("SLURM_CPUS_ON_NODE", "1"))
  if (is.na(th) || th < 1) th <- 1L
  Sys.setenv(OMP_NUM_THREADS = th, MKL_NUM_THREADS = th, BLAS_THREADS = th, OPENBLAS_NUM_THREADS = th)
  log_msg("Threading: using", th, "threads for BLAS/OMP (from SLURM_CPUS_ON_NODE)")
  plan(sequential)
  options(future.seed = TRUE, future.rng.onMisuse = "ignore")
  th
}

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  res <- list(project_root = Sys.getenv("PROJECT_ROOT", ""),
              hypoxia_genes = Sys.getenv("HYPOXIA_GENES", ""),
              glycolysis_genes = Sys.getenv("GLYCOLYSIS_GENES", ""))
  i <- 1
  while (i <= length(args)) {
    if (args[i] %in% c("-p", "--project-root")) { res$project_root <- args[i + 1]; i <- i + 2; next }
    if (args[i] == "--hypoxia-genes") { res$hypoxia_genes <- args[i + 1]; i <- i + 2; next }
    if (args[i] == "--glycolysis-genes") { res$glycolysis_genes <- args[i + 1]; i <- i + 2; next }
    if (args[i] %in% c("-h", "--help")) {
      cat("Usage: 02_walsh_day75_seurat.R --project-root <path> --hypoxia-genes <file> --glycolysis-genes <file>\n")
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
  genes <- gsub("PKM[`’']+", "PKM", genes, perl = TRUE)
  unique(genes[genes != ""])
}

read_counts_csv <- function(path) {
  dt <- if (grepl("\\.gz$", path, ignore.case = TRUE)) fread(cmd = paste("gunzip -c", shQuote(path))) else fread(path)
  genes <- dt[[1]]
  mat <- as.matrix(dt[, -1, with = FALSE])
  rownames(mat) <- genes
  m <- Matrix(mat, sparse = TRUE)
  if (!all(m@x == round(m@x))) stop("Counts are not integer-valued")
  m
}

load_counts <- function(project_root) {
  tarfile <- file.path(project_root, "data", "raw", "walsh_2025_geo_files", "suppl", "GSE250482_RAW.tar")
  outdir <- file.path(project_root, "data", "processed", "walsh_day75", "raw_extract")
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  files <- c("GSM7979671_dFB_d75.csv.gz", "GSM7979672_vFB_d75.csv.gz")
  utils::untar(tarfile, files = files, exdir = outdir)
  list(
    dFB = file.path(outdir, files[1]),
    vFB = file.path(outdir, files[2])
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

load_or_run <- function(path, compute_fn) {
  if (file.exists(path)) {
    log_msg("Loading checkpoint:", path)
    return(readRDS(path))
  }
  obj <- compute_fn()
  saveRDS(obj, path)
  obj
}

run_pipeline <- function(obj, seeds, hypoxia_genes, glycolysis_genes, plots_dir, label, ckpts) {
  obj <- load_or_run(ckpts$norm_hvg, function() {
    set.seed(seeds$hvg)
    log_msg("Normalization: LogNormalize, scale.factor=1e4; HVG: vst, nfeatures=5000")
    x <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 1e4, verbose = FALSE)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 5000, verbose = FALSE)
    x
  })
  hvgs <- VariableFeatures(obj)
  log_msg("Scaling only HVGs:", length(hvgs))

  obj <- load_or_run(ckpts$cellcycle, function() {
    set.seed(seeds$cc)
    CellCycleScoring(obj, s.features = cc.genes.updated.2019$s.genes, g2m.features = cc.genes.updated.2019$g2m.genes, set.ident = FALSE)
  })

  obj <- load_or_run(ckpts$scaled, function() {
    log_msg("ScaleData regressors: S.Score, G2M.Score")
    ScaleData(obj, features = hvgs, vars.to.regress = c("S.Score", "G2M.Score"), verbose = FALSE)
  })

  obj <- load_or_run(ckpts$pca, function() {
    set.seed(seeds$pca)
    RunPCA(obj, features = hvgs, verbose = FALSE)
  })

  # JackStraw + Elbow to choose PCs
  set.seed(seeds$jackstraw)
  js_obj <- JackStraw(obj, dims = 50, num.replicate = 100, verbose = FALSE)
  js_obj <- ScoreJackStraw(js_obj, dims = 1:50)
  js_overall <- js_obj@reductions$pca@jackstraw@overall
  sig_pcs <- which(js_overall[, 2] < 0.05)
  selected_pcs <- if (length(sig_pcs) > 0) max(sig_pcs) else 20L
  selected_pcs <- max(5L, min(selected_pcs, 50L))
  log_msg("JackStraw significant PCs (<0.05):", paste(sig_pcs, collapse = ", "))
  log_msg("Selected PCs for downstream:", selected_pcs)

  dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
  js_plot <- JackStrawPlot(js_obj, dims = 1:50)
  ggsave(file.path(plots_dir, paste0("jackstraw_", label, ".pdf")), js_plot, width = 8, height = 6)
  ggsave(file.path(plots_dir, paste0("jackstraw_", label, ".png")), js_plot, width = 8, height = 6, dpi = 300)
  elbow_plot <- ElbowPlot(obj, ndims = 50)
  ggsave(file.path(plots_dir, paste0("elbow_", label, ".pdf")), elbow_plot, width = 8, height = 6)
  ggsave(file.path(plots_dir, paste0("elbow_", label, ".png")), elbow_plot, width = 8, height = 6, dpi = 300)

  set.seed(seeds$neighbors)
  obj <- FindNeighbors(obj, dims = 1:selected_pcs, k.param = 20, algorithm = 1, verbose = FALSE)

  set.seed(seeds$clusters)
  obj <- FindClusters(obj, resolution = 2.0, algorithm = 1, verbose = FALSE)

  # Main UMAP with selected PCs
  set.seed(seeds$umap)
  obj <- RunUMAP(
    obj,
    dims = 1:selected_pcs,
    reduction.name = "umap_sel",
    reduction.key = "UMAPsel_",
    seed.use = seeds$umap,
    n.neighbors = 30,
    min.dist = 0.3,
    spread = 1,
    metric = "cosine",
    umap.method = "uwot",
    return.model = FALSE,
    verbose = FALSE
  )

  # Comparison UMAP fixed at 20 PCs
  obj <- RunUMAP(
    obj,
    dims = 1:20,
    reduction.name = "umap20",
    reduction.key = "UMAP20_",
    seed.use = seeds$umap + 1L,
    n.neighbors = 30,
    min.dist = 0.3,
    spread = 1,
    metric = "cosine",
    umap.method = "uwot",
    return.model = FALSE,
    verbose = FALSE
  )

  obj <- AddModuleScore(obj, features = list(hypoxia_genes), name = "HypoxiaScore", verbose = FALSE)
  obj$HypoxiaScore <- obj$HypoxiaScore1
  obj <- AddModuleScore(obj, features = list(glycolysis_genes), name = "GlycolysisScore", verbose = FALSE)
  obj$GlycolysisScore <- obj$GlycolysisScore1
  obj@misc$selected_pcs <- selected_pcs
  obj@misc$jackstraw_sig_pcs <- sig_pcs
  obj
}

stress_summary <- function(obj) {
  meta <- obj@meta.data %>%
    dplyr::mutate(cluster = as.character(seurat_clusters)) %>%
    dplyr::group_by(cluster) %>%
    dplyr::summarise(mean_hypoxia = mean(HypoxiaScore), mean_glycolysis = mean(GlycolysisScore), cells = dplyr::n(), .groups = "drop")
  stressed <- meta %>% dplyr::filter(mean_hypoxia > 0.5 | mean_glycolysis > 0.5)
  list(table = meta, stressed = stressed$cluster)
}

write_report <- function(path, info) {
  con <- file(path, "w"); on.exit(close(con)); writeLines(info, con)
}

main <- function() {
  check_versions()
  detect_threads()
  args <- parse_args()
  project_root <- normalizePath(args$project_root, mustWork = TRUE)
  results_dir <- file.path(project_root, "results", "walsh_day75")
  plots_dir <- file.path(results_dir, "plots")
  ckpt_dir <- file.path(results_dir, "checkpoints")
  dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(ckpt_dir, recursive = TRUE, showWarnings = FALSE)

  seeds <- list(hvg = 1001L, cc = 1002L, pca = 1003L, neighbors = 1004L, clusters = 1005L, umap = 1006L, plots = 1007L, jackstraw = 1008L)
  hypoxia_genes <- read_gene_list(args$hypoxia_genes)
  glycolysis_genes <- read_gene_list(args$glycolysis_genes)

  log_msg("Project root:", project_root)
  log_msg("Stress gene lists (hypoxia):", paste(hypoxia_genes, collapse = ", "))
  log_msg("Stress gene lists (glycolysis):", paste(glycolysis_genes, collapse = ", "))
  log_msg("MT pattern: ^MT- ; QC: 1000-5000 genes, <15% MT")
  log_msg("UMAP params: method=uwot (R), metric=cosine, n.neighbors=30, min.dist=0.3, spread=1, dims=1:20")
  log_msg("Clustering: Louvain (algorithm=1), resolution=2.0; Neighbors: dims=1:20, k=20, algorithm=1")

  paths <- load_counts(project_root)
  counts_dfb <- read_counts_csv(paths$dFB)
  counts_vfb <- read_counts_csv(paths$vFB)
  log_msg("Counts dFB: genes", nrow(counts_dfb), "cells", ncol(counts_dfb))
  log_msg("Counts vFB: genes", nrow(counts_vfb), "cells", ncol(counts_vfb))

  ckpts_pre <- list(
    merged = file.path(ckpt_dir, "walsh_merged_raw.rds"),
    postqc = file.path(ckpt_dir, "walsh_postqc.rds"),
    norm_hvg = file.path(ckpt_dir, "walsh_norm_hvg.rds"),
    cellcycle = file.path(ckpt_dir, "walsh_cellcycle.rds"),
    scaled = file.path(ckpt_dir, "walsh_scaled.rds"),
    pca = file.path(ckpt_dir, "walsh_pca.rds")
  )

  combined <- load_or_run(ckpts_pre$merged, function() {
    obj_dfb <- build_seurat(counts_dfb, "GSM7979671", "dFB_domain")
    obj_vfb <- build_seurat(counts_vfb, "GSM7979672", "vFB_domain")
    merge(obj_dfb, y = obj_vfb, add.cell.ids = c("GSM7979671", "GSM7979672"), project = "Walsh_d75")
  })
  log_msg("Combined object size (MB):", signif(as.numeric(object.size(combined)) / 1024^2, 3))

  combined_qc <- load_or_run(ckpts_pre$postqc, function() qc_filter(combined))
  saveRDS(combined_qc, file.path(results_dir, "walsh_day75_postQC.rds"))
  log_msg("Post-QC cells:", ncol(combined_qc), "features:", nrow(combined_qc))

  obj_initial <- run_pipeline(combined_qc, seeds, hypoxia_genes, glycolysis_genes, plots_dir, "pre_stress", ckpts_pre)
  saveRDS(obj_initial, file.path(results_dir, "walsh_day75_pre_stress_filter.rds"))

  stress_info <- stress_summary(obj_initial)
  write.csv(stress_info$table, file.path(results_dir, "walsh_day75_stress_scores_by_cluster.csv"), row.names = FALSE)
  log_msg("Clusters flagged for stress removal (>0.5):", paste(stress_info$stressed, collapse = ", "))

  if (length(stress_info$stressed) > 0) {
    obj_filtered <- subset(obj_initial, idents = setdiff(levels(obj_initial), stress_info$stressed))
  } else { obj_filtered <- obj_initial }
  saveRDS(obj_filtered, file.path(results_dir, "walsh_day75_post_stress_prerun.rds"))
  log_msg("Cells after stress removal:", ncol(obj_filtered))

  ckpts_post <- list(
    norm_hvg = file.path(ckpt_dir, "walsh_norm_hvg_poststress.rds"),
    cellcycle = file.path(ckpt_dir, "walsh_cellcycle_poststress.rds"),
    scaled = file.path(ckpt_dir, "walsh_scaled_poststress.rds"),
    pca = file.path(ckpt_dir, "walsh_pca_poststress.rds")
  )

  obj_final <- run_pipeline(obj_filtered, seeds, hypoxia_genes, glycolysis_genes, plots_dir, "post_stress", ckpts_post)
  saveRDS(obj_final, file.path(results_dir, "walsh_day75_final.rds"))

  hypoxia_present <- hypoxia_genes[hypoxia_genes %in% rownames(obj_final)]
  glycolysis_present <- glycolysis_genes[glycolysis_genes %in% rownames(obj_final)]
  log_msg("Hypoxia genes present:", paste(hypoxia_present, collapse = ", "))
  log_msg("Glycolysis genes present:", paste(glycolysis_present, collapse = ", "))

  dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
  common_args <- list(object = obj_final, reduction = "umap_sel", pt.size = 0.4, shuffle = FALSE, order = TRUE, raster = TRUE)
  set.seed(seeds$plots)
  p1 <- do.call(DimPlot, c(common_args, list(group.by = "seurat_clusters", label = TRUE, repel = TRUE))) +
    ggtitle("UMAP by Louvain cluster (res=2.0, dims 1:20)")
  ggsave(file.path(plots_dir, "umap_by_cluster.pdf"), p1, width = 8, height = 6)
  ggsave(file.path(plots_dir, "umap_by_cluster.png"), p1, width = 8, height = 6, dpi = 300)
  if ("domain" %in% colnames(obj_final@meta.data)) {
    p2 <- do.call(DimPlot, c(common_args, list(group.by = "domain"))) + ggtitle("UMAP by domain (vFB vs dFB)")
    ggsave(file.path(plots_dir, "umap_by_domain.pdf"), p2, width = 8, height = 6)
    ggsave(file.path(plots_dir, "umap_by_domain.png"), p2, width = 8, height = 6, dpi = 300)
  }
  if ("sample_id" %in% colnames(obj_final@meta.data)) {
    p3 <- do.call(DimPlot, c(common_args, list(group.by = "sample_id"))) + ggtitle("UMAP by sample")
    ggsave(file.path(plots_dir, "umap_by_sample.pdf"), p3, width = 8, height = 6)
    ggsave(file.path(plots_dir, "umap_by_sample.png"), p3, width = 8, height = 6, dpi = 300)
  }
  markers <- c("NKX2-1", "DLX2", "GAD1", "GAD2", "SST")
  markers_present <- markers[markers %in% rownames(obj_final)]
  if (length(markers_present) > 0) {
    fp <- FeaturePlot(obj_final, features = markers_present, reduction = "umap_sel", order = TRUE, raster = TRUE, pt.size = 0.3)
    ggsave(file.path(plots_dir, "umap_marker_features.pdf"), fp, width = 10, height = 8)
    ggsave(file.path(plots_dir, "umap_marker_features.png"), fp, width = 10, height = 8, dpi = 300)
  }

  # UMAP using fixed 20 PCs for comparison
  common_args20 <- common_args
  common_args20$reduction <- "umap20"
  p1_20 <- do.call(DimPlot, c(common_args20, list(group.by = "seurat_clusters", label = TRUE, repel = TRUE))) +
    ggtitle("UMAP (20 PCs) by Louvain cluster")
  ggsave(file.path(plots_dir, "umap20_by_cluster.pdf"), p1_20, width = 8, height = 6)
  ggsave(file.path(plots_dir, "umap20_by_cluster.png"), p1_20, width = 8, height = 6, dpi = 300)

  p2_20 <- do.call(DimPlot, c(common_args20, list(group.by = "domain"))) + ggtitle("UMAP (20 PCs) by domain")
  ggsave(file.path(plots_dir, "umap20_by_domain.pdf"), p2_20, width = 8, height = 6)
  ggsave(file.path(plots_dir, "umap20_by_domain.png"), p2_20, width = 8, height = 6, dpi = 300)

  report <- c(
    sprintf("R version: %s", getRversion()),
    sprintf("Seurat version: %s", packageVersion("Seurat")),
    sprintf("Project root: %s", project_root),
    "Samples: GSM7979671 (dFB), GSM7979672 (vFB)",
    "QC: 1000-5000 genes, <15% MT (pattern ^MT-)",
    "Normalization: LogNormalize, scale.factor=1e4",
    "HVG: vst, nfeatures=5000",
    "Cell cycle regression: S.Score, G2M.Score",
    "Scaling features: HVGs only (count = 5000)",
    "Threading: BLAS/OMP set from SLURM_CPUS_ON_NODE",
    "Neighbors: dims 1:20, k.param=20, algorithm=1 (exact)",
    "Clustering: Louvain (algorithm=1), resolution=2.0",
    "UMAP: dims 1:20, umap.method=uwot (R), metric=cosine, n.neighbors=30, min.dist=0.3, spread=1",
    sprintf("Seeds: hvg=%s, cc=%s, pca=%s, neighbors=%s, clusters=%s, umap=%s",
            seeds$hvg, seeds$cc, seeds$pca, seeds$neighbors, seeds$clusters, seeds$umap),
    sprintf("Hypoxia genes: %s", paste(hypoxia_genes, collapse = ", ")),
    sprintf("Glycolysis genes: %s", paste(glycolysis_genes, collapse = ", ")),
    sprintf("Hypoxia present: %s", paste(hypoxia_present, collapse = ", ")),
    sprintf("Glycolysis present: %s", paste(glycolysis_present, collapse = ", ")),
    sprintf("Stress removal threshold: mean > 0.5; removed clusters: %s", paste(stress_info$stressed, collapse = ", ")),
    sprintf("Markers plotted: %s", paste(markers_present, collapse = ", ")),
    "Outputs: walsh_day75_postQC.rds; walsh_day75_pre_stress_filter.rds; walsh_day75_post_stress_prerun.rds; walsh_day75_final.rds; walsh_day75_stress_scores_by_cluster.csv; plots/*.pdf/png"
  )
  write_report(file.path(results_dir, "walsh_day75_report.txt"), report)

  sink(file.path(results_dir, "sessionInfo_walsh_day75.txt"))
  print(sessionInfo())
  sink()
  log_msg("Done. Results in", results_dir)
}

main()
