#!/usr/bin/env Rscript

# Walsh day-75 vFB/dFB domain reconstruction (Fig 2A scope).
# Requires: R 4.1.x, Seurat 4.2.0, Matrix, data.table.
# Input: PROJECT_ROOT env or --project-root; hypoxia/glycolysis gene list files (one gene symbol per line).
# Output: RDS checkpoints and logs under PROJECT_ROOT/results and PROJECT_ROOT/logs.

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(data.table)
  library(dplyr)
})

log_msg <- function(...) {
  msg <- sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste0(..., collapse = " "))
  cat(msg, "\n")
  flush.console()
}

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  res <- list(project_root = Sys.getenv("PROJECT_ROOT", ""), hypoxia_genes = NULL, glycolysis_genes = NULL, seed = 777)
  i <- 1
  while (i <= length(args)) {
    if (args[i] %in% c("-p", "--project-root")) { res$project_root <- args[i + 1]; i <- i + 2; next }
    if (args[i] == "--hypoxia-genes") { res$hypoxia_genes <- args[i + 1]; i <- i + 2; next }
    if (args[i] == "--glycolysis-genes") { res$glycolysis_genes <- args[i + 1]; i <- i + 2; next }
    if (args[i] == "--seed") { res$seed <- as.integer(args[i + 1]); i <- i + 2; next }
    if (args[i] %in% c("-h", "--help")) {
      cat("Usage: 02_walsh_day75_seurat.R --project-root <path> --hypoxia-genes <file> --glycolysis-genes <file> [--seed 777]\n")
      quit(status = 0)
    }
    i <- i + 1
  }
  if (res$project_root == "") stop("PROJECT_ROOT is required via env or --project-root")
  if (is.null(res$hypoxia_genes) || is.null(res$glycolysis_genes)) {
    stop("Both --hypoxia-genes and --glycolysis-genes are required (one gene symbol per line).")
  }
  res
}

read_gene_list <- function(path) {
  if (!file.exists(path)) stop(sprintf("Gene list file not found: %s", path))
  genes <- scan(path, what = character(), quiet = TRUE)
  unique(genes[genes != ""])
}

read_counts_csv <- function(path) {
  if (!file.exists(path)) stop(sprintf("Counts file not found: %s", path))
  log_msg("Reading counts:", path)
  if (grepl("\\.gz$", path, ignore.case = TRUE)) {
    dt <- fread(cmd = paste("gunzip -c", shQuote(path)))
  } else {
    dt <- fread(path)
  }
  gene_col <- dt[[1]]
  if (anyNA(gene_col)) stop("Gene column contains NA values")
  mat <- as.matrix(dt[, -1, with = FALSE])
  rownames(mat) <- gene_col
  Matrix(mat, sparse = TRUE)
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

run_pipeline <- function(obj, hypoxia_genes, glycolysis_genes, seed_val) {
  set.seed(seed_val)
  log_msg("Normalizing data")
  obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 1e4, verbose = FALSE)
  log_msg("Finding variable features (top 5000)")
  obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 5000, verbose = FALSE)
  log_msg("Cell cycle scoring")
  obj <- CellCycleScoring(obj, s.features = cc.genes.updated.2019$s.genes, g2m.features = cc.genes.updated.2019$g2m.genes, set.ident = FALSE)
  log_msg("Scaling data with cell cycle regression")
  obj <- ScaleData(obj, features = rownames(obj), vars.to.regress = c("S.Score", "G2M.Score"), verbose = FALSE)
  log_msg("Running PCA")
  obj <- RunPCA(obj, features = VariableFeatures(obj), verbose = FALSE)
  log_msg("Neighbors/UMAP/Clustering (dims 1:20, resolution 2.0)")
  obj <- FindNeighbors(obj, dims = 1:20, verbose = FALSE)
  obj <- FindClusters(obj, resolution = 2.0, verbose = FALSE)
  obj <- RunUMAP(obj, dims = 1:20, seed.use = seed_val, verbose = FALSE)

  log_msg("Scoring stress signatures")
  obj <- AddModuleScore(obj, features = list(hypoxia_genes), name = "HypoxiaScore", verbose = FALSE)
  obj$HypoxiaScore <- obj$HypoxiaScore1
  obj <- AddModuleScore(obj, features = list(glycolysis_genes), name = "GlycolysisScore", verbose = FALSE)
  obj$GlycolysisScore <- obj$GlycolysisScore1
  obj
}

identify_stress_clusters <- function(obj) {
  meta <- obj@meta.data %>%
    dplyr::mutate(cluster = as.character(seurat_clusters)) %>%
    dplyr::group_by(cluster) %>%
    dplyr::summarise(mean_hypoxia = mean(HypoxiaScore), mean_glycolysis = mean(GlycolysisScore), .groups = "drop")
  stressed <- meta %>% dplyr::filter(mean_hypoxia > 0.5 | mean_glycolysis > 0.5)
  list(table = meta, stressed = stressed$cluster)
}

main <- function() {
  args <- parse_args()
  project_root <- normalizePath(args$project_root, mustWork = TRUE)
  raw_tar <- file.path(project_root, "data", "raw", "walsh_2025_geo_files", "suppl", "GSE250482_RAW.tar")
  extract_dir <- file.path(project_root, "data", "processed", "walsh_day75", "raw_extract")
  results_dir <- file.path(project_root, "results", "walsh_day75")
  dir.create(extract_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

  log_msg("Project root:", project_root)
  log_msg("Using hypoxia genes from:", args$hypoxia_genes)
  log_msg("Using glycolysis genes from:", args$glycolysis_genes)
  hypoxia_genes <- read_gene_list(args$hypoxia_genes)
  glycolysis_genes <- read_gene_list(args$glycolysis_genes)

  if (!file.exists(raw_tar)) stop(sprintf("Raw tar not found at %s", raw_tar))

  files_to_extract <- c("GSM7979671_dFB_d75.csv.gz", "GSM7979672_vFB_d75.csv.gz")
  log_msg("Extracting day-75 files from tar")
  utils::untar(raw_tar, files = files_to_extract, exdir = extract_dir)

  counts_dfb <- read_counts_csv(file.path(extract_dir, "GSM7979671_dFB_d75.csv.gz"))
  counts_vfb <- read_counts_csv(file.path(extract_dir, "GSM7979672_vFB_d75.csv.gz"))

  log_msg("Building Seurat objects")
  obj_dfb <- build_seurat(counts_dfb, "GSM7979671", "dFB_domain")
  obj_vfb <- build_seurat(counts_vfb, "GSM7979672", "vFB_domain")

  log_msg("Merging objects")
  combined <- merge(obj_dfb, y = obj_vfb, add.cell.ids = c("GSM7979671", "GSM7979672"), project = "Walsh_d75")

  log_msg("QC filtering (1,000–5,000 genes, <15% MT)")
  combined_qc <- qc_filter(combined)
  saveRDS(combined_qc, file.path(results_dir, "walsh_day75_postQC.rds"))
  log_msg("Post-QC cells:", ncol(combined_qc))

  log_msg("Running initial pipeline (pre-stress filtering)")
  obj_initial <- run_pipeline(combined_qc, hypoxia_genes, glycolysis_genes, args$seed)
  saveRDS(obj_initial, file.path(results_dir, "walsh_day75_pre_stress_filter.rds"))

  stress_info <- identify_stress_clusters(obj_initial)
  write.csv(stress_info$table, file.path(results_dir, "walsh_day75_stress_scores_by_cluster.csv"), row.names = FALSE)
  log_msg("Clusters flagged for removal (mean score > 0.5):", paste(stress_info$stressed, collapse = ", "))

  if (length(stress_info$stressed) > 0) {
    obj_filtered <- subset(obj_initial, idents = setdiff(levels(obj_initial), stress_info$stressed))
  } else {
    obj_filtered <- obj_initial
  }
  saveRDS(obj_filtered, file.path(results_dir, "walsh_day75_post_stress_prerun.rds"))
  log_msg("Cells after stress-filter cluster removal:", ncol(obj_filtered))

  log_msg("Re-running pipeline after stress cluster removal")
  obj_final <- run_pipeline(obj_filtered, hypoxia_genes, glycolysis_genes, args$seed)
  saveRDS(obj_final, file.path(results_dir, "walsh_day75_final.rds"))

  log_msg("Saving sessionInfo")
  sink(file.path(results_dir, "sessionInfo_walsh_day75.txt"))
  print(sessionInfo())
  sink()
  log_msg("Done.")
}

main()
