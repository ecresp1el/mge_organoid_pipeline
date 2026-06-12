#!/usr/bin/env Rscript

# Build a healthy-only Seurat object and UMAP for Liu et al. 2025 hnbMO
# scRNA-seq (GSE286235). GEO provides Cell Ranger raw_feature_bc_matrix.h5
# files, so this script applies explicit cell-level filtering before merging.
# The primary plot target is a Figure 1C-like UMAP of hnbMOs from H9 day 36,
# H9 day 63, and IMR90-4 day 63, colored by exploratory marker-based cell type.

suppressPackageStartupMessages({
  library(Matrix)
  library(Seurat)
  library(ggplot2)
})

parse_args <- function(args) {
  out <- list(
    `project-root` = Sys.getenv("PROJECT_ROOT", unset = ""),
    `h5-dir` = "",
    `tenx-dir` = "",
    `raw-tar` = "",
    outdir = "",
    `include-ds` = "false",
    `min-features` = "1360",
    `max-percent-mt` = "20",
    `gene-min-cells` = "3",
    nfeatures = "2000",
    dims = "10",
    resolution = "0.5",
    integrate = "false",
    seed = "11",
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
    if (a == "-p") a <- "--project-root"
    if (a == "-o") a <- "--outdir"
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
  cat(
    paste(
      "Usage:",
      "  Rscript scripts/05j_liu_2025_hnbmo_seurat.R --project-root <PROJECT_ROOT> [options]",
      "",
      "Options:",
      "  --h5-dir <path>         Directory containing GSM*_raw_feature_bc_matrix.h5 files",
      "  --tenx-dir <path>       Directory containing converted 10x folders by sample_id",
      "  --raw-tar <path>        Optional GSE286235_RAW.tar to extract if H5s are absent",
      "  --outdir <path>         Output directory (default: PROJECT_ROOT/results/liu_2025_hnbmo)",
      "  --include-ds <bool>     Include DS sample GSM8721443 (default: false)",
      "  --min-features <int>    Lower cell filter on detected genes (default: 1360; gives 14,249 healthy cells in a raw-H5 probe, close to the paper's 14,245)",
      "  --max-percent-mt <num>  Upper mitochondrial percent filter (default: 20)",
      "  --gene-min-cells <int>  Keep genes expressed in at least this many cells (default: 3)",
      "  --nfeatures <int>       HVGs for FindVariableFeatures (default: 2000)",
      "  --dims <int>            PCs used for neighbors/UMAP (default: 10)",
      "  --resolution <num>      Cluster resolution (default: 0.5)",
      "  --integrate <bool>      Use Seurat anchor integration across healthy samples (default: false)",
      "  --seed <int>            Random seed for PCA/UMAP/clustering steps (default: 11)",
      "",
      "Environment:",
      "  PROJECT_ROOT            Alternative way to pass --project-root",
      sep = "\n"
    )
  )
}

as_bool <- function(x) {
  tolower(as.character(x)) %in% c("1", "true", "t", "yes", "y")
}

trim_trailing_slash <- function(x) sub("/+$", "", x)

write_tsv <- function(df, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.table(df, path, sep = "\t", row.names = FALSE, quote = FALSE, na = "NA")
}

extract_accession <- function(path) {
  sub("^(GSM[0-9]+)_.*$", "\\1", basename(path))
}

choose_gene_expression_matrix <- function(x, path) {
  if (is.list(x)) {
    if ("Gene Expression" %in% names(x)) {
      return(x[["Gene Expression"]])
    }
    message("H5 returned multiple feature types for ", basename(path),
            "; using first entry: ", names(x)[[1]])
    return(x[[1]])
  }
  x
}

read_10x_matrix <- function(path, sample_id, tenx_dir) {
  converted_dir <- file.path(tenx_dir, sample_id)
  if (dir.exists(converted_dir) && file.exists(file.path(converted_dir, "matrix.mtx.gz"))) {
    message("Reading converted 10x directory for ", sample_id, ": ", converted_dir)
    x <- Read10X(data.dir = converted_dir, gene.column = 2, unique.features = TRUE)
    return(choose_gene_expression_matrix(x, converted_dir))
  }

  if (requireNamespace("hdf5r", quietly = TRUE)) {
    message("Reading Cell Ranger H5 for ", sample_id, ": ", path)
    x <- Read10X_h5(path, use.names = TRUE, unique.features = TRUE)
    return(choose_gene_expression_matrix(x, path))
  }

  stop(
    "Cannot read ", basename(path), ": hdf5r is not installed and converted 10x directory is absent: ",
    converted_dir,
    "\nRun scripts/05j_convert_gse286235_h5_to_10x.py first."
  )
}

get_expression_data <- function(obj) {
  tryCatch(
    GetAssayData(obj, assay = "RNA", layer = "data"),
    error = function(e) GetAssayData(obj, assay = "RNA", slot = "data")
  )
}

sanitize_label <- function(x) {
  gsub("[^A-Za-z0-9]+", "_", x)
}

cell_type_marker_sets <- list(
  `RG-div` = c("MKI67", "TOP2A", "UBE2C", "CENPF", "HMGB2"),
  RGC = c("SOX2", "VIM", "HES1", "NES", "PAX6"),
  IPC = c("EOMES", "NEUROD1", "NEUROG2", "ASCL1"),
  IM = c("DCX", "STMN2", "TUBB3", "MAP2"),
  GABA = c("GAD1", "GAD2", "DLX1", "DLX2", "DLX5"),
  CHN = c("CHAT", "SLC18A3", "SLC5A7", "LHX8", "ISL1", "NGFR")
)

add_marker_scores <- function(obj, marker_sets) {
  expr <- get_expression_data(obj)
  score_cols <- character()
  marker_presence <- list()

  for (cell_type in names(marker_sets)) {
    markers <- marker_sets[[cell_type]]
    present <- intersect(markers, rownames(expr))
    missing <- setdiff(markers, present)
    score_col <- paste0("liu2025_score_", sanitize_label(cell_type))
    score_cols <- c(score_cols, score_col)
    marker_presence[[cell_type]] <- data.frame(
      cell_type = cell_type,
      marker = markers,
      present = markers %in% present,
      stringsAsFactors = FALSE
    )

    if (length(present) == 0L) {
      obj[[score_col]] <- NA_real_
      next
    }
    obj[[score_col]] <- as.numeric(Matrix::colMeans(expr[present, , drop = FALSE]))
    if (length(missing) > 0L) {
      message("Missing marker(s) for ", cell_type, ": ", paste(missing, collapse = ", "))
    }
  }

  list(
    obj = obj,
    score_cols = score_cols,
    marker_presence = do.call(rbind, marker_presence)
  )
}

annotate_clusters_by_marker_scores <- function(obj, score_cols) {
  md <- obj@meta.data
  cluster_col <- "seurat_clusters"
  if (!(cluster_col %in% colnames(md))) stop("Missing seurat_clusters metadata")

  cluster_levels <- sort(unique(as.character(md[[cluster_col]])))
  summary_rows <- list()
  cluster_to_type <- setNames(rep(NA_character_, length(cluster_levels)), cluster_levels)

  for (cluster_id in cluster_levels) {
    idx <- as.character(md[[cluster_col]]) == cluster_id
    means <- vapply(score_cols, function(col) mean(md[[col]][idx], na.rm = TRUE), numeric(1))
    cell_types <- sub("^liu2025_score_", "", names(means))
    cell_types <- gsub("_", "-", cell_types)
    if (all(is.na(means))) {
      predicted <- "unassigned"
    } else {
      predicted <- cell_types[which.max(means)]
      if (predicted == "RG-div") predicted <- "RG-div"
    }
    cluster_to_type[[cluster_id]] <- predicted
    summary_rows[[cluster_id]] <- data.frame(
      seurat_cluster = cluster_id,
      predicted_cell_type = predicted,
      n_cells = sum(idx),
      score_name = names(means),
      mean_score = as.numeric(means),
      stringsAsFactors = FALSE
    )
  }

  obj$liu2025_exploratory_cell_type <- unname(cluster_to_type[as.character(md[[cluster_col]])])
  obj$liu2025_exploratory_cell_type <- factor(
    obj$liu2025_exploratory_cell_type,
    levels = c("RG-div", "RGC", "IPC", "IM", "GABA", "CHN", "unassigned")
  )
  list(obj = obj, cluster_score_summary = do.call(rbind, summary_rows))
}

sample_info <- data.frame(
  sample_geo_accession = c("GSM8721440", "GSM8721441", "GSM8721442", "GSM8721443"),
  sample_id = c("BF_H9_D36", "BF_H9_D63", "BFCO_IMR_D63", "BF_2DS3_D63"),
  sample_title = c(
    "scRNA-seq of hnbMOs from healthy people: 36 days in vitro",
    "scRNA-seq of hnbMOs from healthy people: 63 days in vitro, H9",
    "scRNA-seq of hnbMOs from healthy people: 63 days in vitro, IMR",
    "scRNA-seq of hnbMOs from DS patients: 63 days in vitro"
  ),
  disease_status = c("healthy", "healthy", "healthy", "DS"),
  day = c("D36", "D63", "D63", "D63"),
  day_numeric = c(36L, 63L, 63L, 63L),
  cell_line = c("H9", "H9", "IMR90-4", "2DS3"),
  stringsAsFactors = FALSE
)

opt <- parse_args(commandArgs(trailingOnly = TRUE))
if (isTRUE(opt$help)) {
  print_usage()
  quit(save = "no", status = 0)
}

if (!nzchar(opt$`project-root`)) {
  stop("PROJECT_ROOT or --project-root is required")
}
project_root <- trim_trailing_slash(opt$`project-root`)

raw_root <- file.path(project_root, "data/raw/liu_2025_hnbmo_geo_files")
h5_dir <- if (nzchar(opt$`h5-dir`)) opt$`h5-dir` else file.path(raw_root, "suppl")
tenx_dir <- if (nzchar(opt$`tenx-dir`)) opt$`tenx-dir` else file.path(raw_root, "suppl", "10x")
raw_tar <- if (nzchar(opt$`raw-tar`)) opt$`raw-tar` else file.path(h5_dir, "GSE286235_RAW.tar")
outdir <- if (nzchar(opt$outdir)) opt$outdir else file.path(project_root, "results/liu_2025_hnbmo")
plot_dir <- file.path(outdir, "plots")
table_dir <- file.path(outdir, "tables")

include_ds <- as_bool(opt$`include-ds`)
min_features <- as.integer(opt$`min-features`)
max_percent_mt <- as.numeric(opt$`max-percent-mt`)
gene_min_cells <- as.integer(opt$`gene-min-cells`)
nfeatures <- as.integer(opt$nfeatures)
dims_n <- as.integer(opt$dims)
resolution <- as.numeric(opt$resolution)
integrate_samples <- as_bool(opt$integrate)
seed <- as.integer(opt$seed)

if (!dir.exists(h5_dir)) dir.create(h5_dir, recursive = TRUE, showWarnings = FALSE)
if (!length(list.files(h5_dir, pattern = "raw_feature_bc_matrix[.]h5$", full.names = TRUE))) {
  if (!file.exists(raw_tar)) {
    stop("No H5 files found in ", h5_dir, " and raw TAR does not exist: ", raw_tar)
  }
  message("Extracting H5 files from: ", raw_tar)
  utils::untar(raw_tar, exdir = h5_dir)
}

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(seed)

h5_files <- list.files(h5_dir, pattern = "^GSM872144[0-3].*raw_feature_bc_matrix[.]h5$",
                       full.names = TRUE)
if (!length(h5_files)) stop("No GSE286235 H5 files found in: ", h5_dir)

file_info <- data.frame(
  path = h5_files,
  sample_geo_accession = vapply(h5_files, extract_accession, character(1)),
  stringsAsFactors = FALSE
)
file_info <- merge(file_info, sample_info, by = "sample_geo_accession", all.x = TRUE)
if (any(is.na(file_info$sample_id))) {
  stop("Unrecognized H5 accession(s): ",
       paste(file_info$sample_geo_accession[is.na(file_info$sample_id)], collapse = ", "))
}
if (!include_ds) file_info <- file_info[file_info$disease_status == "healthy", , drop = FALSE]
file_info <- file_info[order(file_info$sample_geo_accession), , drop = FALSE]

write_tsv(file_info[, c("sample_geo_accession", "sample_id", "sample_title",
                        "disease_status", "day", "day_numeric", "cell_line", "path")],
          file.path(table_dir, "sample_manifest_used.tsv"))

write_tsv(
  data.frame(
    assumption = c(
      "GEO provides raw_feature_bc_matrix.h5 files, not filtered Cell Ranger outputs",
      "Primary Figure 1C-like run uses healthy H9 D36, H9 D63, and IMR90-4 D63 only",
      "Cell calling is explicit because filtered matrices were not deposited",
      "Downstream Seurat workflow follows guided-tutorial-style defaults",
      "Sample handling",
      "Cell-type labels are exploratory marker-score labels mapped at cluster level"
    ),
    value = c(
      "raw H5 converted to 10x matrix.mtx/features/barcodes directories before Seurat",
      "GSM8721440, GSM8721441, GSM8721442",
      paste0("nFeature_RNA >= ", min_features, " and percent.mt < ", max_percent_mt),
      paste0("NormalizeData; FindVariableFeatures nfeatures=", nfeatures,
             "; ScaleData no regression; PCA/Neighbors/Clusters/UMAP dims=1:", dims_n,
             "; resolution=", resolution),
      if (integrate_samples) {
        "Seurat anchor integration across H9 D36, H9 D63, and IMR90-4 D63 before PCA/UMAP/clustering"
      } else {
        "No integration; samples are merged directly before PCA/UMAP/clustering"
      },
      "RG-div, RGC, IPC, IM, GABA, CHN marker sets from paper text/results"
    ),
    stringsAsFactors = FALSE
  ),
  file.path(table_dir, "analysis_assumptions.tsv")
)

message("Building Liu 2025 hnbMO Seurat object")
message("H5 directory: ", h5_dir)
message("Converted 10x directory: ", tenx_dir)
message("Output directory: ", outdir)
message("Include DS sample: ", include_ds)
message("Integrate samples: ", integrate_samples)
message("Samples: ", paste(file_info$sample_id, collapse = ", "))

objects <- list()
qc_rows <- list()

for (i in seq_len(nrow(file_info))) {
  info <- file_info[i, , drop = FALSE]
  path <- info$path[[1]]
  mat <- read_10x_matrix(path, info$sample_id[[1]], tenx_dir)

  n_features_by_cell <- Matrix::colSums(mat > 0)
  n_counts_by_cell <- Matrix::colSums(mat)
  mito_features <- grepl("^MT-", rownames(mat))
  if (any(mito_features)) {
    mt_counts_by_cell <- Matrix::colSums(mat[mito_features, , drop = FALSE])
  } else {
    mt_counts_by_cell <- rep(0, ncol(mat))
  }
  percent_mt_by_cell <- ifelse(n_counts_by_cell > 0, 100 * mt_counts_by_cell / n_counts_by_cell, 0)
  keep_cells <- n_features_by_cell >= min_features & percent_mt_by_cell < max_percent_mt

  qc_rows[[info$sample_id[[1]]]] <- data.frame(
    sample_geo_accession = info$sample_geo_accession,
    sample_id = info$sample_id,
    disease_status = info$disease_status,
    day = info$day,
    cell_line = info$cell_line,
    n_barcodes_raw = ncol(mat),
    n_genes_raw = nrow(mat),
    n_cells_after_filter = sum(keep_cells),
    min_features = min_features,
    max_percent_mt = max_percent_mt,
    median_nFeature_filtered = median(n_features_by_cell[keep_cells]),
    median_nCount_filtered = median(n_counts_by_cell[keep_cells]),
    median_percent_mt_filtered = median(percent_mt_by_cell[keep_cells]),
    stringsAsFactors = FALSE
  )

  mat <- mat[, keep_cells, drop = FALSE]
  colnames(mat) <- paste(info$sample_id[[1]], colnames(mat), sep = "_")

  obj <- CreateSeuratObject(
    counts = mat,
    project = "Liu2025_hnbMO",
    min.cells = gene_min_cells,
    min.features = 0
  )
  obj$sample_geo_accession <- info$sample_geo_accession
  obj$sample_id <- info$sample_id
  obj$orig.ident <- info$sample_id
  obj$sample_title <- info$sample_title
  obj$disease_status <- info$disease_status
  obj$day <- info$day
  obj$day_numeric <- info$day_numeric
  obj$cell_line <- info$cell_line
  obj$nFeature_RNA_prefilter <- as.numeric(n_features_by_cell[keep_cells])
  obj$nCount_RNA_prefilter <- as.numeric(n_counts_by_cell[keep_cells])
  obj$percent.mt.prefilter <- as.numeric(percent_mt_by_cell[keep_cells])
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")

  objects[[info$sample_id[[1]]]] <- obj
  rm(mat, obj, n_features_by_cell, n_counts_by_cell, mt_counts_by_cell,
     percent_mt_by_cell, keep_cells)
  gc()
}

qc_summary <- do.call(rbind, qc_rows)
write_tsv(qc_summary, file.path(table_dir, "qc_summary_by_sample.tsv"))

if (integrate_samples && length(objects) > 1L) {
  message("Running Seurat anchor integration across samples")
  objects <- lapply(objects, function(obj) {
    obj <- NormalizeData(obj)
    obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = nfeatures)
    obj
  })
  anchors <- FindIntegrationAnchors(object.list = objects, dims = seq_len(dims_n))
  seu <- IntegrateData(anchorset = anchors, dims = seq_len(dims_n))
  rm(anchors); gc()
} else {
  if (length(objects) == 1L) {
    seu <- objects[[1L]]
  } else {
    seu <- Reduce(function(x, y) merge(x, y), objects)
  }
}
rm(objects); gc()

if (!integrate_samples && exists("JoinLayers", mode = "function")) {
  message("Joining Seurat v5 assay layers after sample merge")
  seu <- JoinLayers(seu)
}

message("Merged cells: ", ncol(seu), "; genes: ", nrow(seu))

message("Normalization / HVG / scaling / PCA / UMAP / clustering")
if (integrate_samples) {
  DefaultAssay(seu) <- "RNA"
  if (exists("JoinLayers", mode = "function")) {
    message("Joining Seurat v5 RNA layers after integration")
    seu <- JoinLayers(seu, assay = "RNA")
  }
} else {
  seu <- NormalizeData(seu)
}
marker_scoring <- add_marker_scores(seu, cell_type_marker_sets)
seu <- marker_scoring$obj
score_cols <- marker_scoring$score_cols
marker_presence <- marker_scoring$marker_presence
if (integrate_samples) {
  DefaultAssay(seu) <- "integrated"
} else {
  seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = nfeatures)
}
seu <- ScaleData(seu)
seu <- RunPCA(seu, npcs = max(20L, dims_n), seed.use = seed)
seu <- FindNeighbors(seu, dims = seq_len(dims_n))
seu <- FindClusters(seu, resolution = resolution, random.seed = seed)
seu <- RunUMAP(seu, dims = seq_len(dims_n), seed.use = seed)

annotation <- annotate_clusters_by_marker_scores(seu, score_cols)
seu <- annotation$obj
seu$figure1c_sample <- factor(
  paste(seu$cell_line, seu$day, sep = " "),
  levels = c("H9 D36", "H9 D63", "IMR90-4 D63")
)
write_tsv(marker_presence, file.path(table_dir, "cell_type_marker_presence.tsv"))
write_tsv(annotation$cluster_score_summary,
          file.path(table_dir, "cell_type_marker_score_by_cluster.tsv"))

seurat_path <- if (integrate_samples) {
  file.path(outdir, "liu_2025_hnbmo_healthy_integrated_seurat.rds")
} else {
  file.path(outdir, "liu_2025_hnbmo_healthy_seurat.rds")
}
if (include_ds) {
  seurat_path <- if (integrate_samples) {
    file.path(outdir, "liu_2025_hnbmo_with_ds_integrated_seurat.rds")
  } else {
    file.path(outdir, "liu_2025_hnbmo_with_ds_seurat.rds")
  }
}
message("Saving Seurat object: ", seurat_path)
saveRDS(seu, seurat_path)

write_tsv(as.data.frame(table(seu$sample_id, seu$seurat_clusters)),
          file.path(table_dir, "cluster_counts_by_sample.tsv"))
write_tsv(as.data.frame(table(seu$day, seu$seurat_clusters)),
          file.path(table_dir, "cluster_counts_by_day.tsv"))

plot_dim <- function(group_by, filename, title, label = FALSE) {
  p <- DimPlot(seu, reduction = "umap", group.by = group_by, label = label) +
    ggtitle(title)
  ggsave(file.path(plot_dir, paste0(filename, ".png")), p, width = 8, height = 6, dpi = 300)
  ggsave(file.path(plot_dir, paste0(filename, ".pdf")), p, width = 8, height = 6)
}

plot_dim("seurat_clusters", "umap_by_cluster",
         paste0("Liu 2025 hnbMO UMAP by cluster (res=", resolution, ")"), TRUE)
plot_dim("liu2025_exploratory_cell_type", "figure1c_like_umap_by_exploratory_cell_type",
         "Figure 1C-like hnbMO UMAP by exploratory cell type", TRUE)
plot_dim("sample_id", "umap_by_sample", "Liu 2025 hnbMO UMAP by sample")
plot_dim("day", "umap_by_day", "Liu 2025 hnbMO UMAP by day")
plot_dim("cell_line", "umap_by_cell_line", "Liu 2025 hnbMO UMAP by cell line")
plot_dim("figure1c_sample", "figure1c_like_umap_by_h9_imr90_timepoint",
         "Figure 1C-like hnbMO UMAP: H9 D36, H9 D63, IMR90-4 D63")

marker_genes <- unique(c(
  "MKI67", "TOP2A", "SOX2", "VIM", "EOMES", "NEUROD1",
  "DCX", "STMN2", "GAD1", "GAD2", "CHAT", "SLC18A3",
  "FOXG1", "NKX2-1", "LHX8", "ISL1"
))
marker_genes <- marker_genes[marker_genes %in% rownames(seu)]
if (length(marker_genes)) {
  DefaultAssay(seu) <- "RNA"
  p_markers <- FeaturePlot(seu, features = marker_genes, reduction = "umap", ncol = 4)
  ggsave(file.path(plot_dir, "umap_marker_panel.png"), p_markers,
         width = 12, height = max(4, 3 * ceiling(length(marker_genes) / 4)), dpi = 300)
  ggsave(file.path(plot_dir, "umap_marker_panel.pdf"), p_markers,
         width = 12, height = max(4, 3 * ceiling(length(marker_genes) / 4)))
}

p_elbow <- ElbowPlot(seu, ndims = min(50L, max(50L, dims_n))) + ggtitle("PCA elbow")
ggsave(file.path(plot_dir, "pca_elbow.png"), p_elbow, width = 6, height = 4, dpi = 300)
ggsave(file.path(plot_dir, "pca_elbow.pdf"), p_elbow, width = 6, height = 4)

message("Done. Outputs in: ", outdir)
