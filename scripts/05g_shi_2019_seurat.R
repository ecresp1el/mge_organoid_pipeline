#!/usr/bin/env Rscript

# Build a standalone Seurat object + UMAP for Shi et al (GSE135827).
# Input count file is GEO supplementary table with format:
#   row 1: cell barcodes (with week suffix, e.g. -GW09)
#   rows 2+: gene symbol + raw counts per cell

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(Seurat)
  library(ggplot2)
})

parse_args <- function(args) {
  out <- list(
    `project-root` = Sys.getenv("PROJECT_ROOT"),
    counts = NULL,
    `series-matrix` = NULL,
    outdir = NULL,
    `min-features` = "500",
    `max-percent-mt` = "20",
    `nfeatures` = "3000",
    `dims` = "30",
    resolution = "0.8",
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
      "  Rscript scripts/05g_shi_2019_seurat.R --project-root <PROJECT_ROOT> [options]",
      "",
      "Options:",
      "  --counts <path>         Counts table (default: data/raw/shi_2019_geo_files/suppl/GSE135827_GE_mat_raw_count_with_week_info.txt.gz)",
      "  --series-matrix <path>  GEO series matrix for sample metadata (default: data/raw/shi_2019_geo_files/matrix/GSE135827_series_matrix.txt.gz)",
      "  --outdir <path>         Output directory (default: PROJECT_ROOT/results/shi_2019)",
      "  --min-features <int>    Cell filter lower bound on nFeature_RNA (default: 500)",
      "  --max-percent-mt <num>  Cell filter upper bound on percent.mt (default: 20)",
      "  --nfeatures <int>       HVGs for FindVariableFeatures (default: 3000)",
      "  --dims <int>            PCs used for neighbors/UMAP (default: 30)",
      "  --resolution <num>      Cluster resolution (default: 0.8)",
      "",
      "Environment:",
      "  PROJECT_ROOT            Alternative way to pass --project-root",
      sep = "\n"
    )
  )
}

trim_trailing_slash <- function(x) sub("/+$", "", x)

open_text_connection <- function(path) {
  if (grepl("\\.gz$", path, ignore.case = TRUE)) {
    gzfile(path, open = "rt")
  } else {
    file(path, open = "rt")
  }
}

strip_quotes <- function(x) gsub('^"|"$', "", x)

extract_tab_values <- function(line) {
  fields <- strsplit(line, "\t", fixed = TRUE)[[1]]
  values <- fields[-1]
  strip_quotes(values)
}

parse_series_matrix_samples <- function(series_path) {
  if (!file.exists(series_path)) return(data.table())

  con <- open_text_connection(series_path)
  on.exit(close(con), add = TRUE)
  lines <- readLines(con, warn = FALSE)

  extract_first <- function(prefix) {
    idx <- grep(paste0("^", prefix, "\\t"), lines)
    if (length(idx) == 0) return(character())
    extract_tab_values(lines[[idx[1]]])
  }

  extract_all <- function(prefix) {
    idx <- grep(paste0("^", prefix, "\\t"), lines)
    if (length(idx) == 0) return(list())
    lapply(idx, function(i) extract_tab_values(lines[[i]]))
  }

  geo <- extract_first("!Sample_geo_accession")
  if (length(geo) == 0) return(data.table())

  n <- length(geo)
  dt <- data.table(sample_geo_accession = geo)

  add_if_length <- function(colname, values) {
    if (length(values) == n) dt[[colname]] <<- values
  }

  add_if_length("sample_title", extract_first("!Sample_title"))
  add_if_length("sample_source_name_ch1", extract_first("!Sample_source_name_ch1"))
  add_if_length("sample_organism_ch1", extract_first("!Sample_organism_ch1"))

  char_sets <- extract_all("!Sample_characteristics_ch1")
  if (length(char_sets) > 0) {
    for (k in seq_along(char_sets)) {
      vals <- char_sets[[k]]
      if (length(vals) == n) {
        dt[[paste0("sample_characteristics_ch1_", k)]] <- vals
      }
    }
  }

  # Extract week labels from any characteristics column if present.
  week_from_text <- function(x) {
    m <- regexpr("GW[0-9]+(?:-[0-9]+)?", x, perl = TRUE)
    out <- rep(NA_character_, length(x))
    hit <- m > 0
    out[hit] <- regmatches(x, m)
    out
  }

  char_cols <- grep("^sample_characteristics_ch1_", names(dt), value = TRUE)
  if (length(char_cols) > 0) {
    week_mat <- lapply(char_cols, function(cc) week_from_text(dt[[cc]]))
    week_mat <- do.call(cbind, week_mat)
    dt[, sample_week_label := NA_character_]
    for (j in seq_len(ncol(week_mat))) {
      fill <- is.na(dt$sample_week_label) & !is.na(week_mat[, j])
      dt$sample_week_label[fill] <- week_mat[, j][fill]
    }
  }

  dt
}

summarize_metadata <- function(meta_df) {
  cols <- names(meta_df)
  out <- lapply(cols, function(cc) {
    x <- meta_df[[cc]]
    as.character_x <- as.character(x)
    non_na <- sum(!is.na(x))
    uniq <- length(unique(as.character_x[!is.na(as.character_x)]))
    example <- paste(head(unique(as.character_x[!is.na(as.character_x)]), 3), collapse = " | ")
    data.table(
      column = cc,
      class = paste(class(x), collapse = ";"),
      n_non_na = non_na,
      n_unique = uniq,
      example_values = example
    )
  })
  rbindlist(out, fill = TRUE)
}

opt <- parse_args(commandArgs(trailingOnly = TRUE))
if (isTRUE(opt$help)) {
  print_usage()
  quit(save = "no", status = 0)
}

if (!nzchar(opt$`project-root`)) {
  stop("PROJECT_ROOT or --project-root is required")
}
project_root <- trim_trailing_slash(opt$`project-root`)

counts_path <- opt$counts
if (is.null(counts_path) || !nzchar(counts_path)) {
  counts_path <- file.path(project_root, "data/raw/shi_2019_geo_files/suppl/GSE135827_GE_mat_raw_count_with_week_info.txt.gz")
}
if (!file.exists(counts_path)) stop("Counts file not found: ", counts_path)

series_matrix_path <- opt$`series-matrix`
if (is.null(series_matrix_path) || !nzchar(series_matrix_path)) {
  series_matrix_path <- file.path(project_root, "data/raw/shi_2019_geo_files/matrix/GSE135827_series_matrix.txt.gz")
}

outdir <- if (is.null(opt$outdir) || !nzchar(opt$outdir)) {
  file.path(project_root, "results/shi_2019")
} else {
  opt$outdir
}
plot_dir <- file.path(outdir, "plots")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

min_features <- as.integer(opt$`min-features`)
max_percent_mt <- as.numeric(opt$`max-percent-mt`)
nfeatures <- as.integer(opt$nfeatures)
dims_use <- as.integer(opt$dims)
resolution <- as.numeric(opt$resolution)
if (is.na(min_features) || min_features < 0) stop("Invalid --min-features")
if (is.na(max_percent_mt) || max_percent_mt <= 0) stop("Invalid --max-percent-mt")
if (is.na(nfeatures) || nfeatures <= 0) stop("Invalid --nfeatures")
if (is.na(dims_use) || dims_use <= 1) stop("Invalid --dims")
if (is.na(resolution) || resolution <= 0) stop("Invalid --resolution")

message("Reading cell header from counts file...")
con <- open_text_connection(counts_path)
header_line <- readLines(con, n = 1, warn = FALSE)
close(con)
if (length(header_line) != 1 || !nzchar(header_line)) {
  stop("Failed to read counts header from: ", counts_path)
}
cell_ids <- strsplit(header_line, "\t", fixed = TRUE)[[1]]
if (length(cell_ids) < 10) {
  stop("Header parse failed: expected many cell IDs, got ", length(cell_ids))
}
raw_cell_ids <- cell_ids
dup_cell_count <- sum(duplicated(raw_cell_ids))
if (dup_cell_count > 0) {
  message("Detected ", dup_cell_count, " duplicated cell IDs in header; making unique names for Seurat compatibility.")
  cell_ids <- make.unique(raw_cell_ids, sep = "_dup")
}

message("Reading count table body and converting to sparse matrix...")
read_cmd <- if (grepl("\\.gz$", counts_path, ignore.case = TRUE)) {
  paste("gzip -dc", shQuote(counts_path), "| tail -n +2")
} else {
  paste("tail -n +2", shQuote(counts_path))
}

counts_dt <- fread(
  cmd = read_cmd,
  sep = "\t",
  header = FALSE,
  data.table = TRUE,
  showProgress = TRUE
)

if (ncol(counts_dt) < 2) {
  stop("Counts body parse failed: less than 2 columns")
}
if ((ncol(counts_dt) - 1L) != length(cell_ids)) {
  stop(
    "Cell column mismatch: header has ", length(cell_ids),
    " cells but body has ", ncol(counts_dt) - 1L, " count columns"
  )
}

genes <- make.unique(as.character(counts_dt[[1]]))
counts_dt[[1]] <- NULL

# Convert the dense table to matrix then sparse dgCMatrix for Seurat.
dense <- as.matrix(counts_dt)
rm(counts_dt)
gc()

storage.mode(dense) <- "double"
spmat <- Matrix(dense, sparse = TRUE)
rm(dense)
gc()

rownames(spmat) <- genes
colnames(spmat) <- cell_ids

message("Creating Seurat object...")
seu <- CreateSeuratObject(
  counts = spmat,
  project = "Shi2019",
  min.cells = 3,
  min.features = 0
)
rm(spmat)
gc()

# Derive basic cell metadata from barcode suffix.
cell_meta <- data.table(
  cell_id = colnames(seu),
  raw_cell_id = raw_cell_ids
)
cell_meta[, study_id := "shi_2019"]
cell_meta[, has_duplicated_raw_cell_id := duplicated(raw_cell_id) | duplicated(raw_cell_id, fromLast = TRUE)]
cell_meta[, week_label := ifelse(grepl("-GW[0-9]+(?:[-_][0-9A-Za-z]+)?$", raw_cell_id, perl = TRUE),
                                 sub("^.*-(GW[0-9]+(?:[-_][0-9A-Za-z]+)?)$", "\\1", raw_cell_id, perl = TRUE),
                                 NA_character_)]
cell_meta[, week_numeric := suppressWarnings(as.integer(sub("^GW([0-9]+).*$", "\\1", week_label, perl = TRUE)))]
cell_meta[, barcode := sub("-GW[0-9]+(?:[-_][0-9A-Za-z]+)?$", "", raw_cell_id, perl = TRUE)]

rownames(cell_meta) <- cell_meta$cell_id
seu <- AddMetaData(
  seu,
  metadata = as.data.frame(cell_meta[, .(
    raw_cell_id,
    has_duplicated_raw_cell_id,
    study_id,
    week_label,
    week_numeric,
    barcode
  )])
)

# Optional sample-level metadata from GEO series matrix.
sample_meta <- parse_series_matrix_samples(series_matrix_path)
if (nrow(sample_meta) > 0) {
  fwrite(sample_meta, file.path(outdir, "shi_2019_sample_metadata_from_series_matrix.tsv"), sep = "\t")
}

seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")

message("Applying QC filter: nFeature_RNA > ", min_features, ", percent.mt < ", max_percent_mt)
seu <- subset(seu, subset = nFeature_RNA > min_features & percent.mt < max_percent_mt)

message("Normalize / HVG / Scale / PCA / Neighbors / Clusters / UMAP...")
seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = nfeatures)
seu <- ScaleData(seu, vars.to.regress = "percent.mt")
seu <- RunPCA(seu, npcs = max(50L, dims_use))
seu <- FindNeighbors(seu, dims = 1:dims_use)
seu <- FindClusters(seu, resolution = resolution)
seu <- RunUMAP(seu, dims = 1:dims_use)

seu_path <- file.path(outdir, "shi_2019_seurat.rds")
saveRDS(seu, seu_path)

# Save metadata exports for easy downstream inspection.
cell_meta_out <- as.data.table(seu@meta.data, keep.rownames = "cell_id")
fwrite(cell_meta_out, file.path(outdir, "shi_2019_cell_metadata.tsv.gz"), sep = "\t")

meta_summary <- summarize_metadata(seu@meta.data)
fwrite(meta_summary, file.path(outdir, "shi_2019_metadata_columns_summary.tsv"), sep = "\t")

summary_dt <- data.table(
  metric = c("n_cells", "n_features", "default_assay", "reductions", "clusters", "weeks_detected"),
  value = c(
    as.character(ncol(seu)),
    as.character(nrow(seu)),
    DefaultAssay(seu),
    paste(Reductions(seu), collapse = ","),
    as.character(length(unique(Idents(seu)))),
    paste(sort(unique(na.omit(seu$week_label))), collapse = ",")
  )
)
fwrite(summary_dt, file.path(outdir, "shi_2019_summary.tsv"), sep = "\t")

# UMAP and elbow plots.
p_umap <- DimPlot(seu, reduction = "umap", group.by = "seurat_clusters", label = TRUE) +
  ggtitle("Shi et al 2019 (GSE135827) UMAP by cluster")
ggsave(file.path(plot_dir, "umap_by_cluster.png"), p_umap, width = 8, height = 6, dpi = 300)
ggsave(file.path(plot_dir, "umap_by_cluster.pdf"), p_umap, width = 8, height = 6)

p_week <- DimPlot(seu, reduction = "umap", group.by = "week_label") +
  ggtitle("Shi et al 2019 UMAP by week label")
ggsave(file.path(plot_dir, "umap_by_week.png"), p_week, width = 8, height = 6, dpi = 300)
ggsave(file.path(plot_dir, "umap_by_week.pdf"), p_week, width = 8, height = 6)

p_elbow <- ElbowPlot(seu, ndims = max(50L, dims_use)) + ggtitle("Shi 2019 PCA elbow")
ggsave(file.path(plot_dir, "pca_elbow.png"), p_elbow, width = 6, height = 4, dpi = 300)

message("Done. Wrote Seurat object: ", seu_path)
