#!/usr/bin/env Rscript

# Build a minimal Seurat object for He et al from the extracted SCN8A slice.
# This is intentionally scoped for SCN8A cross-study plotting compatibility.

suppressPackageStartupMessages({
  library(Matrix)
  library(Seurat)
  library(data.table)
  library(ggplot2)
})

parse_args <- function(args) {
  out <- list(
    `project-root` = Sys.getenv("PROJECT_ROOT"),
    `slice-dir` = NULL,
    outdir = NULL,
    help = FALSE
  )
  i <- 1L
  while (i <= length(args)) {
    a <- args[[i]]
    if (a == "--help" || a == "-h") {
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
      "  Rscript scripts/05f_he_et_al_scn8a_seurat.R --project-root <PROJECT_ROOT> [--slice-dir <dir>] [--outdir <dir>]",
      "",
      "Defaults:",
      "  --project-root defaults to PROJECT_ROOT env var",
      "  --slice-dir defaults to <PROJECT_ROOT>/data/processed/he_et_al_scn8a_slice",
      "  --outdir defaults to <PROJECT_ROOT>/results/he_et_al",
      sep = "\n"
    )
  )
}

opt <- parse_args(commandArgs(trailingOnly = TRUE))
if (isTRUE(opt$help)) {
  print_usage()
  quit(save = "no", status = 0)
}

trim_trailing_slash <- function(x) sub("/+$", "", x)

fread_any <- function(path, ...) {
  if (grepl("\\.gz$", path, ignore.case = TRUE)) {
    data.table::fread(cmd = paste("gzip -dc", shQuote(path)), ...)
  } else {
    data.table::fread(path, ...)
  }
}

if (!nzchar(opt$`project-root`)) {
  stop("PROJECT_ROOT or --project-root is required")
}
project_root <- trim_trailing_slash(opt$`project-root`)

slice_dir <- opt$`slice-dir`
if (is.null(slice_dir) || !nzchar(slice_dir)) {
  slice_dir <- file.path(project_root, "data/processed/he_et_al_scn8a_slice")
}
outdir <- if (is.null(opt$outdir) || !nzchar(opt$outdir)) {
  file.path(project_root, "results/he_et_al")
} else {
  opt$outdir
}
plot_dir <- file.path(outdir, "plots")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

counts_path <- file.path(slice_dir, "scn8a_counts.mtx")
genes_path <- file.path(slice_dir, "genes.tsv")
barcodes_path <- file.path(slice_dir, "barcodes.tsv")
umap_path <- file.path(slice_dir, "umap.tsv.gz")
obs_meta_path <- file.path(slice_dir, "obs_metadata.tsv.gz")

needed <- c(counts_path, genes_path, barcodes_path, umap_path)
missing <- needed[!file.exists(needed)]
if (length(missing) > 0) {
  stop("Missing slice input files: ", paste(missing, collapse = ", "))
}

message("Reading SCN8A matrix slice...")
counts <- readMM(counts_path)
if (!inherits(counts, "dgCMatrix")) {
  counts <- as(counts, "dgCMatrix")
}

features <- fread_any(genes_path, header = FALSE)$V1
cells <- fread_any(barcodes_path, header = FALSE)$V1
if (length(features) != nrow(counts)) {
  stop("Feature count mismatch: genes.tsv has ", length(features), " rows, matrix has ", nrow(counts))
}
if (length(cells) != ncol(counts)) {
  stop("Cell count mismatch: barcodes.tsv has ", length(cells), " rows, matrix has ", ncol(counts))
}

rownames(counts) <- make.unique(as.character(features))
colnames(counts) <- as.character(cells)

message("Creating Seurat object...")
seu <- CreateSeuratObject(
  counts = counts,
  project = "HeEtAl",
  min.cells = 0,
  min.features = 0
)

# Keep object aligned with existing scripts that may request data slot first.
seu <- NormalizeData(seu, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)

message("Reading UMAP coordinates...")
umap_dt <- fread_any(umap_path)
required_cols <- c("cell_id", "UMAP_1", "UMAP_2")
if (!all(required_cols %in% names(umap_dt))) {
  stop("UMAP file missing required columns: ", paste(setdiff(required_cols, names(umap_dt)), collapse = ","))
}

umap_cells <- as.character(umap_dt$cell_id)
common_cells <- intersect(colnames(seu), umap_cells)
if (length(common_cells) == 0) {
  stop("No overlapping cells between Seurat object and UMAP coordinates")
}

if (length(common_cells) < ncol(seu)) {
  message("Subsetting to cells with UMAP coordinates: ", length(common_cells), " / ", ncol(seu))
  seu <- subset(seu, cells = common_cells)
}

umap_idx <- match(colnames(seu), umap_dt$cell_id)
emb <- as.matrix(umap_dt[umap_idx, .(UMAP_1, UMAP_2)])
rownames(emb) <- colnames(seu)
colnames(emb) <- c("UMAP_1", "UMAP_2")

seu[["umap"]] <- CreateDimReducObject(
  embeddings = emb,
  key = "UMAP_",
  assay = DefaultAssay(seu)
)

# Attach selected metadata when available.
if (file.exists(obs_meta_path)) {
  message("Reading extracted obs metadata...")
  obs_dt <- fread_any(obs_meta_path)
  if ("cell_id" %in% names(obs_dt)) {
    obs_dt <- unique(obs_dt, by = "cell_id")
    meta_idx <- match(colnames(seu), obs_dt$cell_id)
    meta_use <- as.data.frame(obs_dt[meta_idx])
    rownames(meta_use) <- meta_use$cell_id
    keep_cols <- setdiff(colnames(meta_use), "cell_id")
    if (length(keep_cols) > 0) {
      seu <- AddMetaData(seu, metadata = meta_use[, keep_cols, drop = FALSE])
    }
  }
}

seu_path <- file.path(outdir, "he_et_al_scn8a_seurat.rds")
saveRDS(seu, seu_path)

# Lightweight QC table and reference plot for sanity checks.
summary_path <- file.path(outdir, "he_et_al_scn8a_summary.tsv")
summary_dt <- data.table(
  metric = c("n_cells", "n_features", "default_assay", "reductions"),
  value = c(
    as.character(ncol(seu)),
    as.character(nrow(seu)),
    DefaultAssay(seu),
    paste(Reductions(seu), collapse = ",")
  )
)
fwrite(summary_dt, summary_path, sep = "\t")

if ("SCN8A" %in% rownames(seu)) {
  p <- FeaturePlot(seu, features = "SCN8A", reduction = "umap", pt.size = 0.1) +
    ggtitle("He et al (full V2): SCN8A on UMAP")
  ggsave(file.path(plot_dir, "scn8a_umap.png"), p, width = 8, height = 6, dpi = 300)
  ggsave(file.path(plot_dir, "scn8a_umap.pdf"), p, width = 8, height = 6)
}

message("Done. Wrote: ", seu_path)
