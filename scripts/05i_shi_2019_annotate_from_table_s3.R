#!/usr/bin/env Rscript

# Annotate Shi et al clusters with paper major cell types from Table S3 DEGs.
# This step assumes a pre-built Shi Seurat object exists and reuses its UMAP.

suppressPackageStartupMessages({
  library(data.table)
  library(Seurat)
  library(ggplot2)
})

parse_args <- function(args) {
  out <- list(
    `project-root` = Sys.getenv("PROJECT_ROOT"),
    `seurat-rds` = NULL,
    `table-s3-tsv` = NULL,
    outdir = NULL,
    resolution = "0.11",
    `top-n-our` = "200",
    `top-n-paper` = "200",
    `min-pct` = "0.2",
    `logfc-threshold` = "0.25",
    `gene-column` = "",
    `celltype-column` = "",
    `paper-positive-only` = "true",
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
      "  Rscript scripts/05i_shi_2019_annotate_from_table_s3.R --project-root <PROJECT_ROOT> [options]",
      "",
      "Options:",
      "  --seurat-rds <path>      Shi Seurat object (default: PROJECT_ROOT/results/shi_2019_paper_qc/shi_2019_seurat.rds)",
      "  --table-s3-tsv <path>    Table S3 TSV converted from xlsx",
      "                           (default: PROJECT_ROOT/data/raw/shi_2019_geo_files/supplementary/science.abj6641_table_s3.tsv)",
      "  --outdir <path>          Output directory (default: PROJECT_ROOT/results/shi_2019_paper_qc/table_s3_annotation)",
      "  --resolution <num>       Clustering resolution for annotation mapping (default: 0.11)",
      "  --top-n-our <int>        Top marker genes per our cluster (default: 200)",
      "  --top-n-paper <int>      Top marker genes per paper cell type (default: 200)",
      "  --min-pct <num>          FindAllMarkers min.pct (default: 0.2)",
      "  --logfc-threshold <num>  FindAllMarkers logfc.threshold (default: 0.25)",
      "  --gene-column <name>     Optional explicit gene column in Table S3 TSV",
      "  --celltype-column <name> Optional explicit cell-type column in Table S3 TSV",
      "  --paper-positive-only <true|false>  Keep only positive paper logFC rows if logFC column exists (default: true)",
      "",
      "Environment:",
      "  PROJECT_ROOT             Alternative way to pass --project-root",
      sep = "\n"
    )
  )
}

trim_trailing_slash <- function(x) sub("/+$", "", x)

to_bool <- function(x, default = FALSE) {
  if (is.null(x) || !nzchar(x)) return(default)
  y <- tolower(trimws(as.character(x)))
  y %in% c("1", "true", "t", "yes", "y")
}

pick_column <- function(nms, include_patterns, exclude = character()) {
  nms_l <- tolower(nms)
  ex_l <- tolower(exclude)
  for (pat in include_patterns) {
    hits <- grep(pat, nms_l, perl = TRUE)
    if (length(hits) > 0) {
      for (h in hits) {
        if (!(nms[h] %in% exclude) && !(nms_l[h] %in% ex_l)) return(nms[h])
      }
    }
  }
  NULL
}

fc_column_name <- function(nms) {
  pick_column(
    nms,
    include_patterns = c("avg[_\\.]*log2*fc", "log2*[_\\.]*fc", "^logfc$", "avgfc")
  )
}

opt <- parse_args(commandArgs(trailingOnly = TRUE))
if (isTRUE(opt$help)) {
  print_usage()
  quit(save = "no", status = 0)
}

if (!nzchar(opt$`project-root`)) stop("PROJECT_ROOT or --project-root is required")
project_root <- trim_trailing_slash(opt$`project-root`)

seurat_rds <- opt$`seurat-rds`
if (is.null(seurat_rds) || !nzchar(seurat_rds)) {
  seurat_rds <- file.path(project_root, "results/shi_2019_paper_qc/shi_2019_seurat.rds")
}
if (!file.exists(seurat_rds)) stop("Seurat RDS not found: ", seurat_rds)

table_s3_tsv <- opt$`table-s3-tsv`
if (is.null(table_s3_tsv) || !nzchar(table_s3_tsv)) {
  table_s3_tsv <- file.path(project_root, "data/raw/shi_2019_geo_files/supplementary/science.abj6641_table_s3.tsv")
}
if (!file.exists(table_s3_tsv)) stop("Table S3 TSV not found: ", table_s3_tsv)

outdir <- if (is.null(opt$outdir) || !nzchar(opt$outdir)) {
  file.path(project_root, "results/shi_2019_paper_qc/table_s3_annotation")
} else {
  opt$outdir
}
plot_dir <- file.path(outdir, "plots")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

resolution <- as.numeric(opt$resolution)
top_n_our <- as.integer(opt$`top-n-our`)
top_n_paper <- as.integer(opt$`top-n-paper`)
min_pct <- as.numeric(opt$`min-pct`)
logfc_threshold <- as.numeric(opt$`logfc-threshold`)
paper_positive_only <- to_bool(opt$`paper-positive-only`, default = TRUE)

if (is.na(resolution) || resolution <= 0) stop("Invalid --resolution")
if (is.na(top_n_our) || top_n_our < 5) stop("Invalid --top-n-our")
if (is.na(top_n_paper) || top_n_paper < 5) stop("Invalid --top-n-paper")
if (is.na(min_pct) || min_pct <= 0 || min_pct > 1) stop("Invalid --min-pct")
if (is.na(logfc_threshold) || logfc_threshold < 0) stop("Invalid --logfc-threshold")

message("Loading Shi Seurat object: ", seurat_rds)
seu <- readRDS(seurat_rds)
if (!("umap" %in% names(seu@reductions))) {
  stop("UMAP reduction not found in Seurat object. Expected reduction name: 'umap'.")
}

message("Loading Table S3 TSV: ", table_s3_tsv)
paper_dt <- fread(table_s3_tsv, sep = "\t", header = TRUE, fill = TRUE, showProgress = FALSE)
if (nrow(paper_dt) == 0 || ncol(paper_dt) < 2) {
  stop("Table S3 TSV appears empty or malformed: ", table_s3_tsv)
}

raw_cols <- names(paper_dt)
gene_col <- if (nzchar(opt$`gene-column`)) opt$`gene-column` else {
  pick_column(
    raw_cols,
    include_patterns = c("^gene$", "gene[_\\.]*symbol", "gene", "feature", "^symbol$")
  )
}
if (is.null(gene_col) || !(gene_col %in% raw_cols)) {
  stop("Could not detect gene column in Table S3 TSV. Pass --gene-column explicitly. Columns: ",
       paste(raw_cols, collapse = ", "))
}

celltype_col <- if (nzchar(opt$`celltype-column`)) opt$`celltype-column` else {
  pick_column(
    raw_cols,
    include_patterns = c("cell[_\\.]*type", "major[_\\.]*cell", "cluster", "identity", "group"),
    exclude = gene_col
  )
}
if (is.null(celltype_col) || !(celltype_col %in% raw_cols)) {
  stop("Could not detect cell-type column in Table S3 TSV. Pass --celltype-column explicitly. Columns: ",
       paste(raw_cols, collapse = ", "))
}

paper_fc_col <- fc_column_name(raw_cols)

paper_dt[, paper_cell_type := trimws(as.character(get(celltype_col)))]
paper_dt[, gene_symbol := toupper(trimws(as.character(get(gene_col))))]
paper_dt <- paper_dt[nzchar(paper_cell_type) & nzchar(gene_symbol)]

if (!is.null(paper_fc_col) && paper_positive_only) {
  suppressWarnings({
    paper_dt[, .paper_fc := as.numeric(get(paper_fc_col))]
  })
  paper_dt <- paper_dt[is.na(.paper_fc) | .paper_fc > 0]
}

if (!is.null(paper_fc_col)) {
  suppressWarnings({
    paper_dt[, .paper_fc := as.numeric(get(paper_fc_col))]
  })
  setorder(paper_dt, paper_cell_type, -.paper_fc)
  paper_sig <- paper_dt[, head(unique(gene_symbol), top_n_paper), by = paper_cell_type]
} else {
  paper_sig <- paper_dt[, head(unique(gene_symbol), top_n_paper), by = paper_cell_type]
}
setnames(paper_sig, "V1", "gene_symbol")
if (nrow(paper_sig) == 0) {
  stop("No paper signature rows available after filtering Table S3 input.")
}

paper_sig_sizes <- paper_sig[, .(paper_signature_size = uniqueN(gene_symbol)), by = paper_cell_type]

message("Reclustering identities at resolution ", resolution, " (reuses existing graph/UMAP)")
if (!("RNA_snn" %in% names(seu@graphs))) {
  message("RNA_snn graph missing; rebuilding neighbors from PCA dims 1:30")
  if (!("pca" %in% names(seu@reductions))) {
    stop("PCA reduction missing; cannot rebuild neighbors for reclustering.")
  }
  seu <- FindNeighbors(seu, dims = 1:30, verbose = FALSE)
}
seu <- FindClusters(seu, resolution = resolution, verbose = FALSE)
cluster_col <- paste0("shi_cluster_res_", gsub("\\.", "_", formatC(resolution, format = "fg", digits = 6, drop0trailing = TRUE)))
seu[[cluster_col]] <- as.character(Idents(seu))
Idents(seu) <- cluster_col

message("Computing our cluster markers with FindAllMarkers...")
our_markers <- FindAllMarkers(
  seu,
  only.pos = TRUE,
  min.pct = min_pct,
  logfc.threshold = logfc_threshold,
  verbose = FALSE
)
if (nrow(our_markers) == 0) stop("FindAllMarkers returned zero rows")

our_fc_col <- fc_column_name(names(our_markers))
if (is.null(our_fc_col)) {
  stop("Could not detect logFC column in FindAllMarkers output. Columns: ",
       paste(names(our_markers), collapse = ", "))
}

our_dt <- as.data.table(our_markers)
our_dt[, cluster := as.character(cluster)]
our_dt[, gene_symbol := toupper(trimws(as.character(gene)))]
suppressWarnings({
  our_dt[, .our_fc := as.numeric(get(our_fc_col))]
})
setorder(our_dt, cluster, -.our_fc)
our_sig <- our_dt[, head(unique(gene_symbol), top_n_our), by = cluster]
setnames(our_sig, "V1", "gene_symbol")
if (nrow(our_sig) == 0) {
  stop("No cluster signatures available from FindAllMarkers output.")
}

our_sig_sizes <- our_sig[, .(our_signature_size = uniqueN(gene_symbol)), by = cluster]

message("Scoring overlap between our clusters and paper cell-type signatures...")
our_list <- split(our_sig$gene_symbol, our_sig$cluster)
paper_list <- split(paper_sig$gene_symbol, paper_sig$paper_cell_type)

score_rows <- list()
idx <- 1L
for (cl in names(our_list)) {
  g1 <- unique(our_list[[cl]])
  for (ct in names(paper_list)) {
    g2 <- unique(paper_list[[ct]])
    ov <- intersect(g1, g2)
    un <- union(g1, g2)
    j <- if (length(un) > 0) length(ov) / length(un) else 0
    score_rows[[idx]] <- data.table(
      cluster = cl,
      paper_cell_type = ct,
      overlap = length(ov),
      jaccard = j,
      overlap_genes = paste(head(ov, 30), collapse = ",")
    )
    idx <- idx + 1L
  }
}
score_dt <- rbindlist(score_rows, fill = TRUE)

setorder(score_dt, cluster, -jaccard, -overlap, paper_cell_type)
cluster_map <- score_dt[, .SD[1], by = cluster]
cluster_map <- merge(cluster_map, our_sig_sizes, by = "cluster", all.x = TRUE)
cluster_map <- merge(cluster_map, paper_sig_sizes, by = "paper_cell_type", all.x = TRUE)

map_type <- setNames(cluster_map$paper_cell_type, cluster_map$cluster)
map_j <- setNames(cluster_map$jaccard, cluster_map$cluster)
map_ov <- setNames(cluster_map$overlap, cluster_map$cluster)
cluster_ids <- as.character(Idents(seu))
seu$paper_cell_type_table_s3 <- map_type[cluster_ids]
seu$paper_cell_type_table_s3_jaccard <- as.numeric(map_j[cluster_ids])
seu$paper_cell_type_table_s3_overlap <- as.integer(map_ov[cluster_ids])

res_tag <- gsub("\\.", "_", formatC(resolution, format = "fg", digits = 6, drop0trailing = TRUE))

score_path <- file.path(outdir, paste0("cluster_vs_paper_scores_res", res_tag, ".tsv"))
map_path <- file.path(outdir, paste0("cluster_to_paper_celltype_res", res_tag, ".tsv"))
paper_sig_path <- file.path(outdir, "paper_signature_genes_used.tsv")
our_sig_path <- file.path(outdir, paste0("our_signature_genes_used_res", res_tag, ".tsv"))
summary_path <- file.path(outdir, paste0("annotation_summary_res", res_tag, ".tsv"))
rds_path <- file.path(outdir, paste0("shi_2019_seurat_annotated_table_s3_res", res_tag, ".rds"))

fwrite(score_dt, score_path, sep = "\t")
fwrite(cluster_map, map_path, sep = "\t")
fwrite(paper_sig, paper_sig_path, sep = "\t")
fwrite(our_sig, our_sig_path, sep = "\t")

summary_dt <- data.table(
  metric = c(
    "seurat_rds_input",
    "table_s3_tsv_input",
    "resolution",
    "cluster_column",
    "n_cells",
    "n_clusters",
    "paper_celltype_column",
    "paper_gene_column",
    "paper_logfc_column",
    "paper_positive_only",
    "top_n_our",
    "top_n_paper",
    "min_pct",
    "logfc_threshold",
    "n_paper_cell_types",
    "n_cluster_map_rows"
  ),
  value = c(
    seurat_rds,
    table_s3_tsv,
    as.character(resolution),
    cluster_col,
    as.character(ncol(seu)),
    as.character(length(unique(Idents(seu)))),
    celltype_col,
    gene_col,
    ifelse(is.null(paper_fc_col), "NA", paper_fc_col),
    as.character(paper_positive_only),
    as.character(top_n_our),
    as.character(top_n_paper),
    as.character(min_pct),
    as.character(logfc_threshold),
    as.character(uniqueN(paper_sig$paper_cell_type)),
    as.character(nrow(cluster_map))
  )
)
fwrite(summary_dt, summary_path, sep = "\t")

message("Saving annotated object: ", rds_path)
saveRDS(seu, rds_path)

message("Writing UMAP plots...")
p_cluster <- DimPlot(seu, reduction = "umap", group.by = cluster_col, label = TRUE) +
  ggtitle(sprintf("Shi 2019 UMAP by cluster (res=%s)", as.character(resolution)))
p_type <- DimPlot(seu, reduction = "umap", group.by = "paper_cell_type_table_s3", label = TRUE) +
  ggtitle("Shi 2019 UMAP by Table S3 mapped cell type")

ggsave(file.path(plot_dir, paste0("umap_by_cluster_res", res_tag, ".png")), p_cluster, width = 8, height = 6, dpi = 300)
ggsave(file.path(plot_dir, paste0("umap_by_cluster_res", res_tag, ".pdf")), p_cluster, width = 8, height = 6)
ggsave(file.path(plot_dir, paste0("umap_by_table_s3_celltype_res", res_tag, ".png")), p_type, width = 9, height = 6, dpi = 300)
ggsave(file.path(plot_dir, paste0("umap_by_table_s3_celltype_res", res_tag, ".pdf")), p_type, width = 9, height = 6)

message("Done.")
message("Mapping table: ", map_path)
message("Annotated object: ", rds_path)
