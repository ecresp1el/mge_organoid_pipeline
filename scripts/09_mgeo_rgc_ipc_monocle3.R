#!/usr/bin/env Rscript

# Run Monocle3 pseudotime on the focused MGEO RGC/IPC export.

log_msg <- function(...) {
  msg <- paste0(..., collapse = " ")
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), msg))
  flush.console()
}

parse_args <- function(args) {
  out <- list(
    `project-root` = Sys.getenv("PROJECT_ROOT"),
    `input-dir` = NULL,
    outdir = NULL,
    `num-dim` = "50",
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
      "  Rscript scripts/09_mgeo_rgc_ipc_monocle3.R --project-root <PROJECT_ROOT> [--input-dir <dir>] [--outdir <dir>] [--num-dim 50]",
      "",
      "Defaults:",
      "  --input-dir <PROJECT_ROOT>/results/mgeo_rgc_ipc_monocle3/inputs",
      "  --outdir    <PROJECT_ROOT>/results/mgeo_rgc_ipc_monocle3",
      "",
      "Outputs:",
      "  <outdir>/mgeo_rgc_ipc_monocle3_pseudotime.csv",
      "  <outdir>/mgeo_rgc_ipc_monocle3_cds.rds",
      "  <outdir>/mgeo_rgc_ipc_monocle3_summary.tsv",
      sep = "\n"
    )
  )
}

opt <- parse_args(commandArgs(trailingOnly = TRUE))
if (isTRUE(opt$help)) {
  print_usage()
  quit(save = "no", status = 0)
}

required <- c("Matrix", "monocle3", "igraph")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) {
  stop(
    "Missing required R package(s): ", paste(missing, collapse = ", "),
    ". Install/load these packages before running this Slurm stage.",
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(Matrix)
  library(monocle3)
  library(igraph)
})

trim_trailing_slash <- function(x) sub("/+$", "", x)
if (!nzchar(opt$`project-root`)) stop("PROJECT_ROOT or --project-root is required")
project_root <- trim_trailing_slash(opt$`project-root`)
input_dir <- if (is.null(opt$`input-dir`) || !nzchar(opt$`input-dir`)) {
  file.path(project_root, "results/mgeo_rgc_ipc_monocle3/inputs")
} else {
  opt$`input-dir`
}
outdir <- if (is.null(opt$outdir) || !nzchar(opt$outdir)) {
  file.path(project_root, "results/mgeo_rgc_ipc_monocle3")
} else {
  opt$outdir
}
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

num_dim <- as.integer(opt$`num-dim`)
if (is.na(num_dim) || num_dim <= 0L) stop("--num-dim must be a positive integer")

expr_mtx_path <- file.path(input_dir, "mgeo_rgc_ipc_expression_genes_by_cells.mtx")
cell_metadata_path <- file.path(input_dir, "mgeo_rgc_ipc_cell_metadata.csv")
gene_metadata_path <- file.path(input_dir, "mgeo_rgc_ipc_gene_metadata.csv")
pseudotime_path <- file.path(outdir, "mgeo_rgc_ipc_monocle3_pseudotime.csv")
cds_rds_path <- file.path(outdir, "mgeo_rgc_ipc_monocle3_cds.rds")
summary_path <- file.path(outdir, "mgeo_rgc_ipc_monocle3_summary.tsv")

needed <- c(expr_mtx_path, cell_metadata_path, gene_metadata_path)
missing_files <- needed[!file.exists(needed)]
if (length(missing_files) > 0) {
  stop("Missing Monocle3 input file(s): ", paste(missing_files, collapse = ", "))
}

log_msg("Reading matrix: ", expr_mtx_path)
expression_matrix <- Matrix::readMM(expr_mtx_path)
if (!inherits(expression_matrix, "dgCMatrix")) {
  expression_matrix <- as(expression_matrix, "dgCMatrix")
}
log_msg("Reading cell metadata: ", cell_metadata_path)
cell_metadata <- read.csv(cell_metadata_path, row.names = 1, check.names = FALSE)
log_msg("Reading gene metadata: ", gene_metadata_path)
gene_metadata <- read.csv(gene_metadata_path, row.names = 1, check.names = FALSE)

if (nrow(expression_matrix) != nrow(gene_metadata)) {
  stop("Gene metadata rows do not match expression matrix rows")
}
if (ncol(expression_matrix) != nrow(cell_metadata)) {
  stop("Cell metadata rows do not match expression matrix columns")
}
rownames(expression_matrix) <- rownames(gene_metadata)
colnames(expression_matrix) <- rownames(cell_metadata)

log_msg("Creating Monocle3 cell_data_set")
cds <- new_cell_data_set(
  expression_matrix,
  cell_metadata = cell_metadata,
  gene_metadata = gene_metadata
)

set.seed(7)
log_msg("preprocess_cds")
cds <- preprocess_cds(cds, num_dim = num_dim, norm_method = "none")
log_msg("reduce_dimension UMAP")
cds <- reduce_dimension(cds, reduction_method = "UMAP")
log_msg("cluster_cells")
cds <- cluster_cells(cds, reduction_method = "UMAP")
log_msg("learn_graph")
cds <- learn_graph(cds, use_partition = TRUE)

root_cell_ids <- rownames(cell_metadata)[as.logical(cell_metadata$monocle_root_candidate)]
if (length(root_cell_ids) == 0) stop("No monocle_root_candidate cells were provided.")

get_root_pr_nodes <- function(cds, root_cell_ids) {
  closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), , drop = FALSE])
  available_root_cells <- intersect(root_cell_ids, rownames(closest_vertex))
  if (length(available_root_cells) == 0) {
    stop("No root candidate cells are represented in the Monocle3 principal graph mapping.")
  }
  root_vertex_table <- table(closest_vertex[available_root_cells, 1])
  root_vertex <- as.numeric(names(which.max(root_vertex_table)))
  igraph::V(principal_graph(cds)[["UMAP"]])$name[root_vertex]
}

root_pr_nodes <- get_root_pr_nodes(cds, root_cell_ids)
log_msg("order_cells with root_pr_nodes=", paste(root_pr_nodes, collapse = ";"))
cds <- order_cells(cds, reduction_method = "UMAP", root_pr_nodes = root_pr_nodes)

pt <- monocle3::pseudotime(cds)
pt_df <- data.frame(
  cell_id = names(pt),
  rgc_ipc_pseudotime = as.numeric(pt),
  rgc_ipc_pseudotime_method = "monocle3_slurm",
  monocle3_partition = as.character(monocle3::partitions(cds, reduction_method = "UMAP")),
  monocle3_cluster = as.character(monocle3::clusters(cds, reduction_method = "UMAP")),
  monocle3_root_pr_node = paste(root_pr_nodes, collapse = ";"),
  stringsAsFactors = FALSE
)

keep_meta <- intersect(
  c(
    "DIV",
    "original_cell_id",
    "monocle_root_candidate",
    "monocle_root_seed_cell",
    "shi_s5_RGC_score",
    "shi_s5_IPC_score",
    "shi_s5_IPC_minus_RGC_score"
  ),
  colnames(cell_metadata)
)
pt_df <- cbind(pt_df, cell_metadata[pt_df$cell_id, keep_meta, drop = FALSE])

write.csv(pt_df, pseudotime_path, row.names = FALSE)
saveRDS(cds, cds_rds_path)

summary_df <- data.frame(
  key = c(
    "n_cells",
    "n_genes",
    "num_dim",
    "n_root_candidates",
    "root_pr_nodes",
    "pseudotime_csv",
    "cds_rds"
  ),
  value = c(
    ncol(cds),
    nrow(cds),
    num_dim,
    length(root_cell_ids),
    paste(root_pr_nodes, collapse = ";"),
    pseudotime_path,
    cds_rds_path
  )
)
write.table(summary_df, summary_path, sep = "\t", row.names = FALSE, quote = FALSE)

log_msg("Wrote pseudotime CSV: ", pseudotime_path)
log_msg("Saved CDS RDS: ", cds_rds_path)
log_msg("Wrote summary: ", summary_path)
