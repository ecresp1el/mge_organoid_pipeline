#!/usr/bin/env Rscript

# VS Code-friendly interactive Monocle3 workflow for MGEO RGC/IPC pseudotime.
#
# Use this on an allocated compute node:
#   cd /home/elcrespo/Desktop/githubprojects/mge_organoid_pipeline
#   module load Bioinformatics
#   module load Rmonocle3/1.3.7
#   R
#
# Then run sections block-by-block in VS Code's R terminal, or source this file
# after editing the config below. Section markers use "# %%" for editors that
# recognize notebook-style cells.

# %% 0. Setup

options(stringsAsFactors = FALSE)
options(expressions = 500000)

suppressPackageStartupMessages({
  library(Matrix)
  library(monocle3)
  library(igraph)
})

PROJECT_ROOT <- Sys.getenv(
  "PROJECT_ROOT",
  unset = "/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder"
)

RUN_ROOT <- file.path(PROJECT_ROOT, "results/mgeo_rgc_ipc_monocle3")
INPUT_DIR <- file.path(RUN_ROOT, "inputs")
INTERACTIVE_DIR <- file.path(RUN_ROOT, "interactive")
dir.create(INTERACTIVE_DIR, recursive = TRUE, showWarnings = FALSE)

EXPR_MTX_PATH <- file.path(INPUT_DIR, "mgeo_rgc_ipc_expression_genes_by_cells.mtx")
CELL_METADATA_PATH <- file.path(INPUT_DIR, "mgeo_rgc_ipc_cell_metadata.csv")
GENE_METADATA_PATH <- file.path(INPUT_DIR, "mgeo_rgc_ipc_gene_metadata.csv")

CDS_PREPROCESS_RDS <- file.path(INTERACTIVE_DIR, "cds_after_preprocess.rds")
CDS_UMAP_CLUSTER_RDS <- file.path(INTERACTIVE_DIR, "cds_after_umap_cluster_no_graph.rds")
UMAP_CLUSTER_CSV <- file.path(INTERACTIVE_DIR, "mgeo_rgc_ipc_monocle3_umap_clusters.csv")
CDS_UMAP_RDS <- file.path(INTERACTIVE_DIR, "cds_after_umap_cluster_graph.rds")
PSEUDOTIME_CSV <- file.path(INTERACTIVE_DIR, "mgeo_rgc_ipc_monocle3_pseudotime_interactive.csv")
CDS_ORDERED_RDS <- file.path(INTERACTIVE_DIR, "cds_ordered_interactive.rds")

NUM_DIM <- 50

# TRUE is more stable for this full ~79k-cell subset. FALSE tries to force one
# connected graph and can fail with node stack overflow on this dataset.
LEARN_GRAPH_USE_PARTITION <- TRUE

log_step <- function(...) {
  message(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste0(..., collapse = ""))
  flush.console()
}

timed_step <- function(label, expr) {
  log_step("START ", label)
  start_time <- proc.time()[["elapsed"]]
  result <- force(expr)
  elapsed_min <- (proc.time()[["elapsed"]] - start_time) / 60
  log_step("DONE ", label, " (", round(elapsed_min, 2), " min)")
  result
}

safe_save_rds <- function(object, path) {
  tryCatch(
    timed_step(paste0("save RDS ", basename(path)), {
      saveRDS(object, path, compress = FALSE)
    }),
    error = function(e) {
      warning(
        "Could not save RDS checkpoint to ",
        path,
        ": ",
        conditionMessage(e),
        immediate. = TRUE
      )
      invisible(FALSE)
    }
  )
}

get_umap_clusters <- function(cds) {
  tryCatch(
    clusters(cds, reduction_method = "UMAP"),
    error = function(e) clusters(cds)
  )
}

get_umap_partitions <- function(cds) {
  tryCatch(
    partitions(cds, reduction_method = "UMAP"),
    error = function(e) partitions(cds)
  )
}

write_umap_cluster_checkpoint <- function(cds, path) {
  umap <- reducedDims(cds)$UMAP
  cluster_values <- get_umap_clusters(cds)
  checkpoint <- data.frame(
    cell_id = rownames(umap),
    monocle3_umap_1 = umap[, 1],
    monocle3_umap_2 = umap[, 2],
    monocle3_cluster = as.character(cluster_values),
    stringsAsFactors = FALSE
  )
  timed_step(paste0("write CSV ", basename(path)), {
    write.csv(checkpoint, path, row.names = FALSE)
  })
}

has_umap <- function(cds) {
  "UMAP" %in% names(reducedDims(cds)) &&
    !is.null(reducedDims(cds)$UMAP) &&
    nrow(reducedDims(cds)$UMAP) == ncol(cds)
}

has_umap_clusters <- function(cds) {
  cluster_values <- tryCatch(
    get_umap_clusters(cds),
    error = function(e) NULL
  )
  !is.null(cluster_values) && length(cluster_values) == ncol(cds)
}

cat("PROJECT_ROOT:", PROJECT_ROOT, "\n")
cat("INPUT_DIR:", INPUT_DIR, "\n")
cat("INTERACTIVE_DIR:", INTERACTIVE_DIR, "\n")

# %% 0b. Runtime diagnostics

cat("=== Host / Slurm ===\n")
cat("hostname:", Sys.info()[["nodename"]], "\n")
cat("user:", Sys.info()[["user"]], "\n")
cat("SLURM_JOB_ID:", Sys.getenv("SLURM_JOB_ID", unset = "<unset>"), "\n")
cat("SLURM_JOB_NAME:", Sys.getenv("SLURM_JOB_NAME", unset = "<unset>"), "\n")
cat("SLURM_CPUS_PER_TASK:", Sys.getenv("SLURM_CPUS_PER_TASK", unset = "<unset>"), "\n")
cat("SLURM_MEM_PER_NODE:", Sys.getenv("SLURM_MEM_PER_NODE", unset = "<unset>"), "\n")
cat("SLURM_JOB_NODELIST:", Sys.getenv("SLURM_JOB_NODELIST", unset = "<unset>"), "\n")

cat("\n=== R ===\n")
cat("R.home():", R.home(), "\n")
cat("R.version:", R.version.string, "\n")
cat("R executable from shell:", system("which R", intern = TRUE), "\n")
cat("Rscript executable from shell:", system("which Rscript", intern = TRUE), "\n")
cat(".libPaths():\n")
print(.libPaths())

cat("\n=== Modules ===\n")
module_list <- tryCatch(
  system("module list 2>&1", intern = TRUE),
  error = function(e) paste("module list failed:", conditionMessage(e))
)
cat(paste(module_list, collapse = "\n"), "\n")

cat("\n=== Package Availability ===\n")
pkgs <- c("monocle3", "Matrix", "igraph", "HDF5Array", "sf", "spdep")
print(data.frame(package = pkgs, installed = sapply(pkgs, requireNamespace, quietly = TRUE)), row.names = FALSE)

cat("\n=== Parallel / Memory Hints ===\n")
cat("parallel::detectCores():", parallel::detectCores(), "\n")
meminfo <- tryCatch(readLines("/proc/meminfo", n = 3), error = function(e) character(0))
cat(paste(meminfo, collapse = "\n"), "\n")

cat("\n=== Paths ===\n")
cat("getwd():", getwd(), "\n")
cat("PROJECT_ROOT:", PROJECT_ROOT, "\n")
cat("INPUT_DIR:", INPUT_DIR, "\n")
cat("INTERACTIVE_DIR:", INTERACTIVE_DIR, "\n")
cat("Expression matrix exists:", file.exists(EXPR_MTX_PATH), "\n")
cat("Cell metadata exists:", file.exists(CELL_METADATA_PATH), "\n")
cat("Gene metadata exists:", file.exists(GENE_METADATA_PATH), "\n")
if (file.exists(EXPR_MTX_PATH)) {
  cat("Expression matrix GB:", round(file.info(EXPR_MTX_PATH)$size / 1e9, 3), "\n")
}

stopifnot(requireNamespace("monocle3", quietly = TRUE))
stopifnot(grepl("^gl[0-9]+", Sys.info()[["nodename"]]))

if (Sys.getenv("SLURM_JOB_ID", unset = "") == "") {
  warning(
    "SLURM_JOB_ID is unset. This can happen when VS Code launches a fresh ",
    "R process on a compute node instead of using the R terminal that ",
    "inherited salloc. Hostname is still a compute node: ",
    Sys.info()[["nodename"]]
  )
}

# %% 1. Read exported matrix and metadata

stopifnot(file.exists(EXPR_MTX_PATH))
stopifnot(file.exists(CELL_METADATA_PATH))
stopifnot(file.exists(GENE_METADATA_PATH))

log_step("Reading expression matrix from ", EXPR_MTX_PATH)
expression_matrix <- timed_step("Matrix::readMM", {
  Matrix::readMM(EXPR_MTX_PATH)
})
if (!inherits(expression_matrix, "dgCMatrix")) {
  expression_matrix <- timed_step("Convert expression matrix to dgCMatrix", {
    as(expression_matrix, "dgCMatrix")
  })
}

log_step("Reading cell metadata from ", CELL_METADATA_PATH)
cell_metadata <- timed_step("read cell metadata CSV", {
  read.csv(CELL_METADATA_PATH, row.names = 1, check.names = FALSE)
})

log_step("Reading gene metadata from ", GENE_METADATA_PATH)
gene_metadata <- timed_step("read gene metadata CSV", {
  read.csv(GENE_METADATA_PATH, row.names = 1, check.names = FALSE)
})

stopifnot(nrow(expression_matrix) == nrow(gene_metadata))
stopifnot(ncol(expression_matrix) == nrow(cell_metadata))

rownames(expression_matrix) <- rownames(gene_metadata)
colnames(expression_matrix) <- rownames(cell_metadata)

log_step(
  "Loaded matrix with ", nrow(expression_matrix), " genes x ",
  ncol(expression_matrix), " cells; nonzero entries = ", length(expression_matrix@x)
)

dim(expression_matrix)
head(cell_metadata[, c("DIV", "shi_s5_RGC_score", "shi_s5_IPC_score")])

# %% 2. Create CDS and preprocess

cds <- timed_step("new_cell_data_set", {
  new_cell_data_set(
    expression_matrix,
    cell_metadata = cell_metadata,
    gene_metadata = gene_metadata
  )
})

set.seed(7)
cds <- timed_step(paste0("preprocess_cds num_dim=", NUM_DIM), {
  preprocess_cds(cds, num_dim = NUM_DIM, norm_method = "none")
})
safe_save_rds(cds, CDS_PREPROCESS_RDS)
CDS_PREPROCESS_RDS

# %% 3. UMAP, clusters, graph

set.seed(7)
if (!has_umap(cds)) {
  cds <- timed_step("reduce_dimension UMAP", {
    reduce_dimension(cds, reduction_method = "UMAP")
  })
} else {
  log_step("Skipping reduce_dimension; UMAP already exists.")
}
log_step("UMAP coordinates: ", paste(dim(reducedDims(cds)$UMAP), collapse = " x "))

if (!has_umap_clusters(cds)) {
  cds <- timed_step("cluster_cells UMAP", {
    cluster_cells(cds, reduction_method = "UMAP")
  })
} else {
  log_step("Skipping cluster_cells; UMAP clusters already exist.")
}
log_step("Cluster count: ", length(unique(get_umap_clusters(cds))))

write_umap_cluster_checkpoint(cds, UMAP_CLUSTER_CSV)
safe_save_rds(cds, CDS_UMAP_CLUSTER_RDS)

cds <- tryCatch(
  timed_step(paste0("learn_graph use_partition=", LEARN_GRAPH_USE_PARTITION), {
    learn_graph(cds, use_partition = LEARN_GRAPH_USE_PARTITION)
  }),
  error = function(e) {
    write_umap_cluster_checkpoint(cds, UMAP_CLUSTER_CSV)
    safe_save_rds(cds, CDS_UMAP_CLUSTER_RDS)
    stop(
      "learn_graph failed after saving the UMAP/cluster checkpoint to ",
      UMAP_CLUSTER_CSV,
      ". Error: ",
      conditionMessage(e),
      call. = FALSE
    )
  }
)
log_step("Principal graph nodes: ", length(igraph::V(principal_graph(cds)[["UMAP"]])))
log_step("Principal graph edges: ", length(igraph::E(principal_graph(cds)[["UMAP"]])))

safe_save_rds(cds, CDS_UMAP_RDS)
CDS_UMAP_RDS

# Quick visual checkpoint. In VS Code/RStudio this should draw in the plot pane.
plot_cells(
  cds,
  color_cells_by = "DIV",
  label_groups_by_cluster = FALSE,
  label_leaves = FALSE,
  label_branch_points = FALSE
)

plot_cells(
  cds,
  color_cells_by = "shi_s5_RGC_score",
  label_groups_by_cluster = FALSE,
  label_leaves = FALSE,
  label_branch_points = FALSE
)

plot_cells(
  cds,
  color_cells_by = "shi_s5_IPC_score",
  label_groups_by_cluster = FALSE,
  label_leaves = FALSE,
  label_branch_points = FALSE
)

# %% 4. Inspect root candidates

root_cells <- rownames(cell_metadata)[as.logical(cell_metadata$monocle_root_candidate)]
root_seed_cell <- rownames(cell_metadata)[as.logical(cell_metadata$monocle_root_seed_cell)]

length(root_cells)
root_seed_cell

root_summary <- cell_metadata[root_cells, c(
  "DIV",
  "shi_s5_RGC_score",
  "shi_s5_IPC_score",
  "shi_s5_IPC_minus_RGC_score"
)]
summary(root_summary)

root_by_div <- table(root_summary$DIV)
root_by_div

# %% 5. Convert root cells to principal graph nodes

get_root_pr_nodes <- function(cds, root_cell_ids) {
  closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), , drop = FALSE])
  available_root_cells <- intersect(root_cell_ids, rownames(closest_vertex))
  if (length(available_root_cells) == 0) {
    stop("No root cells represented in principal graph mapping.")
  }
  root_vertex_table <- sort(table(closest_vertex[available_root_cells, 1]), decreasing = TRUE)
  root_vertex_ids <- as.numeric(names(root_vertex_table))
  graph_node_names <- igraph::V(principal_graph(cds)[["UMAP"]])$name[root_vertex_ids]
  list(
    table = root_vertex_table,
    primary = graph_node_names[[1]],
    all = graph_node_names
  )
}

root_nodes <- get_root_pr_nodes(cds, root_cells)
root_nodes$table[1:20]
root_nodes$primary

# %% 6. Order cells

# Try root_cells first. This is easier to reason about interactively than only
# passing one graph node. If it fails, use root_pr_nodes below.
cds_ordered <- tryCatch(
  order_cells(cds, reduction_method = "UMAP", root_cells = root_cells),
  error = function(e) {
    message("order_cells(root_cells=...) failed: ", conditionMessage(e))
    NULL
  }
)

if (is.null(cds_ordered)) {
  cds_ordered <- order_cells(
    cds,
    reduction_method = "UMAP",
    root_pr_nodes = root_nodes$primary
  )
}

pt <- monocle3::pseudotime(cds_ordered)
summary(pt)
table(is.finite(pt), useNA = "ifany")

# %% 7. Plot ordered cells

plot_cells(
  cds_ordered,
  color_cells_by = "pseudotime",
  label_groups_by_cluster = FALSE,
  label_leaves = TRUE,
  label_branch_points = TRUE
)

plot_cells(
  cds_ordered,
  color_cells_by = "DIV",
  label_groups_by_cluster = FALSE,
  label_leaves = FALSE,
  label_branch_points = FALSE
)

# %% 8. Export pseudotime table for Python

pt_df <- data.frame(
  cell_id = names(pt),
  rgc_ipc_pseudotime = as.numeric(pt),
  rgc_ipc_pseudotime_method = paste0(
    "monocle3_interactive_use_partition_",
    LEARN_GRAPH_USE_PARTITION
  ),
  monocle3_partition = as.character(get_umap_partitions(cds_ordered)),
  monocle3_cluster = as.character(get_umap_clusters(cds_ordered)),
  monocle3_root_pr_node = root_nodes$primary,
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

write.csv(pt_df, PSEUDOTIME_CSV, row.names = FALSE)
saveRDS(cds_ordered, CDS_ORDERED_RDS)

PSEUDOTIME_CSV
CDS_ORDERED_RDS
