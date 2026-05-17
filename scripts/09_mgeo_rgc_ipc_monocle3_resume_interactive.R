#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
options(expressions = 500000)

suppressPackageStartupMessages({
  library(Matrix)
  library(monocle3)
})

PROJECT_ROOT <- Sys.getenv(
  "PROJECT_ROOT",
  unset = "/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder"
)

RUN_ROOT <- file.path(PROJECT_ROOT, "results/mgeo_rgc_ipc_monocle3")
INTERACTIVE_DIR <- file.path(RUN_ROOT, "interactive")
PLOT_DIR <- file.path(INTERACTIVE_DIR, "plots")
dir.create(INTERACTIVE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)

CDS_PREPROCESS_RDS <- file.path(INTERACTIVE_DIR, "cds_after_preprocess.rds")
CDS_UMAP_CLUSTER_RDS <- file.path(INTERACTIVE_DIR, "cds_after_umap_cluster_no_graph.rds")
UMAP_CLUSTER_CSV <- file.path(INTERACTIVE_DIR, "mgeo_rgc_ipc_monocle3_umap_clusters.csv")
CDS_UMAP_RDS <- file.path(INTERACTIVE_DIR, "cds_after_umap_cluster_graph.rds")
PSEUDOTIME_CSV <- file.path(INTERACTIVE_DIR, "mgeo_rgc_ipc_monocle3_pseudotime_interactive.csv")
CDS_ORDERED_RDS <- file.path(INTERACTIVE_DIR, "cds_ordered_interactive.rds")

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
  monocle3::clusters(cds, reduction_method = "UMAP")
}

get_umap_partitions <- function(cds) {
  monocle3::partitions(cds, reduction_method = "UMAP")
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

save_plot <- function(filename, expr, width = 8, height = 6) {
  path <- file.path(PLOT_DIR, filename)
  png(path, width = width, height = height, units = "in", res = 180)
  print(expr)
  dev.off()
  log_step("Saved plot: ", path)
  path
}

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

log_step("hostname=", Sys.info()[["nodename"]])
log_step("SLURM_JOB_ID=", Sys.getenv("SLURM_JOB_ID", unset = "<unset>"))
log_step("PROJECT_ROOT=", PROJECT_ROOT)
stopifnot(file.exists(CDS_PREPROCESS_RDS))

cds <- timed_step("read preprocessed CDS RDS", {
  readRDS(CDS_PREPROCESS_RDS)
})

cell_metadata <- as.data.frame(colData(cds))
gene_metadata <- as.data.frame(rowData(cds))
log_step("Loaded preprocessed CDS: ", nrow(cds), " genes x ", ncol(cds), " cells")

set.seed(7)
cds <- timed_step("reduce_dimension UMAP", {
  reduce_dimension(cds, reduction_method = "UMAP")
})
log_step("UMAP coordinates: ", paste(dim(reducedDims(cds)$UMAP), collapse = " x "))

cds <- timed_step("cluster_cells UMAP", {
  cluster_cells(cds, reduction_method = "UMAP")
})
log_step("Cluster count: ", length(unique(get_umap_clusters(cds))))
write_umap_cluster_checkpoint(cds, UMAP_CLUSTER_CSV)
safe_save_rds(cds, CDS_UMAP_CLUSTER_RDS)

cds <- timed_step(paste0("learn_graph use_partition=", LEARN_GRAPH_USE_PARTITION), {
  learn_graph(cds, use_partition = LEARN_GRAPH_USE_PARTITION)
})
log_step("Principal graph nodes: ", length(igraph::V(principal_graph(cds)[["UMAP"]])))
log_step("Principal graph edges: ", length(igraph::E(principal_graph(cds)[["UMAP"]])))
safe_save_rds(cds, CDS_UMAP_RDS)

save_plot(
  "monocle3_div.png",
  plot_cells(
    cds,
    color_cells_by = "DIV",
    label_groups_by_cluster = FALSE,
    label_leaves = FALSE,
    label_branch_points = FALSE
  )
)
save_plot(
  "monocle3_rgc_score.png",
  plot_cells(
    cds,
    color_cells_by = "shi_s5_RGC_score",
    label_groups_by_cluster = FALSE,
    label_leaves = FALSE,
    label_branch_points = FALSE
  )
)
save_plot(
  "monocle3_ipc_score.png",
  plot_cells(
    cds,
    color_cells_by = "shi_s5_IPC_score",
    label_groups_by_cluster = FALSE,
    label_leaves = FALSE,
    label_branch_points = FALSE
  )
)

root_cells <- rownames(cell_metadata)[as.logical(cell_metadata$monocle_root_candidate)]
root_seed_cell <- rownames(cell_metadata)[as.logical(cell_metadata$monocle_root_seed_cell)]
if (length(root_cells) == 0) stop("No monocle_root_candidate cells were provided.")
log_step("Root candidate cells: ", length(root_cells))
log_step("Root seed cell: ", paste(root_seed_cell, collapse = ";"))

root_summary <- cell_metadata[root_cells, c(
  "DIV",
  "shi_s5_RGC_score",
  "shi_s5_IPC_score",
  "shi_s5_IPC_minus_RGC_score"
)]
print(summary(root_summary))
print(table(root_summary$DIV))

root_nodes <- get_root_pr_nodes(cds, root_cells)
log_step("Primary root principal node: ", root_nodes$primary)
print(head(root_nodes$table, 20))

cds_ordered <- tryCatch(
  timed_step("order_cells root_cells", {
    order_cells(cds, reduction_method = "UMAP", root_cells = root_cells)
  }),
  error = function(e) {
    log_step("order_cells(root_cells=...) failed: ", conditionMessage(e))
    NULL
  }
)

if (is.null(cds_ordered)) {
  cds_ordered <- timed_step("order_cells root_pr_nodes primary", {
    order_cells(cds, reduction_method = "UMAP", root_pr_nodes = root_nodes$primary)
  })
}

pt <- monocle3::pseudotime(cds_ordered)
print(summary(pt))
print(table(is.finite(pt), useNA = "ifany"))

save_plot(
  "monocle3_pseudotime.png",
  plot_cells(
    cds_ordered,
    color_cells_by = "pseudotime",
    label_groups_by_cluster = FALSE,
    label_leaves = TRUE,
    label_branch_points = TRUE
  )
)

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

timed_step(paste0("write CSV ", basename(PSEUDOTIME_CSV)), {
  write.csv(pt_df, PSEUDOTIME_CSV, row.names = FALSE)
})
safe_save_rds(cds_ordered, CDS_ORDERED_RDS)

log_step("Pseudotime CSV: ", PSEUDOTIME_CSV)
log_step("Done")
