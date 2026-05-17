#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
options(expressions = 500000)

suppressPackageStartupMessages({
  library(Matrix)
  library(monocle3)
  library(ggplot2)
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

parse_bool <- function(x) {
  x <- tolower(trimws(as.character(x)))
  if (x %in% c("true", "t", "1", "yes", "y")) return(TRUE)
  if (x %in% c("false", "f", "0", "no", "n")) return(FALSE)
  stop("Cannot parse boolean value: ", x)
}

LEARN_GRAPH_USE_PARTITION <- parse_bool(Sys.getenv("MONOCLE_USE_PARTITION", unset = "true"))
PARTITION_LABEL <- paste0("partition_", tolower(as.character(LEARN_GRAPH_USE_PARTITION)))

with_partition_suffix <- function(prefix, ext) {
  paste0(prefix, "_", PARTITION_LABEL, ext)
}

CDS_UMAP_CLUSTER_RDS <- file.path(INTERACTIVE_DIR, with_partition_suffix("cds_after_umap_cluster_no_graph", ".rds"))
UMAP_CLUSTER_CSV <- file.path(INTERACTIVE_DIR, with_partition_suffix("mgeo_rgc_ipc_monocle3_umap_clusters", ".csv"))
CDS_UMAP_RDS <- file.path(INTERACTIVE_DIR, with_partition_suffix("cds_after_umap_cluster_graph", ".rds"))
PSEUDOTIME_CSV <- file.path(INTERACTIVE_DIR, with_partition_suffix("mgeo_rgc_ipc_monocle3_pseudotime_interactive", ".csv"))
CDS_ORDERED_RDS <- file.path(INTERACTIVE_DIR, with_partition_suffix("cds_ordered_interactive", ".rds"))
STATUS_TSV <- file.path(INTERACTIVE_DIR, with_partition_suffix("mgeo_rgc_ipc_monocle3_run_status", ".tsv"))

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

write_status <- function(status, message_text = "") {
  status_df <- data.frame(
    timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    use_partition = LEARN_GRAPH_USE_PARTITION,
    status = status,
    message = message_text,
    stringsAsFactors = FALSE
  )
  write.table(status_df, STATUS_TSV, sep = "\t", quote = FALSE, row.names = FALSE)
  log_step("Status TSV: ", STATUS_TSV)
}

scale_to_unit <- function(x) {
  finite_x <- x[is.finite(x)]
  if (length(finite_x) == 0) {
    return(rep(NA_real_, length(x)))
  }
  lo <- min(finite_x, na.rm = TRUE)
  hi <- stats::quantile(finite_x, 0.99, na.rm = TRUE, names = FALSE)
  if (!is.finite(hi) || hi <= lo) {
    return(ifelse(is.finite(x), 0, NA_real_))
  }
  pmin(pmax((x - lo) / (hi - lo), 0), 1)
}

get_principal_graph_edges <- function(cds) {
  graph <- principal_graph(cds)[["UMAP"]]
  coords <- cds@principal_graph_aux[["UMAP"]]$dp_mst
  if (is.null(graph) || is.null(coords)) {
    return(data.frame())
  }
  if (nrow(coords) != 2 && ncol(coords) == 2) {
    coords <- t(coords)
  }
  edges <- igraph::as_data_frame(graph, what = "edges")
  node_names <- igraph::V(graph)$name
  if (is.null(colnames(coords))) {
    colnames(coords) <- node_names
  }
  edges <- edges[edges$from %in% colnames(coords) & edges$to %in% colnames(coords), ]
  data.frame(
    x = coords[1, edges$from],
    y = coords[2, edges$from],
    xend = coords[1, edges$to],
    yend = coords[2, edges$to],
    stringsAsFactors = FALSE
  )
}

save_marker_validation_figures <- function(cds_ordered) {
  marker_genes <- c(
    "NES",
    "VIM",
    "HES1",
    "ASCL1",
    "DLX2",
    "DCX",
    "STMN2",
    "FABP7",
    "HES5",
    "SOX2",
    "DLX1",
    "NKX2-1"
  )

  available_marker_genes <- intersect(marker_genes, rownames(cds_ordered))
  missing_marker_genes <- setdiff(marker_genes, available_marker_genes)
  if (length(missing_marker_genes) > 0) {
    warning("Missing marker genes: ", paste(missing_marker_genes, collapse = ", "))
  }
  if (length(available_marker_genes) == 0) {
    stop("None of the requested marker genes were found in cds_ordered.")
  }

  umap <- reducedDims(cds_ordered)$UMAP
  umap_df <- data.frame(
    cell_id = rownames(umap),
    UMAP_1 = umap[, 1],
    UMAP_2 = umap[, 2],
    pseudotime = as.numeric(monocle3::pseudotime(cds_ordered)),
    stringsAsFactors = FALSE
  )

  expr_mat <- exprs(cds_ordered)[available_marker_genes, umap_df$cell_id, drop = FALSE]
  if (inherits(expr_mat, "sparseMatrix")) {
    expr_mat <- as.matrix(expr_mat)
  }

  umap_panel_df <- data.frame(
    cell_id = umap_df$cell_id,
    UMAP_1 = umap_df$UMAP_1,
    UMAP_2 = umap_df$UMAP_2,
    panel = "pseudotime",
    value = umap_df$pseudotime,
    value_scaled = scale_to_unit(umap_df$pseudotime),
    stringsAsFactors = FALSE
  )

  for (gene in available_marker_genes) {
    expr_values <- as.numeric(expr_mat[gene, ])
    umap_panel_df <- rbind(
      umap_panel_df,
      data.frame(
        cell_id = umap_df$cell_id,
        UMAP_1 = umap_df$UMAP_1,
        UMAP_2 = umap_df$UMAP_2,
        panel = gene,
        value = expr_values,
        value_scaled = scale_to_unit(expr_values),
        stringsAsFactors = FALSE
      )
    )
  }

  panel_levels <- c("pseudotime", marker_genes[marker_genes %in% available_marker_genes])
  umap_panel_df$panel <- factor(umap_panel_df$panel, levels = panel_levels)
  graph_edges <- get_principal_graph_edges(cds_ordered)

  umap_grid <- ggplot(umap_panel_df, aes(UMAP_1, UMAP_2)) +
    geom_point(aes(color = value_scaled, alpha = value_scaled), size = 0.08, stroke = 0) +
    geom_segment(
      data = graph_edges,
      aes(x = x, y = y, xend = xend, yend = yend),
      inherit.aes = FALSE,
      color = "grey25",
      linewidth = 0.18,
      alpha = 0.45
    ) +
    facet_wrap(~panel, ncol = 4) +
    scale_color_gradientn(
      colors = c("grey92", "gold", "orange", "firebrick"),
      limits = c(0, 1),
      na.value = "grey92",
      name = "scaled value"
    ) +
    scale_alpha(range = c(0.08, 0.95), limits = c(0, 1), guide = "none") +
    coord_equal() +
    labs(x = NULL, y = NULL) +
    theme_void(base_size = 10) +
    theme(
      strip.text = element_text(face = "bold", size = 11),
      legend.position = "right",
      panel.spacing = unit(0.6, "lines")
    )

  umap_png <- file.path(PLOT_DIR, with_partition_suffix("monocle3_marker_umap_grid", ".png"))
  umap_pdf <- file.path(PLOT_DIR, with_partition_suffix("monocle3_marker_umap_grid", ".pdf"))
  timed_step("save marker UMAP PNG/PDF", {
    ggsave(umap_png, umap_grid, width = 12, height = 10, dpi = 220)
    ggsave(umap_pdf, umap_grid, width = 12, height = 10, device = cairo_pdf)
  })

  finite_pt <- is.finite(umap_df$pseudotime)
  pseudotime_gene_df <- data.frame()
  for (gene in available_marker_genes) {
    pseudotime_gene_df <- rbind(
      pseudotime_gene_df,
      data.frame(
        cell_id = umap_df$cell_id[finite_pt],
        pseudotime = umap_df$pseudotime[finite_pt],
        gene = gene,
        expression = as.numeric(expr_mat[gene, finite_pt]),
        stringsAsFactors = FALSE
      )
    )
  }
  pseudotime_gene_df$gene <- factor(
    pseudotime_gene_df$gene,
    levels = marker_genes[marker_genes %in% available_marker_genes]
  )

  pseudotime_grid <- ggplot(pseudotime_gene_df, aes(pseudotime, expression)) +
    geom_point(color = "grey72", alpha = 0.06, size = 0.08, stroke = 0) +
    geom_smooth(
      method = "gam",
      formula = y ~ s(x, k = 10),
      color = "#1f78b4",
      fill = "#a6cee3",
      linewidth = 0.55,
      se = TRUE
    ) +
    facet_wrap(~gene, ncol = 4, scales = "free_y") +
    labs(x = "Monocle3 pseudotime", y = "log-normalized expression") +
    theme_bw(base_size = 10) +
    theme(
      strip.text = element_text(face = "bold", size = 11),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey90", linewidth = 0.2)
    )

  pseudotime_png <- file.path(PLOT_DIR, with_partition_suffix("monocle3_marker_pseudotime_grid", ".png"))
  pseudotime_pdf <- file.path(PLOT_DIR, with_partition_suffix("monocle3_marker_pseudotime_grid", ".pdf"))
  timed_step("save marker pseudotime PNG/PDF", {
    ggsave(pseudotime_png, pseudotime_grid, width = 12, height = 9, dpi = 220)
    ggsave(pseudotime_pdf, pseudotime_grid, width = 12, height = 9, device = cairo_pdf)
  })

  log_step("Marker UMAP PNG: ", umap_png)
  log_step("Marker UMAP PDF: ", umap_pdf)
  log_step("Marker pseudotime PNG: ", pseudotime_png)
  log_step("Marker pseudotime PDF: ", pseudotime_pdf)
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
log_step("LEARN_GRAPH_USE_PARTITION=", LEARN_GRAPH_USE_PARTITION)
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

cds <- tryCatch(
  timed_step(paste0("learn_graph use_partition=", LEARN_GRAPH_USE_PARTITION), {
    learn_graph(cds, use_partition = LEARN_GRAPH_USE_PARTITION)
  }),
  error = function(e) {
    write_status("learn_graph_failed", conditionMessage(e))
    stop(e)
  }
)
log_step("Principal graph nodes: ", length(igraph::V(principal_graph(cds)[["UMAP"]])))
log_step("Principal graph edges: ", length(igraph::E(principal_graph(cds)[["UMAP"]])))
safe_save_rds(cds, CDS_UMAP_RDS)

save_plot(
  with_partition_suffix("monocle3_div", ".png"),
  plot_cells(
    cds,
    color_cells_by = "DIV",
    label_groups_by_cluster = FALSE,
    label_leaves = FALSE,
    label_branch_points = FALSE
  )
)
save_plot(
  with_partition_suffix("monocle3_rgc_score", ".png"),
  plot_cells(
    cds,
    color_cells_by = "shi_s5_RGC_score",
    label_groups_by_cluster = FALSE,
    label_leaves = FALSE,
    label_branch_points = FALSE
  )
)
save_plot(
  with_partition_suffix("monocle3_ipc_score", ".png"),
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
  with_partition_suffix("monocle3_pseudotime", ".png"),
  plot_cells(
    cds_ordered,
    color_cells_by = "pseudotime",
    label_groups_by_cluster = FALSE,
    label_leaves = TRUE,
    label_branch_points = TRUE
  )
)

save_marker_validation_figures(cds_ordered)

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
write_status("completed", "Run completed successfully.")
log_step("Done")
