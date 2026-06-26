#!/usr/bin/env Rscript

# Save true Seurat anchor/projection assets for Varela DIV30/DIV90 to Shi 2019.
#
# This targeted helper reruns the same whole-Shi FindTransferAnchors setup used
# by scripts/13_run_cross_study_shi_seurat_label_transfer.R, but it does not
# overwrite finalized prediction outputs. It saves the TransferAnchorSet objects,
# full anchor-pair tables, top anchor per query cell, and static 3D validation
# plots showing query cells and Shi reference cells on separated planes.

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
})

PROJECT_ROOT_DEFAULT <- "/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder"
RUN_LABEL_DEFAULT <- "cross_study_shi_seurat_anchor_projection_v1"
TRANSFER_HELPER_SCRIPT_DEFAULT <- "scripts/13_run_cross_study_shi_seurat_label_transfer.R"
REFERENCE_DEFAULT <- file.path(PROJECT_ROOT_DEFAULT, "results/shi_2019_paper_qc/shi_2019_seurat.rds")

timestamp <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")
log_msg <- function(...) message("[", timestamp(), "] ", paste0(..., collapse = ""))

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  out <- list(
    project_root = Sys.getenv("PROJECT_ROOT", PROJECT_ROOT_DEFAULT),
    run_label = Sys.getenv("SHI_ANCHOR_PROJECTION_RUN_LABEL", RUN_LABEL_DEFAULT),
    outdir = "",
    reference = Sys.getenv("SHI_REFERENCE_RDS", REFERENCE_DEFAULT),
    reference_labels_tsv = Sys.getenv(
      "SHI_REFERENCE_LABELS_TSV",
      file.path(
        Sys.getenv("PROJECT_ROOT", PROJECT_ROOT_DEFAULT),
        "results/shi_reference_div30_seurat_label_transfer/shi_reference_div30_seurat_label_transfer_v1/seurat/shi_reference_labels_for_seurat.tsv"
      )
    ),
    helper_script = Sys.getenv("SHI_TRANSFER_HELPER_SCRIPT", TRANSFER_HELPER_SCRIPT_DEFAULT),
    study_id = c("varela_div30", "varela_div90"),
    dims = "50",
    min_shared_features = "500",
    normalization_method = "LogNormalize",
    seed = "17",
    max_query_points = "6000",
    max_reference_points = "6000",
    max_anchor_lines = "3500",
    save_anchor_rds = "true"
  )
  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    value <- if (i < length(args)) args[[i + 1L]] else ""
    if (key == "--project-root") {
      out$project_root <- value; i <- i + 2L
    } else if (key == "--run-label") {
      out$run_label <- value; i <- i + 2L
    } else if (key == "--outdir") {
      out$outdir <- value; i <- i + 2L
    } else if (key == "--reference") {
      out$reference <- value; i <- i + 2L
    } else if (key == "--reference-labels-tsv") {
      out$reference_labels_tsv <- value; i <- i + 2L
    } else if (key == "--helper-script") {
      out$helper_script <- value; i <- i + 2L
    } else if (key == "--study-id") {
      out$study_id <- unlist(strsplit(value, "[,;[:space:]]+")); i <- i + 2L
    } else if (key == "--dims") {
      out$dims <- value; i <- i + 2L
    } else if (key == "--min-shared-features") {
      out$min_shared_features <- value; i <- i + 2L
    } else if (key == "--normalization-method") {
      out$normalization_method <- value; i <- i + 2L
    } else if (key == "--seed") {
      out$seed <- value; i <- i + 2L
    } else if (key == "--max-query-points") {
      out$max_query_points <- value; i <- i + 2L
    } else if (key == "--max-reference-points") {
      out$max_reference_points <- value; i <- i + 2L
    } else if (key == "--max-anchor-lines") {
      out$max_anchor_lines <- value; i <- i + 2L
    } else if (key == "--save-anchor-rds") {
      out$save_anchor_rds <- value; i <- i + 2L
    } else {
      stop("Unknown argument: ", key, call. = FALSE)
    }
  }
  out$study_id <- out$study_id[nzchar(out$study_id)]
  if (!nzchar(out$outdir)) {
    out$outdir <- file.path(out$project_root, "results/cross_study_shi_seurat_label_transfer", out$run_label)
  }
  out
}

to_bool <- function(x, default = FALSE) {
  if (is.null(x) || !nzchar(x)) return(default)
  tolower(trimws(as.character(x))) %in% c("1", "true", "t", "yes", "y")
}

load_transfer_helpers <- function(path) {
  if (!file.exists(path)) stop("Missing transfer helper script: ", path, call. = FALSE)
  lines <- readLines(path, warn = FALSE)
  keep <- !grepl("^\\s*main\\(\\)\\s*$", lines)
  env <- new.env(parent = globalenv())
  eval(parse(text = lines[keep]), envir = env)
  env
}

write_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
}

write_tsv_gz <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  con <- gzfile(path, open = "wt")
  on.exit(close(con), add = TRUE)
  utils::write.table(x, con, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
}

first_existing_reduction <- function(obj, preferred) {
  for (red in preferred) {
    if (red %in% names(obj@reductions)) return(red)
  }
  ""
}

coordinate_table <- function(obj, preferred, label = "object") {
  reduction <- first_existing_reduction(obj, preferred)
  if (!nzchar(reduction)) {
    stop(label, " has none of the requested reductions: ", paste(preferred, collapse = ", "), call. = FALSE)
  }
  emb <- Seurat::Embeddings(obj, reduction = reduction)
  if (ncol(emb) < 2L) stop(label, " reduction has fewer than 2 dimensions: ", reduction, call. = FALSE)
  data.frame(
    cell_id = rownames(emb),
    coord_1 = as.numeric(emb[, 1]),
    coord_2 = as.numeric(emb[, 2]),
    coordinate_reduction = reduction,
    stringsAsFactors = FALSE
  )
}

resolve_anchor_cells <- function(values, candidates) {
  if (is.numeric(values) || is.integer(values)) {
    idx <- as.integer(values)
    out <- rep(NA_character_, length(idx))
    ok <- !is.na(idx) & idx >= 1L & idx <= length(candidates)
    out[ok] <- candidates[idx[ok]]
    return(out)
  }
  as.character(values)
}

anchor_pair_table <- function(anchors, reference, query, label_col, week_col) {
  anchor_df <- as.data.frame(anchors@anchors)
  if (!all(c("cell1", "cell2") %in% colnames(anchor_df))) {
    stop("Anchor table does not contain expected cell1/cell2 columns", call. = FALSE)
  }
  ref_candidates <- colnames(reference)
  query_candidates <- colnames(query)
  if ("object.list" %in% slotNames(anchors)) {
    ref_candidates <- tryCatch(colnames(anchors@object.list[[1]]), error = function(e) ref_candidates)
    query_candidates <- tryCatch(colnames(anchors@object.list[[2]]), error = function(e) query_candidates)
  }
  anchor_df$reference_cell_id <- resolve_anchor_cells(anchor_df$cell1, ref_candidates)
  anchor_df$query_cell_id <- resolve_anchor_cells(anchor_df$cell2, query_candidates)
  ref_meta <- reference@meta.data
  query_meta <- query@meta.data
  anchor_df$reference_shi_label <- as.character(ref_meta[anchor_df$reference_cell_id, label_col])
  anchor_df$reference_shi_week_label <- as.character(ref_meta[anchor_df$reference_cell_id, week_col])
  anchor_df$query_cluster <- if ("seurat_clusters" %in% colnames(query_meta)) as.character(query_meta[anchor_df$query_cell_id, "seurat_clusters"]) else ""
  anchor_df$query_sample <- if ("orig.ident" %in% colnames(query_meta)) as.character(query_meta[anchor_df$query_cell_id, "orig.ident"]) else ""
  anchor_df
}

top_anchor_per_query <- function(anchor_df) {
  score_col <- if ("score" %in% colnames(anchor_df)) "score" else ""
  if (nzchar(score_col)) {
    anchor_df <- anchor_df[order(anchor_df$query_cell_id, -as.numeric(anchor_df[[score_col]])), , drop = FALSE]
  } else {
    anchor_df <- anchor_df[order(anchor_df$query_cell_id), , drop = FALSE]
  }
  anchor_df[!duplicated(anchor_df$query_cell_id), , drop = FALSE]
}

scale_to_unit <- function(x) {
  x <- as.numeric(x)
  rng <- range(x, finite = TRUE)
  if (!all(is.finite(rng)) || abs(diff(rng)) < 1e-12) return(rep(0, length(x)))
  ((x - rng[[1]]) / diff(rng) - 0.5) * 2
}

sample_rows <- function(df, n, seed) {
  if (nrow(df) <= n) return(df)
  set.seed(seed)
  df[sort(sample.int(nrow(df), n)), , drop = FALSE]
}

plot_anchor_projection <- function(study_id, ref_coords, query_coords, links, out_prefix, opt) {
  if (!requireNamespace("scatterplot3d", quietly = TRUE)) {
    log_msg("Package scatterplot3d is unavailable; skipping static anchor projection plot for ", study_id)
    return(FALSE)
  }
  max_query <- as.integer(opt$max_query_points)
  max_reference <- as.integer(opt$max_reference_points)
  max_lines <- as.integer(opt$max_anchor_lines)
  seed <- as.integer(opt$seed)
  ref_plot <- sample_rows(ref_coords, max_reference, seed + 1L)
  query_plot <- sample_rows(query_coords, max_query, seed + 2L)
  line_plot <- sample_rows(links[links$plot_include, , drop = FALSE], max_lines, seed + 3L)

  ref_plot$x <- scale_to_unit(ref_plot$coord_1) * 10
  ref_plot$y <- scale_to_unit(ref_plot$coord_2) * 10
  query_plot$x <- scale_to_unit(query_plot$coord_1) * 10
  query_plot$y <- scale_to_unit(query_plot$coord_2) * 10
  ref_xy <- ref_plot[, c("cell_id", "x", "y")]
  query_xy <- query_plot[, c("cell_id", "x", "y")]
  names(ref_xy) <- c("reference_cell_id", "reference_x", "reference_y")
  names(query_xy) <- c("query_cell_id", "query_x", "query_y")
  line_plot <- merge(line_plot, ref_xy, by = "reference_cell_id", all.x = FALSE, sort = FALSE)
  line_plot <- merge(line_plot, query_xy, by = "query_cell_id", all.x = FALSE, sort = FALSE)
  if (nrow(line_plot) == 0L) {
    log_msg("No sampled line links survived coordinate merge for ", study_id)
    return(FALSE)
  }
  label_cols <- c(
    MGE = "#16697a", LGE = "#4895ef", CGE = "#52b788", progenitor = "#f4a261",
    `Excitatory IPC` = "#e76f51", `Excitatory neuron` = "#d62828",
    `Thalamic neurons` = "#7b2cbf", Microglia = "#6c757d",
    OPC = "#a7c957", Endothelial = "#2a9d8f"
  )
  ref_col <- unname(label_cols[as.character(ref_plot$shi_label)])
  ref_col[is.na(ref_col)] <- "#888888"
  line_col <- unname(label_cols[as.character(line_plot$reference_shi_label)])
  line_col[is.na(line_col)] <- "#666666"
  score <- if ("score" %in% colnames(line_plot)) as.numeric(line_plot$score) else rep(0.65, nrow(line_plot))
  line_lwd <- pmax(0.2, 0.25 + 2.8 * (score - min(score, na.rm = TRUE)) / max(diff(range(score, finite = TRUE)), 1e-6))

  draw_one <- function(path, device) {
    if (device == "png") grDevices::png(path, width = 5200, height = 3800, res = 450)
    if (device == "pdf") grDevices::pdf(path, width = 11.6, height = 8.4, useDingbats = FALSE)
    if (device == "svg") grDevices::svg(path, width = 11.6, height = 8.4)
    on.exit(grDevices::dev.off(), add = TRUE)
    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar), add = TRUE)
    par(mar = c(2.5, 2.5, 3.2, 1.0), family = "sans")
    s3d <- scatterplot3d::scatterplot3d(
      x = c(query_plot$x, ref_plot$x),
      y = c(query_plot$y, ref_plot$y),
      z = c(rep(0, nrow(query_plot)), rep(8, nrow(ref_plot))),
      type = "n",
      angle = 45,
      scale.y = 0.65,
      xlab = "scaled UMAP / reference coordinate 1",
      ylab = "scaled UMAP / reference coordinate 2",
      zlab = "projection plane",
      main = paste0(study_id, ": query cells to Shi reference anchors")
    )
    s3d$points3d(query_plot$x, query_plot$y, rep(0, nrow(query_plot)), pch = 16, cex = 0.22, col = grDevices::adjustcolor("#bdbdbd", alpha.f = 0.45))
    s3d$points3d(ref_plot$x, ref_plot$y, rep(8, nrow(ref_plot)), pch = 16, cex = 0.28, col = grDevices::adjustcolor(ref_col, alpha.f = 0.65))
    for (i in seq_len(nrow(line_plot))) {
      pts <- s3d$xyz.convert(
        c(line_plot$query_x[[i]], line_plot$reference_x[[i]]),
        c(line_plot$query_y[[i]], line_plot$reference_y[[i]]),
        c(0, 8)
      )
      graphics::segments(
        pts$x[1], pts$y[1], pts$x[2], pts$y[2],
        col = grDevices::adjustcolor(line_col[[i]], alpha.f = 0.18),
        lwd = line_lwd[[i]]
      )
    }
    s3d$points3d(query_plot$x, query_plot$y, rep(0, nrow(query_plot)), pch = 16, cex = 0.18, col = grDevices::adjustcolor("#4d4d4d", alpha.f = 0.55))
    s3d$points3d(ref_plot$x, ref_plot$y, rep(8, nrow(ref_plot)), pch = 16, cex = 0.28, col = grDevices::adjustcolor(ref_col, alpha.f = 0.85))
    graphics::legend(
      "bottom",
      legend = c("query cells", names(label_cols)),
      col = c("#4d4d4d", unname(label_cols)),
      pch = 16,
      horiz = TRUE,
      inset = -0.02,
      xpd = TRUE,
      bty = "n",
      cex = 0.72
    )
    graphics::mtext(
      paste0(
        "Lines are sampled top Seurat anchor links per query cell; reference plane uses ",
        unique(ref_coords$coordinate_reduction),
        ", query plane uses ", unique(query_coords$coordinate_reduction), "."
      ),
      side = 1,
      line = 1.2,
      cex = 0.75
    )
  }
  draw_one(paste0(out_prefix, ".png"), "png")
  draw_one(paste0(out_prefix, ".pdf"), "pdf")
  draw_one(paste0(out_prefix, ".svg"), "svg")
  TRUE
}

main <- function() {
  opt <- parse_args()
  set.seed(as.integer(opt$seed))
  helper_path <- if (grepl("^/", opt$helper_script)) opt$helper_script else file.path(getwd(), opt$helper_script)
  helper <- load_transfer_helpers(helper_path)
  for (fn in c("default_studies", "resolve_path", "join_layers_if_needed", "ensure_log_normalized", "attach_reference_labels_from_tsv", "find_existing_column", "parse_gw_numeric")) {
    if (!exists(fn, envir = helper, inherits = FALSE)) stop("Missing helper function after loading script 13: ", fn, call. = FALSE)
  }

  anchors_dir <- file.path(opt$outdir, "seurat", "anchors")
  tables_dir <- file.path(opt$outdir, "tables")
  plots_dir <- file.path(opt$outdir, "plots")
  diagnostics_dir <- file.path(opt$outdir, "diagnostics")
  for (path in c(anchors_dir, tables_dir, plots_dir, diagnostics_dir)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }

  if (!file.exists(opt$reference)) stop("Missing Shi reference RDS: ", opt$reference, call. = FALSE)
  log_msg("Loading Shi reference: ", opt$reference)
  reference <- readRDS(opt$reference)
  if (!inherits(reference, "Seurat")) stop("Reference is not a Seurat object", call. = FALSE)
  if (!("RNA" %in% Seurat::Assays(reference))) stop("Shi reference lacks RNA assay", call. = FALSE)
  label_col <- helper$find_existing_column(
    reference@meta.data,
    "",
    c("shi_label", "shi_transfer_label", "Major types", "major_type", "major_cell_type", "cell_type", "celltype"),
    "label"
  )
  if (!nzchar(label_col)) {
    reference <- helper$attach_reference_labels_from_tsv(reference, opt$reference_labels_tsv)
    label_col <- helper$find_existing_column(
      reference@meta.data,
      "",
      c("shi_label", "shi_transfer_label", "Major types", "major_type", "major_cell_type", "cell_type", "celltype"),
      "label"
    )
  }
  week_col <- helper$find_existing_column(
    reference@meta.data,
    "",
    c("week_label", "shi_week_label", "shi_transfer_week_label", "sample_week_label", "GW", "gw", "week"),
    "week"
  )
  if (!nzchar(label_col) || !nzchar(week_col)) stop("Could not identify Shi reference label/week columns", call. = FALSE)
  keep <- !is.na(reference@meta.data[[label_col]]) & nzchar(as.character(reference@meta.data[[label_col]])) &
    !is.na(reference@meta.data[[week_col]]) & nzchar(as.character(reference@meta.data[[week_col]]))
  reference <- subset(reference, cells = colnames(reference)[keep])
  reference@meta.data[[label_col]] <- as.character(reference@meta.data[[label_col]])
  reference@meta.data[[week_col]] <- as.character(reference@meta.data[[week_col]])
  reference$shi_label <- reference@meta.data[[label_col]]
  reference$shi_week_label <- reference@meta.data[[week_col]]
  write_tsv(as.data.frame(table(reference$shi_label)), file.path(diagnostics_dir, "shi_reference_labels_used_by_anchor_projection.tsv"))
  write_tsv(as.data.frame(table(reference$shi_week_label)), file.path(diagnostics_dir, "shi_reference_weeks_used_by_anchor_projection.tsv"))

  studies <- helper$default_studies(opt$project_root)
  studies <- studies[studies$study_id %in% opt$study_id, , drop = FALSE]
  if (nrow(studies) == 0L) stop("No selected studies matched default_studies()", call. = FALSE)
  write_tsv(studies, file.path(tables_dir, "shi_anchor_projection_selected_studies.tsv"))

  diagnostics <- list()
  for (idx in seq_len(nrow(studies))) {
    study <- studies[idx, , drop = FALSE]
    study_id <- study$study_id[[1]]
    object_path <- helper$resolve_path(study$object_path[[1]], opt$project_root)
    if (!file.exists(object_path)) stop(study_id, " query object not found: ", object_path, call. = FALSE)
    log_msg("Loading query ", study_id, ": ", object_path)
    query <- readRDS(object_path)
    if (!inherits(query, "Seurat")) stop(study_id, " query is not a Seurat object", call. = FALSE)
    if (!("RNA" %in% Seurat::Assays(query))) stop(study_id, " query lacks RNA assay", call. = FALSE)
    if (!(study$reduction[[1]] %in% names(query@reductions))) {
      stop(study_id, " query missing UMAP reduction: ", study$reduction[[1]], call. = FALSE)
    }

    query <- helper$join_layers_if_needed(query, "RNA")
    reference <- helper$join_layers_if_needed(reference, "RNA")
    query <- helper$ensure_log_normalized(query, "RNA")
    reference <- helper$ensure_log_normalized(reference, "RNA")

    shared_features <- intersect(rownames(reference[["RNA"]]), rownames(query[["RNA"]]))
    min_shared <- as.integer(opt$min_shared_features)
    dims_n <- as.integer(opt$dims)
    if (length(shared_features) < min_shared) {
      stop(study_id, " has too few shared features with Shi reference: ", length(shared_features), " < ", min_shared, call. = FALSE)
    }
    if (length(shared_features) <= dims_n) {
      stop(study_id, " has only ", length(shared_features), " shared features; dims=", dims_n, " requires more", call. = FALSE)
    }
    dims_use <- seq_len(dims_n)

    log_msg("ScaleData/RunPCA reference for ", study_id, " using ", length(shared_features), " shared features")
    reference_for_transfer <- Seurat::ScaleData(reference, assay = "RNA", features = shared_features, verbose = FALSE)
    reference_for_transfer <- Seurat::RunPCA(reference_for_transfer, assay = "RNA", features = shared_features, npcs = dims_n, verbose = FALSE)

    log_msg("FindTransferAnchors whole-Shi for ", study_id)
    anchors <- Seurat::FindTransferAnchors(
      reference = reference_for_transfer,
      query = query,
      normalization.method = opt$normalization_method,
      reference.assay = "RNA",
      query.assay = "RNA",
      features = shared_features,
      reference.reduction = "pca",
      reduction = "pcaproject",
      dims = dims_use,
      verbose = TRUE
    )
    anchor_df <- anchor_pair_table(anchors, reference_for_transfer, query, label_col = "shi_label", week_col = "shi_week_label")
    top_df <- top_anchor_per_query(anchor_df)
    top_df$plot_include <- TRUE
    if (study_id == "varela_div90" && "query_cluster" %in% colnames(top_df)) {
      top_df$plot_include <- !(as.character(top_df$query_cluster) %in% c("6", "7", "6.0", "7.0"))
    }

    ref_coords <- coordinate_table(reference_for_transfer, c("umap", "UMAP", "tsne", "pca"), label = "Shi reference")
    ref_coords$shi_label <- as.character(reference_for_transfer@meta.data[ref_coords$cell_id, "shi_label"])
    ref_coords$shi_week_label <- as.character(reference_for_transfer@meta.data[ref_coords$cell_id, "shi_week_label"])
    query_coords <- coordinate_table(query, c(study$reduction[[1]], "umap", "UMAP", "pca"), label = paste0(study_id, " query"))
    query_coords$query_cluster <- if ("seurat_clusters" %in% colnames(query@meta.data)) as.character(query@meta.data[query_coords$cell_id, "seurat_clusters"]) else ""
    query_coords$plot_include <- TRUE
    if (study_id == "varela_div90") {
      query_coords$plot_include <- !(as.character(query_coords$query_cluster) %in% c("6", "7", "6.0", "7.0"))
    }

    if (to_bool(opt$save_anchor_rds, default = TRUE)) {
      anchor_path <- file.path(anchors_dir, paste0(study_id, "_shi_full_transfer_anchors.rds"))
      log_msg("Saving anchor object: ", anchor_path)
      saveRDS(anchors, anchor_path, compress = "gzip")
    }
    write_tsv_gz(anchor_df, file.path(tables_dir, paste0(study_id, "_shi_full_anchor_pairs.tsv.gz")))
    write_tsv_gz(top_df, file.path(tables_dir, paste0(study_id, "_shi_full_top_anchor_per_query.tsv.gz")))
    write_tsv_gz(ref_coords, file.path(tables_dir, paste0(study_id, "_shi_reference_coordinates_for_anchor_plot.tsv.gz")))
    write_tsv_gz(query_coords, file.path(tables_dir, paste0(study_id, "_query_coordinates_for_anchor_plot.tsv.gz")))

    plot_ok <- plot_anchor_projection(
      study_id,
      ref_coords = ref_coords,
      query_coords = query_coords,
      links = top_df,
      out_prefix = file.path(plots_dir, paste0(study_id, "_shi_full_anchor_projection_static")),
      opt = opt
    )

    diagnostics[[length(diagnostics) + 1L]] <- data.frame(
      study_id = study_id,
      status = "ok",
      query_object_path = object_path,
      reference_object_path = opt$reference,
      n_query_cells = ncol(query),
      n_reference_cells = ncol(reference_for_transfer),
      n_shared_features = length(shared_features),
      dims = paste0("1:", dims_n),
      n_anchors = nrow(anchor_df),
      n_unique_query_anchor_cells = length(unique(anchor_df$query_cell_id)),
      n_unique_reference_anchor_cells = length(unique(anchor_df$reference_cell_id)),
      reference_coordinate_reduction = unique(ref_coords$coordinate_reduction)[[1]],
      query_coordinate_reduction = unique(query_coords$coordinate_reduction)[[1]],
      div90_plot_filter = ifelse(study_id == "varela_div90", "static_plot_excludes_query_clusters_6_7", "none"),
      anchor_rds_saved = to_bool(opt$save_anchor_rds, default = TRUE),
      static_plot_rendered = plot_ok,
      stringsAsFactors = FALSE
    )
    rm(query, reference_for_transfer, anchors, anchor_df, top_df, ref_coords, query_coords)
    gc()
  }
  diagnostics <- do.call(rbind, diagnostics)
  write_tsv(diagnostics, file.path(diagnostics_dir, "shi_anchor_projection_diagnostics.tsv"))
  write_tsv(diagnostics, file.path(tables_dir, "shi_anchor_projection_diagnostics.tsv"))
  readme <- c(
    "# Shi Query Anchor Projection Assets",
    "",
    "This run saves true Seurat whole-Shi TransferAnchorSet assets for Varela DIV30 and DIV90.",
    "It is separate from the finalized plot-ready prediction tables and does not overwrite them.",
    "",
    "Important:",
    "- Anchor objects and full anchor-pair tables are generated from FindTransferAnchors.",
    "- Static plots sample top anchor links for readability.",
    "- DIV90 static plots exclude current query clusters 6/7 for visualization consistency, but the saved anchor object and full anchor table include all query cells.",
    "- These are anchor links, not TransferData score-routing summaries.",
    "",
    paste0("Reference RDS: ", opt$reference),
    paste0("Run label: ", opt$run_label)
  )
  writeLines(readme, file.path(opt$outdir, "README.md"))
  print(diagnostics)
  log_msg("Finished Shi anchor projection assets: ", opt$outdir)
}

main()
