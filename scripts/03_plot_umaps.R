#!/usr/bin/env Rscript

# Plot UMAPs for Walsh day-75 final object and Bershteyn Seurat object.
# Saves PDF and PNG plots under PROJECT_ROOT/results/<dataset>/plots.

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})

log_msg <- function(...) {
  msg <- sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste0(..., collapse = " "))
  cat(msg, "\n")
  flush.console()
}

save_plot <- function(p, out_base) {
  ggsave(paste0(out_base, ".pdf"), plot = p, width = 8, height = 6, dpi = 300)
  ggsave(paste0(out_base, ".png"), plot = p, width = 8, height = 6, dpi = 300)
}

plot_umaps <- function(obj, out_dir, dataset_label) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  # cluster
  p1 <- DimPlot(obj, reduction = "umap", group.by = "seurat_clusters", label = TRUE, repel = TRUE) +
    ggtitle(paste0(dataset_label, " UMAP by cluster"))
  save_plot(p1, file.path(out_dir, "umap_by_cluster"))
  # sample_id if present
  if ("sample_id" %in% colnames(obj@meta.data)) {
    p2 <- DimPlot(obj, reduction = "umap", group.by = "sample_id", label = FALSE) +
      ggtitle(paste0(dataset_label, " UMAP by sample_id"))
    save_plot(p2, file.path(out_dir, "umap_by_sample"))
  }
  # domain if present
  if ("domain" %in% colnames(obj@meta.data)) {
    p3 <- DimPlot(obj, reduction = "umap", group.by = "domain", label = FALSE) +
      ggtitle(paste0(dataset_label, " UMAP by domain"))
    save_plot(p3, file.path(out_dir, "umap_by_domain"))
  }
}

load_object <- function(path) {
  if (grepl("\\.gz$", path, ignore.case = TRUE)) {
    con <- gzcon(file(path, "rb"))
    out <- tryCatch(readRDS(con), error = function(e) NULL)
    close(con)
    if (!is.null(out)) return(out)
    tmp <- tempfile(fileext = ".rds")
    system2("gunzip", args = c("-c", shQuote(path)), stdout = tmp)
    on.exit(unlink(tmp), add = TRUE)
    return(readRDS(tmp))
  }
  readRDS(path)
}

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  project_root <- Sys.getenv("PROJECT_ROOT", "")
  if (project_root == "") stop("PROJECT_ROOT env must be set")

  datasets <- c("walsh", "bershteyn_2025")
  for (ds in datasets) {
    if (ds == "walsh") {
      in_path <- file.path(project_root, "results", "walsh_day75", "walsh_day75_final.rds")
      out_dir <- file.path(project_root, "results", "walsh_day75", "plots")
      label <- "Walsh day75"
    } else {
      canonical_path <- file.path(project_root, "results", "bershteyn_2025", "bershteyn_2025_seurat.rds")
      raw_path <- file.path(project_root, "data", "raw", "bershteyn_2025_geo_files", "suppl", "GSE283775_Seurat_scRNA_seq.rds.gz")
      in_path <- if (file.exists(canonical_path)) canonical_path else raw_path
      out_dir <- file.path(project_root, "results", "bershteyn_2025", "plots")
      label <- "Bershteyn 2025"
    }
    if (!file.exists(in_path)) {
      log_msg("Skipping", ds, "- missing file:", in_path)
      next
    }
    log_msg("Loading", ds, "object from", in_path)
    obj <- load_object(in_path)
    if (!"umap" %in% names(obj@reductions)) {
      log_msg("No UMAP found for", ds, "- skipping plots")
      next
    }
    log_msg("Plotting", ds)
    plot_umaps(obj, out_dir, label)
    log_msg("Done", ds, "->", out_dir)
  }
}

main()
