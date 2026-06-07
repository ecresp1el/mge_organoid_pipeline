#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  out <- list(
    study_table = "",
    gene_table = "",
    outdir = "",
    study_id = character(),
    project_root = Sys.getenv("PROJECT_ROOT", "/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder")
  )
  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    value <- if (i < length(args)) args[[i + 1L]] else ""
    if (key == "--study-table") {
      out$study_table <- value; i <- i + 2L
    } else if (key == "--gene-table") {
      out$gene_table <- value; i <- i + 2L
    } else if (key == "--outdir") {
      out$outdir <- value; i <- i + 2L
    } else if (key == "--project-root") {
      out$project_root <- value; i <- i + 2L
    } else if (key == "--study-id") {
      out$study_id <- c(out$study_id, unlist(strsplit(value, "[,;[:space:]]+"))); i <- i + 2L
    } else {
      stop("Unknown argument: ", key)
    }
  }
  if (!nzchar(out$study_table)) stop("--study-table is required")
  if (!nzchar(out$gene_table)) stop("--gene-table is required")
  if (!nzchar(out$outdir)) stop("--outdir is required")
  out$study_id <- out$study_id[nzchar(out$study_id)]
  out
}

resolve_path <- function(path, project_root) {
  if (!nzchar(path)) return("")
  if (grepl("^/", path)) return(path)
  file.path(project_root, path)
}

log_msg <- function(...) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste0(..., collapse = "")))
  flush.console()
}

gene_aliases <- list(
  "NKX2-1" = c("NKX2-1", "NKX2.1"),
  "NKX6-2" = c("NKX6-2", "NKX6.2")
)

match_genes <- function(genes, features) {
  feature_set <- unique(as.character(features))
  upper <- toupper(feature_set)
  upper_tab <- table(upper)
  upper_lookup <- setNames(feature_set[!duplicated(upper)], upper[!duplicated(upper)])
  rows <- list()
  matches <- character()
  for (gene in genes) {
    candidates <- gene_aliases[[gene]]
    if (is.null(candidates)) candidates <- gene
    matched <- ""
    match_type <- "missing"
    for (candidate in candidates) {
      key <- toupper(candidate)
      if (candidate %in% feature_set) {
        matched <- candidate
        match_type <- if (identical(candidate, gene)) "exact" else "alias_exact"
        break
      }
      upper_count <- if (key %in% names(upper_tab)) as.integer(upper_tab[[key]]) else 0L
      if (upper_count > 1L) {
        match_type <- "ambiguous_case_insensitive"
      } else if (key %in% names(upper_lookup)) {
        matched <- upper_lookup[[key]]
        match_type <- if (identical(candidate, gene)) "case_insensitive" else "alias_case_insensitive"
        break
      }
    }
    if (nzchar(matched)) matches[[gene]] <- matched
    rows[[length(rows) + 1L]] <- data.frame(
      gene = gene,
      matched_feature = matched,
      matched = nzchar(matched),
      match_type = match_type,
      aliases_considered = paste(candidates, collapse = ","),
      stringsAsFactors = FALSE
    )
  }
  list(matches = matches, table = do.call(rbind, rows))
}

get_matrix_layer <- function(obj, assay, layer) {
  if (!assay %in% names(obj@assays)) {
    stop("Assay not found: ", assay)
  }
  assay_obj <- obj[[assay]]
  if (inherits(assay_obj, "Assay5")) {
    layers <- Layers(assay_obj)
    layer_hits <- layers[layers == layer]
    if (length(layer_hits) == 0L) {
      layer_hits <- layers[grepl(paste0("^", layer, "(\\.|$)"), layers)]
    }
    if (length(layer_hits) != 1L) {
      stop("Expected exactly one Assay5 layer for ", assay, "/", layer, "; found: ", paste(layer_hits, collapse = ","))
    }
    return(LayerData(assay_obj, layer = layer_hits[[1L]]))
  }
  GetAssayData(obj, assay = assay, slot = layer)
}

write_tsv_gz <- function(df, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  con <- gzfile(path, open = "wt")
  on.exit(close(con), add = TRUE)
  write.table(df, con, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
}

export_one_study <- function(row, genes, outdir, project_root) {
  study_id <- row$study_id[[1]]
  label <- row$study_label[[1]]
  object_path <- resolve_path(row$seurat_path[[1]], project_root)
  assay <- row$assay[[1]]
  layer <- row$expression_layer[[1]]
  reduction <- row$reduction[[1]]
  sample_col <- row$sample_col[[1]]
  cluster_col <- row$cluster_col[[1]]
  out_path <- file.path(outdir, paste0(study_id, "_marker_expression.tsv.gz"))
  match_path <- file.path(outdir, paste0(study_id, "_gene_matches.tsv"))

  log_msg("Loading ", study_id, ": ", object_path)
  obj <- readRDS(object_path)
  expr <- get_matrix_layer(obj, assay = assay, layer = layer)
  features <- rownames(expr)
  gene_match <- match_genes(genes, features)
  missing <- setdiff(genes, names(gene_match$matches))
  if (length(missing) > 0L) {
    stop(study_id, " missing marker genes: ", paste(missing, collapse = ","))
  }
  emb <- Embeddings(obj, reduction = reduction)
  meta <- obj@meta.data
  cells <- Reduce(intersect, list(colnames(expr), rownames(emb), rownames(meta)))
  if (length(cells) == 0L) stop(study_id, " has no common cells across expression, UMAP, and metadata")

  expr_sub <- expr[unname(gene_match$matches[genes]), cells, drop = FALSE]
  expr_dense <- as.matrix(t(expr_sub))
  colnames(expr_dense) <- genes
  emb_sub <- emb[cells, , drop = FALSE]
  meta_sub <- meta[cells, , drop = FALSE]
  sample <- if (nzchar(sample_col) && sample_col %in% colnames(meta_sub)) as.character(meta_sub[[sample_col]]) else rep("", length(cells))
  cluster <- if (nzchar(cluster_col) && cluster_col %in% colnames(meta_sub)) as.character(meta_sub[[cluster_col]]) else rep("", length(cells))

  out <- data.frame(
    cell_id = cells,
    study_id = study_id,
    study_label = label,
    sample = sample,
    cluster = cluster,
    umap_1 = emb_sub[, 1],
    umap_2 = emb_sub[, 2],
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  out <- cbind(out, as.data.frame(expr_dense, check.names = FALSE))
  write_tsv_gz(out, out_path)
  write.table(gene_match$table, match_path, sep = "\t", quote = FALSE, row.names = FALSE)
  data.frame(
    study_id = study_id,
    status = "ok",
    output_path = out_path,
    n_cells = nrow(out),
    n_genes = length(genes),
    stringsAsFactors = FALSE
  )
}

main <- function() {
  opt <- parse_args()
  studies <- read.delim(opt$study_table, stringsAsFactors = FALSE, check.names = FALSE)
  genes <- read.delim(opt$gene_table, stringsAsFactors = FALSE, check.names = FALSE)
  genes <- genes$gene
  studies <- studies[studies$include_in_first_plot %in% c(TRUE, "TRUE", "true", "1"), , drop = FALSE]
  if (length(opt$study_id) > 0L) {
    studies <- studies[studies$study_id %in% opt$study_id, , drop = FALSE]
  }
  if (nrow(studies) == 0L) stop("No studies selected for export")
  dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)
  summaries <- list()
  for (idx in seq_len(nrow(studies))) {
    summaries[[length(summaries) + 1L]] <- export_one_study(studies[idx, , drop = FALSE], genes, opt$outdir, opt$project_root)
    gc()
  }
  summary <- do.call(rbind, summaries)
  write.table(
    summary,
    file.path(dirname(opt$outdir), "cross_study_marker_expression_seurat_export_summary.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  print(summary)
}

main()
