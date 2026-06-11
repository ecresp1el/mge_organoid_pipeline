#!/usr/bin/env Rscript

# Build the first URD lineage tree from an existing DIV30 URD pseudotime object.
#
# Programmatic flow:
#   1. Read a completed `div30_first_urd_object.rds` containing diffusion maps
#      and flood pseudotime.
#   2. Define the root cells from `paper_cluster_annotation == "Radial glia"`.
#   3. Define candidate terminal tips from paper/manual annotations:
#        1 = SST+ cIN
#        2 = PV neuron precursor
#        3 = MGE subpallial neurons
#      Numeric tip IDs are intentional because URD's tree internals expect
#      numeric-like segment names.
#   4. Fit the pseudotime transition logistic curve.
#   5. Weight the transition matrix by pseudotime.
#   6. Simulate random walks from each candidate tip back to the radial glia root.
#   7. Process random walks into tip visitation frequencies.
#   8. Run `buildTree()` to ask whether neuronal populations separate into
#      distinct branches in the diffusion manifold.
#   9. Save the tree object, parameter tables, tree status tables, branch-gene
#      tables, and tree plots.

log_msg <- function(...) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste0(..., collapse = " ")))
  flush.console()
}

parse_args <- function(args) {
  out <- list(
    `urd-rds` = NULL,
    outdir = NULL,
    `annotation-col` = "paper_cluster_annotation",
    `pseudotime-name` = "",
    `root-label` = "Radial glia",
    `tip-labels` = "SST+ cIN,PV neuron precursor,MGE subpallial neurons",
    `seed` = "7",
    `optimal-cells-forward` = "100",
    `max-cells-back` = "40",
    `pseudotime-direction` = "<",
    `n-per-tip` = "5000",
    `root-visits` = "1",
    `max-steps` = "",
    `process-n-subsample` = "10",
    `cells-per-pseudotime-bin` = "80",
    `bins-per-pseudotime-window` = "5",
    `minimum-visits` = "10",
    `visit-threshold` = "0.7",
    `divergence-method` = "ks",
    `p-thresh` = "0.01",
    `min-cells-per-segment` = "1",
    `min-pseudotime-per-segment` = "0.01",
    `dendro-node-size` = "100",
    `top-n-branch-genes` = "50",
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
  cat(paste(
    "Usage:",
    "  Rscript scripts/16_div30_urd_build_lineage_tree.R --urd-rds <div30_first_urd_object.rds> --outdir <tree_dir>",
    "",
    "Required inputs:",
    "  --urd-rds     Existing URD object with diffusion map and flood pseudotime.",
    "  --outdir      Directory for tree outputs.",
    "",
    "Main outputs:",
    "  div30_urd_lineage_tree_object.rds",
    "  tables/lineage_tree_parameters.tsv",
    "  tables/tree_tip_mapping.tsv",
    "  tables/tree_status.tsv",
    "  tables/tree_segment_joins.tsv",
    "  tables/branch_specific_genes.tsv",
    "  plots/pseudotime_logistic.png",
    "  plots/urd_tree_annotation.png",
    "  plots/urd_tree_pseudotime.png",
    "  urd_lineage_tree_report.md",
    sep = "\n"
  ))
}

required <- c("URD", "Matrix")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) stop("Missing required R package(s): ", paste(missing, collapse = ", "), call. = FALSE)
suppressPackageStartupMessages({
  library(URD)
  library(Matrix)
})

as_int <- function(x, name) {
  value <- suppressWarnings(as.integer(x))
  if (is.na(value)) stop(name, " must be an integer; got ", x, call. = FALSE)
  value
}

as_num <- function(x, name) {
  value <- suppressWarnings(as.numeric(x))
  if (is.na(value)) stop(name, " must be numeric; got ", x, call. = FALSE)
  value
}

write_tsv <- function(x, path) {
  write.table(x, path, sep = "\t", row.names = FALSE, quote = FALSE, na = "")
}

extract_pseudotime_name <- function(object, requested) {
  if (!is.data.frame(object@pseudotime) && !is.matrix(object@pseudotime)) {
    stop("Expected object@pseudotime to be a data.frame or matrix.", call. = FALSE)
  }
  if (ncol(object@pseudotime) == 0) stop("URD object has no pseudotime columns.", call. = FALSE)
  if (!nzchar(requested)) return(colnames(object@pseudotime)[[1]])
  if (!(requested %in% colnames(object@pseudotime))) {
    stop("Pseudotime column not found: ", requested, ". Available: ", paste(colnames(object@pseudotime), collapse = ", "), call. = FALSE)
  }
  requested
}

make_tip_mapping <- function(tip_labels) {
  labels <- trimws(strsplit(tip_labels, ",", fixed = TRUE)[[1]])
  labels <- labels[nzchar(labels)]
  if (length(labels) < 2) stop("At least two tip labels are required to build a tree.", call. = FALSE)
  data.frame(
    tip_id = as.character(seq_along(labels)),
    paper_cluster_annotation = labels,
    stringsAsFactors = FALSE
  )
}

add_tree_groups <- function(object, annotation_col, root_label, tip_mapping) {
  if (!(annotation_col %in% colnames(object@meta))) stop("Missing annotation column: ", annotation_col, call. = FALSE)
  meta <- object@meta
  root_cells <- rownames(meta)[meta[[annotation_col]] == root_label]
  if (length(root_cells) == 0) stop("No root cells found for label: ", root_label, call. = FALSE)

  tip_id_by_label <- setNames(tip_mapping$tip_id, tip_mapping$paper_cluster_annotation)
  tip_ids <- unname(tip_id_by_label[as.character(meta[[annotation_col]])])
  tip_ids[is.na(tip_ids)] <- NA_character_
  names(tip_ids) <- rownames(meta)

  if (!all(tip_mapping$paper_cluster_annotation %in% meta[[annotation_col]])) {
    missing_labels <- setdiff(tip_mapping$paper_cluster_annotation, unique(meta[[annotation_col]]))
    stop("No cells found for tip label(s): ", paste(missing_labels, collapse = ", "), call. = FALSE)
  }

  object@group.ids[rownames(meta), "paper_tree_tip_id"] <- tip_ids[rownames(meta)]
  object@meta[rownames(meta), "paper_tree_tip_id"] <- tip_ids[rownames(meta)]
  object@meta[rownames(meta), "paper_tree_root_candidate"] <- rownames(meta) %in% root_cells
  list(object = object, root_cells = root_cells, tip_ids = tip_mapping$tip_id)
}

plot_logistic <- function(object, pseudotime_name, optimal_cells_forward, max_cells_back, direction, path) {
  png(path, width = 1800, height = 1400, res = 220)
  on.exit(dev.off(), add = TRUE)
  params <- pseudotimeDetermineLogistic(
    object,
    pseudotime = pseudotime_name,
    optimal.cells.forward = optimal_cells_forward,
    max.cells.back = max_cells_back,
    pseudotime.direction = direction,
    do.plot = TRUE,
    print.values = TRUE
  )
  title(main = "URD pseudotime logistic transition weighting")
  params
}

plot_tree_safe <- function(object, label, path, title) {
  ok <- FALSE
  png(path, width = 2200, height = 1600, res = 220)
  tryCatch(
    {
      plotTree(object, label = label, title = title)
      ok <<- TRUE
    },
    error = function(e) {
      plot.new()
      text(0.5, 0.5, paste("plotTree failed:", conditionMessage(e)), cex = 0.8)
    },
    finally = dev.off()
  )
  ok
}

tree_status <- function(object, tip_mapping) {
  joins <- object@tree$segment.joins
  n_joins <- if (is.null(joins)) 0L else nrow(joins)
  segments <- object@tree$segments
  data.frame(
    tree_slot_length = length(object@tree),
    n_requested_tips = nrow(tip_mapping),
    n_segment_joins = n_joins,
    n_segments = if (is.null(segments)) 0L else length(segments),
    has_distinct_branching = n_joins > 0,
    tips = paste(tip_mapping$tip_id, tip_mapping$paper_cluster_annotation, sep = "=", collapse = "; "),
    stringsAsFactors = FALSE
  )
}

segment_join_table <- function(object) {
  joins <- object@tree$segment.joins
  if (is.null(joins) || nrow(joins) == 0) {
    return(data.frame(child_1 = character(), child_2 = character(), parent = character(), pseudotime = numeric()))
  }
  out <- as.data.frame(joins, stringsAsFactors = FALSE)
  colnames(out) <- sub("\\.", "_", colnames(out))
  out
}

tip_composition <- function(object, annotation_col, tip_mapping) {
  rows <- lapply(seq_len(nrow(tip_mapping)), function(i) {
    id <- tip_mapping$tip_id[[i]]
    cells <- rownames(object@meta)[object@meta$paper_tree_tip_id %in% id]
    data.frame(
      tip_id = id,
      expected_annotation = tip_mapping$paper_cluster_annotation[[i]],
      n_tip_cells = length(cells),
      median_pseudotime = median(object@tree$pseudotime[cells], na.rm = TRUE),
      annotation_composition = paste(capture.output(print(table(object@meta[cells, annotation_col]), quote = FALSE)), collapse = "; "),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

branch_genes <- function(object, tip_mapping, top_n) {
  genes <- intersect(object@var.genes, rownames(object@logupx.data))
  if (length(genes) == 0) genes <- rownames(object@logupx.data)
  rows <- lapply(seq_len(nrow(tip_mapping)), function(i) {
    tip_id <- tip_mapping$tip_id[[i]]
    cells_1 <- intersect(rownames(object@meta)[object@meta$paper_tree_tip_id == tip_id], colnames(object@logupx.data))
    cells_2 <- intersect(
      rownames(object@meta)[!is.na(object@meta$paper_tree_tip_id) & object@meta$paper_tree_tip_id != tip_id],
      colnames(object@logupx.data)
    )
    if (length(cells_1) == 0 || length(cells_2) == 0) return(NULL)
    mat_1 <- object@logupx.data[genes, cells_1, drop = FALSE]
    mat_2 <- object@logupx.data[genes, cells_2, drop = FALSE]

    # Direct branch contrast keeps the table independent of URD marker helpers,
    # which can be brittle after buildTree() subsets internal visitation slots.
    mean_1 <- Matrix::rowMeans(mat_1)
    mean_2 <- Matrix::rowMeans(mat_2)
    frac_1 <- Matrix::rowSums(mat_1 > 0) / ncol(mat_1)
    frac_2 <- Matrix::rowSums(mat_2 > 0) / ncol(mat_2)
    log2_fc <- log2((2^mean_1 - 1 + 1e-6) / (2^mean_2 - 1 + 1e-6))
    out <- data.frame(
      tip_id = tip_id,
      tip_annotation = tip_mapping$paper_cluster_annotation[[i]],
      gene = genes,
      mean_logupx_tip = as.numeric(mean_1),
      mean_logupx_other_tips = as.numeric(mean_2),
      frac_expressed_tip = as.numeric(frac_1),
      frac_expressed_other_tips = as.numeric(frac_2),
      log2_fold_change_tip_vs_other_tips = as.numeric(log2_fc),
      expression_fraction_delta = as.numeric(frac_1 - frac_2),
      n_tip_cells = length(cells_1),
      n_other_tip_cells = length(cells_2),
      stringsAsFactors = FALSE
    )
    out$branch_specificity_score <- out$log2_fold_change_tip_vs_other_tips * pmax(out$expression_fraction_delta, 0)
    out <- out[is.finite(out$branch_specificity_score), , drop = FALSE]
    out <- out[order(out$branch_specificity_score, out$log2_fold_change_tip_vs_other_tips, decreasing = TRUE), , drop = FALSE]
    head(out, top_n)
  })
  out <- do.call(rbind, rows)
  if (is.null(out)) {
    return(data.frame(status = "not_available", note = "No branch marker genes could be computed.", stringsAsFactors = FALSE))
  }
  rownames(out) <- NULL
  out
}

write_report <- function(path, cfg, status, tip_mapping, tip_comp, joins, branch_gene_path) {
  lines <- c(
    "# DIV30 URD Lineage Tree Report",
    "",
    paste0("- Input URD object: `", cfg$urd_rds, "`"),
    paste0("- Output tree object: `", file.path(cfg$outdir, "div30_urd_lineage_tree_object.rds"), "`"),
    paste0("- Root annotation: `", cfg$root_label, "`"),
    paste0("- Tip group column written to URD: `paper_tree_tip_id`"),
    paste0("- Pseudotime column: `", cfg$pseudotime_name, "`"),
    "",
    "## Tip Mapping",
    "",
    paste(capture.output(print(tip_mapping, row.names = FALSE)), collapse = "\n"),
    "",
    "## Tree Status",
    "",
    paste(capture.output(print(status, row.names = FALSE)), collapse = "\n"),
    "",
    "## Tip Composition",
    "",
    paste(capture.output(print(tip_comp, row.names = FALSE)), collapse = "\n"),
    "",
    "## Segment Joins",
    "",
    if (nrow(joins) > 0) paste(capture.output(print(joins, row.names = FALSE)), collapse = "\n") else "No segment joins were produced.",
    "",
    "## Branch Genes",
    "",
    paste0("Branch-specific genes were written to `", branch_gene_path, "`. They are computed after the tree stage by comparing each requested tip population against the other requested tip populations over variable genes."),
    "",
    "## Core URD Calls",
    "",
    "- `pseudotimeDetermineLogistic()`",
    "- `pseudotimeWeightTransitionMatrix()`",
    "- `simulateRandomWalksFromTips()`",
    "- `processRandomWalksFromTips()`",
    "- `buildTree()`"
  )
  writeLines(lines, path)
}

opt <- parse_args(commandArgs(trailingOnly = TRUE))
if (isTRUE(opt$help)) {
  print_usage()
  quit(save = "no", status = 0)
}
if (is.null(opt$`urd-rds`) || !nzchar(opt$`urd-rds`)) stop("--urd-rds is required", call. = FALSE)
if (is.null(opt$outdir) || !nzchar(opt$outdir)) stop("--outdir is required", call. = FALSE)

cfg <- list(
  urd_rds = opt$`urd-rds`,
  outdir = opt$outdir,
  annotation_col = opt$`annotation-col`,
  pseudotime_name = opt$`pseudotime-name`,
  root_label = opt$`root-label`,
  tip_labels = opt$`tip-labels`,
  seed = as_int(opt$seed, "seed"),
  optimal_cells_forward = as_int(opt$`optimal-cells-forward`, "optimal-cells-forward"),
  max_cells_back = as_int(opt$`max-cells-back`, "max-cells-back"),
  pseudotime_direction = opt$`pseudotime-direction`,
  n_per_tip = as_int(opt$`n-per-tip`, "n-per-tip"),
  root_visits = as_int(opt$`root-visits`, "root-visits"),
  max_steps = opt$`max-steps`,
  process_n_subsample = as_int(opt$`process-n-subsample`, "process-n-subsample"),
  cells_per_pseudotime_bin = as_int(opt$`cells-per-pseudotime-bin`, "cells-per-pseudotime-bin"),
  bins_per_pseudotime_window = as_int(opt$`bins-per-pseudotime-window`, "bins-per-pseudotime-window"),
  minimum_visits = as_int(opt$`minimum-visits`, "minimum-visits"),
  visit_threshold = as_num(opt$`visit-threshold`, "visit-threshold"),
  divergence_method = opt$`divergence-method`,
  p_thresh = as_num(opt$`p-thresh`, "p-thresh"),
  min_cells_per_segment = as_int(opt$`min-cells-per-segment`, "min-cells-per-segment"),
  min_pseudotime_per_segment = as_num(opt$`min-pseudotime-per-segment`, "min-pseudotime-per-segment"),
  dendro_node_size = as_int(opt$`dendro-node-size`, "dendro-node-size"),
  top_n_branch_genes = as_int(opt$`top-n-branch-genes`, "top-n-branch-genes")
)

table_dir <- file.path(cfg$outdir, "tables")
plot_dir <- file.path(cfg$outdir, "plots")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

set.seed(cfg$seed)
log_msg("Reading URD object: ", cfg$urd_rds)
urd <- readRDS(cfg$urd_rds)
if (!inherits(urd, "URD")) stop("Input is not an URD object: ", cfg$urd_rds, call. = FALSE)
cfg$pseudotime_name <- extract_pseudotime_name(urd, cfg$pseudotime_name)
if (!all(rownames(urd@meta) %in% rownames(urd@pseudotime))) stop("Metadata and pseudotime cell names do not align.", call. = FALSE)

tip_mapping <- make_tip_mapping(cfg$tip_labels)
grouped <- add_tree_groups(urd, cfg$annotation_col, cfg$root_label, tip_mapping)
urd <- grouped$object
root_cells <- grouped$root_cells
tip_ids <- grouped$tip_ids
cfg$max_steps <- if (nzchar(cfg$max_steps)) as_int(cfg$max_steps, "max-steps") else ncol(urd@logupx.data)

parameters <- data.frame(
  parameter = names(cfg),
  value = vapply(cfg, as.character, character(1)),
  stringsAsFactors = FALSE
)
parameters <- rbind(
  parameters,
  data.frame(parameter = "n_root_cells", value = as.character(length(root_cells)), stringsAsFactors = FALSE),
  data.frame(parameter = "n_cells", value = as.character(nrow(urd@meta)), stringsAsFactors = FALSE),
  data.frame(parameter = "n_genes", value = as.character(nrow(urd@logupx.data)), stringsAsFactors = FALSE),
  data.frame(parameter = "n_variable_genes", value = as.character(length(urd@var.genes)), stringsAsFactors = FALSE)
)

write_tsv(parameters, file.path(table_dir, "lineage_tree_parameters.tsv"))
write_tsv(tip_mapping, file.path(table_dir, "tree_tip_mapping.tsv"))

log_msg("Loading tip cells into URD tree slot")
urd <- loadTipCells(urd, tips = "paper_tree_tip_id")

log_msg("Fitting pseudotime logistic transition model")
logistic_params <- plot_logistic(
  urd,
  cfg$pseudotime_name,
  cfg$optimal_cells_forward,
  cfg$max_cells_back,
  cfg$pseudotime_direction,
  file.path(plot_dir, "pseudotime_logistic.png")
)
write_tsv(data.frame(parameter = names(logistic_params), value = unlist(logistic_params), stringsAsFactors = FALSE), file.path(table_dir, "pseudotime_logistic_parameters.tsv"))

log_msg("Weighting transition matrix by pseudotime")
transition_matrix <- pseudotimeWeightTransitionMatrix(
  urd,
  pseudotime = cfg$pseudotime_name,
  logistic.params = logistic_params,
  pseudotime.direction = cfg$pseudotime_direction,
  verbose = TRUE
)

log_msg("Simulating random walks from tips: ", paste(tip_ids, collapse = ", "))
walks <- simulateRandomWalksFromTips(
  urd,
  tip.group.id = "paper_tree_tip_id",
  root.cells = root_cells,
  transition.matrix = transition_matrix,
  n.per.tip = cfg$n_per_tip,
  root.visits = cfg$root_visits,
  max.steps = cfg$max_steps,
  verbose = TRUE
)
saveRDS(walks, file.path(cfg$outdir, "div30_urd_tip_random_walks.rds"))

log_msg("Processing random walks into tip visitation frequencies")
urd <- processRandomWalksFromTips(
  urd,
  walks.list = walks,
  n.subsample = cfg$process_n_subsample,
  verbose = TRUE
)
rm(walks, transition_matrix)
gc()

log_msg("Building URD lineage tree")
breakpoint_dir <- file.path(plot_dir, "breakpoint_decisions")
dir.create(breakpoint_dir, recursive = TRUE, showWarnings = FALSE)
urd <- buildTree(
  urd,
  pseudotime = cfg$pseudotime_name,
  tips.use = tip_ids,
  divergence.method = cfg$divergence_method,
  weighted.fusion = TRUE,
  use.only.original.tips = TRUE,
  cells.per.pseudotime.bin = cfg$cells_per_pseudotime_bin,
  bins.per.pseudotime.window = cfg$bins_per_pseudotime_window,
  minimum.visits = cfg$minimum_visits,
  visit.threshold = cfg$visit_threshold,
  save.breakpoint.plots = breakpoint_dir,
  save.all.breakpoint.info = TRUE,
  p.thresh = cfg$p_thresh,
  min.cells.per.segment = cfg$min_cells_per_segment,
  min.pseudotime.per.segment = cfg$min_pseudotime_per_segment,
  dendro.node.size = cfg$dendro_node_size,
  verbose = TRUE
)

tree_rds <- file.path(cfg$outdir, "div30_urd_lineage_tree_object.rds")
saveRDS(urd, tree_rds)
log_msg("Saved tree object: ", tree_rds)

finalize_command <- paste(
  "Rscript scripts/17_div30_urd_finalize_lineage_tree_report.R",
  "--tree-rds", shQuote(tree_rds),
  "--outdir", shQuote(cfg$outdir),
  "--annotation-col", shQuote(cfg$annotation_col),
  "--pseudotime-name", shQuote(cfg$pseudotime_name),
  "--tip-mapping", shQuote(file.path(table_dir, "tree_tip_mapping.tsv")),
  "--top-n-branch-genes", shQuote(as.character(cfg$top_n_branch_genes))
)
writeLines(finalize_command, file.path(cfg$outdir, "finalize_lineage_tree_report.command.txt"))
log_msg("Tree build stage complete. Finalize tables/plots in a fresh R process with:")
log_msg(finalize_command)
quit(save = "no", status = 0)

status <- tree_status(urd, tip_mapping)
joins <- segment_join_table(urd)
tip_comp <- tip_composition(urd, cfg$annotation_col, tip_mapping)
branches <- branch_genes(urd, tip_mapping, cfg$top_n_branch_genes)

write_tsv(status, file.path(table_dir, "tree_status.tsv"))
write_tsv(joins, file.path(table_dir, "tree_segment_joins.tsv"))
write_tsv(tip_comp, file.path(table_dir, "tree_tip_composition.tsv"))
branch_gene_path <- file.path(table_dir, "branch_specific_genes.tsv")
write_tsv(branches, branch_gene_path)

plot_tree_safe(urd, cfg$annotation_col, file.path(plot_dir, "urd_tree_annotation.png"), "URD tree colored by paper/manual annotation")
plot_tree_safe(urd, cfg$pseudotime_name, file.path(plot_dir, "urd_tree_pseudotime.png"), "URD tree colored by flood pseudotime")

write_report(
  file.path(cfg$outdir, "urd_lineage_tree_report.md"),
  cfg,
  status,
  tip_mapping,
  tip_comp,
  joins,
  branch_gene_path
)

log_msg("Done. Tree report: ", file.path(cfg$outdir, "urd_lineage_tree_report.md"))
