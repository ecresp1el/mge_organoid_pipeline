#!/usr/bin/env Rscript

# Build a clean DIV90 UMAP audit figure showing cluster numbers on the UMAP
# and the exact `cluster_number_name` metadata labels in numeric order.
#
# This is upstream of URD root/tip selection. It is intentionally a metadata
# audit figure, not a lineage reconstruction output.

log_msg <- function(...) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste0(..., collapse = " ")))
  flush.console()
}

parse_args <- function(args) {
  out <- list(
    obs = NULL,
    umap = NULL,
    outdir = NULL,
    `label-col` = "cluster_number_name",
    `cluster-col` = "cluster_id",
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

required <- c("ggplot2", "ggrepel", "cowplot")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) stop("Missing required R package(s): ", paste(missing, collapse = ", "), call. = FALSE)
suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(cowplot)
})

write_tsv <- function(x, path) {
  write.table(x, path, sep = "\t", row.names = FALSE, quote = FALSE, na = "")
}

read_tsv <- function(path) {
  read.delim(path, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
}

clean_biology_name <- function(exact_name) {
  sub("^[0-9]+ - ", "", exact_name)
}

cluster_role_guess <- function(cluster_id, biology_name) {
  if (grepl("Stressed", biology_name, ignore.case = TRUE)) return("exclude/track separately: stressed")
  if (grepl("Astrocytes", biology_name, ignore.case = TRUE)) return("exclude from neuronal URD: astrocyte lineage")
  if (grepl("OPC", biology_name, ignore.case = TRUE)) return("exclude from neuronal URD: OPC lineage")
  if (grepl("Dividing", biology_name, ignore.case = TRUE)) return("candidate root/proliferative progenitor")
  if (grepl("SST|NPY", biology_name, ignore.case = TRUE)) return("candidate SST/cortical interneuron tip")
  if (grepl("CRABP1", biology_name, ignore.case = TRUE)) return("candidate CRABP1/PV-associated branch")
  if (grepl("PV", biology_name, ignore.case = TRUE)) return("candidate PV lineage tip")
  if (grepl("LHX8|vMGE|Striatal|GP|MGE", biology_name, ignore.case = TRUE)) return("candidate MGE/LHX8/subpallial tip")
  "review"
}

palette_for_clusters <- function(n) {
  base <- c(
    "#1f78b4", "#33a02c", "#6a3d9a", "#b15928", "#a6cee3",
    "#fb9a99", "#999999", "#bdbdbd", "#e31a1c", "#cab2d6",
    "#fdbf6f", "#ff7f00", "#ffff99", "#8dd3c7", "#bebada"
  )
  rep(base, length.out = n)
}

print_usage <- function() {
  cat("Usage: Rscript scripts/26_div90_umap_cluster_label_audit.R --obs <obs.tsv> --umap <umap.tsv> --outdir <dir>\n")
}

opt <- parse_args(commandArgs(trailingOnly = TRUE))
if (isTRUE(opt$help)) {
  print_usage()
  quit(save = "no", status = 0)
}
if (is.null(opt$obs) || !nzchar(opt$obs)) stop("--obs is required", call. = FALSE)
if (is.null(opt$umap) || !nzchar(opt$umap)) stop("--umap is required", call. = FALSE)
if (is.null(opt$outdir) || !nzchar(opt$outdir)) stop("--outdir is required", call. = FALSE)

plot_dir <- file.path(opt$outdir, "plots")
table_dir <- file.path(opt$outdir, "tables")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

log_msg("Reading obs: ", opt$obs)
obs <- read_tsv(opt$obs)
log_msg("Reading UMAP: ", opt$umap)
umap <- read_tsv(opt$umap)

needed <- c("cell_id", opt$`label-col`, opt$`cluster-col`)
missing_cols <- setdiff(needed, colnames(obs))
if (length(missing_cols) > 0) stop("Missing obs columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
missing_umap <- setdiff(c("cell_id", "UMAP_1", "UMAP_2"), colnames(umap))
if (length(missing_umap) > 0) stop("Missing UMAP columns: ", paste(missing_umap, collapse = ", "), call. = FALSE)

df <- merge(umap, obs[, needed, drop = FALSE], by = "cell_id")
df$cluster_id_numeric <- suppressWarnings(as.integer(df[[opt$`cluster-col`]]))
df$exact_metadata_name <- df[[opt$`label-col`]]
df$biology_name <- clean_biology_name(df$exact_metadata_name)

mapping <- aggregate(cell_id ~ cluster_id_numeric + exact_metadata_name + biology_name, df, length)
colnames(mapping)[colnames(mapping) == "cell_id"] <- "n_cells"
mapping <- mapping[order(mapping$cluster_id_numeric), , drop = FALSE]
mapping$role_guess_for_urd_smoke <- mapply(cluster_role_guess, mapping$cluster_id_numeric, mapping$biology_name)
mapping$paper_cluster_reconciliation_note <- ifelse(
  mapping$cluster_id_numeric %in% c(6, 7),
  "two Seurat clusters share the same exact metadata biology name: Stressed Cells",
  ifelse(
    mapping$cluster_id_numeric %in% c(4, 10),
    "two Seurat clusters are separate Pre-Astrocytes/Astrocytes states",
    ifelse(
      mapping$cluster_id_numeric %in% c(5, 8),
      "two Seurat clusters split the LHX8+ vMGE GABAergic Striatal/GP-fated biology",
      "one Seurat cluster label"
    )
  )
)
write_tsv(mapping, file.path(table_dir, "div90_cluster_number_name_to_biology_mapping.tsv"))

centers <- aggregate(cbind(UMAP_1, UMAP_2) ~ cluster_id_numeric + exact_metadata_name, df, median)
centers <- centers[order(centers$cluster_id_numeric), , drop = FALSE]

df$cluster_id_factor <- factor(df$cluster_id_numeric, levels = mapping$cluster_id_numeric)
colors <- setNames(palette_for_clusters(nrow(mapping)), as.character(mapping$cluster_id_numeric))

umap_plot <- ggplot(df, aes(UMAP_1, UMAP_2, color = cluster_id_factor)) +
  geom_point(size = 0.18, alpha = 0.68) +
  ggrepel::geom_label_repel(
    data = centers,
    aes(UMAP_1, UMAP_2, label = cluster_id_numeric),
    inherit.aes = FALSE,
    size = 4.8,
    fontface = "bold",
    color = "black",
    fill = "white",
    label.size = 0.25,
    min.segment.length = 0,
    max.overlaps = Inf
  ) +
  scale_color_manual(values = colors, guide = "none") +
  coord_equal() +
  theme_bw(base_size = 11) +
  theme(panel.grid = element_blank()) +
  labs(
    title = "DIV90 UMAP by exact metadata cluster_number_name",
    subtitle = "Numbers are Seurat cluster IDs; exact metadata names are listed at right.",
    x = "UMAP 1",
    y = "UMAP 2"
  )

mapping$legend_text <- sprintf("%2d = %s (n=%s)", mapping$cluster_id_numeric, mapping$biology_name, format(mapping$n_cells, big.mark = ","))
legend_text <- paste(mapping$legend_text, collapse = "\n")
legend_plot <- cowplot::ggdraw() +
  cowplot::draw_label("Exact metadata name key", x = 0, y = 0.98, hjust = 0, vjust = 1, fontface = "bold", size = 12) +
  cowplot::draw_label(legend_text, x = 0, y = 0.91, hjust = 0, vjust = 1, size = 9, lineheight = 1.05)

combined <- cowplot::plot_grid(umap_plot, legend_plot, ncol = 2, rel_widths = c(1.25, 1))

highlight_rows <- lapply(mapping$cluster_id_numeric, function(cluster_id) {
  tmp <- df
  tmp$highlight_cluster <- cluster_id
  tmp$highlight <- tmp$cluster_id_numeric == cluster_id
  tmp$facet_label <- sprintf(
    "%d - %s\nn=%s",
    cluster_id,
    mapping$biology_name[mapping$cluster_id_numeric == cluster_id],
    format(mapping$n_cells[mapping$cluster_id_numeric == cluster_id], big.mark = ",")
  )
  tmp
})
highlight_df <- do.call(rbind, highlight_rows)
highlight_df$facet_label <- factor(
  highlight_df$facet_label,
  levels = vapply(mapping$cluster_id_numeric, function(cluster_id) {
    sprintf(
      "%d - %s\nn=%s",
      cluster_id,
      mapping$biology_name[mapping$cluster_id_numeric == cluster_id],
      format(mapping$n_cells[mapping$cluster_id_numeric == cluster_id], big.mark = ",")
    )
  }, character(1))
)
highlight_plot <- ggplot() +
  geom_point(
    data = highlight_df[!highlight_df$highlight, , drop = FALSE],
    aes(UMAP_1, UMAP_2),
    color = "grey82",
    size = 0.08,
    alpha = 0.28
  ) +
  geom_point(
    data = highlight_df[highlight_df$highlight, , drop = FALSE],
    aes(UMAP_1, UMAP_2, color = cluster_id_factor),
    size = 0.17,
    alpha = 0.9
  ) +
  facet_wrap(~facet_label, ncol = ceiling(sqrt(nrow(mapping)))) +
  scale_color_manual(values = colors, guide = "none") +
  coord_equal() +
  theme_bw(base_size = 7) +
  theme(
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    strip.text = element_text(size = 7.2, face = "bold"),
    plot.title = element_text(size = 13, face = "bold")
  ) +
  labs(title = "DIV90 UMAP cluster overlays: one cluster highlighted per panel")

ggsave(file.path(plot_dir, "div90_umap_cluster_number_name_labeled.png"), combined, width = 14, height = 7, dpi = 300, bg = "white")
ggsave(file.path(plot_dir, "div90_umap_cluster_number_name_labeled.pdf"), combined, width = 14, height = 7, bg = "white")
ggsave(file.path(plot_dir, "div90_umap_cluster_numbers_only.png"), umap_plot, width = 7, height = 6, dpi = 300, bg = "white")
ggsave(file.path(plot_dir, "div90_umap_cluster_numbers_only.pdf"), umap_plot, width = 7, height = 6, bg = "white")
ggsave(file.path(plot_dir, "div90_umap_cluster_overlay_grid.png"), highlight_plot, width = 12, height = 12, dpi = 300, bg = "white")
ggsave(file.path(plot_dir, "div90_umap_cluster_overlay_grid.pdf"), highlight_plot, width = 12, height = 12, bg = "white")

report <- c(
  "# DIV90 UMAP Cluster Label Audit",
  "",
  paste0("- obs: `", opt$obs, "`"),
  paste0("- umap: `", opt$umap, "`"),
  paste0("- label column: `", opt$`label-col`, "`"),
  paste0("- cluster column: `", opt$`cluster-col`, "`"),
  "",
  "## Outputs",
  "",
  "- `plots/div90_umap_cluster_number_name_labeled.png`",
  "- `plots/div90_umap_cluster_number_name_labeled.pdf`",
  "- `plots/div90_umap_cluster_numbers_only.png`",
  "- `plots/div90_umap_cluster_overlay_grid.png`",
  "- `plots/div90_umap_cluster_overlay_grid.pdf`",
  "- `tables/div90_cluster_number_name_to_biology_mapping.tsv`",
  "",
  "## Mapping",
  "",
  paste(capture.output(print(mapping, row.names = FALSE)), collapse = "\n")
)
writeLines(report, file.path(opt$outdir, "div90_umap_cluster_label_audit_report.md"))
log_msg("Done: ", file.path(opt$outdir, "div90_umap_cluster_label_audit_report.md"))
