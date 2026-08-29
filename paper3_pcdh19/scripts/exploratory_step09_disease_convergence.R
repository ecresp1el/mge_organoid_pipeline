#!/usr/bin/env Rscript

# PRELIMINARY / EXPLORATORY follow-up to frozen Step 09.
# This is deliberately not a numbered pipeline step. It reads the complete edgeR
# tables without changing Steps 03-09 and writes an isolated ad hoc result bundle.

suppressPackageStartupMessages({
  library(data.table)
  library(fgsea)
  library(ggplot2)
  library(msigdbr)
})

options(stringsAsFactors = FALSE)
set.seed(20260828)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: exploratory_step09_disease_convergence.R STEP09_DIR OUTPUT_DIR")
}
step09_dir <- normalizePath(args[[1]], mustWork = TRUE)
output_dir <- args[[2]]
if (dir.exists(output_dir)) {
  stop("Refusing to overwrite existing exploratory output directory: ", output_dir)
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "figures"), showWarnings = FALSE)
dir.create(file.path(output_dir, "tables"), showWarnings = FALSE)

status_label <- "PRELIMINARY_EXPLORATORY"
contrasts <- c("het_f_vs_wt_f", "ko_m_vs_wt_m")
strata <- c("all_cells", "progenitor", "immature_neuron")
contrast_labels <- c(
  het_f_vs_wt_f = "HET-F vs WT-F",
  ko_m_vs_wt_m = "KO-M vs WT-M"
)
stratum_labels <- c(
  all_cells = "All cells",
  progenitor = "Progenitors",
  immature_neuron = "Immature neurons"
)

read_de <- function(contrast, stratum) {
  path <- file.path(
    step09_dir, "differential_expression",
    sprintf("step_09_preliminary_%s__%s__edger.tsv.gz", contrast, stratum)
  )
  if (!file.exists(path)) stop("Missing Step 09 table: ", path)
  x <- fread(path)
  required <- c("ensembl_gene_id", "gene_symbol", "logFC", "F", "PValue", "FDR")
  missing <- setdiff(required, names(x))
  if (length(missing)) stop("Missing columns in ", path, ": ", paste(missing, collapse = ", "))
  x[, `:=`(
    contrast_id = contrast,
    stratum_id = stratum,
    signed_stat = sign(logFC) * sqrt(pmax(F, 0))
  )]
  x
}

de_list <- list()
for (s in strata) for (c in contrasts) {
  de_list[[paste(c, s, sep = "__")]] <- read_de(c, s)
}

# Native mouse MSigDB (2026.1.Mm at analysis time).
message("Loading native mouse MSigDB collections...")
msig_specs <- list(
  Hallmark = list(collection = "MH", subcollection = NULL),
  GO_BP = list(collection = "M5", subcollection = "GO:BP"),
  Reactome = list(collection = "M2", subcollection = "CP:REACTOME")
)
msig_tables <- lapply(names(msig_specs), function(label) {
  spec <- msig_specs[[label]]
  d <- msigdbr(
    db_species = "MM", species = "Mus musculus",
    collection = spec$collection, subcollection = spec$subcollection
  )
  d <- as.data.table(d)[, .(gene_symbol, gs_name, db_version)]
  d[, collection := label]
  unique(d)
})
names(msig_tables) <- names(msig_specs)
msig <- rbindlist(msig_tables, use.names = TRUE)
msig_version <- paste(unique(msig$db_version), collapse = ";")
pathways <- split(msig$gene_symbol, paste(msig$collection, msig$gs_name, sep = "::"))
pathways <- lapply(pathways, unique)

deduplicate_ranks <- function(x) {
  y <- x[!is.na(gene_symbol) & nzchar(gene_symbol) & is.finite(signed_stat)]
  y[, abs_signed_stat := abs(signed_stat)]
  setorder(y, gene_symbol, -abs_signed_stat, PValue, ensembl_gene_id)
  y <- y[!duplicated(gene_symbol)]
  ranks <- y$signed_stat
  names(ranks) <- y$gene_symbol
  sort(ranks, decreasing = TRUE)
}

message("Running six preranked GSEA analyses...")
gsea_all <- vector("list", length(de_list))
names(gsea_all) <- names(de_list)
for (key in names(de_list)) {
  ranks <- deduplicate_ranks(de_list[[key]])
  fg <- suppressWarnings(fgseaMultilevel(
    pathways = pathways, stats = ranks, minSize = 15, maxSize = 500,
    eps = 0, scoreType = "std"
  ))
  fg <- as.data.table(fg)
  fg[, c("collection", "pathway") := tstrsplit(pathway, "::", fixed = TRUE)]
  ids <- tstrsplit(key, "__", fixed = TRUE)
  fg[, `:=`(
    analysis_status = status_label,
    contrast_id = ids[[1]],
    stratum_id = ids[[2]],
    rank_metric = "sign(logFC)*sqrt(edgeR_QLF_F)",
    leading_edge_genes = vapply(leadingEdge, paste, collapse = ";", FUN.VALUE = character(1))
  )]
  fg[, leadingEdge := NULL]
  setcolorder(fg, c(
    "analysis_status", "contrast_id", "stratum_id", "collection", "pathway",
    "size", "ES", "NES", "pval", "padj", "log2err", "rank_metric",
    "leading_edge_genes"
  ))
  gsea_all[[key]] <- fg
}
gsea <- rbindlist(gsea_all, use.names = TRUE, fill = TRUE)
fwrite(gsea, file.path(output_dir, "tables", "preliminary_exploratory_all_pathway_enrichment.tsv.gz"), sep = "\t")

top_pathways <- gsea[, {
  pos <- .SD[NES > 0][order(pval, padj, -NES)][seq_len(min(10L, .N))]
  neg <- .SD[NES < 0][order(pval, padj, NES)][seq_len(min(10L, .N))]
  pos[, direction := "positive"]
  neg[, direction := "negative"]
  rbind(pos, neg, use.names = TRUE)
}, by = .(contrast_id, stratum_id, collection)]
fwrite(top_pathways, file.path(output_dir, "tables", "preliminary_exploratory_strongest_pathways.tsv"), sep = "\t")

safe_cor <- function(x, y, method) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3L) return(NA_real_)
  cor(x[ok], y[ok], method = method)
}

top_set <- function(d, n, direction = c("positive", "negative", "absolute")) {
  direction <- match.arg(direction)
  d <- d[is.finite(signed_stat) & !is.na(gene_symbol) & nzchar(gene_symbol)]
  if (direction == "positive") setorder(d, -signed_stat, gene_symbol)
  if (direction == "negative") setorder(d, signed_stat, gene_symbol)
  if (direction == "absolute") {
    d[, abs_signed_stat := abs(signed_stat)]
    setorder(d, -abs_signed_stat, gene_symbol)
  }
  unique(head(d$gene_symbol, n))
}

gene_convergence <- list()
gene_overlap <- list()
paired_genes <- list()
for (s in strata) {
  h <- copy(de_list[[paste("het_f_vs_wt_f", s, sep = "__")]])
  k <- copy(de_list[[paste("ko_m_vs_wt_m", s, sep = "__")]])
  h <- h[!duplicated(ensembl_gene_id)]
  k <- k[!duplicated(ensembl_gene_id)]
  m <- merge(
    h[, .(ensembl_gene_id, gene_symbol_het = gene_symbol, het_logFC = logFC, het_signed_stat = signed_stat)],
    k[, .(ensembl_gene_id, gene_symbol_ko = gene_symbol, ko_logFC = logFC, ko_signed_stat = signed_stat)],
    by = "ensembl_gene_id"
  )
  m[, gene_symbol := fifelse(!is.na(gene_symbol_het) & nzchar(gene_symbol_het), gene_symbol_het, gene_symbol_ko)]
  m[, `:=`(
    stratum_id = s,
    logFC_direction = fifelse(het_logFC * ko_logFC > 0, "concordant", fifelse(het_logFC * ko_logFC < 0, "opposite", "zero"))
  )]
  paired_genes[[s]] <- m
  gene_convergence[[s]] <- data.table(
    analysis_status = status_label,
    stratum_id = s,
    n_common_genes = nrow(m),
    logFC_pearson = safe_cor(m$het_logFC, m$ko_logFC, "pearson"),
    logFC_spearman = safe_cor(m$het_logFC, m$ko_logFC, "spearman"),
    signed_stat_pearson = safe_cor(m$het_signed_stat, m$ko_signed_stat, "pearson"),
    signed_stat_spearman = safe_cor(m$het_signed_stat, m$ko_signed_stat, "spearman"),
    fraction_logFC_concordant = mean(m$logFC_direction == "concordant", na.rm = TRUE)
  )
  hp <- top_set(h, 200, "positive"); hn <- top_set(h, 200, "negative")
  kp <- top_set(k, 200, "positive"); kn <- top_set(k, 200, "negative")
  ha <- top_set(h, 500, "absolute"); ka <- top_set(k, 500, "absolute")
  abs_overlap <- intersect(ha, ka)
  abs_m <- m[gene_symbol %in% abs_overlap]
  gene_overlap[[s]] <- data.table(
    analysis_status = status_label,
    stratum_id = s,
    top_positive_n = 200L,
    top_negative_n = 200L,
    top_absolute_n = 500L,
    positive_overlap_n = length(intersect(hp, kp)),
    negative_overlap_n = length(intersect(hn, kn)),
    het_positive_ko_negative_n = length(intersect(hp, kn)),
    het_negative_ko_positive_n = length(intersect(hn, kp)),
    top_absolute_overlap_n = length(abs_overlap),
    top_absolute_overlap_concordant_n = sum(abs_m$logFC_direction == "concordant"),
    top_absolute_overlap_opposite_n = sum(abs_m$logFC_direction == "opposite"),
    positive_overlap_genes = paste(intersect(hp, kp), collapse = ";"),
    negative_overlap_genes = paste(intersect(hn, kn), collapse = ";"),
    opposite_overlap_genes = paste(c(intersect(hp, kn), intersect(hn, kp)), collapse = ";")
  )
}
fwrite(rbindlist(gene_convergence), file.path(output_dir, "tables", "preliminary_exploratory_genomewide_convergence.tsv"), sep = "\t")
fwrite(rbindlist(gene_overlap), file.path(output_dir, "tables", "preliminary_exploratory_top_gene_overlap.tsv"), sep = "\t")

pathway_convergence <- list()
pathway_overlap <- list()
for (s in strata) {
  h <- gsea[contrast_id == "het_f_vs_wt_f" & stratum_id == s]
  k <- gsea[contrast_id == "ko_m_vs_wt_m" & stratum_id == s]
  m <- merge(
    h[, .(collection, pathway, het_NES = NES, het_pval = pval, het_padj = padj)],
    k[, .(collection, pathway, ko_NES = NES, ko_pval = pval, ko_padj = padj)],
    by = c("collection", "pathway")
  )
  pathway_convergence[[s]] <- data.table(
    analysis_status = status_label,
    stratum_id = s,
    n_common_pathways = nrow(m),
    NES_pearson = safe_cor(m$het_NES, m$ko_NES, "pearson"),
    NES_spearman = safe_cor(m$het_NES, m$ko_NES, "spearman"),
    fraction_NES_concordant = mean(m$het_NES * m$ko_NES > 0, na.rm = TRUE),
    shared_FDR_0_25_n = sum(m$het_padj < 0.25 & m$ko_padj < 0.25, na.rm = TRUE),
    shared_FDR_0_25_concordant_n = sum(m$het_padj < 0.25 & m$ko_padj < 0.25 & m$het_NES * m$ko_NES > 0, na.rm = TRUE)
  )
  get_top <- function(d, direction, n = 30L) {
    z <- if (direction == "positive") d[NES > 0] else d[NES < 0]
    z <- z[order(pval, padj, -abs(NES))]
    head(paste(z$collection, z$pathway, sep = "::"), n)
  }
  hp <- get_top(h, "positive"); hn <- get_top(h, "negative")
  kp <- get_top(k, "positive"); kn <- get_top(k, "negative")
  pathway_overlap[[s]] <- data.table(
    analysis_status = status_label,
    stratum_id = s,
    top_per_direction_n = 30L,
    positive_overlap_n = length(intersect(hp, kp)),
    negative_overlap_n = length(intersect(hn, kn)),
    het_positive_ko_negative_n = length(intersect(hp, kn)),
    het_negative_ko_positive_n = length(intersect(hn, kp)),
    positive_overlap_pathways = paste(intersect(hp, kp), collapse = ";"),
    negative_overlap_pathways = paste(intersect(hn, kn), collapse = ";"),
    opposite_overlap_pathways = paste(c(intersect(hp, kn), intersect(hn, kp)), collapse = ";")
  )
}
fwrite(rbindlist(pathway_convergence), file.path(output_dir, "tables", "preliminary_exploratory_pathway_convergence.tsv"), sep = "\t")
fwrite(rbindlist(pathway_overlap), file.path(output_dir, "tables", "preliminary_exploratory_top_pathway_overlap.tsv"), sep = "\t")

# Prespecified, balanced target panel; every requested category is retained whether
# or not its members support an anticipated biological direction.
target_genes <- list(
  MGE_progenitor_neurogenesis = c("Nkx2-1", "Ascl1", "Dlx1", "Dlx2", "Sox2", "Hes1", "Hes5", "Neurog1", "Neurog2", "Olig2"),
  interneuron_differentiation = c("Lhx6", "Lhx8", "Dlx5", "Dlx6", "Sox6", "Maf", "Mafb", "Arx", "Satb1", "Gsx2"),
  adhesion_protocadherin_cadherin = c("Pcdh19", "Pcdh10", "Pcdh17", "Pcdh20", "Cdh2", "Cdh11", "Cdh13", "Ctnna1", "Ctnna2", "Ctnnb1"),
  cytoskeleton_migration = c("Dcx", "Tubb3", "Map1b", "Flna", "Dbn1", "Cfl1", "Pfn1", "Rac1", "Cdc42", "Rhoa"),
  WNT = c("Ctnnb1", "Axin1", "Axin2", "Apc", "Lef1", "Tcf7l2", "Wnt5a", "Fzd3", "Dvl1", "Gsk3b"),
  Notch = c("Notch1", "Notch2", "Hes1", "Hes5", "Dll1", "Jag1", "Rbpj", "Hey1", "Hey2", "Numb"),
  cell_cycle = c("Mki67", "Top2a", "Ccnb1", "Ccnb2", "Cdk1", "Ccnd1", "Ccnd2", "Pcna", "Mcm2", "Mcm5"),
  synapse_neuronal_maturation = c("Stmn2", "Map2", "Rbfox3", "Syp", "Syn1", "Dlg4", "Snap25", "Gria1", "Grin1", "Tubb3"),
  GABAergic_differentiation = c("Gad1", "Gad2", "Slc6a1", "Slc32a1", "Slc12a5", "Gabra1", "Dlx5", "Dlx6", "Lhx6", "Sox6")
)
target_map <- unique(rbindlist(lapply(names(target_genes), function(category) {
  data.table(category = category, gene_symbol = target_genes[[category]])
})))
target_results <- rbindlist(lapply(names(de_list), function(key) {
  d <- de_list[[key]]
  m <- merge(target_map, d[, .(gene_symbol, logFC, signed_stat, PValue, FDR)], by = "gene_symbol", all.x = TRUE)
  ids <- tstrsplit(key, "__", fixed = TRUE)
  m[, `:=`(analysis_status = status_label, contrast_id = ids[[1]], stratum_id = ids[[2]], tested = !is.na(logFC))]
  m
}), use.names = TRUE)
fwrite(target_results, file.path(output_dir, "tables", "preliminary_exploratory_target_gene_panel.tsv"), sep = "\t")
target_summary <- target_results[tested == TRUE, .(
  n_tested = .N,
  median_logFC = median(logFC, na.rm = TRUE),
  median_signed_stat = median(signed_stat, na.rm = TRUE),
  fraction_logFC_positive = mean(logFC > 0, na.rm = TRUE),
  minimum_nominal_p = min(PValue, na.rm = TRUE),
  minimum_FDR = min(FDR, na.rm = TRUE)
), by = .(analysis_status, contrast_id, stratum_id, category)]
fwrite(target_summary, file.path(output_dir, "tables", "preliminary_exploratory_target_category_summary.tsv"), sep = "\t")

plot_scatter <- function(s) {
  d <- copy(paired_genes[[s]])
  d[, joint_stat := sqrt(het_signed_stat^2 + ko_signed_stat^2)]
  d[, point_class := "Other tested genes"]
  htop <- head(d[order(-abs(het_signed_stat))]$gene_symbol, 500)
  ktop <- head(d[order(-abs(ko_signed_stat))]$gene_symbol, 500)
  shared <- intersect(htop, ktop)
  d[gene_symbol %in% shared & logFC_direction == "concordant", point_class := "Shared top-500, concordant"]
  d[gene_symbol %in% shared & logFC_direction == "opposite", point_class := "Shared top-500, opposite"]
  label_d <- d[!is.na(gene_symbol) & nzchar(gene_symbol)][order(-joint_stat)][1:min(.N, 14)]
  cor_p <- safe_cor(d$het_logFC, d$ko_logFC, "pearson")
  cor_s <- safe_cor(d$het_logFC, d$ko_logFC, "spearman")
  ggplot(d, aes(het_logFC, ko_logFC)) +
    geom_hline(yintercept = 0, colour = "grey75", linewidth = 0.35) +
    geom_vline(xintercept = 0, colour = "grey75", linewidth = 0.35) +
    geom_abline(slope = 1, intercept = 0, linetype = 3, colour = "grey55", linewidth = 0.4) +
    geom_point(aes(colour = point_class), alpha = 0.52, size = 1.05) +
    geom_text(data = label_d, aes(label = gene_symbol), size = 2.7, check_overlap = TRUE, vjust = -0.55) +
    scale_colour_manual(values = c(
      "Other tested genes" = "grey70",
      "Shared top-500, concordant" = "#6A3D9A",
      "Shared top-500, opposite" = "#E66101"
    )) +
    labs(
      title = paste("PRELIMINARY / EXPLORATORY —", stratum_labels[[s]]),
      subtitle = sprintf("Genome-wide Step 09 log2FC; Pearson r = %.3f, Spearman rho = %.3f", cor_p, cor_s),
      x = "HET-F vs WT-F log2FC", y = "KO-M vs WT-M log2FC", colour = NULL,
      caption = "No gene-level FDR filter. Highlighting uses signed-statistic top ranks in both comparisons."
    ) +
    theme_bw(base_size = 10) +
    theme(legend.position = "bottom", plot.title = element_text(face = "bold"))
}

ggsave(
  file.path(output_dir, "figures", "preliminary_exploratory_logfc_convergence_progenitor.png"),
  plot_scatter("progenitor"), width = 7.5, height = 6.2, dpi = 300
)
ggsave(
  file.path(output_dir, "figures", "preliminary_exploratory_logfc_convergence_immature_neuron.png"),
  plot_scatter("immature_neuron"), width = 7.5, height = 6.2, dpi = 300
)

# Heatmap selection is objective and collection-balanced: for each collection,
# retain six pathways having the smallest nominal P in any of the six analyses.
heat_candidates <- gsea[, .(
  minimum_nominal_p = min(pval, na.rm = TRUE),
  minimum_FDR = min(padj, na.rm = TRUE),
  maximum_abs_NES = max(abs(NES), na.rm = TRUE)
), by = .(collection, pathway)]
setorder(heat_candidates, collection, minimum_nominal_p, minimum_FDR, -maximum_abs_NES, pathway)
heat_selected <- heat_candidates[, head(.SD, 6L), by = collection]
heat <- merge(gsea, heat_selected[, .(collection, pathway)], by = c("collection", "pathway"))
heat[, column_id := factor(
  paste(contrast_id, stratum_id, sep = "__"),
  levels = c(
    "het_f_vs_wt_f__all_cells", "ko_m_vs_wt_m__all_cells",
    "het_f_vs_wt_f__progenitor", "ko_m_vs_wt_m__progenitor",
    "het_f_vs_wt_f__immature_neuron", "ko_m_vs_wt_m__immature_neuron"
  ),
  labels = c("HET | All", "KO | All", "HET | Prog", "KO | Prog", "HET | Imm", "KO | Imm")
)]
heat[, pathway_label := gsub("_", " ", sub("^(HALLMARK_|GOBP_|REACTOME_)", "", pathway))]
heat[, row_id := paste(collection, pathway_label, sep = " | ")]
row_order <- heat_selected[order(collection, minimum_nominal_p)]
row_labels <- msig[match(paste(row_order$collection, row_order$pathway), paste(collection, gs_name)),
                   paste(collection, gsub("_", " ", sub("^(HALLMARK_|GOBP_|REACTOME_)", "", gs_name)), sep = " | ")]
heat[, row_id := factor(row_id, levels = rev(unique(row_labels)))]
p_heat <- ggplot(heat, aes(column_id, row_id, fill = NES)) +
  geom_tile(colour = "white", linewidth = 0.35) +
  geom_text(aes(label = sprintf("%.2f", NES)), size = 2.4) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, name = "Signed NES") +
  labs(
    title = "PRELIMINARY / EXPLORATORY — ranked pathway enrichment",
    subtitle = "Six strongest pathways per collection by minimum nominal P across disease/state analyses",
    x = NULL, y = NULL,
    caption = "Native mouse MSigDB; all tested genes ranked by sign(logFC) × sqrt(edgeR QL F)."
  ) +
  theme_minimal(base_size = 9) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1, face = "bold"),
    axis.text.y = element_text(size = 7.6),
    panel.grid = element_blank(), plot.title = element_text(face = "bold")
  )
ggsave(
  file.path(output_dir, "figures", "preliminary_exploratory_pathway_signed_nes_heatmap.png"),
  p_heat, width = 10.8, height = 8.4, dpi = 300
)

input_files <- file.path(
  step09_dir, "differential_expression",
  as.vector(outer(
    contrasts, strata,
    function(c, s) sprintf("step_09_preliminary_%s__%s__edger.tsv.gz", c, s)
  ))
)
input_manifest <- data.table(
  path = input_files,
  md5 = vapply(input_files, function(p) unname(tools::md5sum(p)), character(1)),
  note = "MD5 recorded for compact ad hoc provenance"
)
fwrite(input_manifest, file.path(output_dir, "input_files_and_checksums.tsv"), sep = "\t")

metadata <- data.table(
  metadata_field = c(
    "analysis_status", "formal_pipeline_step", "source_step", "rank_metric",
    "gene_filter_for_GSEA", "GSEA_method", "GSEA_min_size", "GSEA_max_size",
    "MSigDB_database", "MSigDB_version", "collections", "WT_sex_comparison_included",
    "figure_count", "created_at"
  ),
  value = c(
    status_label, "no", "frozen Step 09 complete edgeR tables",
    "sign(logFC)*sqrt(edgeR_QLF_F)", "all tested genes; no FDR prerequisite",
    "fgseaMultilevel", "15", "500", "native mouse", msig_version,
    "Hallmark;GO Biological Process;Reactome", "no", "3",
    format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
  )
)
fwrite(metadata, file.path(output_dir, "analysis_metadata.tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(output_dir, "R_session_info.txt"))

files <- list.files(output_dir, recursive = TRUE, full.names = TRUE)
files <- files[basename(files) != "output_manifest.tsv"]
manifest <- data.table(
  relative_path = substring(files, nchar(output_dir) + 2L),
  bytes = file.info(files)$size,
  md5 = unname(tools::md5sum(files))
)
fwrite(manifest, file.path(output_dir, "output_manifest.tsv"), sep = "\t")

message("Completed PRELIMINARY / EXPLORATORY analysis: ", output_dir)
