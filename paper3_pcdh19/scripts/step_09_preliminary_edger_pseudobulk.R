#!/usr/bin/env Rscript

# PRELIMINARY Step 09 edgeR module. Input counts are sample/state-summed raw
# integer gene-level Flex UMI/probe-ligation evidence from filtered cells.

args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 4L) stop("Expected counts, metadata, contrasts, output root")
counts_path <- args[[1]]
metadata_path <- args[[2]]
contrasts_path <- args[[3]]
output_root <- args[[4]]

suppressPackageStartupMessages(library(edgeR))

STATUS <- "PRELIMINARY"
strata <- c("all_cells", "progenitor", "immature_neuron")
group_order <- c("WT_M", "WT_F", "HET_F", "KO_M")
group_colors <- c(WT_M="#2166ac", WT_F="#67a9cf", HET_F="#fdae61", KO_M="#b2182b")

counts_frame <- read.delim(gzfile(counts_path), check.names=FALSE, stringsAsFactors=FALSE)
metadata <- read.delim(metadata_path, check.names=FALSE, stringsAsFactors=FALSE)
contrasts <- read.delim(contrasts_path, check.names=FALSE, stringsAsFactors=FALSE)
if (!all(metadata$analysis_status == STATUS) || !all(contrasts$analysis_status == STATUS)) stop("PRELIMINARY status missing")
if (anyDuplicated(counts_frame$ensembl_gene_id)) stop("Duplicate Ensembl gene IDs")
gene_ids <- counts_frame$ensembl_gene_id
gene_symbols <- counts_frame$gene_symbol
count_matrix <- as.matrix(counts_frame[, metadata$pseudobulk_id, drop=FALSE])
storage.mode(count_matrix) <- "integer"
rownames(count_matrix) <- gene_ids
if (any(count_matrix < 0L) || any(abs(count_matrix - round(count_matrix)) > 0)) stop("Counts are not nonnegative integers")
if (!identical(colnames(count_matrix), metadata$pseudobulk_id)) stop("Count/metadata order mismatch")
if (!all(colSums(count_matrix) == metadata$raw_gene_umi_ligation_count)) stop("Library sums do not match metadata")

figures <- file.path(output_root, "figures")
de_root <- file.path(output_root, "differential_expression")
dir.create(figures, showWarnings=FALSE, recursive=TRUE)
dir.create(de_root, showWarnings=FALSE, recursive=TRUE)

safe_png <- function(path, width=1800, height=1400, res=180) {
  png(path, width=width, height=height, res=res)
}

sample_label <- function(x) sub("15662-", "", x, fixed=TRUE)

normalization_rows <- list()
summary_rows <- list()
marker_rows <- list()
test_index <- 0L

# Raw pseudobulk library-size overview.
safe_png(file.path(figures, "step_09_preliminary_raw_pseudobulk_library_sizes.png"), 2400, 1200, 180)
par(mar=c(9,5,4,1))
bar_cols <- group_colors[metadata$design_group]
barplot(metadata$raw_gene_umi_ligation_count / 1e6, names.arg=paste(sample_label(metadata$biological_sample_id), metadata$stratum, sep="\n"),
        las=2, cex.names=0.65, col=bar_cols, border=NA, ylab="Raw pseudobulk gene-level UMI/ligation evidence (millions)",
        main="PRELIMINARY Step 09: raw pseudobulk library sizes")
legend("topright", legend=names(group_colors), fill=group_colors, bty="n", ncol=4)
dev.off()

# Per-stratum normalization, MDS, and marker context using all 12 samples.
marker_symbols <- c("Pcdh19", "Sox2", "Hes1", "Hes5", "Mki67", "Top2a", "Lhx6", "Gad1", "Gad2", "Dlx5", "Dcx", "Stmn2")
marker_matrices <- list()
for (stratum in strata) {
  meta_s <- metadata[metadata$stratum == stratum, , drop=FALSE]
  y_all <- DGEList(counts=count_matrix[, meta_s$pseudobulk_id, drop=FALSE], genes=data.frame(ensembl_gene_id=gene_ids, gene_symbol=gene_symbols, stringsAsFactors=FALSE))
  keep_all <- filterByExpr(y_all, group=factor(meta_s$design_group, levels=group_order))
  y_all <- calcNormFactors(y_all[keep_all, , keep.lib.sizes=FALSE], method="TMM")
  for (i in seq_len(nrow(meta_s))) {
    normalization_rows[[length(normalization_rows)+1L]] <- data.frame(
      analysis_status=STATUS, stratum=stratum, pseudobulk_id=meta_s$pseudobulk_id[i],
      biological_sample_id=meta_s$biological_sample_id[i], design_group=meta_s$design_group[i],
      cells_aggregated=meta_s$cells_aggregated[i], raw_library_size=y_all$samples$lib.size[i],
      TMM_normalization_factor=y_all$samples$norm.factors[i],
      effective_library_size=y_all$samples$lib.size[i] * y_all$samples$norm.factors[i], stringsAsFactors=FALSE)
  }
  safe_png(file.path(figures, paste0("step_09_preliminary_mds_", stratum, ".png")), 1600, 1350, 180)
  plotMDS(y_all, col=group_colors[meta_s$design_group], pch=19, cex=1.4,
          labels=sample_label(meta_s$biological_sample_id), main=paste("PRELIMINARY Step 09 MDS:", stratum))
  legend("topright", legend=names(group_colors), col=group_colors, pch=19, bty="n")
  dev.off()
  logcpm_all <- cpm(y_all, log=TRUE, prior.count=2)
  available <- marker_symbols[marker_symbols %in% y_all$genes$gene_symbol]
  marker_matrix <- matrix(NA_real_, nrow=length(available), ncol=length(group_order), dimnames=list(available, group_order))
  for (marker in available) {
    idx <- which(y_all$genes$gene_symbol == marker)[1]
    for (group in group_order) {
      marker_matrix[marker, group] <- mean(logcpm_all[idx, meta_s$design_group == group])
      marker_rows[[length(marker_rows)+1L]] <- data.frame(analysis_status=STATUS, stratum=stratum, gene_symbol=marker,
        design_group=group, mean_TMM_logCPM=marker_matrix[marker, group], stringsAsFactors=FALSE)
    }
  }
  marker_matrices[[stratum]] <- marker_matrix
}

safe_png(file.path(figures, "step_09_preliminary_marker_gene_group_mean_logcpm.png"), 2200, 1700, 180)
par(mfrow=c(1,3), mar=c(7,7,5,2))
global_marker_range <- range(unlist(marker_matrices), finite=TRUE)
for (stratum in strata) {
  mat <- marker_matrices[[stratum]]
  image(seq_len(ncol(mat)), seq_len(nrow(mat)), t(mat[nrow(mat):1,,drop=FALSE]), col=hcl.colors(64, "YlOrRd"),
        zlim=global_marker_range, axes=FALSE, xlab="", ylab="", main=stratum)
  axis(1, at=seq_len(ncol(mat)), labels=colnames(mat), las=2)
  axis(2, at=seq_len(nrow(mat)), labels=rev(rownames(mat)), las=2)
  for (i in seq_len(nrow(mat))) for (j in seq_len(ncol(mat)))
    text(j, nrow(mat)-i+1, sprintf("%.1f", mat[i,j]), cex=0.68)
}
mtext("PRELIMINARY Step 09: group-mean TMM logCPM marker context", outer=TRUE, line=-1.5, cex=1.2)
dev.off()

# Nine prespecified contrast-by-stratum edgeR QL tests and diagnostic plots.
for (contrast_i in seq_len(nrow(contrasts))) {
  contrast <- contrasts[contrast_i,]
  for (stratum in strata) {
    test_index <- test_index + 1L
    groups <- c(contrast$reference_group, contrast$case_group)
    meta_test <- metadata[metadata$stratum == stratum & metadata$design_group %in% groups, , drop=FALSE]
    meta_test$group <- factor(meta_test$design_group, levels=groups)
    if (!identical(as.integer(table(meta_test$group)), c(3L,3L))) stop("Each test requires 3+3 samples")
    y <- DGEList(counts=count_matrix[, meta_test$pseudobulk_id, drop=FALSE],
                 genes=data.frame(ensembl_gene_id=gene_ids, gene_symbol=gene_symbols, stringsAsFactors=FALSE),
                 group=meta_test$group)
    keep <- filterByExpr(y, group=meta_test$group)
    y <- calcNormFactors(y[keep, , keep.lib.sizes=FALSE], method="TMM")
    design <- model.matrix(~0 + group, data=meta_test)
    colnames(design) <- levels(meta_test$group)
    y <- estimateDisp(y, design, robust=TRUE)
    fit <- glmQLFit(y, design, robust=TRUE)
    contrast_vector <- rep(0, ncol(design)); names(contrast_vector) <- colnames(design)
    contrast_vector[contrast$case_group] <- 1
    contrast_vector[contrast$reference_group] <- -1
    qlf <- glmQLFTest(fit, contrast=contrast_vector)
    result <- topTags(qlf, n=Inf, sort.by="PValue")$table
    logcpm <- cpm(y, log=TRUE, prior.count=2)
    ref_cols <- meta_test$design_group == contrast$reference_group
    case_cols <- meta_test$design_group == contrast$case_group
    result$reference_mean_TMM_logCPM <- rowMeans(logcpm[rownames(result), ref_cols, drop=FALSE])
    result$case_mean_TMM_logCPM <- rowMeans(logcpm[rownames(result), case_cols, drop=FALSE])
    result$analysis_status <- STATUS
    result$contrast_id <- contrast$contrast_id
    result$stratum <- stratum
    result$case_group <- contrast$case_group
    result$reference_group <- contrast$reference_group
    result$ensembl_gene_id <- rownames(result)
    result$gene_symbol <- y$genes[rownames(result), "gene_symbol"]
    result <- result[, c("analysis_status", "contrast_id", "stratum", "case_group", "reference_group", "ensembl_gene_id", "gene_symbol",
                         "logFC", "logCPM", "F", "PValue", "FDR", "reference_mean_TMM_logCPM", "case_mean_TMM_logCPM")]
    stem <- paste(contrast$contrast_id, stratum, sep="__")
    write.table(result, gzfile(file.path(de_root, paste0("step_09_preliminary_", stem, "__edger.tsv.gz"))), sep="\t", quote=FALSE, row.names=FALSE)
    fdr05 <- result$FDR <= 0.05
    summary_rows[[test_index]] <- data.frame(
      analysis_status=STATUS, contrast_id=contrast$contrast_id, interpretation=contrast$interpretation,
      case_group=contrast$case_group, reference_group=contrast$reference_group, stratum=stratum,
      case_samples=3L, reference_samples=3L, genes_input=nrow(count_matrix), genes_tested=nrow(result),
      genes_FDR_le_0.05=sum(fdr05), genes_FDR_le_0.05_up_in_case=sum(fdr05 & result$logFC > 0),
      genes_FDR_le_0.05_down_in_case=sum(fdr05 & result$logFC < 0),
      genes_FDR_le_0.05_abs_logFC_ge_1=sum(fdr05 & abs(result$logFC) >= 1),
      top_gene_symbol=if(nrow(result)) result$gene_symbol[1] else NA,
      top_gene_logFC=if(nrow(result)) result$logFC[1] else NA,
      top_gene_FDR=if(nrow(result)) result$FDR[1] else NA, stringsAsFactors=FALSE)

    # Volcano.
    safe_png(file.path(figures, paste0("step_09_preliminary_volcano_", stem, ".png")), 1500, 1350, 180)
    sig <- result$FDR <= 0.05
    cols <- ifelse(sig & result$logFC > 0, "#b2182b", ifelse(sig & result$logFC < 0, "#2166ac", "#bdbdbd"))
    plot(result$logFC, -log10(pmax(result$FDR, .Machine$double.xmin)), pch=16, cex=0.45, col=cols,
         xlab=paste0("log2 fold change (", contrast$case_group, " / ", contrast$reference_group, ")"), ylab="-log10 FDR",
         main=paste("PRELIMINARY", contrast$contrast_id, stratum, sep="\n"))
    abline(h=-log10(0.05), lty=2, col="#555555"); abline(v=c(-1,1), lty=3, col="#777777")
    top_label <- head(which(!is.na(result$gene_symbol) & result$gene_symbol != ""), 10)
    text(result$logFC[top_label], -log10(pmax(result$FDR[top_label], .Machine$double.xmin)), labels=result$gene_symbol[top_label], cex=0.65, pos=3)
    legend("topright", legend=c("FDR<=0.05, up", "FDR<=0.05, down", "not FDR<=0.05"), col=c("#b2182b", "#2166ac", "#bdbdbd"), pch=16, bty="n")
    dev.off()

    # MA/MD plot.
    safe_png(file.path(figures, paste0("step_09_preliminary_ma_", stem, ".png")), 1500, 1350, 180)
    plot(result$logCPM, result$logFC, pch=16, cex=0.45, col=cols, xlab="Average logCPM", ylab="log2 fold change",
         main=paste("PRELIMINARY MA plot", contrast$contrast_id, stratum, sep="\n"))
    abline(h=0, col="black"); abline(h=c(-1,1), lty=3, col="#777777")
    text(result$logCPM[top_label], result$logFC[top_label], labels=result$gene_symbol[top_label], cex=0.65, pos=3)
    dev.off()

    # Top-30 gene sample heatmap on TMM logCPM.
    top_ids <- head(rownames(result)[apply(logcpm[rownames(result),,drop=FALSE], 1, sd) > 0], 30)
    heat <- logcpm[top_ids,,drop=FALSE]
    rownames(heat) <- make.unique(ifelse(y$genes[top_ids,"gene_symbol"] == "", top_ids, y$genes[top_ids,"gene_symbol"]))
    safe_png(file.path(figures, paste0("step_09_preliminary_heatmap_top30_", stem, ".png")), 1700, 1800, 180)
    heatmap(heat, scale="row", Colv=NA, margins=c(10,10), col=hcl.colors(64, "Blue-Red 3"),
            labCol=paste(sample_label(meta_test$biological_sample_id), meta_test$design_group, sep="\n"),
            ColSideColors=group_colors[meta_test$design_group], cexRow=0.7, cexCol=0.75,
            main=paste("PRELIMINARY top-30", contrast$contrast_id, stratum, sep="\n"))
    legend("topright", legend=groups, fill=group_colors[groups], bty="n", cex=0.8)
    dev.off()
  }
}

summary_table <- do.call(rbind, summary_rows)
write.table(summary_table, file.path(output_root, "step_09_preliminary_differential_expression_summary.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
write.table(do.call(rbind, normalization_rows), file.path(output_root, "step_09_preliminary_TMM_normalization_factors.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
write.table(do.call(rbind, marker_rows), file.path(output_root, "step_09_preliminary_marker_group_mean_logcpm.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

safe_png(file.path(figures, "step_09_preliminary_deg_counts_summary.png"), 1900, 1350, 180)
par(mar=c(10,5,4,1))
labels <- paste(summary_table$contrast_id, summary_table$stratum, sep="\n")
heights <- rbind(summary_table$genes_FDR_le_0.05_up_in_case, summary_table$genes_FDR_le_0.05_down_in_case)
barplot(heights, names.arg=labels, las=2, col=c("#b2182b", "#2166ac"), border=NA,
        ylab="Genes with FDR <= 0.05", main="PRELIMINARY Step 09: pseudobulk DE summary")
legend("topright", legend=c("Up in case", "Down in case"), fill=c("#b2182b", "#2166ac"), bty="n")
dev.off()

writeLines(capture.output(sessionInfo()), file.path(output_root, "step_09_preliminary_R_session_info.txt"))
