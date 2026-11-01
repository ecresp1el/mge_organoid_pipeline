args <- commandArgs(trailingOnly=TRUE)
out <- args[1]
suppressPackageStartupMessages(library(edgeR))
counts <- read.delim(gzfile(file.path(out,'raw_pseudobulk_counts.tsv.gz')),check.names=FALSE)
meta <- read.delim(file.path(out,'sample_metadata.tsv'),stringsAsFactors=FALSE)
stopifnot(identical(names(counts)[-(1:2)],meta$library))
mat <- as.matrix(counts[,-(1:2)])
rownames(mat) <- counts$gene_id
genes <- counts[,1:2]
rownames(genes) <- genes$gene_id
summary <- list()
pdf(file.path(out,'top10_sample_expression.pdf'),width=8,height=6)
for (case in c('HET_WT_like','HET_KO_like')) {
  use <- meta$group %in% c('WT_F_all',case)
  md <- meta[use,]
  group <- factor(md$group,levels=c('WT_F_all',case))
  design <- model.matrix(~group)
  y <- DGEList(mat[,use],genes=genes)
  keep <- filterByExpr(y,design=design)
  stopifnot(sum(keep)>=10)
  y <- calcNormFactors(y[keep,,keep.lib.sizes=FALSE],method='TMM')
  y <- estimateDisp(y,design,robust=TRUE)
  fit <- glmQLFit(y,design,robust=TRUE)
  test <- glmQLFTest(fit,coef=2)
  tab <- topTags(test,n=Inf,sort.by='PValue')$table
  # Positive logFC means higher in the HET subset than in unfiltered WT females.
  logcpm <- cpm(y,log=TRUE,prior.count=2)
  tab$direction <- ifelse(tab$logFC>0,'higher_in_HET_subset','lower_in_HET_subset')
  tab$probe_gate_gene <- tab$gene_symbol=='Pcdh19'
  tab$all_9_sample_pairs_same_direction <- vapply(rownames(tab),function(g) {
    a <- logcpm[g,group=='WT_F_all']; b <- logcpm[g,group==case]
    if(tab[g,'logFC']>0) min(b)>max(a) else max(b)<min(a)
  },logical(1))
  write.table(tab,file.path(out,paste0(case,'_vs_WT_F_all_genes.tsv')),sep='\t',quote=FALSE,row.names=FALSE)
  # Pcdh19 is reported in the full table but excluded from the discovery top ten because it defines the gate.
  ranked <- tab[!tab$probe_gate_gene,]
  top <- head(ranked,10)
  write.table(top,file.path(out,paste0(case,'_vs_WT_F_top10.tsv')),sep='\t',quote=FALSE,row.names=FALSE)
  write.table(data.frame(gene_id=rownames(logcpm),logcpm,check.names=FALSE),file.path(out,paste0(case,'_sample_logCPM.tsv')),sep='\t',quote=FALSE,row.names=FALSE)
  write.table(data.frame(md,normalization_factor=y$samples$norm.factors,effective_library_size=y$samples$lib.size*y$samples$norm.factors),file.path(out,paste0(case,'_library_QC.tsv')),sep='\t',quote=FALSE,row.names=FALSE)
  vals <- logcpm[top$gene_id,,drop=FALSE]
  rownames(vals) <- paste0(top$gene_symbol,' (',sprintf('%+.2f',top$logFC),')')
  heatmap(vals,Rowv=NA,Colv=NA,scale='row',margins=c(10,12),col=colorRampPalette(c('#337EAC','white','#B54D60'))(101),main=paste(case,'versus WT females\nTop 10 by sample-level edgeR QL P value'),labCol=paste(md$group,md$sample,sep='\n'))
  summary[[case]] <- data.frame(comparison=case,genes_tested=nrow(tab),FDR_005=sum(tab$FDR<=.05),FDR_005_non_gate_gene=sum(tab$FDR<=.05 & !tab$probe_gate_gene),smallest_FDR=min(tab$FDR),minimum_cells=min(md$cells),maximum_cells=max(md$cells))
}
dev.off()
write.table(do.call(rbind,summary),file.path(out,'DE_summary.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
capture.output(sessionInfo(),file=file.path(out,'R_session_info.txt'))
