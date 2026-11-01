# Exploratory sample-level edgeR contrasts, including paired and gate-interaction controls.
args <- commandArgs(trailingOnly=TRUE)
out <- args[1]; stratum <- args[2]
suppressPackageStartupMessages(library(edgeR))
dest <- file.path(out,paste0('DE_',stratum));dir.create(dest,showWarnings=FALSE)
meta <- read.delim(file.path(out,'pseudobulk_metadata.tsv'),stringsAsFactors=FALSE)
data <- read.delim(gzfile(file.path(out,paste0(stratum,'_counts.tsv.gz'))),check.names=FALSE)
genes <- data[,1:2];rownames(genes)<-genes$gene_id
counts<-as.matrix(data[,-(1:2)]);rownames(counts)<-genes$gene_id
meta<-meta[match(colnames(counts),meta$library),];stopifnot(!anyNA(meta$library))
spec <- list()
add <- function(id,ref,case,refgeno='WT_F',casegeno='HET_F',paired=FALSE,interaction=FALSE) {
  spec[[length(spec)+1L]] <<- list(id=id,ref=ref,case=case,refgeno=refgeno,casegeno=casegeno,paired=paired,interaction=interaction)
}
for(g in c('all','A_positive','Aneg_BC1','Aneg_BC2','Aneg_BC3','frozen_KO','KOtail10','KOtail20')) add(paste0('HET_',g,'_vs_WTall'),'all',g)
if(stratum=='all_cells') {
  for(g in c('frozen_WT','WTtail05','WTtail10','WTtail20','KOtail05')) add(paste0('HET_',g,'_vs_WTall'),'all',g)
}
for(g in c('A_positive','Aneg_BC1')) add(paste0('matched_',g),g,g)
for(g in c('Aneg_BC1','Aneg_BC2')) add(paste0('within_HET_',g,'_vs_Apositive'),'A_positive',g,'HET_F','HET_F',TRUE)
add('gate_interaction_BC1_vs_Apositive','A_positive','Aneg_BC1',interaction=TRUE)
if(stratum=='all_cells') {
  for(g in c('Aneg_BC2','Aneg_BC3')) add(paste0('matched_',g),g,g)
  add('within_HET_BC3_vs_Apositive','A_positive','Aneg_BC3','HET_F','HET_F',TRUE)
  add('within_HET_frozenKO_vs_frozenWT','frozen_WT','frozen_KO','HET_F','HET_F',TRUE)
  for(pct in c('05','10','20')) add(paste0('within_HET_tail',pct),paste0('WTtail',pct),paste0('KOtail',pct),'HET_F','HET_F',TRUE)
  for(g in c('Aneg_BC1','Aneg_BC2','Aneg_BC3')) add(paste0('within_WT_',g,'_vs_Apositive'),'A_positive',g,'WT_F','WT_F',TRUE)
  for(g in c('Aneg_BC2','Aneg_BC3')) add(paste0('gate_interaction_',g,'_vs_Apositive'),'A_positive',g,interaction=TRUE)
}
specframe<-do.call(rbind,lapply(spec,as.data.frame));write.table(specframe,file.path(dest,'contrast_definitions.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
summaries<-list();toprows<-list();scores<-list()
programs<-read.delim(file.path(out,'marker_program_gene_coverage.tsv'),stringsAsFactors=FALSE)
pdf(file.path(dest,'selected_diagnostics.pdf'),width=8,height=6)
for(sp in spec) {
  cat('FITTING',stratum,sp$id,'\n');flush.console()
  if(sp$interaction) {
    use<-meta$gate %in% c(sp$ref,sp$case)
  } else {
    use<-(meta$design_group==sp$refgeno & meta$gate==sp$ref) | (meta$design_group==sp$casegeno & meta$gate==sp$case)
  }
  md<-meta[use,];xx<-counts[,use,drop=FALSE]
  required<-if(sp$interaction) 12L else 6L
  reason<-''
  if(nrow(md)!=required) reason<-'unexpected_library_number'
  if(any(md$cells<3 | md$UMI==0)) reason<-'fewer_than_3_cells_in_a_library'
  if(reason!='') {
    summaries[[length(summaries)+1]]<-data.frame(stratum=stratum,contrast=sp$id,status='SKIPPED',reason=reason,min_cells=min(md$cells),genes_tested=0,FDR05=0,FDR10=0,nominal01_concordant_FC05=0,minimum_FDR=NA)
    next
  }
  if(sp$interaction) {
    gate<-as.integer(md$gate==sp$case);het<-as.integer(md$design_group=='HET_F')
    design<-cbind(model.matrix(~0+factor(md$sample)),gate_effect=gate,HET_gate_interaction=gate*het)
    case_mask<-md$design_group=='HET_F' & md$gate==sp$case
    ref_mask<-md$design_group=='HET_F' & md$gate==sp$ref
  } else {
    case_mask<-md$design_group==sp$casegeno & md$gate==sp$case
    ref_mask<-!case_mask
    if(sp$paired) design<-cbind(model.matrix(~factor(md$sample)),case=as.integer(case_mask))
    else design<-model.matrix(~as.integer(case_mask))
  }
  result<-tryCatch({
    stopifnot(qr(design)$rank==ncol(design))
    y<-DGEList(xx,genes=genes)
    keep<-filterByExpr(y,design=design,min.count=5,min.total.count=10)
    if(sum(keep)<100) stop('fewer_than_100_expressed_genes')
    y<-calcNormFactors(y[keep,,keep.lib.sizes=FALSE],method='TMM')
    y<-estimateDisp(y,design,robust=TRUE)
    fit<-glmQLFit(y,design,robust=TRUE)
    test<-glmQLFTest(fit,coef=ncol(design))
    tab<-topTags(test,n=Inf,sort.by='PValue')$table
    lc<-cpm(y,log=TRUE,prior.count=2)
    if(sp$paired || sp$interaction) {
      ids<-intersect(md$sample[case_mask],md$sample[ref_mask])
      diff<-sapply(ids,function(s) lc[,md$sample==s & case_mask]-lc[,md$sample==s & ref_mask])
      if(sp$interaction) {
        wcase<-md$design_group=='WT_F' & md$gate==sp$case;wref<-md$design_group=='WT_F' & md$gate==sp$ref
        wids<-intersect(md$sample[wcase],md$sample[wref])
        wd<-sapply(wids,function(s) lc[,md$sample==s & wcase]-lc[,md$sample==s & wref])
        # Direction concordance here asks whether each HET donor's gate contrast
        # lies on the interaction's side of the average WT gate contrast.
        diff<-diff-rowMeans(wd)
      }
      agreement<-rowSums(diff>0)==ncol(diff) | rowSums(diff<0)==ncol(diff)
    } else agreement<-apply(lc[,case_mask,drop=FALSE],1,min)>apply(lc[,ref_mask,drop=FALSE],1,max) |
                       apply(lc[,case_mask,drop=FALSE],1,max)<apply(lc[,ref_mask,drop=FALSE],1,min)
    tab$direction_concordant<-agreement[match(tab$gene_id,rownames(lc))]
    tab$gate_gene<-tab$gene_symbol=='Pcdh19'
    gz<-gzfile(file.path(dest,paste0(sp$id,'.tsv.gz')),'wt');write.table(tab,gz,sep='\t',quote=FALSE,row.names=FALSE);close(gz)
    top<-head(tab[!tab$gate_gene,],30);top$stratum<-stratum;top$contrast<-sp$id
    toprows[[length(toprows)+1]]<-top
    for(p in unique(programs$program)) {
      ix<-tolower(tab$gene_symbol) %in% tolower(programs$gene[programs$program==p]) & !tab$gate_gene
      if(sum(ix)>=3) scores[[length(scores)+1]]<-data.frame(stratum=stratum,contrast=sp$id,program=p,genes=sum(ix),median_logFC=median(tab$logFC[ix]),fraction_positive=mean(tab$logFC[ix]>0),nominal01=sum(tab$PValue[ix]<.01))
    }
    if(sp$id %in% c('HET_all_vs_WTall','HET_Aneg_BC1_vs_WTall','matched_A_positive','within_HET_Aneg_BC1_vs_Apositive','gate_interaction_BC1_vs_Apositive')) {
      plotMDS(y,labels=paste(md$design_group,md$gate,md$sample,sep='\n'),col=ifelse(case_mask,'#A678AC','#549E89'),main=paste(stratum,sp$id,sep='\n'),cex=.6)
      g<-head(tab$gene_id[!tab$gate_gene],15)
      v<-lc[g,,drop=FALSE];rownames(v)<-tab$gene_symbol[match(g,tab$gene_id)]
      heatmap(v,Rowv=NA,Colv=NA,scale='row',margins=c(11,9),labCol=paste(md$design_group,md$gate,md$sample,sep='\n'),main=paste(stratum,sp$id,sep='\n'))
    }
    sel<-!tab$gate_gene
    data.frame(stratum=stratum,contrast=sp$id,status='COMPLETE',reason='',min_cells=min(md$cells),genes_tested=nrow(tab),FDR05=sum(tab$FDR[sel]<.05),FDR10=sum(tab$FDR[sel]<.1),nominal01_concordant_FC05=sum(tab$PValue[sel]<.01 & abs(tab$logFC[sel])>=.5 & tab$direction_concordant[sel]),minimum_FDR=min(tab$FDR[sel]))
  },error=function(e) data.frame(stratum=stratum,contrast=sp$id,status='FAILED',reason=conditionMessage(e),min_cells=min(md$cells),genes_tested=0,FDR05=0,FDR10=0,nominal01_concordant_FC05=0,minimum_FDR=NA))
  summaries[[length(summaries)+1]]<-result
}
dev.off()
write.table(do.call(rbind,summaries),file.path(dest,'summary.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
if(length(toprows)) write.table(do.call(rbind,toprows),file.path(dest,'top30_all_definitions.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
if(length(scores)) write.table(do.call(rbind,scores),file.path(dest,'program_effects.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
capture.output(sessionInfo(),file=file.path(dest,'R_session_info.txt'))
if(any(vapply(summaries,function(x) x$status=='FAILED',logical(1)))) quit(status=2)
writeLines('COMPLETE',file.path(dest,'COMPLETE.txt'))
