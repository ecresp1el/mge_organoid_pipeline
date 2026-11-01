args<-commandArgs(trailingOnly=TRUE);out<-args[1];suppressPackageStartupMessages(library(edgeR))
meta<-read.delim(file.path(out,'male_metadata.tsv'),stringsAsFactors=FALSE);summary<-list()
for(st in unique(meta$stratum)) {
 d<-read.delim(file.path(out,paste0('male_',st,'_counts.tsv')),check.names=FALSE);xx<-as.matrix(d[,-(1:2)]);md<-meta[match(colnames(xx),meta$library),];case<-as.integer(md$group=='KO male');design<-model.matrix(~case)
 y<-DGEList(xx,genes=d[,1:2]);keep<-filterByExpr(y,design,min.count=5,min.total.count=10);y<-calcNormFactors(y[keep,,keep.lib.sizes=FALSE]);y<-estimateDisp(y,design,robust=TRUE);fit<-glmQLFit(y,design,robust=TRUE);tab<-topTags(glmQLFTest(fit,coef=2),n=Inf)$table
 write.table(tab,file.path(out,paste0('male_',st,'_DE.tsv')),sep='\t',quote=FALSE,row.names=FALSE)
 cp<-cpm(calcNormFactors(DGEList(xx)),log=FALSE);write.table(data.frame(gene_id=d$gene_id,gene_symbol=d$gene_symbol,cp,check.names=FALSE),file.path(out,paste0('male_',st,'_CPM.tsv')),sep='\t',quote=FALSE,row.names=FALSE)
 summary[[st]]<-data.frame(stratum=st,genes=nrow(tab),FDR05=sum(tab$FDR<.05 & tab$gene_symbol!='Pcdh19'),FDR10=sum(tab$FDR<.1 & tab$gene_symbol!='Pcdh19'),minimum_FDR=min(tab$FDR[tab$gene_symbol!='Pcdh19']))
}
write.table(do.call(rbind,summary),file.path(out,'male_DE_summary.tsv'),sep='\t',quote=FALSE,row.names=FALSE);writeLines('COMPLETE',file.path(out,'MALE_DE_COMPLETE.txt'))
