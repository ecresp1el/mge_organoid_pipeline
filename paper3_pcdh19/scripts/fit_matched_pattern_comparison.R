args<-commandArgs(trailingOnly=TRUE);base<-args[1];out<-args[2]
suppressPackageStartupMessages(library(edgeR))
d<-read.delim(gzfile(file.path(base,'C0006_counts.tsv.gz')),check.names=FALSE)
md<-read.delim(file.path(base,'pseudobulk_metadata.tsv'),stringsAsFactors=FALSE);md<-md[match(names(d)[-(1:2)],md$library),]
use<-(md$design_group=='WT_F' & md$gate=='A_positive')|(md$design_group=='HET_F' & md$gate=='Aneg_BC1')
xx<-as.matrix(d[,-(1:2)])[,use];m<-md[use,];case<-as.integer(m$design_group=='HET_F');design<-model.matrix(~case)
y<-DGEList(xx,genes=d[,1:2]);keep<-filterByExpr(y,design,min.count=5,min.total.count=10);y<-calcNormFactors(y[keep,,keep.lib.sizes=FALSE]);y<-estimateDisp(y,design,robust=TRUE);fit<-glmQLFit(y,design,robust=TRUE)
tab<-topTags(glmQLFTest(fit,coef=2),n=Inf)$table;write.table(tab,file.path(out,'KOlike_vs_matched_WT_DE.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
