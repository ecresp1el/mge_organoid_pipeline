args<-commandArgs(trailingOnly=TRUE);base<-args[1];out<-args[2]
suppressPackageStartupMessages(library(edgeR))
meta<-read.delim(file.path(base,'pseudobulk_metadata.tsv'),stringsAsFactors=FALSE)
d<-read.delim(gzfile(file.path(base,'C0006_counts.tsv.gz')),check.names=FALSE)
md<-meta[match(names(d)[-(1:2)],meta$library),]
use<-(md$design_group=='WT_F' & md$gate=='A_positive') | (md$design_group=='HET_F' & md$gate %in% c('A_positive','Aneg_BC1'))
md<-md[use,];xx<-as.matrix(d[,-(1:2)])[,use,drop=FALSE];rownames(xx)<-d$gene_id
y<-calcNormFactors(DGEList(xx),method='TMM');lc<-cpm(y,log=TRUE,prior.count=2);cp<-cpm(y)
write.table(data.frame(gene_id=d$gene_id,gene_symbol=d$gene_symbol,lc,check.names=FALSE),file.path(out,'nine_sample_log2CPM.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
write.table(md,file.path(out,'nine_sample_metadata.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
keep<-rowSums(cp>=1)>=3 & rowSums(xx)>=10 & d$gene_symbol!='Pcdh19'
v<-apply(lc[keep,,drop=FALSE],1,var);hvg<-names(sort(v,decreasing=TRUE))[seq_len(min(60,length(v)))];
write.table(data.frame(gene_id=hvg,variance=v[hvg]),file.path(out,'top60_variable_genes.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
write.table(data.frame(gene_id=rownames(lc)[keep]),file.path(out,'expressed_gene_universe.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
capture.output(sessionInfo(),file=file.path(out,'R_session_info.txt'));writeLines('COMPLETE',file.path(out,'EXPRESSION_COMPLETE.txt'))
