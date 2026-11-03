args<-commandArgs(trailingOnly=TRUE);base<-args[1];out<-args[2];dir.create(out,recursive=TRUE,showWarnings=FALSE)
suppressPackageStartupMessages(library(edgeR))
meta<-read.delim(file.path(base,'pseudobulk_metadata.tsv'),stringsAsFactors=FALSE)
d<-read.delim(gzfile(file.path(base,'C0006_counts.tsv.gz')),check.names=FALSE)
meta<-meta[match(names(d)[-(1:2)],meta$library),];use<-meta$design_group=='HET_F' & meta$gate %in% c('A_positive','Aneg_BC1')
md<-meta[use,];xx<-as.matrix(d[,-(1:2)])[,use,drop=FALSE];rownames(xx)<-d$gene_id
stopifnot(nrow(md)==6,all(md$cells>0));design<-cbind(model.matrix(~factor(md$sample)),case=as.integer(md$gate=='Aneg_BC1'))
y<-DGEList(xx);keep<-filterByExpr(y,design=design,min.count=5,min.total.count=10);filtered<-calcNormFactors(y[keep,,keep.lib.sizes=FALSE],method='TMM')
# Use the identical effective library sizes as the existing paired-DE normalization;
# also display requested low-count genes that were excluded from that DE fit.
y$samples$lib.size<-filtered$samples$lib.size;y$samples$norm.factors<-filtered$samples$norm.factors
cp<-cpm(y);lc<-cpm(y,log=TRUE,prior.count=2);result<-list()
for(gene in c('Egr1','Npas4','Nr4a1')) {
 i<-which(d$gene_symbol==gene);stopifnot(length(i)==1)
 for(j in seq_len(nrow(md)))result[[length(result)+1]]<-data.frame(gene=gene,sample=md$sample[j],group=ifelse(md$gate[j]=='A_positive','WT cells','KO cells'),cells=md$cells[j],raw_pseudobulk_counts=xx[i,j],CPM=cp[i,j],log2CPM=lc[i,j],included_in_existing_DE=keep[i],effective_library_size=y$samples$lib.size[j]*y$samples$norm.factors[j])
}
write.table(do.call(rbind,result),file.path(out,'paired_IEG_expression.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
capture.output(sessionInfo(),file=file.path(out,'R_session_info.txt'))
