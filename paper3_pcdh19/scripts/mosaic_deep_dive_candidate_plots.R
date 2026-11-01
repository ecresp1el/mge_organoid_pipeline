args<-commandArgs(trailingOnly=TRUE)
out<-args[1]
suppressPackageStartupMessages(library(edgeR))
meta<-read.delim(file.path(out,'pseudobulk_metadata.tsv'),stringsAsFactors=FALSE)
candidates<-data.frame(
 gene=c('Gabrr2','Kcnc2','Pcdhga12','Cacng8','Gria1','Gadd45g','Unc5b','Cxcl12','Htra1'),
 stratum=c('neural_vascular_low','C0006','C0003','C0003','neural_depth_trimmed','C0001_neural_vascular_low','vascular_marked','vascular_marked','vascular_marked'),
 label=c('Neural; vascular-low','Immature MGE neurons (C0006)','MGE IPCs (C0003)','MGE IPCs (C0003)','Neural; depth-trimmed','MGE progenitors (C0001)','Vascular-marker positive','Vascular-marker positive','Vascular-marker positive'),stringsAsFactors=FALSE)
comparisons<-c('HET_A_positive_vs_WTall','HET_Aneg_BC1_vs_WTall','matched_A_positive','matched_Aneg_BC1','within_HET_Aneg_BC1_vs_Apositive','gate_interaction_BC1_vs_Apositive')
points<-list();stats<-list()
for(st in unique(candidates$stratum)) {
 d<-read.delim(gzfile(file.path(out,paste0(st,'_counts.tsv.gz'))),check.names=FALSE)
 md<-meta[match(names(d)[-(1:2)],meta$library),]
 use<-(md$design_group=='WT_F' & md$gate=='all') | (md$design_group=='HET_F' & md$gate %in% c('A_positive','Aneg_BC1'))
 md<-md[use,];xx<-as.matrix(d[,-(1:2)])[,use,drop=FALSE];rownames(xx)<-d$gene_id
 y<-calcNormFactors(DGEList(xx),method='TMM');lc<-cpm(y,log=TRUE,prior.count=2)
 for(g in candidates$gene[candidates$stratum==st]) {
  id<-d$gene_id[match(g,d$gene_symbol)];stopifnot(length(id)==1)
  rr<-data.frame(gene=g,stratum=st,library=md$library,sample=md$sample,group=ifelse(md$design_group=='WT_F','WT_F_all',md$gate),cells=md$cells,raw_gene_UMI=xx[id,],TMM_log2CPM=lc[id,],row.names=NULL)
  points[[length(points)+1]]<-rr
  for(co in comparisons) {
   tab<-read.delim(gzfile(file.path(out,paste0('DE_',st),paste0(co,'.tsv.gz'))),stringsAsFactors=FALSE)
   r<-tab[tab$gene_symbol==g,,drop=FALSE]
   r$gene<-g;r$stratum<-st;r$contrast<-co
   stats[[length(stats)+1]]<-r
  }
 }
}
points<-do.call(rbind,points);stats<-do.call(rbind,stats)
write.table(points,file.path(out,'selected_biological_candidates_sample_values.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
write.table(stats,file.path(out,'selected_biological_candidates_statistics.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
colors<-c('#549E89','#337EAC','#A678AC')
draw<-function(g,st,title,letter) {
 r<-points[points$gene==g & points$stratum==st,]
 groups<-c('WT_F_all','A_positive','Aneg_BC1')
 yrange<-range(r$TMM_log2CPM);pad<-max(diff(yrange)*.18,.25)
 plot(NA,xlim=c(.6,3.4),ylim=yrange+c(-pad,pad),xaxt='n',xlab='',ylab='log2 normalized CPM',bty='l',main=paste(g,title,sep='\n'),cex.main=.9)
 axis(1,at=1:3,labels=c('WT F','HET F\nWT-like','HET F\nKO-enriched'),tick=FALSE,cex.axis=.75)
 for(s in unique(r$sample[r$group=='A_positive'])) {
  a<-r$TMM_log2CPM[r$sample==s & r$group=='A_positive'];b<-r$TMM_log2CPM[r$sample==s & r$group=='Aneg_BC1']
  lines(c(2,3),c(a,b),col='#bbbbbb',lwd=.8)
 }
 for(i in 1:3) {
  vals<-r$TMM_log2CPM[r$group==groups[i]]
  points(i+c(-.1,0,.1),vals,pch=21,bg=colors[i],col='white',cex=1.05,lwd=.5)
  segments(i-.22,mean(vals),i+.22,mean(vals),col=colors[i],lwd=2)
 }
 s<-stats[stats$gene==g & stats$stratum==st & stats$contrast=='within_HET_Aneg_BC1_vs_Apositive',]
 mtext(sprintf('Within HET: log2FC %+.2f; P %.2g; q %.2g',s$logFC,s$PValue,s$FDR),side=3,line=.2,cex=.65)
 mtext(letter,side=3,line=2.4,adj=-.12,font=2,cex=1.2)
}
pdf(file.path(out,'exploratory_neuronal_mosaic_candidates.pdf'),width=10,height=7.3)
par(mfrow=c(2,3),mar=c(4.2,4,4.4,1),oma=c(2.2,0,1,0),mgp=c(2.2,.7,0),cex=.9)
for(i in 1:6)draw(candidates$gene[i],candidates$stratum[i],candidates$label[i],letters[i])
mtext('Exploratory candidates selected after screening. Each dot: one sample; gray lines pair the two HET subsets.',side=1,outer=TRUE,line=.4,cex=.75)
mtext('WT-like = A-positive; KO-enriched = A-negative, B/C-positive. Small q values are not required for this exploratory display.',side=1,outer=TRUE,line=1.4,cex=.7)
dev.off()
pdf(file.path(out,'exploratory_vascular_mosaic_candidates.pdf'),width=10,height=3.8)
par(mfrow=c(1,3),mar=c(4.2,4,4.4,1),oma=c(2.2,0,1,0),mgp=c(2.2,.7,0),cex=.9)
for(i in 7:9)draw(candidates$gene[i],candidates$stratum[i],candidates$label[i],letters[i-6])
mtext('Broad probe gates within vascular-marker-positive cells; three samples/group. Gray lines pair HET subsets.',side=1,outer=TRUE,line=.4,cex=.75)
mtext('WT/HET gate-interaction controls and all comparison statistics are supplied separately.',side=1,outer=TRUE,line=1.4,cex=.7)
dev.off()
