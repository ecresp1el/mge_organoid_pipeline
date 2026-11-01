# Paired HET DEG-first screen, curated signaling annotation, and cell-type expression support.
args<-commandArgs(trailingOnly=TRUE);base<-args[1];out<-args[2]
suppressPackageStartupMessages(library(edgeR))
dir.create(out,showWarnings=FALSE,recursive=TRUE)
e<-new.env();load(file.path(out,'CellChatDB.mouse.rda'),envir=e);db<-e$CellChatDB.mouse
write.table(db$interaction,file.path(out,'CellChat_mouse_interactions.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
write.table(data.frame(complex=rownames(db$complex),db$complex),file.path(out,'CellChat_mouse_complexes.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
expand<-function(x) {
 if(x %in% rownames(db$complex)) {
  v<-unlist(db$complex[x,,drop=FALSE],use.names=FALSE);return(unique(v[!is.na(v) & v!='']))
 }
 x
}
interactions<-db$interaction
if('annotation' %in% names(interactions)) interactions<-interactions[interactions$annotation!='Non-protein Signaling',]
rolemap<-list();pairs<-list()
for(i in seq_len(nrow(interactions))) {
 r<-interactions[i,];lg<-expand(r$ligand);rc<-expand(r$receptor)
 pairs[[i]]<-data.frame(interaction=r$interaction_name,pathway=r$pathway_name,ligand=r$ligand,receptor=r$receptor,
                         ligand_genes=paste(lg,collapse=';'),receptor_genes=paste(rc,collapse=';'),annotation=r$annotation)
 for(role in c('ligand','receptor'))for(g in if(role=='ligand')lg else rc) rolemap[[length(rolemap)+1]]<-data.frame(gene=g,role=role,pathway=r$pathway_name,interaction=r$interaction_name)
}
pairs<-do.call(rbind,pairs);roles<-unique(do.call(rbind,rolemap));write.table(pairs,file.path(out,'expanded_signaling_pairs.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
write.table(roles,file.path(out,'signaling_gene_annotations.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
meta<-read.delim(file.path(base,'pseudobulk_metadata.tsv'),stringsAsFactors=FALSE)
plan<-read.delim(file.path(base,'all_definition_DE_summary.tsv'),stringsAsFactors=FALSE)
strata<-unique(plan$stratum);allrows<-list();enrichment<-list();expression<-list();plotted<-list();ptstats<-list()
blood<-c('Ppbp','Pf4','Gp9','Gp1ba','Gp1bb','Itga2b','Itgb3','Mpl','Treml1','Hba-a1','Hba-a2','Hbb-bs','Hbb-bt','Hbb-bh1','Hbb-y','Slc4a1','Alas2','Gata1','Klf1')
sex<-c('Uty','Kdm5d','Ddx3y','Eif2s3y','Zfy1','Zfy2','Sry','Usp9y')
for(st in strata) {
 path<-file.path(base,paste0('DE_',st),'within_HET_Aneg_BC1_vs_Apositive.tsv.gz')
 if(!file.exists(path))next
 tab<-read.delim(gzfile(path),stringsAsFactors=FALSE)
 tab$stratum<-st
 tab$roles<-vapply(tab$gene_symbol,function(g)paste(unique(roles$role[roles$gene==g]),collapse=';'),character(1))
 tab$pathways<-vapply(tab$gene_symbol,function(g)paste(unique(roles$pathway[roles$gene==g]),collapse=';'),character(1))
 tab$flag<-ifelse(tab$gene_symbol=='Pcdh19','gate_gene',ifelse(tab$gene_symbol %in% blood,'blood_associated',ifelse(tab$gene_symbol %in% sex,'Y_linked','')))
 allrows[[st]]<-tab
 universe<-unique(tab$gene_symbol[tab$flag==''])
 hit<-unique(tab$gene_symbol[tab$FDR<.05 & tab$flag==''])
 for(p in unique(roles$pathway)) {
  gs<-intersect(universe,roles$gene[roles$pathway==p]);k<-length(intersect(hit,gs))
  if(length(gs)<3)next
  enrichment[[length(enrichment)+1]]<-data.frame(stratum=st,pathway=p,tested_genes=length(gs),DE_genes=k,total_DE=length(hit),universe=length(universe),PValue=phyper(k-1,length(gs),length(universe)-length(gs),length(hit),lower.tail=FALSE),genes=paste(intersect(hit,gs),collapse=';'))
 }
 d<-read.delim(gzfile(file.path(base,paste0(st,'_counts.tsv.gz'))),check.names=FALSE)
 md<-meta[match(names(d)[-(1:2)],meta$library),]
 keep<-md$design_group=='HET_F' & md$gate %in% c('A_positive','Aneg_BC1')
 md6<-md[keep,];xx<-as.matrix(d[,-(1:2)])[,keep,drop=FALSE];rownames(xx)<-d$gene_symbol
 if(any(md6$cells<3))next
 y<-calcNormFactors(DGEList(xx));cp<-cpm(y,log=FALSE)
 signaling<-intersect(unique(roles$gene),rownames(cp))
 for(g in signaling)for(gate in c('A_positive','Aneg_BC1')) {
  v<-cp[g,md6$gate==gate]
  expression[[length(expression)+1]]<-data.frame(stratum=st,gene=g,gate=gate,mean_CPM=mean(v),minimum_CPM=min(v),samples_CPM_ge1=sum(v>=1),sample_CPM=paste(signif(v,5),collapse=';'))
 }
 # Pcdh19 all-cell panel: common normalization for all three displayed groups.
 if(st=='all_cells') {
  use<-(md$design_group=='WT_F' & md$gate=='all') | (md$design_group=='HET_F' & md$gate %in% c('A_positive','Aneg_BC1'))
  m<-md[use,];x<-as.matrix(d[,-(1:2)])[,use,drop=FALSE];rownames(x)<-d$gene_symbol
  lc<-cpm(calcNormFactors(DGEList(x)),log=TRUE,prior.count=2)
  plotted[[1]]<-data.frame(gene='Pcdh19',stratum='all_cells',library=m$library,sample=m$sample,group=ifelse(m$design_group=='WT_F','WT_F_all',m$gate),cells=m$cells,raw_gene_UMI=x['Pcdh19',],TMM_log2CPM=lc['Pcdh19',])
  z<-tab[tab$gene_symbol=='Pcdh19',];z$gene<-'Pcdh19';z$contrast<-'within_HET_Aneg_BC1_vs_Apositive';ptstats[[1]]<-z
 }
}
annotated<-do.call(rbind,allrows);write.table(annotated,file.path(out,'paired_DE_all_celltypes_with_signaling_annotations.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
sig<-annotated[annotated$roles!='' & annotated$flag=='',];sig<-sig[order(sig$stratum,sig$PValue),]
write.table(sig,file.path(out,'paired_signaling_DE_all_celltypes.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
tops<-do.call(rbind,lapply(split(sig,sig$stratum),function(d)head(d,20)))
write.table(tops,file.path(out,'top20_signaling_genes_by_celltype.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
topgenes<-unique(annotated$gene_symbol[annotated$stratum=='all_cells' & annotated$FDR<.05 & annotated$flag==''])
follow<-annotated[annotated$gene_symbol %in% topgenes,];write.table(follow,file.path(out,'all_cell_DEGs_followed_by_celltype.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
if(length(enrichment)) {
 enr<-do.call(rbind,enrichment);enr$FDR<-ave(enr$PValue,enr$stratum,FUN=function(v)p.adjust(v,'BH'))
 enr<-enr[order(enr$stratum,enr$PValue),];write.table(enr,file.path(out,'signaling_pathway_overrepresentation.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
}
expression<-do.call(rbind,expression);write.table(expression,file.path(out,'signaling_gene_expression_support.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
write.table(plotted[[1]],file.path(out,'Pcdh19_three_group_sample_values.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
write.table(ptstats[[1]],file.path(out,'Pcdh19_paired_statistics.tsv'),sep='\t',quote=FALSE,row.names=FALSE)

# Curated pair support: expression availability, not inferred communication probability.
supported<-list()
for(st in unique(expression$stratum)) {
 ex<-expression[expression$stratum==st,]
 get<-function(gs,gate) {
  rr<-ex[ex$gene %in% gs & ex$gate==gate,]
  if(nrow(rr)!=length(gs) || any(rr$samples_CPM_ge1<2))return(NA_real_)
  exp(mean(log(pmax(rr$mean_CPM,1e-12))))
 }
 de<-annotated[annotated$stratum==st,]
 for(i in seq_len(nrow(pairs))) {
  pair<-pairs[i,];lg<-strsplit(pair$ligand_genes,';',fixed=TRUE)[[1]];rg<-strsplit(pair$receptor_genes,';',fixed=TRUE)[[1]]
  hits<-de[de$gene_symbol %in% c(lg,rg) & de$flag=='' & (de$FDR<.1 | de$PValue<.01),]
  if(!nrow(hits))next
  for(direction in c('WT_like_to_KO_enriched','KO_enriched_to_WT_like')) {
   a<-get(lg,if(direction=='WT_like_to_KO_enriched')'A_positive' else 'Aneg_BC1')
   b<-get(rg,if(direction=='WT_like_to_KO_enriched')'Aneg_BC1' else 'A_positive')
   if(is.na(a) || is.na(b))next
   supported[[length(supported)+1]]<-data.frame(stratum=st,interaction=pair$interaction,pathway=pair$pathway,ligand=pair$ligand,receptor=pair$receptor,direction=direction,
     ligand_CPM_geomean=a,receptor_CPM_geomean=b,expression_support=sqrt(a*b),DE_genes=paste(hits$gene_symbol,collapse=';'),best_P=min(hits$PValue),best_FDR=min(hits$FDR))
  }
 }
}
if(length(supported)) {
 sp<-do.call(rbind,supported);sp<-sp[order(sp$stratum,sp$best_FDR,sp$best_P),]
 write.table(sp,file.path(out,'curated_signaling_pair_expression_support.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
}
overview<-do.call(rbind,lapply(split(annotated,annotated$stratum),function(d)data.frame(stratum=d$stratum[1],genes_tested=nrow(d),DE_FDR05=sum(d$FDR<.05 & d$flag!='gate_gene'),signaling_DE_FDR05=sum(d$FDR<.05 & d$roles!='' & d$flag==''),signaling_nominal01=sum(d$PValue<.01 & d$roles!='' & d$flag==''))))
write.table(overview,file.path(out,'paired_DE_overview.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
capture.output(sessionInfo(),file=file.path(out,'R_session_info.txt'));writeLines('COMPLETE',file.path(out,'ANALYSIS_COMPLETE.txt'))
print(overview)
