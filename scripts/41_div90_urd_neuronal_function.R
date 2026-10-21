#!/usr/bin/env Rscript
# Read-only extension: existing URD statistics are reused byte-for-byte.
suppressPackageStartupMessages({library(URD);library(Matrix);library(data.table);library(ggplot2)})
data.table::setDTthreads(1L)
args <- commandArgs(TRUE); stopifnot(length(args)==3L)
run <- normalizePath(args[1]); previous <- normalizePath(args[2]); source <- normalizePath(args[3])
log <- function(...) {cat(format(Sys.time(),'%FT%T%z'),..., '\n');flush.console()}
out <- function(x,n) fwrite(x,file.path(run,'tables',n),na='NA')
stamp <- format(Sys.time(),'%FT%T%z'); run_id <- basename(run)
log('START annotation/statistics reuse; no differential-expression rerun')
ann <- fread(file.path(run,'config/functional_gene_sets.tsv'))
stat <- fread(file.path(previous,'DIV90_URD_gene_results.csv'))
stopifnot(setequal(unique(stat$branchpoint),c('BP11','BP9')), uniqueN(stat$comparison_id)==6L)
stopifnot(all(stat$candidate == (stat$log2fc>=.25 & stat$global_fdr<=.05)))
annotation <- ann[,.(functional_categories=paste(sort(unique(unlist(strsplit(categories,'|',fixed=TRUE)))),collapse='|'),
  functional_sets=paste(sort(unique(set_id)),collapse='|'),annotation_sources=paste(sort(unique(source)),collapse='|')),by=gene]
# Avoid masking the protected source path with the annotation column named source.
master <- merge(stat,annotation,by='gene',all.x=TRUE,sort=FALSE)
master[,channel_family:=fcase(grepl('^SCN',gene),'Sodium',grepl('^KCN',gene),'Potassium',
 grepl('^CACN',gene),'Calcium',grepl('^HCN',gene),'HCN',grepl('^CLCN|^CLIC|^ANO|^BEST',gene),'Chloride/anion-associated',default='Other/see ontology')]
master[,opposing_daughter:=ifelse(reported_daughter==daughter_A,daughter_B,daughter_A)]
flags <- master[candidate==TRUE,.(siblings_beaten=uniqueN(opposing_daughter)),by=.(branchpoint,reported_daughter,gene)]
if('siblings_beaten' %in% names(master))master[,siblings_beaten:=NULL]
master <- merge(master,flags,by=c('branchpoint','reported_daughter','gene'),all.x=TRUE,sort=FALSE)
eligibility<-master[,.(sibling_tests_available=uniqueN(comparison_id)),by=.(branchpoint,reported_daughter,gene)]
master<-merge(master,eligibility,by=c('branchpoint','reported_daughter','gene'),all.x=TRUE,sort=FALSE)
master[is.na(siblings_beaten),siblings_beaten:=0L]
master[,sibling_status:=fifelse(siblings_beaten==2L,'higher_than_both',fifelse(siblings_beaten==1L,'higher_than_one','not_significant'))]
master[,`:=`(percent_expressing_A=100*detection_A,percent_expressing_B=100*detection_B,
  n_expressing_A=as.integer(round(detection_A*n_A)),n_expressing_B=as.integer(round(detection_B*n_B)),
  expression_scope='original_matched_DE_window',functional_run_id=run_id)]
out(master,'DIV90_neuronal_function_master.csv')
for(cat in c('ion_channel','synaptic','activity','maturation')) {
  g <- unique(ann[grepl(cat,categories,fixed=TRUE),gene])
  out(master[gene %in% g],paste0(cat,'_all_tested_genes.csv'))
}
out(ann,'functional_gene_annotations.csv')
log('DE_REUSE_DONE',nrow(master),'direction rows; FDR and candidate values preserved')
log('READ_TREE',source)
obj <- readRDS(source)
stopifnot(inherits(obj,'URD'))
expr <- obj@logupx.data; genes <- rownames(expr); cells <- colnames(expr)
pt <- setNames(obj@pseudotime[cells,'div90_jia_rootscore_root'],cells)
stopifnot(max(abs(obj@tree$pseudotime[cells]-pt))<1e-10)
members <- lapply(obj@tree$cells.in.segment,function(x)intersect(as.character(x),cells))
branches <- c('3','5','9','1','2','4'); segments <- c('12','11',branches)
labels <- c('3'='CRABP1/ANGPT2','5'='NR2F1/NR2F2','9'='Shared neuronal',
 '1'='LHX8/ISL1 state 1','2'='LHX8/ISL1 state 2','4'='LHX6/MEF2C cortical',
 '11'='Shared ancestor BP11','12'='Earlier shared ancestor')
palette <- c('3'='#D55E00','5'='#CC79A7','9'='#009E73','1'='#0072B2','2'='#56B4E9','4'='#E69F00','11'='#666666','12'='#999999')
selected <- intersect(unique(ann$gene),genes)
coverage <- data.table(gene=unique(ann$gene));coverage[,present_in_tree:=gene %in% genes]
coverage[,tested_in_existing_DE:=gene %in% stat$gene];out(coverage,'annotation_coverage.csv')
log('FUNCTIONAL_GENES',length(selected),'of',length(genes),'object features')
samples <- setNames(as.character(obj@meta[cells,'orig.ident']),cells)
summaries <- list(); trajectories <- list(); coexpression <- list(); score_rows <- list()
summarize <- function(cc,gs) {
  x <- expr[gs,cc,drop=FALSE]; lin <- x;lin@x <- 2^lin@x-1
  data.table(gene=gs,n_cells=length(cc),n_expressing=as.integer(Matrix::rowSums(x>0)),
    percent_expressing=100*Matrix::rowMeans(x>0),mean_normalized=Matrix::rowMeans(lin),mean_log2upx=Matrix::rowMeans(x))
}
# Scores use one fixed all-neuronal reference scale, never separately scale each branch.
neuron_cells <- unique(unlist(members[c('11',branches)]))
module_defs <- split(ann[grepl('^PANEL:',set_id),gene],ann[grepl('^PANEL:',set_id),set_name])
module_defs$curated_ion_channels <- unique(ann[grepl('ion_channel',categories,fixed=TRUE),gene])
for(name in names(module_defs)) {
  gs <- intersect(unique(module_defs[[name]]),genes)
  x <- expr[gs,neuron_cells,drop=FALSE];mu<-Matrix::rowMeans(x)
  vari<-pmax(0,Matrix::rowMeans(x^2)-mu^2);keep<-vari>1e-8
  gs<-gs[keep];mu<-mu[keep];sdv<-sqrt(vari[keep]);stopifnot(length(gs)>1)
  for(s in branches) {
    cc<-members[[s]];xx<-expr[gs,cc,drop=FALSE]
    z<-as.matrix(xx);z<-(z-mu)/sdv
    score_rows[[length(score_rows)+1L]]<-data.table(cell=cc,segment=s,sample=samples[cc],pseudotime=pt[cc],
      module=name,score=colMeans(z),genes_detected=as.integer(Matrix::colSums(xx>0)),n_module_genes=length(gs),
      library_counts=as.numeric(Matrix::colSums(obj@count.data[,cc,drop=FALSE])))
  }
}
scores<-rbindlist(score_rows);out(scores,'cell_module_scores.csv')
score_qc<-scores[,.(n_cells=.N,mean_score=mean(score),median_score=median(score),
  score_library_spearman=if(uniqueN(library_counts)>1&&uniqueN(score)>1)cor(score,library_counts,method='spearman')else NA_real_),by=.(segment,sample,module)]
out(score_qc,'module_sample_depth_summary.csv')
for(s in segments) {
  cc<-members[[s]];if(length(cc)<5)next
  log('SEGMENT_START',s,'cells',length(cc),'genes',length(selected))
  bins<-cut(pt[cc],seq(min(pt[cc]),max(pt[cc]),length.out=13),include.lowest=TRUE,labels=FALSE)
  for(sample in c('ALL',sort(unique(samples[cc])))) {
    use<-if(sample=='ALL')cc else cc[samples[cc]==sample];if(length(use)<2)next
    d<-summarize(use,selected);d[,`:=`(segment=s,sample=sample,scope='whole_native_segment')]
    summaries[[length(summaries)+1L]]<-d
    for(k in sort(unique(bins))) {
      bc<-intersect(use,cc[bins==k]);if(length(bc)<5)next
      d<-summarize(bc,selected);d[,`:=`(segment=s,sample=sample,bin=k,pseudotime=mean(pt[bc]),
        time_low=min(pt[bc]),time_high=max(pt[bc]))]
      trajectories[[length(trajectories)+1L]]<-d
    }
  }
  log('SEGMENT_DONE',s)
}
summary<-rbindlist(summaries);traj<-rbindlist(trajectories)
setorder(traj,segment,sample,gene,pseudotime)
trends<-traj[,.(n_bins=.N,first_pseudotime=first(pseudotime),last_pseudotime=last(pseudotime),
  first_mean_log2upx=first(mean_log2upx),last_mean_log2upx=last(mean_log2upx),
  delta_mean_log2upx=last(mean_log2upx)-first(mean_log2upx),
  first_percent_expressing=first(percent_expressing),last_percent_expressing=last(percent_expressing),
  pseudotime_expression_spearman=if(.N>2&&sd(mean_log2upx)>0)cor(pseudotime,mean_log2upx,method='spearman')else NA_real_),by=.(segment,sample,gene)]
out(trends,'descriptive_pseudotime_trends.csv')
out(summary,'branch_expression_summary.csv');fwrite(traj,file.path(run,'tables/pseudotime_expression_trajectories.tsv.gz'),sep='\t')
whole<-summary[sample=='ALL' & segment %in% branches,.(reported_daughter=segment,gene,
  whole_segment_n_cells=n_cells,whole_segment_percent_expressing=percent_expressing,
  whole_segment_mean_normalized=mean_normalized,whole_segment_mean_log2upx=mean_log2upx)]
whole[,reported_daughter:=as.character(reported_daughter)];master[,reported_daughter:=as.character(reported_daughter)]
master<-merge(master,whole,by=c('reported_daughter','gene'),all.x=TRUE,sort=FALSE)
out(master,'DIV90_neuronal_function_master.csv')
core<-intersect(c('FOS','EGR1','NPAS4','NR4A1','ARC'),genes)
for(s in branches) {
  cc<-members[[s]];n<-as.integer(Matrix::colSums(expr[core,cc,drop=FALSE]>0))
  coexpression[[s]]<-data.table(segment=s,cell=cc,sample=samples[cc],pseudotime=pt[cc],
    core_IEGs_detected=n,n_available_core_IEGs=length(core))
}
core_cells<-rbindlist(coexpression);out(core_cells,'core_IEG_cell_coexpression.csv')
core_distribution<-core_cells[,.(n_cells=.N),by=.(segment,sample,core_IEGs_detected)]
core_distribution[,percent_cells:=100*n_cells/sum(n_cells),by=.(segment,sample)]
out(core_distribution,'core_IEG_coexpression_distribution.csv')
log('ENRICHMENT_START')
sets<-split(ann$gene,ann$set_id);sets<-lapply(sets,unique)
setmeta<-unique(ann[,.(set_id,set_name,source)])
enrich<-list()
for(bp in c('BP11','BP9'))for(s in unique(master[branchpoint==bp,reported_daughter])) {
  d<-master[branchpoint==bp & reported_daughter==s]
  eligible<-d[,.(comparisons=uniqueN(comparison_id)),by=gene]
  for(mode in c('both','one_only')) {
    bg<-if(mode=='both')eligible[comparisons==2,gene]else eligible$gene
    required<-if(mode=='both')2L else 1L
    hits<-unique(d[gene %in% bg & siblings_beaten==required,gene])
    for(id in names(sets)) {
      gs<-intersect(sets[[id]],bg);if(length(gs)<5||length(gs)>2000)next
      overlap<-intersect(hits,gs);N<-length(bg);K<-length(gs);n<-length(hits);k<-length(overlap)
      enrich[[length(enrich)+1L]]<-data.table(branchpoint=bp,segment=s,sibling_category=mode,set_id=id,
        background_genes=N,hit_genes=n,set_background_genes=K,overlap=k,overlap_genes=paste(sort(overlap),collapse='|'),
        fold_enrichment=if(n>0)k/(n*K/N)else NA_real_,p_value=phyper(k-1,K,N-K,n,lower.tail=FALSE))
    }
  }
}
enrichment<-rbindlist(enrich);enrichment[,adjusted_p:=p.adjust(p_value,'BH')]
enrichment<-merge(enrichment,setmeta,by='set_id',all.x=TRUE);out(enrichment,'functional_gene_set_enrichment.csv')
log('ENRICHMENT_DONE',nrow(enrichment),'tests; one global BH family')
caption<-paste(run_id,'|',stamp,'| Original URD pseudotime; descriptive expression; no functional activity measured.')
save_plot<-function(p,name,w=13,h=10) {
  log('PLOT',name)
  for(ext in c('pdf','png','svg')) {
    f<-file.path(run,'plots',paste0(name,'.',ext))
    if(ext=='pdf')ggsave(f,p,width=w,height=h,device=grDevices::cairo_pdf,limitsize=FALSE)
    if(ext=='png')ggsave(f,p,width=w,height=h,dpi=300,bg='white',limitsize=FALSE)
    if(ext=='svg')ggsave(f,p,width=w,height=h,device=svglite::svglite,limitsize=FALSE)
  }
}
avg<-summary[sample=='ALL' & segment %in% branches]
avg[,segment:=factor(segment,levels=branches)]
dot_pages<-function(gs,prefix,title,terminal=FALSE) {
  gs<-sort(intersect(gs,selected));if(!length(gs)){log('NO_GENES',prefix);return(invisible(NULL))}
  pages<-split(gs,ceiling(seq_along(gs)/45))
  for(i in seq_along(pages)) {
    d<-avg[gene %in% pages[[i]]];if(terminal)d<-d[segment!='9']
    d[,z:=if(sd(mean_log2upx)>0)(mean_log2upx-mean(mean_log2upx))/sd(mean_log2upx)else 0,by=gene]
    p<-ggplot(d,aes(segment,factor(gene,levels=rev(pages[[i]]))))+
      geom_point(aes(size=percent_expressing,color=mean_log2upx))+
      scale_size_area(max_size=5,limits=c(0,100))+scale_color_viridis_c()+
      scale_x_discrete(labels=labels)+theme_bw(base_size=10)+theme(axis.text.x=element_text(angle=35,hjust=1))+
      labs(title=title,subtitle='Whole native segments; branch 9 is shared internal context',
        x=NULL,y=NULL,size='% expressing',color='Mean log2(UPX+1)',caption=caption)
    save_plot(p,sprintf('%s_page%02d',prefix,i),13,max(7,length(pages[[i]])*.18+3))
  }
}
ranked<-copy(master[candidate==TRUE & grepl('ion_channel|synaptic|activity',functional_categories)])
setorder(ranked,global_fdr,-log2fc)
# Bounded figure panels; complete statistics and expression tables remain exhaustive.
plot_selection<-ranked[,head(unique(gene),15),by=.(branchpoint,reported_daughter)]
setnames(plot_selection,'V1','gene');out(plot_selection,'figure1_gene_selection.csv')
sig<-unique(plot_selection$gene)
for(i in seq_along(split(sort(sig),ceiling(seq_along(sort(sig))/50)))) {
  gs<-split(sort(sig),ceiling(seq_along(sort(sig))/50))[[i]];d<-avg[gene %in% gs]
  d[,z:=if(sd(mean_log2upx)>0)(mean_log2upx-mean(mean_log2upx))/sd(mean_log2upx)else 0,by=gene]
  p<-ggplot(d,aes(segment,factor(gene,levels=rev(gs)),fill=z))+geom_tile()+
    scale_fill_gradient2(low='#2166AC',mid='white',high='#B2182B')+scale_x_discrete(labels=labels)+
    theme_minimal(base_size=10)+theme(axis.text.x=element_text(angle=35,hjust=1))+
    labs(title='Differential functional genes across neuronal segments',subtitle='Gene-wise standardized whole-segment mean expression; significance comes from matched windows',x=NULL,y=NULL,fill='Gene z score',caption=caption)
  save_plot(p,sprintf('Figure1_functional_heatmap_page%02d',i),13,max(8,length(gs)*.17+3))
}
dot_pages(unique(ann[grepl('ion_channel',categories,fixed=TRUE),gene]),'Figure2_channels','Ion-channel expression across terminal neuronal branches',TRUE)
synaptic_genes<-unique(ann[grepl('synaptic',categories,fixed=TRUE),gene])
synaptic_rank<-master[gene %in% synaptic_genes & candidate==TRUE][order(global_fdr,-log2fc)]
synaptic_selection<-synaptic_rank[,head(unique(gene),10),by=.(branchpoint,reported_daughter)]
setnames(synaptic_selection,'V1','gene');out(synaptic_selection,'figure3_gene_selection.csv')
dot_pages(unique(c(synaptic_selection$gene,intersect(c('GAD1','GAD2','SLC6A1','SLC32A1','GABRA1','GABRA2','GABRB2','GABRG2','GRIA1','GRIA2','GRIN1','GRIN2A','GRIN2B','SYT1','SYN1','SNAP25','STX1A','VAMP2','SYP','NRXN1','NRXN2','NLGN2','NLGN3','DLG4','SHANK2','SHANK3'),selected))), 'Figure3_synaptic','Synaptic transmission: ranked DE genes and requested landmarks')
dot_pages(unique(ann[grepl('^PANEL:',set_id) & grepl('activity',categories,fixed=TRUE),gene]),'Figure4_activity','Immediate-early and broader activity-associated genes')
module_plot<-ggplot(scores,aes(factor(segment,levels=branches),score,fill=segment))+geom_violin(scale='width',linewidth=.25)+
  facet_wrap(~module,scales='free_y',ncol=2)+scale_fill_manual(values=palette)+scale_x_discrete(labels=labels)+
  theme_bw(base_size=10)+theme(axis.text.x=element_text(angle=35,hjust=1),legend.position='none')+
  labs(title='Maturation and activity-response scores shown separately',subtitle='Descriptive scores using a fixed neuronal reference; inspect individual genes and sample summaries',x=NULL,y='Mean standardized gene expression',caption=caption)
save_plot(module_plot,'Figure4b_maturation_vs_activity_modules',14,12)
representatives<-intersect(c('SCN1A','KCNC1','KCNC2','GAD1','SLC32A1','ANK3','ATP1A3','SNAP25','SYT1',
 'FOS','EGR1','NPAS4','NR4A1','ARC','MEF2C','CRABP1'),genes)
# CRABP1 is a requested lineage landmark even if outside functional sets.
missing_rep<-setdiff(representatives,selected)
if(length(missing_rep)) {
  extra<-list()
  for(s in segments){cc<-members[[s]];if(length(cc)<5)next
    bins<-cut(pt[cc],seq(min(pt[cc]),max(pt[cc]),length.out=13),include.lowest=TRUE,labels=FALSE)
    for(k in sort(unique(bins))){bc<-cc[bins==k];if(length(bc)<5)next
      d<-summarize(bc,missing_rep);d[,`:=`(segment=s,sample='ALL',bin=k,pseudotime=mean(pt[bc]),time_low=min(pt[bc]),time_high=max(pt[bc]))];extra[[length(extra)+1L]]<-d}}
  traj<-rbindlist(c(list(traj),extra));fwrite(traj,file.path(run,'tables/pseudotime_expression_trajectories.tsv.gz'),sep='\t')
}
for(mode in c('all_neuronal','CRABP1_vs_MEF2C')) {
  seg<-if(mode=='all_neuronal')segments else c('12','11','9','3','4')
  d<-traj[sample=='ALL' & gene %in% representatives & segment %in% seg]
  p<-ggplot(d,aes(pseudotime,mean_log2upx,color=segment,group=segment))+geom_line(linewidth=.6)+
    geom_point(aes(size=percent_expressing),alpha=.8)+
    geom_vline(xintercept=c(.235455231461149,.272363130443725),linetype=2,color='grey45')+
    facet_wrap(~gene,scales='free_y',ncol=4)+scale_color_manual(values=palette,labels=labels)+
    scale_size_area(max_size=2.8,limits=c(0,100))+theme_bw(base_size=10)+
    labs(title=paste('Pseudotime dynamics:',mode),subtitle='Dashed lines: BP11 and BP9; shared ancestors retain their own identities',x='Original URD pseudotime',y='Mean log2(UPX+1)',size='% expressing',caption=caption)
  save_plot(p,paste0('Figure5_trajectories_',mode),15,13)
}
methods<-c('DIV90 neuronal function downstream analysis',paste('Generated',stamp),paste('Run',run_id),
 paste('Existing DE source',previous),paste('Read-only URD object',source),
 'Six original native URD detection-frequency comparisons; no DE rerun or new trajectory fit.',
 'Native detection-test global BH FDR <=0.05 and directional log2FC >=0.25; no AUCPR or detection-difference gate.',
 'All tested genes retained in master; ontology-only genes below prior prefilter marked not tested in coverage table.',
 'Expression: percentage with nonzero normalized expression; linear mean UPX and mean log2(UPX+1), both include zeros.',
 'Figure 1: up to 15 functional DE genes per daughter by global FDR then effect; Figure 3: up to 10 synaptic DE genes per daughter plus requested landmarks. Full tables are exhaustive.',
 'Descriptive trajectories: 12 equal-width bins within each original segment, >=5 cells per reported bin/sample.',
 'Shared ancestors are context and are not assigned future terminal identities. Internal branch 9 is not a sixth terminal lineage.',
 'Targeted activity panels are labelled hypotheses, separate from GO/SynGO curated memberships.',
 'Module scores: mean gene-wise z scores on log2(UPX+1), fixed pooled neuronal reference; constant genes omitted.',
 'Core IEG coexpression is a count distribution, not a threshold declaring activation.',
 'Enrichment: one-sided hypergeometric overrepresentation; 5-2000 background members per set; global BH across all reported enrichment tests.',
 'Both-sibling enrichment background requires testing against both siblings; one-only uses union of eligible tests.',
 'No new biological-replicate inference: cell-level original statistics and module summaries can reflect sample/depth differences.',
 'Ion-channel and synaptic expression indicate molecular competence; coordinated IEGs support activity-responsive transcription but may reflect stress or handling.',
 'Neither establishes firing, synaptic activity, recent activation timing, fast-spiking physiology, or causal lineage commitment. Pseudotime is not elapsed developmental time.')
writeLines(methods,file.path(run,'METHODS_AND_INTERPRETATION.txt'))
counts<-unique(master[candidate==TRUE & !is.na(functional_categories),
  .(branchpoint,reported_daughter,gene,sibling_status,functional_categories)])
out(counts[,.(functional_DE_genes=uniqueN(gene)),by=.(branchpoint,reported_daughter,sibling_status)],'functional_candidate_counts.csv')
report<-c(paste('DIV90 neuronal function summary',stamp),paste('Run',run_id),
  paste('Original direction rows retained:',nrow(master)),paste('Annotated genes in object:',length(selected)),
  paste('Functional gene/branch associations:',nrow(counts)),
  'See functional_candidate_counts.csv for branch counts and core_IEG_cell_coexpression.csv for responding-subset distributions.',
  'A one-sibling result with sibling_tests_available=1 leaves the other comparison untested, rather than demonstrating nonsignificance.',
  'See module_sample_depth_summary.csv before interpreting branch-level maturation/activity differences.',
  'Canonical IEGs, broader activity-associated regulators and maturation programs are separate outputs; no activation label is assigned.',
  'No firing, synaptic activity or activation timing is established by this transcriptional analysis.')
writeLines(report,file.path(run,'RESULTS_SUMMARY.txt'))
capture.output(sessionInfo(),file=file.path(run,'provenance/sessionInfo.txt'))
stopifnot(nrow(master)==nrow(stat),all(master$candidate==(master$log2fc>=.25 & master$global_fdr<=.05)))
out(data.table(check=c('six_original_comparisons','all_direction_rows_retained','candidate_rule_preserved','no_DE_rerun'),passed=TRUE),'validation.csv')
log('ANALYSIS_DONE',nrow(master),'DE rows;',nrow(traj),'trajectory rows;',length(selected),'functional genes')
