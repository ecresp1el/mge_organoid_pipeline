#!/usr/bin/env Rscript
# Individual genes and native segment identities are the primary evidence.
suppressPackageStartupMessages({library(URD);library(Matrix);library(data.table);library(ggplot2)})
setDTthreads(1L)
args<-commandArgs(TRUE);stopifnot(length(args)==4)
run<-args[1];early_run<-args[2];function_run<-args[3];tree<-args[4]
log<-function(...){cat(format(Sys.time(),'%FT%T%z'),...,'\n');flush.console()}
out<-function(x,n)fwrite(x,file.path(run,'tables',n),na='NA')
# Import only the reviewed selection and native detection helpers, never the old main.
for(x in parse(file.path(run,'code/40_div90_urd_branchpoint_genes.R')))
 if(is.call(x)&&identical(x[[1]],as.name('<-'))&&as.character(x[[2]])[1] %in% c('select_window','native_detection'))eval(x)
ann<-fread(file.path(function_run,'config/functional_gene_sets.tsv'))
log('EXPAND_GO_PROGRAMS')
obo<-readLines(file.path(function_run,'config/go-basic.obo'),warn=FALSE)
ids<-which(obo=='[Term]');ends<-c(ids[-1]-1,length(obo));names_go<-parents<-list()
for(i in seq_along(ids)){
 b<-obo[ids[i]:ends[i]];id<-sub('id: ','',b[grepl('^id: GO:',b)])[1]
 if(is.na(id)||any(b=='is_obsolete: true'))next
 names_go[[id]]<-sub('name: ','',b[grepl('^name: ',b)])[1]
 parents[[id]]<-c(sub('is_a: (GO:[0-9]+).*','\\1',b[grepl('^is_a: ',b)]),sub('relationship: part_of (GO:[0-9]+).*','\\1',b[grepl('^relationship: part_of ',b)]))
}
cache<-new.env(parent=emptyenv())
anc<-function(id){if(exists(id,cache,inherits=FALSE))return(get(id,cache));a<-unique(c(id,unlist(lapply(parents[[id]],anc))));assign(id,a,cache);a}
program_patterns<-c(calcium_signaling='calcium-mediated signaling|calmodulin|calcineurin|calcium ion transport|cellular response to calcium ion$',
 activity_transcription='CREB|cAMP-mediated signaling|response to electrical stimulus$|response to membrane depolarization',
 functional_maturation='synapse maturation|axon initial segment|action potential|regulation of membrane potential')
sets<-lapply(program_patterns,function(pattern)names(names_go)[grepl(pattern,unlist(names_go),ignore.case=TRUE)])
term_program<-lapply(names(names_go),function(id)names(sets)[vapply(sets,function(ss)any(anc(id)%in%ss),logical(1))]);names(term_program)<-names(names_go)
con<-gzfile(file.path(function_run,'config/goa_human.gaf.gz'),'rt');extra<-list();nr<-0L
repeat{
 lines<-readLines(con,n=50000,warn=FALSE);if(!length(lines))break
 for(line in lines[!startsWith(lines,'!')]){a<-strsplit(line,'\t',fixed=TRUE)[[1]]
  if(length(a)<15||'NOT'%in%strsplit(a[4],'|',fixed=TRUE)[[1]]||!grepl('taxon:9606',a[13],fixed=TRUE))next
  pp<-term_program[[a[5]]];if(!length(pp))next
  nr<-nr+1L;extra[[nr]]<-list(gene=a[3],set_id=a[5],set_name=names_go[[a[5]]],categories=paste(pp,collapse='|'),source='GO_human_expanded',evidence='Frozen GO annotations; calcium/activity-associated context, not proof of activation')
 };log('GO_CHUNK',length(extra),'relevant annotations')
};close(con)
if(length(extra))ann<-rbindlist(list(ann,rbindlist(extra)),use.names=TRUE)
panel<-strsplit('CAMK2A CAMK2B CAMK2D CAMK2G CAMK4 CALM1 CALM2 CALM3 PPP3CA PPP3CB PPP3CC PPP3R1 PPP3R2 CRTC1 CRTC2 CRTC3 CREBBP EP300 NFATC1 NFATC2 NFATC3 NFATC4 FOSL1 FOSL2 JUND ATF3 EGR4 BDNF NTRK2 NPTX1 NPTX2 HOMER1 IER2 IER3 IER5 DUSP1 DUSP4 DUSP5 DUSP6 RGS2 RGS4 GADD45A GADD45B GADD45G PER1 SIK1 RHEB RASGRF1 RASGRF2 MAPK1 MAPK3 MAP2K1 MAP2K2 RBFOX1 SATB1',' ')[[1]]
ann<-unique(rbindlist(list(ann,data.table(gene=panel,set_id='PANEL:expanded_signaling',set_name='Expanded signaling and transcription landmarks',categories='activity_signaling_landmarks',source='targeted_hypotheses_v2',evidence='Context-dependent signaling/development genes; not an activation signature')),use.names=TRUE))
out(ann,'expanded_gene_program_annotations.csv')
log('READ_PROTECTED_TREE');obj<-readRDS(tree);expr<-obj@logupx.data;genes<-rownames(expr);cells<-colnames(expr)
pt<-setNames(obj@pseudotime[cells,'div90_jia_rootscore_root'],cells)
stopifnot(max(abs(obj@tree$pseudotime[cells]-pt))<1e-10)
members<-lapply(obj@tree$cells.in.segment,function(cc)intersect(cc,cells));samples<-setNames(as.character(obj@meta[cells,'orig.ident']),cells)
segments<-c('12','11','3','5','9','1','2','4');daughters<-c('3','5','9','1','2','4');gs<-intersect(unique(ann$gene),genes)
out(data.table(gene=unique(ann$gene),present_in_object=unique(ann$gene)%in%genes),'gene_coverage.csv')
state_cols<-grep('cluster_number_name|cell.?type|annotation|cluster|label',colnames(obj@meta),ignore.case=TRUE,value=TRUE)
writeLines(c('Existing state metadata columns:',state_cols),file.path(run,'provenance/state_metadata_columns.txt'))
metric<-function(cc,g){x<-expr[g,cc,drop=FALSE];y<-x;y@x<-2^y@x-1
 data.table(gene=g,n_cells=length(cc),n_expressing=as.integer(rowSums(x>0)),percent_expressing=100*rowMeans(x>0),mean_normalized=rowMeans(y),mean_log2upx=rowMeans(x))}
summary<-traj<-list()
for(s in segments){cc<-members[[s]];if(length(cc)<5)next;log('GENE_PROFILES',s,length(cc),'cells',length(gs),'genes')
 bins<-cut(pt[cc],seq(min(pt[cc]),max(pt[cc]),length.out=13),include.lowest=TRUE,labels=FALSE)
 for(sample in c('ALL',unique(samples[cc]))){use<-if(sample=='ALL')cc else cc[samples[cc]==sample];if(length(use)<2)next
  z<-metric(use,gs);z[,`:=`(segment=s,sample=sample)];summary[[length(summary)+1L]]<-z
  for(k in sort(unique(bins))){bc<-intersect(use,cc[bins==k]);if(length(bc)<5)next
   z<-metric(bc,gs);z[,`:=`(segment=s,sample=sample,bin=k,pseudotime=mean(pt[bc]),time_low=min(pt[bc]),time_high=max(pt[bc]))];traj[[length(traj)+1L]]<-z}
 }
}
summary<-rbindlist(summary);traj<-rbindlist(traj);setorder(traj,segment,sample,gene,pseudotime)
out(summary,'individual_gene_branch_sample_expression.csv');fwrite(traj,file.path(run,'tables/individual_gene_pseudotime_profiles.tsv.gz'),sep='\t')
trend<-traj[,.(bins=.N,delta_mean_log2upx=last(mean_log2upx)-first(mean_log2upx),delta_detection_pp=last(percent_expressing)-first(percent_expressing),
 rho=if(.N>2&&sd(mean_log2upx)>0)cor(pseudotime,mean_log2upx,method='spearman')else NA_real_),by=.(segment,sample,gene)]
out(trend,'individual_gene_sample_trends.csv')
consistency<-trend[sample!='ALL',.(samples_with_support=.N,samples_increasing=sum(delta_mean_log2upx>0),samples_decreasing=sum(delta_mean_log2upx<0),median_delta=median(delta_mean_log2upx),min_delta=min(delta_mean_log2upx),max_delta=max(delta_mean_log2upx)),by=.(segment,gene)]
out(consistency,'individual_gene_sample_consistency.csv')
log('FOS_LOCALIZATION')
foscells<-rbindlist(lapply(c('3','4'),function(s){cc<-members[[s]];d<-data.table(cell=cc,segment=s,sample=samples[cc],pseudotime=pt[cc],FOS_log2upx=as.numeric(expr['FOS',cc]),library_counts=as.numeric(colSums(obj@count.data[,cc,drop=FALSE])))
 for(col in state_cols)d[,(col):=as.character(obj@meta[cc,col])];d}),fill=TRUE)
edges<-seq(min(foscells$pseudotime),max(foscells$pseudotime),length.out=17);foscells[,bin:=cut(pseudotime,edges,include.lowest=TRUE,labels=FALSE)]
out(foscells,'FOS_cell_sample_state_localization.csv')
fos<-foscells[,.(n_cells=.N,n_FOS_positive=sum(FOS_log2upx>0),percent_FOS_positive=100*mean(FOS_log2upx>0),mean_log2upx=mean(FOS_log2upx),mean_library_counts=mean(library_counts)),by=.(segment,sample,bin)]
fos[,`:=`(time_low=edges[bin],time_high=edges[bin+1])];out(fos,'FOS_sample_pseudotime_localization.csv')
fos_states<-rbindlist(lapply(state_cols,function(col)foscells[,.(n_cells=.N,n_FOS_positive=sum(FOS_log2upx>0),percent_FOS_positive=100*mean(FOS_log2upx>0),mean_log2upx=mean(FOS_log2upx)),by=.(segment,sample,state=get(col))][,state_column:=col]),fill=TRUE)
out(fos_states,'FOS_existing_cell_state_localization.csv')
# Predeclared adjacent later window: [split+0.05R,split+0.10R], unchanged support rules.
cfg<-list(minimum_cells=30,minimum_bin_cells=5,minimum_joint_bins=3,matching_bins=5,seed=7)
w<-.05*diff(range(pt));support<-late<-list();selected<-list()
for(bp in c('11','9')){split<-if(bp=='11').235455231461149 else .272363130443725;kids<-if(bp=='11')c('5','9','3')else c('1','2','4')
 for(pair in combn(kids,2,simplify=FALSE)){sel<-select_window(pt,members[[pair[1]]],members[[pair[2]]],split+w,split+2*w,cfg)
  id<-paste(bp,paste(pair,collapse='_'),sep='_');selected[[id]]<-sel
  support[[id]]<-data.table(id=id,branchpoint=paste0('BP',bp),daughter_A=pair[1],daughter_B=pair[2],time_low=split+w,time_high=split+2*w,n_A=length(sel$a),n_B=length(sel$b),joint_bins=sel$joint_bins,status=sel$status)
 }}
support<-rbindlist(support);out(support,'later_comparison_support.csv')
for(i in which(support$status=='supported')){r<-support[i];a<-selected[[r$id]]$a;b<-selected[[r$id]]$b;log('LATER_NATIVE_COMPARISON',r$id,length(a),length(b))
 fa<-rowMeans(expr[,a,drop=FALSE]>0);fb<-rowMeans(expr[,b,drop=FALSE]>0);gg<-genes[pmax(fa,fb)>=.10]
 z<-native_detection(obj,a,b,gg,fa[gg],fb[gg],'div90_jia_rootscore_root');ma<-metric(a,gg);mb<-metric(b,gg)
 z[,`:=`(log2fc=log2((ma$mean_normalized+1e-6)/(mb$mean_normalized+1e-6)),detection_A=fa[gg],detection_B=fb[gg],mean_normalized_A=ma$mean_normalized,mean_normalized_B=mb$mean_normalized,
  id=r$id,branchpoint=r$branchpoint,daughter_A=r$daughter_A,daughter_B=r$daughter_B)]
 z[,p_direction_corrected:=fifelse(is.na(native_p_enriched),1,pmin(1,2*native_p_enriched))];late[[r$id]]<-z
}
late<-rbindlist(late);late[,later_global_fdr:=p.adjust(p_direction_corrected,'BH')];late[,candidate:=abs(log2fc)>=.25 & later_global_fdr<=.05]
out(late,'later_native_URD_all_tested_genes.csv')
out(merge(late,unique(ann[,.(gene,program=categories)]),by='gene',all.x=TRUE,allow.cartesian=TRUE),'later_gene_program_evidence.csv')
early<-fread(file.path(early_run,'DIV90_URD_gene_results.csv'))
early[,`:=`(daughter_A=as.character(daughter_A),daughter_B=as.character(daughter_B))]
early_one<-early[,.(early_effect=first(log2fc_A_vs_B),early_fdr=first(global_fdr),early_significant=any(candidate)),by=.(branchpoint,daughter_A,daughter_B,gene)]
timing<-merge(late,early_one,by=c('branchpoint','daughter_A','daughter_B','gene'),all=TRUE)
timing[,timing_evidence:=fcase(!is.na(candidate)&candidate & !is.na(early_significant)&early_significant,'early_and_later',
 !is.na(candidate)&candidate & !is.na(early_significant)&!early_significant,'later_detected_early_not_significant',
 !is.na(candidate)&candidate & is.na(early_significant),'later_detected_early_untested',
 !is.na(early_significant)&early_significant & is.na(candidate),'early_supported_later_unsupported_or_untested',
 !is.na(early_significant)&early_significant,'early_detected_later_not_significant',default='no_significant_difference_in_available_tests')]
out(timing,'near_split_vs_later_evidence.csv')
annotation<-ann[,.(programs=paste(sort(unique(categories)),collapse='|'),sources=paste(sort(unique(source)),collapse='|')),by=gene]
evidence<-merge(early,annotation,by='gene',all.x=TRUE,sort=FALSE);out(evidence,'consolidated_gene_evidence_all_early_tests.csv')
palette<-c('3'='#D55E00','5'='#CC79A7','9'='#009E73','1'='#0072B2','2'='#56B4E9','4'='#E69F00','11'='#666666','12'='#999999')
saveplot<-function(p,n,w=15,h=11){log('PLOT',n);for(ext in c('pdf','png','svg'))ggsave(file.path(run,'plots',paste0(n,'.',ext)),p,width=w,height=h,dpi=300,device=if(ext=='pdf')grDevices::cairo_pdf else if(ext=='svg')svglite::svglite else 'png',limitsize=FALSE)}
representatives<-intersect(c('FOS','FOSB','EGR1','EGR2','NPAS4','NR4A1','ARC','CREB1','CREM','MEF2C','CAMK4','PPP3CA','CRTC1','BDNF','NPTX2','HOMER1','GAD1','SLC32A1','SCN1A','KCNC1','KCNC2','CACNA1A','ANK3','SNAP25','SYT1','GRIN1','GRIA2','DLG4'),gs)
for(mode in c('all_branches','CRABP1_vs_MEF2C')){
 ss<-if(mode=='all_branches')segments else c('12','11','9','3','4')
 d<-traj[gene%in%representatives & segment%in%ss & sample!='ALL']
 p<-ggplot(d,aes(pseudotime,mean_log2upx,color=segment,group=interaction(segment,sample)))+geom_line(alpha=.45,linewidth=.4)+
  geom_vline(xintercept=c(.235455231461149,.272363130443725),linetype=2,color='grey50')+facet_wrap(~gene,scales='free_y',ncol=4)+scale_color_manual(values=palette)+theme_bw(base_size=9)+labs(title=paste('Individual gene trajectories:',mode),subtitle='Each line is one existing segment/sample; no requirement for IEG coexpression',x='Original URD pseudotime',y='Mean log2(UPX+1)',caption=basename(run))
 saveplot(p,paste0('individual_gene_trajectories_',mode),15,17)
}
p<-ggplot(fos,aes((time_low+time_high)/2,percent_FOS_positive,color=segment,group=segment))+geom_line()+geom_point(aes(size=n_cells))+facet_wrap(~sample)+scale_color_manual(values=palette)+theme_bw()+labs(title='FOS localization by sample and common pseudotime bins',x='Original URD pseudotime',y='% FOS expressing',caption='Small bins retained and labelled by cell count; descriptive only')
saveplot(p,'FOS_sample_pseudotime_localization')
capture.output(sessionInfo(),file=file.path(run,'provenance/sessionInfo.txt'))
writeLines(c('Individual genes, expression detection and sample consistency are primary evidence; module scores and IEG coexpression are not activation criteria.',
 'Three interpretation axes: electrophysiological machinery, activity-associated transcription, lineage specification.',
 'Later window is the adjacent 5% full-pseudotime-range interval; same 30-cell/5-bin/3-supported-bin rules as the original analysis.',
 'Only two of six later comparisons pass support; unsupported comparisons receive no test. Later BH family includes every eligible gene in both supported comparisons.',
 'Early FDR is unchanged. Later detection does not prove biological onset; early/later differences are not an interaction test.',
 'Per-sample trends are descriptive, with >=5 cells/bin. Sample support, depth, and state composition can affect apparent differences.',
 'FOS cell states use existing metadata; no reclustering or reassignment. Shared ancestors remain shared context.',
 'GO calcium/signaling terms and targeted landmarks include context-dependent genes; they are not specific evidence of neuronal activation.',
 'Pseudotime is not elapsed time. Expression does not establish firing, synaptic function, activation timing or causal lineage specification.'),file.path(run,'METHODS_AND_LIMITATIONS.txt'))
log('COMPLETED',length(gs),'profiled genes;',nrow(late),'later tests;',nrow(traj),'trajectory rows')
