from pathlib import Path
from pcdh19_exploratory_paths import RUNTIME_ROOT
import csv,gzip,json,math,copy
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.cluster.hierarchy import linkage,dendrogram
ROOT=Path(__file__).resolve().parents[1];OUT=RUNTIME_ROOT/'results/female_male_signaling_20261031';FEMALE=RUNTIME_ROOT/'results/neuronal_matched_selection_20261031';DEEP=RUNTIME_ROOT/'results/mosaic_deep_dive_20261031';SIG=RUNTIME_ROOT/'results/paired_HET_signaling_20261031'
def read(p):
 with (gzip.open(str(p),'rt') if str(p).endswith('.gz') else open(str(p))) as f:return list(csv.DictReader(f,delimiter='\t'))
def write(p,rows):
 with open(str(p),'w') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':7.5,'axes.titlesize':9,'axes.labelsize':7.5,'xtick.labelsize':7,'ytick.labelsize':7,'pdf.fonttype':42,'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':.65,'xtick.major.width':.65,'ytick.major.width':.65})
BLUE='#326C9D';RED='#B64B50';INK='#2F3238';GRAY='#C4C7CC'
fig=plt.figure(figsize=(13.8,7.6));gs=fig.add_gridspec(2,8,left=.055,right=.975,top=.865,bottom=.17,wspace=.95,hspace=.80,height_ratios=[1,1])
fig.text(.055,.963,'Pcdh19 mosaic states and neuronal expression',fontsize=15,weight='bold',color=INK)
fig.text(.055,.915,'FEMALES   WT reference and HET mosaic states',fontsize=9,weight='bold',color=INK)
by=gs[1,0].get_position(fig).y1;fig.text(.055,by+.050,'MALES   Known WT and KO controls',fontsize=9,weight='bold',color=INK)
fig.text(.78,.955,'Lower in KO-like / KO',fontsize=7.5,color=BLUE,ha='right');fig.text(.975,.955,'Higher in KO-like / KO',fontsize=7.5,color=RED,ha='right')
def title(slot,letter,text):
 b=slot.get_position(fig);fig.text(b.x0-.022,b.y1+.018,letter,fontsize=11,weight='bold',color=INK);fig.text(b.x0,b.y1+.018,text,fontsize=9,color=INK)
def polish(ax):
 ax.tick_params(length=2.5,pad=2);ax.spines['left'].set_color('#555');ax.spines['bottom'].set_color('#555')
# A: supported broad annotations; unresolved clusters remain explicitly unresolved.
broad=[('RG-like progenitors',['C0002'],'#7186A6'),('MGE IPCs',['C0003','C0008'],'#CE9667'),('Neurogenic progenitors',['C0001','C0007','C0012','C0024'],'#B8A363'),('Cycling neuroblasts',['C0011'],'#9B789D'),('Immature inhibitory neurons',['C0004','C0005','C0006','C0020'],'#4D928A'),('Mixed neural / vascular',['C0015'],'#B77D86'),('Unresolved populations',[],GRAY)]
lookup={c:i for i,(_,cc,_) in enumerate(broad) for c in cc};cache=[r for r in read(RUNTIME_ROOT/'results/mosaic_nature_figure_20261031/displayed_umap_cells.tsv.gz') if r['panel']=='background'];xy=np.array([[float(r['UMAP1']),float(r['UMAP2'])] for r in cache]);labels=np.array([lookup.get(r['consensus_cluster'],6) for r in cache]);sub=gs[0,:2].subgridspec(2,1,height_ratios=[1,.36],hspace=.02);ax=fig.add_subplot(sub[0]);
for i in [6,5,0,3,2,1,4]:
 m=labels==i;ax.scatter(xy[m,0],xy[m,1],s=.65,color=broad[i][2],linewidths=0,rasterized=True)
ax.set_aspect('equal');ax.axis('off');la=fig.add_subplot(sub[1]);la.axis('off');handles=[Line2D([0],[0],marker='o',ls='',color=color,markersize=4,label=name) for name,cc,color in broad];la.legend(handles=handles,loc='upper center',ncol=2,frameon=False,fontsize=6.3,columnspacing=.8,handletextpad=.4,handlelength=.6,labelspacing=.45,borderaxespad=0)
title(gs[0,:2],'a','Female MGE atlas')
write(OUT/'UMAP_broad_annotation_mapping.tsv',[dict(cluster=c,broad_annotation=name,color=color) for name,cc,color in broad for c in cc]+[dict(cluster=r['consensus_cluster'],broad_annotation=broad[6][0],color=GRAY) for r in cache if r['consensus_cluster'] not in lookup][:0])
# B: same defining exon-1 readout; no unfiltered WT reference.
e=np.load(str(FEMALE/'neuronal_exon1_signal.npz'));dist=[e['exon1_CP10K'][e['group']==k] for k in range(3)];ax=fig.add_subplot(gs[0,2]);nz=[k for k,d in enumerate(dist) if np.ptp(d)>0];v=ax.violinplot([dist[k] for k in nz],positions=nz,widths=.75,showextrema=False,points=150)
for b in v['bodies']:b.set_facecolor('#B3B5B8');b.set_edgecolor('#74767A');b.set_alpha(1);b.set_linewidth(.65)
for k in range(3):
 if k not in nz:ax.plot([k-.28,k+.28],[0,0],color=INK,lw=1.7,clip_on=False)
upper=max(np.percentile(d,95) for d in dist)*1.08;ax.set_ylim(0,upper);ax.plot(range(3),[d.mean() for d in dist],'-D',color=INK,ms=2.5,lw=.65);ax.set_xticks(range(3));ax.set_xticklabels(['WT female','HET WT-like','HET KO-like'],rotation=50,ha='right',fontsize=6.5);ax.set_ylabel('Exon-1 UMI / 10,000',fontsize=7);polish(ax);title(gs[0,2],'b','Pcdh19 exon 1')
ax.text(.5,1.035,'≤5% above range',transform=ax.transAxes,fontsize=5.5,ha='center',color='#777')
# C: one aligned tree / expression / effect panel replaces two repetitive displays.
md=read(FEMALE/'nine_sample_metadata.tsv');grouporder=[('WT_F','A_positive'),('HET_F','A_positive'),('HET_F','Aneg_BC1')];md=sorted(md,key=lambda r:(grouporder.index((r['design_group'],r['gate'])),r['sample']));libs=[r['library'] for r in md];vals={r['gene_symbol']:np.array([float(r[l]) for l in libs]) for r in read(FEMALE/'nine_sample_log2CPM.tsv')};fd={r['gene_symbol']:r for r in read(DEEP/'DE_C0006/within_HET_Aneg_BC1_vs_Apositive.tsv.gz')}
genes=['Ank3','Kcnc2','Gnao1','Ptprd','Epha5','Snap47','Pcdha5','Pcdha10','Pcdhga5','Lhx8','Sptbn2','Mctp1','Dlx1','Robo2','Erbb4'];means=np.array([[vals[g][k:k+3].mean() for k in [0,3,6]] for g in genes]);z=(means-means.mean(1)[:,None])/means.std(1)[:,None]
sg=gs[0,3:6].subgridspec(1,4,width_ratios=[.35,.55,.9,.7],wspace=.10);ta=fig.add_subplot(sg[0]);ro=dendrogram(linkage(z,method='average'),orientation='left',no_labels=True,color_threshold=0,above_threshold_color='#555',ax=ta)['leaves'];ta.invert_yaxis();ta.axis('off');na=fig.add_subplot(sg[1]);na.set_ylim(len(genes)-.5,-.5);na.set_xlim(0,1);na.axis('off')
for i,k in enumerate(ro):na.text(.98,i,genes[k],ha='right',va='center',fontsize=6.6)
ha=fig.add_subplot(sg[2]);ha.imshow(z[ro],aspect='auto',cmap='Purples',vmin=-1.5,vmax=1.5);ha.set_yticks([]);ha.set_xticks(range(3));ha.set_xticklabels(['WT female','HET WT-like','HET KO-like'],rotation=50,ha='right',fontsize=6.2);ha.tick_params(length=2,pad=2);ha.set_title('Low → High',fontsize=7,pad=5)
ba=fig.add_subplot(sg[3]);fc=np.array([float(fd[genes[k]]['logFC']) for k in ro]);ba.barh(range(len(ro)),fc,height=.7,color=[BLUE if q<0 else RED for q in fc]);ba.set_ylim(len(ro)-.5,-.5);ba.axvline(0,color='#aaa',lw=.6);ba.set_yticks([]);ba.set_xticks([-.5,0,1]);ba.tick_params(axis='x',labelsize=6,length=2);ba.set_title('KO-like / WT-like',fontsize=7,pad=5);ba.set_xlabel('log2FC',fontsize=6.5,labelpad=3);ba.spines['left'].set_visible(False)
title(gs[0,3:6],'c','Neuronal expression and HET differences')
# D: exact ligand/receptor components, not an inferred communication cartoon.
cts=['C0002','C0001','C0003','C0011','C0006'];names=['RG-like progenitors','Neurogenic progenitors','MGE IPCs','Cycling neuroblasts','Immature neurons'];de={(r['stratum'],r['gene_symbol']):r for r in read(SIG/'paired_DE_all_celltypes_with_signaling_annotations.tsv')};arr=np.array([[float(de[(ct,g)]['logFC']) if (ct,g) in de else np.nan for g in ['Nrg1','Erbb4']] for ct in cts]);ds=gs[0,6:].subgridspec(1,3,width_ratios=[1.5,1,.065],wspace=.08);nm=fig.add_subplot(ds[0]);nm.set_ylim(4.5,-.5);nm.set_xlim(0,1);nm.axis('off')
for i,n in enumerate(names):nm.text(.98,i,n,ha='right',va='center',fontsize=6.7)
ax=fig.add_subplot(ds[1]);im=ax.imshow(arr,aspect='auto',cmap='RdBu_r',vmin=-1.2,vmax=1.2);ax.set_yticks([]);ax.set_xticks([0,1]);ax.set_xticklabels(['Nrg1\nligand','Erbb4\nreceptor'],fontsize=7);ax.tick_params(length=0,pad=4)
for i in range(5):
 for j in range(2):
  q=arr[i,j];ax.text(j,i,'{:+.2f}'.format(q),ha='center',va='center',fontsize=7,color='white' if abs(q)>.7 else INK)
ca=fig.add_subplot(ds[2]);cb=fig.colorbar(im,cax=ca,ticks=[-1,0,1]);cb.ax.tick_params(labelsize=6,length=2);cb.set_label('KO-like / WT-like log2FC',fontsize=6.5,labelpad=4);title(gs[0,6:],'d','Candidate Nrg1–Erbb4 axis')
write(OUT/'signaling_component_panel_values.tsv',[dict(cell_type=names[i],gene=g,log2FC=float(arr[i,j]),PValue=de[(ct,g)]['PValue'],FDR=de[(ct,g)]['FDR']) for i,ct in enumerate(cts) for j,g in enumerate(['Nrg1','Erbb4'])])
# E: independent genotype control, with the same neuronal compartment.
ma=read(OUT/'male_exon1_detection.tsv');ax=fig.add_subplot(gs[1,:2]);mv=[]
for k,g in enumerate(['WT males','KO males']):
 vv=np.array([100*float(r['fraction']) for r in ma if r['group']==g]);mv.append(float(vv.mean()));ax.bar(k,vv.mean(),width=.5,color='#A5A8AE' if k==0 else BLUE);ax.errorbar(k,vv.mean(),yerr=vv.std(ddof=1)/np.sqrt(3),fmt='none',ecolor=INK,lw=.85,capsize=3);ax.text(k,vv.mean()+.75,'{:.2f}%'.format(vv.mean()),ha='center',fontsize=8,color=INK)
ax.set_ylim(0,7);ax.set_xlim(-.65,1.65);ax.set_xticks([0,1]);ax.set_xticklabels(['WT males','KO males']);ax.set_ylabel('Exon-1-positive neurons (%)');polish(ax);title(gs[1,:2],'e','Pcdh19 exon-1 detection')
# F: same male neuronal candidates and estimates, across biological populations.
old=json.loads((OUT/'figure_validation.json').read_text());selected=old['male_selected_genes'];male={ct:{r['gene_symbol']:r for r in read(OUT/('male_'+ct+'_DE.tsv'))} for ct in cts};mat=np.array([[float(male[ct][g]['logFC']) if g in male[ct] else np.nan for ct in cts] for g in selected]);sg=gs[1,2:6].subgridspec(1,4,width_ratios=[.4,.6,1.7,.07],wspace=.06);ta=fig.add_subplot(sg[0]);ro=dendrogram(linkage(np.nan_to_num(mat,nan=0),method='average'),orientation='left',no_labels=True,color_threshold=0,above_threshold_color='#555',ax=ta)['leaves'];ta.invert_yaxis();ta.axis('off');na=fig.add_subplot(sg[1]);na.set_ylim(len(selected)-.5,-.5);na.set_xlim(0,1);na.axis('off')
for i,k in enumerate(ro):na.text(.98,i,selected[k],ha='right',va='center',fontsize=6.5)
ax=fig.add_subplot(sg[2]);cm=copy.copy(plt.get_cmap('RdBu_r'));cm.set_bad('#ccc');im=ax.imshow(np.ma.masked_invalid(mat[ro]),aspect='auto',cmap=cm,vmin=-1.5,vmax=1.5);ax.set_yticks([]);ax.set_xticks(range(5));ax.set_xticklabels(['RG-like\nprogenitors','Neurogenic\nprogenitors','MGE IPCs','Cycling\nneuroblasts','Immature\nneurons'],rotation=35,ha='right',fontsize=6.5);ax.tick_params(length=2,pad=3);ca=fig.add_subplot(sg[3]);cb=fig.colorbar(im,cax=ca,ticks=[-1,0,1]);cb.ax.tick_params(labelsize=6,length=2);cb.set_label('KO / WT log2FC',fontsize=6.5,labelpad=3);title(gs[1,2:6],'f','Exploratory neuronal changes across cell types')
# G: compact, intuitive percent changes in immature neurons.
m6=male['C0006'];leaders=[g for g in selected if abs(float(m6[g]['logFC']))>=.5][:9];leaders=sorted(leaders,key=lambda g:float(m6[g]['logFC']));pc=np.array([100*(2**float(m6[g]['logFC'])-1) for g in leaders]);ax=fig.add_subplot(gs[1,6:]);ax.barh(range(len(leaders)),pc,height=.65,color=[BLUE if q<0 else RED for q in pc]);ax.set_yticks(range(len(leaders)));ax.set_yticklabels(leaders,fontsize=7);ax.invert_yaxis();ax.axvline(0,color='#aaa',lw=.6);ax.set_xlabel('KO change relative to WT (%)',fontsize=7);polish(ax);title(gs[1,6:],'g','Immature-neuron candidates')
for ext in ['png','pdf','svg']:fig.savefig(str(OUT/('cohesive_female_male_figure.'+ext)),dpi=350,facecolor='white')
# Complete annotation provenance retains every cluster and its broad category.
write(OUT/'UMAP_broad_annotation_mapping.tsv',[dict(cluster=c,broad_annotation=broad[lookup.get(c,6)][0],color=broad[lookup.get(c,6)][2]) for c in sorted(set(r['consensus_cluster'] for r in cache))])
(OUT/'cohesive_figure_validation.json').write_text(json.dumps({'UMAP_cells':len(cache),'UMAP_named_broad_groups':[r[0] for r in broad],'broad_annotation_source':'composition_review_20261031/full_consensus_working_names.tsv; unresolved groups not assigned invented identities','female_gene_list':'restored mixed-direction list; 15 genes','female_exon1_values_unchanged':True,'signaling_display':'actual paired ligand and receptor gene log2FC; no arrows or communication claim','male_expression_statistics_unchanged':True,'no_gene_wide_FDR_significant_male_DEGs':True,'candidate_signaling_not_FDR_significant':True,'sample_ids_displayed':False,'cluster_codes_displayed':False},indent=2))
(OUT/'COHESIVE_FIGURE_LEGEND.txt').write_text('Top: female MGE atlas (A), matched exon-1 selection readout (B), provisional immature inhibitory-neuron expression and paired HET effects (C), candidate ligand/receptor gene expression by neural population (D). A uses seven broad categories consolidated from saved marker-supported working annotations; unresolved populations remain gray. C combines gene-wise standardized equal-sample mean expression with the corresponding paired KO-like/WT-like edgeR effect estimates. D gives fitted paired gene log2FC for Nrg1 and Erbb4, not measured ligand-receptor interaction or pathway activity. Bottom: known WT/KO male controls without probe selection. E gives exon-1 detection in all cells of the same neuronal cluster; bars are equal-sample means +/- SEM, three samples per group. F/G retain the previous exploratory male neuronal candidate selection and model estimates. Gene differences and signaling candidates do not pass gene-wide FDR correction; neither WT/HET similarity nor functional communication is established. B zero signal in KO-like cells is imposed by the selection rule; WT and HET WT-like both use A-positive cells. Violin upper limits and above-range fractions are shown. Full test universes and previous network figures are retained separately. Color and type changes do not alter expression measurements or statistics.\n')
print('COHESIVE_FIGURE_COMPLETE',flush=True)
