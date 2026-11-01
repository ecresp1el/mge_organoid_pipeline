"""Saved-embedding, three-group exploratory mosaic figure with sample-level expression."""
import csv
import gzip
import hashlib
import json
from pathlib import Path
from pcdh19_exploratory_paths import RUNTIME_ROOT
import numpy as np
import h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib import patheffects

ROOT=Path(__file__).resolve().parents[1]
OUT=RUNTIME_ROOT/'results/mosaic_nature_figure_20261031'
SOURCE=RUNTIME_ROOT/'results/mosaic_deep_dive_20261031'
BASE=Path('/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/paper3_pcdh19/results')
H5=BASE/'primary_processing/06_technical_sample_batch_diagnostics/06_technical_sample_batch_diagnostics_20260903_114110_6087eeb/objects/pcdh19_step06_unintegrated_diagnostics.h5ad'
ASSIGN=BASE/'hicat/03_full_data_consensus_benchmark/hicat_consensus_review_20260909_142348_efd9e13d/outputs/tables/all_cell_assignments.tsv.gz'
CALLS=BASE/'step_07_pcdh19_het_female_wt_ko_like_classification/step_07_het_female_cell_classifications.tsv'
def read(path):
    with (gzip.open(path,'rt') if path.suffix=='.gz' else path.open()) as f:return list(csv.DictReader(f,delimiter='\t'))
OUT.mkdir(exist_ok=True,parents=True)
female_samples=['15662-JZ-'+str(i) for i in range(4,10)]
cache=OUT/'displayed_umap_cells.tsv.gz'
assert cache.exists(),'Reuse the saved, validated UMAP display cells from the initial rendering.'
cached=read(cache)
ids=np.array([r['cell_id'] for r in cached]);samples=np.array([r['sample'] for r in cached])
clusters=np.array([r['consensus_cluster'] for r in cached]);xy=np.array([[float(r['UMAP1']),float(r['UMAP2'])] for r in cached])
tags=np.array([r['panel'] for r in cached])
background=np.flatnonzero(tags=='background');wt_display=np.flatnonzero(tags=='WT_like_HET');ko_display=np.flatnonzero(tags=='KO_enriched_HET')
assert len(background)==48000 and len(wt_display)==len(ko_display)==1350
previous=json.loads((OUT/'provenance.json').read_text())
wt_population=previous['WT_like_cells'];ko_population=previous['KO_enriched_cells']
assert wt_population==3949 and ko_population==16172
# Read coordinates and categorical sample codes only; no count matrix or large string index is reloaded.
with h5py.File(H5,'r') as f:
    full_xy=f['obsm/X_umap'][:]
    node=f['obs/technical_sample_id']
    if isinstance(node,h5py.Group):
        cats=node['categories'].asstr()[:];full_samples=cats[node['codes'][:]]
    else:full_samples=node.asstr()[:]
female=np.isin(full_samples,female_samples)
limits=[(full_xy[female,k].min(),full_xy[female,k].max()) for k in [0,1]]

plt.rcParams.update({'font.family':'DejaVu Sans','font.size':7.5,'axes.linewidth':.65,
    'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none','xtick.major.width':.65,'ytick.major.width':.65})
green,blue,purple='#549E89','#337EAC','#A678AC'
classcolors={'C0001':'#C5A35D','C0003':'#D27A50','C0006':'#438D82'}
fig=plt.figure(figsize=(8.4,6.2))
gs=fig.add_gridspec(2,1,left=.065,right=.985,top=.92,bottom=.20,hspace=.80)
top=gs[0].subgridspec(1,3,wspace=.30);bottom=gs[1].subgridspec(1,4,wspace=.62)
axes=[fig.add_subplot(top[0,j]) for j in range(3)]+[fig.add_subplot(bottom[0,j]) for j in range(4)]
for i,ax in enumerate(axes[:3]):
    ax.scatter(xy[background,0],xy[background,1],s=.5,c='#d8dadd',linewidths=0,rasterized=True,zorder=1)
    if i==0:
        for cl,col in classcolors.items():
            ix=background[clusters[background]==cl]
            ax.scatter(xy[ix,0],xy[ix,1],s=.8,c=col,linewidths=0,rasterized=True,zorder=2)
        ax.set_title('Female MGE cell landscape',fontsize=8,pad=10)
    else:
        ix,col=(wt_display,blue) if i==1 else (ko_display,purple)
        ax.scatter(xy[ix,0],xy[ix,1],s=2,c=col,linewidths=0,alpha=.8,rasterized=True,zorder=3)
        ax.set_title('HET: WT-like (A+)' if i==1 else 'HET: KO-enriched (A−, B/C+)',fontsize=8,pad=10)
        ax.text(.02,.02,'n = {:,} cells'.format(wt_population if i==1 else ko_population),transform=ax.transAxes,fontsize=7)
    for k in [0,1]:
        lo,hi=limits[k];pad=.05*(hi-lo)
        (ax.set_xlim if k==0 else ax.set_ylim)(lo-pad,hi+pad)
    ax.set_aspect('equal',adjustable='box');ax.set_xticks([]);ax.set_yticks([])
    for sp in ax.spines.values():sp.set_visible(False)
    if i==0:
        # Small coordinate arrows convey the shared embedding without arbitrary tick values.
        ax.annotate('',xy=(.22,.04),xytext=(.04,.04),xycoords='axes fraction',arrowprops=dict(arrowstyle='->',lw=.65))
        ax.annotate('',xy=(.04,.22),xytext=(.04,.04),xycoords='axes fraction',arrowprops=dict(arrowstyle='->',lw=.65))
        ax.text(.13,-.04,'UMAP 1',transform=ax.transAxes,ha='center',fontsize=6)
        ax.text(-.05,.13,'UMAP 2',transform=ax.transAxes,rotation=90,va='center',fontsize=6)
    ax.text(-.12,1.09,'abc'[i],transform=ax.transAxes,fontweight='bold',fontsize=11)

# Context legend is outside the embedding; no new biological taxonomy is inferred.
legend=[Line2D([0],[0],marker='o',linestyle='',markerfacecolor=col,markeredgecolor='none',markersize=4,label=label)
 for cl,col,label in [('C0001',classcolors['C0001'],'MGE neurogenic progenitors (C0001)'),
                     ('C0003',classcolors['C0003'],'MGE IPCs (C0003)'),('C0006',classcolors['C0006'],'Immature MGE neurons (C0006)')]]
axes[0].legend(handles=legend,loc='upper left',bbox_to_anchor=(-.04,-.08),frameon=False,fontsize=6,handletextpad=.35,borderaxespad=0,labelspacing=.3)

values=read(SOURCE/'selected_biological_candidates_sample_values.tsv')
statistics=read(SOURCE/'selected_biological_candidates_statistics.tsv')
signal=RUNTIME_ROOT/'results/paired_HET_signaling_20261031'
values+=read(signal/'Pcdh19_three_group_sample_values.tsv')
statistics+=read(signal/'Pcdh19_paired_statistics.tsv')
panels=[('Kcnc2','C0006','Immature MGE neurons'),('Pcdhga12','C0003','MGE IPCs'),('Gabrr2','neural_vascular_low','Neural; vascular-low'),('Pcdh19','all_cells','Gene-level signal')]
for i,(ax,(gene,st,title)) in enumerate(zip(axes[3:],panels)):
    rows=[r for r in values if r['gene']==gene and r['stratum']==st];assert len(rows)==9
    for s in female_samples[3:]:
        y=[float(next(r for r in rows if r['sample']==s and r['group']==g)['TMM_log2CPM']) for g in ['A_positive','Aneg_BC1']]
        ax.plot([1,2],y,c='#bfc2c5',lw=.65,zorder=1)
    for x,(group,color) in enumerate(zip(['WT_F_all','A_positive','Aneg_BC1'],[green,blue,purple])):
        y=[float(r['TMM_log2CPM']) for r in rows if r['group']==group]
        ax.scatter(x+np.array([-.07,0,.07]),y,s=21,c=color,edgecolors='white',linewidths=.4,zorder=3)
        ax.plot([x-.18,x+.18],[np.mean(y),np.mean(y)],c=color,lw=1.6,zorder=4)
    r=next(r for r in statistics if r['gene']==gene and r['stratum']==st and r['contrast']=='within_HET_Aneg_BC1_vs_Apositive')
    ax.set_title(gene+'\n'+title,fontsize=8,pad=24)
    ax.text(.5,1.035,'Paired HET P = {:.2g}'.format(float(r['PValue'])),transform=ax.transAxes,ha='center',fontsize=6.5)
    ax.set_xticks([0,1,2]);ax.set_xticklabels(['WT F','HET F\nWT-like','HET F\nKO-enr.'],fontsize=6.5)
    ax.set_xlim(-.4,2.4);ax.set_ylabel('log2 normalized CPM',fontsize=7)
    ax.tick_params(labelsize=6.5)
    for sp in ['top','right']:ax.spines[sp].set_visible(False)
    ax.text(-.24,1.28,'defg'[i],transform=ax.transAxes,fontweight='bold',fontsize=11)

fig.text(.075,.065,'Each expression point is one sample (n = 3/group); lines pair the two HET subsets.',fontsize=7)
fig.text(.075,.035,'Exploratory comparisons. UMAP overlays show 450 cells per HET sample per gate; expression uses all selected cells.',fontsize=6.5)
fig.text(.075,.01,'Pcdh19 totals combine probe evidence; both HET subsets were selected for detection, while WT females were unfiltered.',fontsize=6.5)
for ext in ['pdf','svg','png']:
    fig.savefig(OUT/('mosaic_umap_and_expression.'+ext),dpi=500,bbox_inches='tight',facecolor='white')
plt.close(fig)
legend_text='''Probe-defined heterogeneity in female E14.5 mouse MGE. (a) Saved joint, unintegrated Step 06 UMAP, displaying 8,000 cells from each of six female samples. Selected full-data consensus clusters are colored using provisional broad marker-based names; other clusters are gray. This figure does not recompute UMAP or clustering. (b,c) Locations of A-positive WT-like and A-negative B/C-positive KO-enriched cells from HET females on identical coordinates and axis limits. Each overlay contains a deterministic balanced display of 450 cells from each HET sample (1,350 points per overlay), so point density does not encode the unequal gate sizes. Printed n values give complete retained gate populations (3,949 and 16,172). (d–f) Gene-level TMM log2 CPM for Kcnc2 in C0006, Pcdhga12 in C0003, and Gabrr2 in neural/vascular-low cells. Each dot is one registered sample; horizontal lines are unweighted sample means and gray lines pair subsets from the same HET sample. Unfiltered WT females are the reference group. Expression panels use all selected cells, not the UMAP display subset. Printed P values are paired sample-level edgeR quasi-likelihood tests. Corresponding within-comparison FDRs are 0.758, 0.736 and 0.280, respectively. These genes were selected after an exploratory definition screen; nominal signals are not independent confirmation or proof of cellular interference. Broad KO-enriched gates are probe-evidence groups, not definitive DNA genotypes. Source sample independence and preparation covariates remain unverified. Raw plotted values and all comparison statistics accompany the source analysis.
'''
legend_text=legend_text.replace('(d–f)','(d–g)').replace('and Gabrr2 in neural/vascular-low cells.','Gabrr2 in neural/vascular-low cells, and total Pcdh19 counts across all cells. Pcdh19 is selection-dependent: both HET subsets require probe detection whereas WT females are unfiltered.').replace('0.758, 0.736 and 0.280, respectively.','0.758, 0.736, 0.280 and 0.00000332, respectively.')
(OUT/'FIGURE_LEGEND.txt').write_text(legend_text)
provenance=dict(embedding=str(H5),embedding_key='obsm/X_umap',embedding_refit=False,total_cells=len(full_xy),female_cells=int(sum(female)),
    seed=20261031,context_display_per_sample=8000,HET_overlay_display_per_sample=450,
    WT_like_cells=wt_population,KO_enriched_cells=ko_population,
    cell_labels='provisional broad names of consensus C0001,C0003,C0006; not pilot cluster IDs',
    source_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [CALLS,ASSIGN,SOURCE/'selected_biological_candidates_sample_values.tsv',SOURCE/'selected_biological_candidates_statistics.tsv']})
(OUT/'provenance.json').write_text(json.dumps(provenance,indent=2))
print('FIGURE_COMPLETE',json.dumps({k:provenance[k] for k in ['female_cells','WT_like_cells','KO_enriched_cells']}))
