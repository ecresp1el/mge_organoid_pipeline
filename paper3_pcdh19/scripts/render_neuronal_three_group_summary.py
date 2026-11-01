from pathlib import Path
from pcdh19_exploratory_paths import RUNTIME_ROOT
import csv, gzip, json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from scipy.cluster.hierarchy import linkage, dendrogram

ROOT=Path(__file__).resolve().parents[1]
OUT=RUNTIME_ROOT/'results/neuronal_mosaic_focus_20261031'
BASE=RUNTIME_ROOT/'results/mosaic_deep_dive_20261031'
def read(p):
    f=gzip.open(str(p),'rt') if str(p).endswith('.gz') else open(str(p))
    with f:return list(csv.DictReader(f,delimiter='\t'))
data=read(OUT/'nine_sample_log2CPM.tsv');md=read(OUT/'nine_sample_metadata.tsv')
groups=[('WT_F','all'),('HET_F','A_positive'),('HET_F','Aneg_BC1')]
md=sorted(md,key=lambda x:(groups.index((x['design_group'],x['gate'])),x['sample']))
libs=[x['library'] for x in md]
names=['WT {}'.format(x['sample'].split('-')[-1]) if x['design_group']=='WT_F' else 'HET {} {}'.format(x['sample'].split('-')[-1],'WT-like' if x['gate']=='A_positive' else 'KO-enriched') for x in md]
values={x['gene_symbol']:np.array([float(x[l]) for l in libs]) for x in data}
de={x['gene_symbol']:x for x in read(BASE/'DE_C0006/within_HET_Aneg_BC1_vs_Apositive.tsv.gz')}
# Neuronal candidates chosen for nominal paired differences, not FDR discoveries.
genes=['Ank3','Kcnc2','Gnao1','Ptprd','Epha5','Snap47','Pcdha5','Pcdha10','Pcdhga5']
assert all(g in values and g in de for g in genes)
x=np.stack([values[g] for g in genes]);z=(x-x.mean(1)[:,None])/x.std(1)[:,None]
diff=x[:,6:9]-x[:,3:6]
colors=['#549E89']*3+['#337EAC']*3+['#A678AC']*3
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'pdf.fonttype':42,'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False})
fig=plt.figure(figsize=(10,5.4))
gs=fig.add_gridspec(2,3,left=.075,right=.97,top=.88,bottom=.16,hspace=.65,wspace=.46,height_ratios=[2.5,1],width_ratios=[1,1,1.1])
labels=['WT females','HET WT-like','HET KO-like']
fig.text(.075,.955,'Female MGE · Lhx6-positive immature inhibitory neurons',fontsize=11,weight='bold')
means=np.stack([x[:,i:i+3].mean(1) for i in [0,3,6]],axis=1)

def panel(ax,letter,title):
    ax.set_title(title,loc='left',fontsize=9,pad=8)
    ax.text(-.16,1.06,letter,transform=ax.transAxes,fontsize=11,weight='bold')

def standardize(a):return (a-a.mean(1)[:,None])/np.maximum(a.std(1)[:,None],1e-9)
orders={}
cache=read(RUNTIME_ROOT/'results/mosaic_nature_figure_20261031/displayed_umap_cells.tsv.gz')
r=[q for q in cache if q['panel']=='background'];xy=np.array([[float(q['UMAP1']),float(q['UMAP2'])] for q in r]);focus=np.array([q['consensus_cluster']=='C0006' for q in r])
ax=fig.add_subplot(gs[0,0]);ax.scatter(xy[:,0],xy[:,1],s=.45,color='#ddd',linewidths=0,rasterized=True);ax.scatter(xy[focus,0],xy[focus,1],s=.7,color='#444',linewidths=0,rasterized=True)
ax.set_aspect('equal');ax.axis('off');panel(ax,'a','Neuronal population')
ax.text(.02,-.03,'Dark: inhibitory neurons · Gray: other cells',transform=ax.transAxes,fontsize=7)

# Separate axes reserve space for the tree and gene names; neither overlaps the other.
outer=gs[0,1].subgridspec(1,4,width_ratios=[.55,.48,1,.07],wspace=.07)
tree=fig.add_subplot(outer[0,0]);mat=standardize(means)
order=dendrogram(linkage(mat,method='average'),orientation='left',no_labels=True,color_threshold=0,above_threshold_color='#333',ax=tree)['leaves'];tree.invert_yaxis();tree.axis('off')
for c in tree.collections:c.set_linewidth(1.2)
nameax=fig.add_subplot(outer[0,1]);nameax.set_xlim(0,1);nameax.set_ylim(len(genes)-.5,-.5);nameax.axis('off')
for i,k in enumerate(order):nameax.text(.98,i,genes[k],ha='right',va='center',fontsize=7)
ax=fig.add_subplot(outer[0,2]);im=ax.imshow(mat[order],aspect='auto',cmap='Purples',vmin=-1.5,vmax=1.5);ax.set_yticks([]);ax.set_xticks(range(3));ax.set_xticklabels(labels,rotation=30,ha='right',fontsize=6.5)
cax=fig.add_subplot(outer[0,3]);cb=fig.colorbar(im,cax=cax);cb.set_ticks([-1.5,1.5]);cb.set_ticklabels(['Low','High']);cb.ax.tick_params(labelsize=7)
# Shift the complete block left and give branches extra horizontal room.
for aa in [tree,nameax,ax,cax]:
    bb=aa.get_position();aa.set_position([bb.x0-.035,bb.y0,bb.width,bb.height])
bb=tree.get_position();tree.set_position([bb.x0-.025,bb.y0,bb.width+.025,bb.height])
# Put the panel title over the complete tree/name/matrix block.
box=gs[0,1].get_position(fig);fig.text(box.x0-.05,box.y1+.018,'b',fontsize=11,weight='bold');fig.text(box.x0-.015,box.y1+.018,'Clustered neuronal expression',fontsize=9)
orders['candidates']=[genes[k] for k in order]

expanded_genes=genes+['Lhx8','Sptbn2','Mctp1','Dlx1','Robo2','Erbb4']
assert all(g in de for g in expanded_genes)
ax=fig.add_subplot(gs[:,2]);effectgenes=sorted(expanded_genes,key=lambda g:float(de[g]['logFC']))
fc=np.array([float(de[g]['logFC']) for g in effectgenes]);ax.barh(range(len(fc)),fc,color=['#2166AC' if v<0 else '#B2182B' for v in fc])
ax.axvline(0,color='#888',lw=.6);ax.set_yticks(range(len(fc)));ax.set_yticklabels(effectgenes,fontsize=7);ax.invert_yaxis();ax.set_xlabel('KO-like / WT-like (log2 fold change)',fontsize=7)
panel(ax,'c','Within HET females')
ax.text(0,-.12,'Lower in KO-like',color='#2166AC',transform=ax.transAxes,fontsize=7,ha='left');ax.text(1,-.12,'Higher in KO-like',color='#B2182B',transform=ax.transAxes,fontsize=7,ha='right')
small=gs[1,:2].subgridspec(1,3,wspace=.6)
for j,g in enumerate(['Pcdh19','Ank3','Kcnc2']):
    ax=fig.add_subplot(small[0,j]);v=values[g]
    cell=np.load(str(OUT/'neuronal_violin_cells.npz'));column=list(cell['genes']).index(g)
    distributions=[cell['CP10K'][cell['group']==k,column] for k in range(3)]
    vp=ax.violinplot(distributions,positions=[0,1,2],widths=.8,showmeans=False,showmedians=True,showextrema=False,points=150)
    for body in vp['bodies']:body.set_facecolor('#888888');body.set_edgecolor('#444444');body.set_alpha(.8);body.set_linewidth(.6)
    vp['cmedians'].set_color('#222222');vp['cmedians'].set_linewidth(1)
    upper=max(float(np.percentile(d,95)) for d in distributions)
    upper=max(upper,1.)*1.08
    ax.set_ylim(0,upper)
    avg=np.array([d.mean() for d in distributions]);ax.plot([0,1,2],avg,color='#222',lw=.7,marker='D',ms=3,zorder=5)
    for k,d in enumerate(distributions):
        frac=float((d>upper).mean())*100
        ax.plot(k,upper*.98,marker='^',ms=3,color='#777',clip_on=False)
        ax.text(k,upper*1.025,'{:.1f}% above'.format(frac),ha='center',fontsize=5,color='#666',clip_on=False)
    if j==0:ax.text(.02,.88,'◆ mean',transform=ax.transAxes,fontsize=6,color='#222')
    ax.set_xlim(-.4,2.4);ax.set_xticks(range(3));ax.set_xticklabels(labels,fontsize=6,rotation=30,ha='right');ax.tick_params(axis='y',labelsize=7)
    ax.set_ylabel('UMI per 10,000',fontsize=6.5);panel(ax,chr(ord('d')+j),g)
for ext in ['pdf','png','svg']:fig.savefig(str(OUT/('neuronal_three_group_summary.'+ext)),dpi=350,facecolor='white')
with open(str(OUT/'neuronal_three_group_means.tsv'),'w') as f:
    w=csv.writer(f,delimiter='\t');w.writerow(['Gene','WT females','WT-like cells in HET females','KO-like cells in HET females'])
    for g in ['Pcdh19']+expanded_genes:w.writerow([g]+[values[g][k:k+3].mean() for k in [0,3,6]])
(OUT/'group_summary_validation.json').write_text(json.dumps({'displayed_groups':labels,'individual_samples_displayed':False,'cluster_codes_displayed':False,'equal_sample_weighting':True,'gene_order':orders,'panel_C_genes':effectgenes,'panel_C_gene_count':len(effectgenes),'expression_display':'full-data individual-cell violin; linear UMI per 10000; zeros included; displayed range at 1.08 times maximum group 95th percentile; above-range fractions marked; diamonds are cell means','bottom_text_removed':True,'Gria1_displayed':False,'UMAP_added':True,'Pcdh19_expression_displayed':True,'tree_has_dedicated_axis':True,'direction_reference':'KO-like relative to WT-like within HET females','direction_colors':{'lower':'blue','higher':'red'}},indent=2))
with open(str(OUT/'expanded_neuronal_effects.tsv'),'w') as f:
    w=csv.writer(f,delimiter='\t');w.writerow(['Gene','log2FC','PValue','FDR'])
    for g in effectgenes:w.writerow([g,de[g]['logFC'],de[g]['PValue'],de[g]['FDR']])
print('THREE_GROUP_FIGURE_COMPLETE')
