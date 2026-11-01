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
genes=['Gria1','Ank3','Kcnc2','Gnao1','Ptprd','Epha5','Snap47','Pcdha5','Pcdha10','Pcdhga5']
assert all(g in values and g in de for g in genes)
x=np.stack([values[g] for g in genes]);z=(x-x.mean(1)[:,None])/x.std(1)[:,None]
diff=x[:,6:9]-x[:,3:6]
colors=['#549E89']*3+['#337EAC']*3+['#A678AC']*3
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'pdf.fonttype':42,'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False})
fig=plt.figure(figsize=(13,9))
gs=fig.add_gridspec(2,3,left=.075,right=.97,top=.82,bottom=.19,hspace=.85,wspace=.55)
fig.text(.075,.965,'Neuronal expression within the same MGE cluster',fontsize=17,weight='bold')
fig.text(.075,.93,'C0006: Lhx6 / Gad1 / Gad2-positive cells. Three female samples per group; HET subsets paired within embryo.',fontsize=10)
fig.text(.075,.898,'Exploratory candidates: none of the neuronal genes shown passes gene-wide FDR < 0.05 in this cluster.',fontsize=9,color='#6f4242')

def label(ax,l,title):
    ax.set_title(title,loc='left',fontsize=10,pad=12)
    ax.text(-.18,1.12,l,transform=ax.transAxes,fontsize=14,weight='bold')

# Identity evidence: gene-wise scaling across all consensus clusters, not across chosen rows.
markers=read(RUNTIME_ROOT/'results/composition_review_20261031/full_consensus_marker_means_source.tsv')
mg=['Nkx2-1','Lhx6','Dlx2','Gad1','Gad2','Mki67','Sox2','Sst','Pvalb']
allm=np.array([[float(r[g]) for g in mg] for r in markers]);scaled=allm/np.maximum(allm.max(0),1e-9)
ids=['C0002','C0003','C0001','C0006'];rr=[next(i for i,r in enumerate(markers) if r['']==c) for c in ids]
ax=fig.add_subplot(gs[0,0]);im=ax.imshow(allm[rr],aspect='auto',cmap='Greys',vmin=0,vmax=8)
ax.set_xticks(range(len(mg)));ax.set_xticklabels(mg,rotation=60,ha='right');ax.set_yticks(range(4));ax.set_yticklabels(ids)
label(ax,'a','Why focus on C0006?\nMarker expression across clusters')
cb=fig.colorbar(im,ax=ax,fraction=.04,pad=.03);cb.set_label('Mean expression (saved marker matrix)',fontsize=7)

def heat(slot,mat,cols,title,letter,clustercols):
    outer=gs[slot].subgridspec(2,2,height_ratios=[.19,1],width_ratios=[.16,1],hspace=.02,wspace=.02)
    lr=linkage(mat,method='average');rt=fig.add_subplot(outer[1,0]);r=dendrogram(lr,orientation='left',no_labels=True,color_threshold=0,above_threshold_color='#555555',ax=rt)['leaves'];rt.invert_yaxis();rt.axis('off')
    rb=rt.get_position();rt.set_position([rb.x0-.035,rb.y0,rb.width*.55,rb.height])
    top=fig.add_subplot(outer[0,1])
    if clustercols:
        co=dendrogram(linkage(mat.T,method='average'),no_labels=True,color_threshold=0,above_threshold_color='#555555',ax=top)['leaves']
    else:co=list(range(mat.shape[1]))
    top.axis('off');top.set_title(title,fontsize=10,pad=10,loc='left');top.text(-.35,1.5,letter,transform=top.transAxes,fontsize=14,weight='bold')
    ax=fig.add_subplot(outer[1,1]);im=ax.imshow(mat[np.ix_(r,co)],aspect='auto',cmap='RdBu_r',vmin=-2,vmax=2)
    ax.set_yticks(range(len(genes)));ax.set_yticklabels([genes[i] for i in r],fontsize=7)
    ax.set_xticks(range(len(cols)));ax.set_xticklabels([cols[i] for i in co],rotation=65,ha='right',fontsize=6.5)
    fig.colorbar(im,ax=ax,fraction=.035,pad=.03,extend='both')
    return [genes[i] for i in r],[cols[i] for i in co]
orders={}
orders['expression']=heat((0,1),z,names,'Neuronal candidates\nGene-wise expression z-score','b',True)
orders['paired']=heat((0,2),diff,['Embryo 7','Embryo 8','Embryo 9'],'Paired HET changes\nKO-enriched minus WT-like (log2 CPM)','c',False)
for j,g in enumerate(['Gria1','Ank3','Kcnc2']):
    ax=fig.add_subplot(gs[1,j]);v=values[g]
    for k in range(3):ax.plot([1+(k-1)*.09,2+(k-1)*.09],[v[3+k],v[6+k]],color='#aaa',lw=.9,zorder=1)
    for k in range(3):
        a=v[k*3:k*3+3];pos=k+np.array([-.09,0,.09]);ax.scatter(pos,a,s=32,color=colors[3*k],edgecolors='white',linewidth=.5,zorder=3)
        ax.plot([k-.18,k+.18],[a.mean(),a.mean()],color=colors[3*k],lw=2)
        for t in range(3):ax.annotate(str(t+4 if k==0 else t+7),(pos[t],a[t]),xytext=(4,3),textcoords='offset points',fontsize=6,color='#555')
    ax.set_xticks([0,1,2]);ax.set_xticklabels(['WT females\nall C0006 cells','HET females\nWT-like cells','HET females\nKO-enriched cells'],fontsize=7)
    ax.set_ylabel('TMM log2 CPM');d=de[g]
    label(ax,chr(ord('d')+j),g+'\nPaired log2FC {:+.2f}; P={:.2g}; FDR={:.2g}'.format(float(d['logFC']),float(d['PValue']),float(d['FDR'])))
fig.legend(handles=[Patch(color=colors[i*3],label=l) for i,l in enumerate(['WT female controls','WT-like cells within HET females','KO-enriched cells within HET females'])],loc='lower left',bbox_to_anchor=(.07,.092),frameon=False,ncol=3,fontsize=8)
fig.text(.075,.074,'WT-like = A-probe positive. KO-enriched = A-negative with B/C-probe evidence. These are probe-defined subsets, not confirmed genotypes.',fontsize=8)
fig.text(.075,.047,'b–c: nominally selected neuronal candidates, not an unbiased test of group separation. Lines pair HET embryos; numbers identify samples.',fontsize=8)
fig.text(.075,.02,'C0006 is a provisional broad neuronal annotation; no PV/SST subtype claim. Marker intensities do not establish cell purity.',fontsize=8)
for ext in ['pdf','png','svg']:fig.savefig(str(OUT/('neuronal_mosaic_focus.'+ext)),dpi=300,facecolor='white')
with open(str(OUT/'neuronal_candidates_expression_and_statistics.tsv'),'w') as f:
    fields=['gene','logFC','PValue','FDR']+names;w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader()
    for g in genes:
        row={'gene':g};row.update({k:de[g][k] for k in ['logFC','PValue','FDR']});row.update(dict(zip(names,values[g])));w.writerow(row)
(OUT/'validation.json').write_text(json.dumps({'cluster':'C0006','samples':9,'neuronal_candidates':len(genes),'candidate_FDR05':sum(float(de[g]['FDR'])<.05 for g in genes),'clustering':orders},indent=2))
print('NEURONAL_FOCUS_COMPLETE')
