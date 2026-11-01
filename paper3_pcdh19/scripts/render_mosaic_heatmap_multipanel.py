"""Heatmap-centered multi-panel exploration with genuine hierarchical clustering."""
import csv
import gzip
import copy
import json
from pathlib import Path
from pcdh19_exploratory_paths import RUNTIME_ROOT
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgba
from matplotlib.lines import Line2D
from scipy.cluster.hierarchy import linkage,dendrogram,leaves_list

ROOT=Path(__file__).resolve().parents[1]
OUT=RUNTIME_ROOT/'results/mosaic_heatmap_multipanel_20261031'
BASE=RUNTIME_ROOT/'results/mosaic_deep_dive_20261031'
SIG=RUNTIME_ROOT/'results/paired_HET_signaling_20261031'
def read(p):
    with (gzip.open(p,'rt') if p.suffix=='.gz' else p.open()) as f:return list(csv.DictReader(f,delimiter='\t'))
def write(name,header,rows):
    with (OUT/name).open('w') as f:
        w=csv.writer(f,delimiter='\t');w.writerow(header);w.writerows(rows)
meta=read(OUT/'nine_sample_metadata.tsv');raw=read(OUT/'nine_sample_log2CPM.tsv')
libs=[r['library'] for r in meta];ids=[r['gene_id'] for r in raw];symbols=[r['gene_symbol'] for r in raw]
expr=np.array([[float(r[l]) for l in libs] for r in raw]);index={g:i for i,g in enumerate(ids)}
symbol_index={g:i for i,g in enumerate(symbols)}
groups=['WT_F_all' if r['design_group']=='WT_F' else r['gate'] for r in meta]
colors={'WT_F_all':'#549E89','A_positive':'#337EAC','Aneg_BC1':'#A678AC'}
column_colors=[colors[g] for g in groups]
labels=[('W'+r['sample'].split('-')[-1]) if r['design_group']=='WT_F' else ('H'+r['sample'].split('-')[-1]+(' A+' if r['gate']=='A_positive' else ' BC+')) for r in meta]
ann=read(SIG/'paired_DE_all_celltypes_with_signaling_annotations.tsv')
lookup={(r['stratum'],r['gene_symbol']):r for r in ann}
all_de=[r for r in ann if r['stratum']=='all_cells' and r['gene_symbol']!='Pcdh19' and float(r['FDR'])<.05]
all_de.sort(key=lambda r:float(r['PValue']))
top=all_de[:40]
rows=[index[r['gene_id']] for r in top]
def zscore(x):
    return (x-x.mean(axis=1,keepdims=True))/np.maximum(x.std(axis=1,keepdims=True),1e-8)
z=zscore(expr[rows]);names=[r['gene_symbol'] for r in top]
paired=[]
for sample in ['15662-JZ-7','15662-JZ-8','15662-JZ-9']:
    a=next(i for i,r in enumerate(meta) if r['sample']==sample and r['gate']=='A_positive')
    b=next(i for i,r in enumerate(meta) if r['sample']==sample and r['gate']=='Aneg_BC1')
    paired.append(expr[rows,b]-expr[rows,a])
paired=np.column_stack(paired)
write('top40_DEG_row_zscores.tsv',['gene']+labels,[[n]+list(v) for n,v in zip(names,z)])
write('top40_paired_HET_log2CPM_differences.tsv',['gene','H7','H8','H9'],[[n]+list(v) for n,v in zip(names,paired)])

plt.rcParams.update({'font.family':'DejaVu Sans','font.size':7,'axes.linewidth':.65,'pdf.fonttype':42,'svg.fonttype':'none'})
clustering=[]
def heat(fig,ax,matrix,rowlabels,collabels,title,letter,cmap='RdBu_r',vmin=-2,vmax=2,cc=None,fontsize=5.2):
    box=ax.get_position();w,h=box.width,box.height
    clean=np.nan_to_num(matrix,nan=0.)
    rl=linkage(clean,method='average',metric='euclidean') if len(rowlabels)>1 else None
    cl=linkage(clean.T,method='average',metric='euclidean') if len(collabels)>1 else None
    ro=leaves_list(rl) if rl is not None else np.array([0]);co=leaves_list(cl) if cl is not None else np.array([0])
    ax.set_position([box.x0+.28*w,box.y0+.06*h,.62*w,.74*h])
    cm=copy.copy(plt.get_cmap(cmap));cm.set_bad('#dddddd')
    im=ax.imshow(np.ma.masked_invalid(matrix[np.ix_(ro,co)]),aspect='auto',cmap=cm,vmin=vmin,vmax=vmax,interpolation='nearest')
    ax.set_yticks(range(len(ro)));ax.set_yticklabels([rowlabels[i] for i in ro],fontsize=fontsize)
    ax.set_xticks(range(len(co)));ax.set_xticklabels([collabels[i] for i in co],rotation=60,ha='right',fontsize=5.5)
    ax.tick_params(length=0,pad=2)
    for sp in ax.spines.values():sp.set_visible(False)
    rd=fig.add_axes([box.x0,box.y0+.06*h,.095*w,.74*h]);dendrogram(rl,orientation='left',no_labels=True,ax=rd,color_threshold=0,above_threshold_color='#555555') if rl is not None else None
    rd.invert_yaxis();rd.axis('off')
    cd=fig.add_axes([box.x0+.28*w,box.y0+.84*h,.62*w,.10*h]);dendrogram(cl,no_labels=True,ax=cd,color_threshold=0,above_threshold_color='#555555') if cl is not None else None
    cd.axis('off')
    if cc is not None:
        strip=fig.add_axes([box.x0+.28*w,box.y0+.815*h,.62*w,.015*h]);strip.imshow(np.array([[to_rgba(cc[i]) for i in co]]),aspect='auto');strip.axis('off')
    cbax=fig.add_axes([box.x0+.94*w,box.y0+.23*h,.025*w,.40*h]);bar=fig.colorbar(im,cax=cbax,extend='both');bar.ax.tick_params(labelsize=5,length=2)
    fig.text(box.x0+.5*w,box.y1+.012,title,ha='center',fontsize=8,fontweight='medium')
    fig.text(box.x0-.018,box.y1+.012,letter,fontweight='bold',fontsize=12)
    clustering.append(dict(panel=letter,title=title,row_order=[rowlabels[i] for i in ro],column_order=[collabels[i] for i in co],method='average linkage; Euclidean distance; missing effects imputed zero only for clustering'))
    return ro,co

fig=plt.figure(figsize=(12,13.5));gs=fig.add_gridspec(3,3,left=.055,right=.97,top=.95,bottom=.10,wspace=.48,hspace=.70,height_ratios=[1.25,1.1,1])
axes=[fig.add_subplot(gs[i,j]) for i in range(3) for j in range(3)]

# a: one compact context panel; the added information is in the heatmaps.
cache=read(RUNTIME_ROOT/'results/mosaic_nature_figure_20261031/displayed_umap_cells.tsv.gz')
r=[r for r in cache if r['panel']=='background'];xy=np.array([[float(x['UMAP1']),float(x['UMAP2'])] for x in r]);ct=np.array([x['consensus_cluster'] for x in r])
ax=axes[0];ax.scatter(xy[:,0],xy[:,1],s=.4,c='#d6d8db',linewidths=0,rasterized=True)
for c,col in [('C0001','#C5A35D'),('C0003','#D27A50'),('C0006','#438D82')]:
    mask=ct==c;ax.scatter(xy[mask,0],xy[mask,1],s=.7,c=col,linewidths=0,rasterized=True)
ax.set_aspect('equal');ax.set_xticks([]);ax.set_yticks([]);ax.axis('off');ax.set_title('Cell landscape',fontsize=9,pad=15)
ax.text(-.08,1.06,'a',transform=ax.transAxes,fontweight='bold',fontsize=12)
ax.legend(handles=[Line2D([0],[0],marker='o',ls='',color=c,label=n,markersize=4) for c,n in [('#C5A35D','MGE progenitors'),('#D27A50','MGE IPCs'),('#438D82','Immature MGE neurons')]],frameon=False,loc='upper left',bbox_to_anchor=(0,-.05),fontsize=6)

# b,c: row and column dendrograms are actually computed, not imposed group ordering.
heat(fig,axes[1],z,names,labels,'Top 40 paired HET DEGs\nGene-wise z-score','b',cc=column_colors)
heat(fig,axes[2],paired,names,['H7','H8','H9'],'Paired KO-enriched − WT-like\nlog2 CPM differences','c',vmin=-2,vmax=2)

# d: gene effect localization across biological strata.
signaling=[r for r in all_de if r['roles']!='' and r['flag']==''][:25]
sg=[r['gene_symbol'] for r in signaling]
strata=['all_cells','vascular_marked','neural_vascular_low','C0001','C0002','C0003','C0006','C0011']
sn=['All cells','Vascular+','Neural/vasc-low','MGE prog.','RG-like','IPC','MGE neuron','Cycling NB']
m=np.full((len(sg),len(strata)),np.nan)
for i,g in enumerate(sg):
    for j,st in enumerate(strata):
        r=lookup.get((st,g))
        if r:m[i,j]=float(r['logFC'])
heat(fig,axes[3],m,sg,sn,'Signaling DEG localization\nKO-enriched / WT-like log2FC','d',fontsize=5.6)

# e: pathway overrepresentation, not pathway activation.
enr=read(SIG/'signaling_pathway_overrepresentation.tsv')
paths=sorted({r['pathway'] for r in enr if r['stratum']=='all_cells' and float(r['PValue'])<.1},key=lambda p:min(float(r['PValue']) for r in enr if r['pathway']==p))[:15]
em=np.zeros((len(paths),len(strata)))
ed={(r['pathway'],r['stratum']):r for r in enr}
for i,p in enumerate(paths):
    for j,st in enumerate(strata):
        r=ed.get((p,st))
        if r:em[i,j]=min(5.,-np.log10(max(float(r['FDR']),1e-300)))
heat(fig,axes[4],em,paths,sn,'Signaling gene-set enrichment\n−log10 BH FDR','e',cmap='YlGnBu',vmin=0,vmax=5,fontsize=6)

# f: complete corrected and nominal counts, not only selected positive examples.
overview=read(SIG/'paired_DE_overview.tsv');ov={r['stratum']:r for r in overview}
ax=axes[5];ii=np.arange(len(strata));q=np.array([int(ov[st]['signaling_DE_FDR05']) for st in strata]);nom=np.array([int(ov[st]['signaling_nominal01']) for st in strata])
ax.barh(ii-.17,nom,height=.32,color='#d5d7db',label='Nominal P <0.01');ax.barh(ii+.17,q,height=.32,color='#A678AC',label='FDR <0.05')
ax.set_yticks(ii);ax.set_yticklabels(sn,fontsize=6);ax.invert_yaxis();ax.set_xlabel('Signaling-associated genes',fontsize=7);ax.set_title('Where signaling DEGs occur',fontsize=8)
ax.spines['top'].set_visible(False);ax.spines['right'].set_visible(False);ax.legend(frameon=False,fontsize=6,loc='lower right');ax.text(-.20,1.06,'f',transform=ax.transAxes,fontweight='bold',fontsize=12)

# g: annotated ligand/receptor expression availability. Dot size is sample support, not cell fraction.
exp=read(SIG/'signaling_gene_expression_support.tsv');ex={(r['stratum'],r['gene'],r['gate']):r for r in exp}
lrgenes=['Cxcl12','Cxcr4','Ackr3','Vegfa','Vegfc','Flt1','Kdr','Ntn1','Unc5b','Sema3g','Nrp2','Plxna1']
contexts=[(st,g) for st in ['vascular_marked','C0001','C0003','C0006'] for g in ['A_positive','Aneg_BC1']]
ax=axes[6]
for i,g in enumerate(lrgenes):
    for j,(st,gate) in enumerate(contexts):
        r=ex.get((st,g,gate))
        if r:
            n=int(r['samples_CPM_ge1']);val=np.log2(1+float(r['mean_CPM']))
            dot=ax.scatter(j,i,s=6+n*12,c=[val],vmin=0,vmax=8,cmap='viridis',edgecolors='none')
ax.set_yticks(range(len(lrgenes)));ax.set_yticklabels(lrgenes,fontsize=6);ax.invert_yaxis()
ax.set_xticks(range(len(contexts)));ax.set_xticklabels([n+' '+g for n in ['Vasc','Prog','IPC','Neuron'] for g in ['A+','BC+']],rotation=60,ha='right',fontsize=5.5)
ax.set_title('Ligand–receptor expression support',fontsize=8);ax.tick_params(length=0)
fig.colorbar(dot,ax=ax,fraction=.045,pad=.02,label='log2(1 + mean CPM)');ax.text(-.20,1.06,'g',transform=ax.transAxes,fontweight='bold',fontsize=12)

# h: sample-level composition under the same broad gates.
comp=read(BASE/'sample_gate_stratum_composition.tsv');ci={(r['sample'],r['gate'],r['stratum']):r for r in comp}
ax=axes[7];cats=['C0001','C0002','C0003','C0006','C0011'];catcols=['#C5A35D','#7C94B8','#D27A50','#438D82','#B78CA3','#d6d8db'];bottom=np.zeros(9)
fixed=[i for group in ['WT_F_all','A_positive','Aneg_BC1'] for i,g in enumerate(groups) if g==group]
fractions=[]
for cat in cats:
    vals=np.array([float(ci[meta[i]['sample'],meta[i]['gate'],cat]['fraction']) for i in fixed]);fractions.append(vals)
fractions.append(1-np.sum(fractions,axis=0))
for name,col,v in zip(['MGE prog.','RG-like','IPC','MGE neuron','Cycling NB','Other'],catcols,fractions):
    ax.bar(np.arange(9),v,bottom=bottom,color=col,width=.82,label=name);bottom+=v
assert np.allclose(bottom,1)
ax.set_xticks(range(9));ax.set_xticklabels([labels[i] for i in fixed],rotation=60,ha='right',fontsize=5.5)
ax.set_ylabel('Fraction within sample/gate',fontsize=7);ax.set_ylim(0,1);ax.set_title('Cell composition by sample',fontsize=8);ax.spines['top'].set_visible(False);ax.spines['right'].set_visible(False)
ax.legend(frameon=False,fontsize=5.5,ncol=2,loc='upper center',bbox_to_anchor=(.5,-.34));ax.text(-.2,1.06,'h',transform=ax.transAxes,fontweight='bold',fontsize=12)

# i: sample relationships from variable genes selected without a DE test.
hvgids=[r['gene_id'] for r in read(OUT/'top60_variable_genes.tsv')]
hvgrows=[index[g] for g in hvgids];hz=zscore(expr[hvgrows]);corr=np.corrcoef(hz.T)
heat(fig,axes[8],corr,labels,labels,'Unsupervised sample relationships\nHVG-centered expression correlation','i',cc=column_colors,vmin=-1,vmax=1,fontsize=5.5)

fig.text(.055,-.015,'WT F: green; HET A+: blue; HET A−/BC+: purple. Each expression column is one registered sample/subset; H7/H8/H9 are paired.',fontsize=7)
fig.text(.055,-.035,'b–c use genes selected by paired HET DE and are not independent validation. Gray heatmap entries are not tested/expressed.',fontsize=7)
fig.text(.055,-.055,'Average-linkage dendrograms are data-derived. Signaling panels show expression/enrichment, not demonstrated communication. Exploratory analysis.',fontsize=7)
for ext in ['pdf','png','svg']:fig.savefig(OUT/('expanded_heatmap_multipanel.'+ext),dpi=350,bbox_inches='tight',facecolor='white')
plt.close(fig)

# Larger standalone heatmaps preserve all gene labels and provide an unselected counterpart.
for name,ix,title in [('all161_paired_DEGs',[index[r['gene_id']] for r in all_de],'All 161 non-gate paired HET DEGs'),
                      ('unsupervised_top60_variable_genes',hvgrows,'Top 60 variable genes; no DE-based selection')]:
    fig,ax=plt.subplots(figsize=(7,18 if len(ix)>100 else 10))
    fig.subplots_adjust(left=.08,right=.91,top=.90,bottom=.13)
    heat(fig,ax,zscore(expr[ix]),[symbols[i] for i in ix],labels,title,'',cc=column_colors,fontsize=5.5 if len(ix)>100 else 6)
    fig.text(.08,.025,'Gene-wise z-score of TMM log2CPM; each column is one sample/subset. Green = WT F, blue = HET A+, purple = HET A−/BC+.',fontsize=7)
    for ext in ['pdf','png']:fig.savefig(OUT/(name+'.'+ext),dpi=300,bbox_inches='tight',facecolor='white')
    plt.close(fig)
(OUT/'clustering_orders_and_methods.json').write_text(json.dumps(clustering,indent=2))
legend='''Expanded exploratory HET mosaic-expression figure. (a) Saved Step 06 UMAP context; no refitting. (b) Forty genes with the smallest paired-HET P values among FDR<0.05 non-Pcdh19 genes, clustered across nine sample/subset pseudobulks. Values are per-gene z-scores of joint TMM log2CPM; both rows and columns use average-linkage Euclidean clustering. Group strips show known WT females (green), A-positive HET cells (blue), and A-negative B/C-positive HET cells (purple). (c) The same genes shown as paired KO-enriched minus WT-like log2CPM differences in H7,H8,H9, with data-derived dendrograms. These two panels use DE-selected genes and are not independent evidence of group separation. (d) The first 25 curated signaling-associated all-cell DEGs followed across cell strata; color is fitted log2FC, gray is not tested/expressed. Missing values are imputed to neutral zero for clustering only and remain gray in the display. (e) Curated pathway overrepresentation, expressed as minus log10 BH FDR, not pathway activity. (f) Corrected and nominal signaling-gene counts in every displayed stratum; thresholds are not nested categories. (g) Ligand/receptor component expression availability: color is log2(1+mean CPM), dot size reflects zero to three samples with CPM>=1, not the fraction of cells. Optional cofactors and spatial adjacency are not assessed. (h) Relative consensus-cluster composition within each sample/gate; these selected populations do not estimate unbiased tissue proportions. (i) Sample expression correlation using the top60 variable genes chosen by expression variance alone, after gene-wise centering/scaling; this is not based on the DE gene list. Larger standalone figures display all161 DEGs and the independent top60 variable-gene set. Three HET samples are the paired biological units; overlapping gates, marker strata and reused embryos are not independent replication. Original independent embryo/litter structure and preparation covariates remain unverified. Pcdh19 is excluded from DE/HVG heatmap selection because it defines the gate and is retained in the separately provided UMAP/gene-panel figure.
'''
(OUT/'FIGURE_LEGEND.txt').write_text(legend)
(OUT/'VALIDATION.json').write_text(json.dumps(dict(samples=9,paired_HET_samples=3,selected_main_DEGs=len(top),all_DEGs=len(all_de),HVGs=len(hvgrows),composition_sums_to_one=True,clustering='average linkage / Euclidean',DE_heatmap_selection_explicit=True,raw_expression_unchanged=True),indent=2))
print('EXPANDED_FIGURE_COMPLETE',len(top),len(all_de),len(hvgrows),flush=True)
