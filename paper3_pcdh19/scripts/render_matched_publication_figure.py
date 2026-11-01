from pathlib import Path
from pcdh19_exploratory_paths import RUNTIME_ROOT
import csv,gzip,json,copy
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.cluster.hierarchy import linkage,dendrogram
ROOT=Path(__file__).resolve().parents[1];OUT=RUNTIME_ROOT/'results/female_male_signaling_20261031';FEMALE=RUNTIME_ROOT/'results/neuronal_matched_selection_20261031';DEEP=RUNTIME_ROOT/'results/mosaic_deep_dive_20261031';SIG=RUNTIME_ROOT/'results/paired_HET_signaling_20261031'
def read(p):
 with (gzip.open(str(p),'rt') if str(p).endswith('.gz') else open(str(p))) as f:return list(csv.DictReader(f,delimiter='\t'))
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':6,'axes.titlesize':6.5,'axes.labelsize':6,'xtick.labelsize':5.5,'ytick.labelsize':5.5,'pdf.fonttype':42,'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':.5,'xtick.major.width':.5,'ytick.major.width':.5})
BLUE='#326C9D';RED='#B64B50';INK='#26292E';GRAY='#BFC4CC'
fig=plt.figure(figsize=(7.2,6.5))
rects={'a':[.055,.59,.265,.33],'b':[.385,.59,.125,.33],'c':[.58,.59,.395,.33],'d':[.055,.13,.265,.33],'e':[.385,.13,.125,.33],'f':[.58,.13,.395,.33]}
fig.text(.055,.978,'Pcdh19 mosaic states and male-control comparisons',fontsize=7,weight='bold',color=INK)
def title(k,text):
 x,y,w,h=rects[k];fig.text(x-.027,y+h+.020,k,fontsize=8,weight='bold',color=INK);fig.text(x,y+h+.020,text,fontsize=6.5,color=INK)
def clean(ax):ax.tick_params(length=2,pad=2)
# a: contextual atlas with broad, provenance-supported labels.
broad=[('RG-like progenitors',['C0002'],'#7186A6'),('MGE IPCs',['C0003','C0008'],'#CE9667'),('Neurogenic progenitors',['C0001','C0007','C0012','C0024'],'#B8A363'),('Cycling neuroblasts',['C0011'],'#9B789D'),('Immature inhibitory neurons',['C0004','C0005','C0006','C0020'],'#4D928A'),('Mixed neural / vascular',['C0015'],'#B77D86'),('Unresolved populations',[],GRAY)]
lookup={c:i for i,(_,cc,_) in enumerate(broad) for c in cc};cache=[r for r in read(RUNTIME_ROOT/'results/mosaic_nature_figure_20261031/displayed_umap_cells.tsv.gz') if r['panel']=='background'];xy=np.array([[float(r['UMAP1']),float(r['UMAP2'])] for r in cache]);category=np.array([lookup.get(r['consensus_cluster'],6) for r in cache]);x,y,w,h=rects['a'];ax=fig.add_axes([x,y+.095,w,h-.095])
for i in [6,5,0,3,2,1,4]:
 m=category==i;ax.scatter(xy[m,0],xy[m,1],s=.23,color=broad[i][2],linewidths=0,rasterized=True)
ax.set_aspect('equal');ax.axis('off');lg=fig.add_axes([x-.008,y-.002,w+.03,.094]);lg.axis('off');handles=[Line2D([0],[0],marker='o',ls='',color=color,label=name,markersize=2.8) for name,cc,color in broad];lg.legend(handles=handles,ncol=1,frameon=False,fontsize=5.1,loc='upper left',handletextpad=.5,handlelength=.6,labelspacing=.22,borderaxespad=0)
title('a','Female MGE atlas')
# b/e: identical cell-level metric, axis range, violin settings and mean markers.
fcache=np.load(str(FEMALE/'neuronal_exon1_signal.npz'));mcache=np.load(str(OUT/'male_exon1_display.npz'));display_fractions={}
def violin(k,data,group,labels,text):
 ax=fig.add_axes(rects[k]);ds=[data[group==i] for i in range(len(labels))];non=[i for i,d in enumerate(ds) if np.ptp(d)>0];vp=ax.violinplot([ds[i] for i in non],positions=non,widths=.74,showextrema=False,points=180,bw_method=.22)
 for body in vp['bodies']:body.set_facecolor('#B1B5BB');body.set_edgecolor('#656B74');body.set_alpha(1);body.set_linewidth(.55)
 for i in range(len(labels)):
  if i not in non:ax.plot([i-.27,i+.27],[0,0],color=INK,lw=1,clip_on=False)
 ax.plot(range(len(ds)),[d.mean() for d in ds],'-D',color=INK,ms=2,lw=.6);ax.set_ylim(-.04,8);ax.set_yticks([0,2,4,6,8]);ax.set_xlim(-.46,len(ds)-.54);ax.set_xticks(range(len(ds)));ax.set_xticklabels(labels,fontsize=5.2,rotation=52,ha='right');ax.set_ylabel('Exon-1 UMI per 10,000',fontsize=5.5,labelpad=2);clean(ax);title(k,text);display_fractions[k]=[float((d>8).mean()) for d in ds]
 ax.text(.5,1.012,'0–8 shown',transform=ax.transAxes,ha='center',fontsize=5,color='#6B7077')
violin('b',fcache['exon1_CP10K'],fcache['group'],['WT female','HET WT-like','HET KO-like'],'Female probe states')
violin('e',mcache['CP10K'],mcache['group'],['WT males','KO males'],'Male genotype groups')
# A common pre-existing candidate panel makes sex comparisons direct.
fgenes=['Ank3','Kcnc2','Gnao1','Ptprd','Epha5','Snap47','Pcdha5','Pcdha10','Pcdhga5'];mgenes=['Eps8l1','Serpine2','Syt9','Ecel1','Lrfn5','Grm7','Cntn1','Sphkap','Nrp2'];genes=fgenes+mgenes
rows=read(OUT/'matched_display_logCPM.tsv');values={r['gene_symbol']:r for r in rows};md=read(OUT/'matched_female_display_metadata.tsv');go=[('WT_F','A_positive'),('HET_F','A_positive'),('HET_F','Aneg_BC1')];md=sorted(md,key=lambda r:(go.index((r['design_group'],r['gate'])),r['sample']));flibs=[r['library'] for r in md];mlibs=['15662-JZ-'+str(i)+'__C0006' for i in [1,2,3,10,11,12]]
fm=np.array([[np.mean([float(values[g][l]) for l in flibs[k:k+3]]) for k in [0,3,6]] for g in genes]);mm=np.array([[np.mean([float(values[g][l]) for l in mlibs[k:k+3]]) for k in [0,3]] for g in genes])
def standardize(a):return (a-a.mean(1)[:,None])/np.maximum(a.std(1)[:,None],1e-9)
fx=np.array([[float(values[g][l]) for l in flibs] for g in genes]);mx=np.array([[float(values[g][l]) for l in mlibs] for g in genes])
fz=(fm-fx.mean(1)[:,None])/np.maximum(fx.std(1)[:,None],1e-9);mz=(mm-mx.mean(1)[:,None])/np.maximum(mx.std(1)[:,None],1e-9);link=linkage(np.column_stack([fz,mz]),method='average');fd={r['gene_symbol']:r for r in read(DEEP/'DE_C0006/within_HET_Aneg_BC1_vs_Apositive.tsv.gz')};mdict={r['gene_symbol']:r for r in read(OUT/'male_C0006_DE.tsv')};fef=np.array([float(fd[g]['logFC']) if g in fd else np.nan for g in genes]);mef=np.array([float(mdict[g]['logFC']) if g in mdict else np.nan for g in genes]);orders={}
def matched_gene_block(k,z,effect,labels,text,contrast):
 x,y,w,h=rects[k];tw=.075*w;nw=.22*w;hw=.37*w;bw=.245*w;gap=.025*w
 ta=fig.add_axes([x,y,tw,h]);order=dendrogram(link,no_labels=True,orientation='left',color_threshold=0,above_threshold_color='#555',ax=ta)['leaves'];ta.invert_yaxis();ta.axis('off')
 for col in ta.collections:col.set_linewidth(.6)
 na=fig.add_axes([x+tw+gap,y,nw,h]);na.set_xlim(0,1);na.set_ylim(len(genes)-.5,-.5);na.axis('off')
 for i,j in enumerate(order):na.text(.99,i,genes[j],ha='right',va='center',fontsize=5.6)
 ha=fig.add_axes([x+tw+nw+2*gap,y,hw,h]);ha.imshow(z[order],aspect='auto',cmap='Purples',vmin=-1.5,vmax=1.5);ha.set_yticks([]);ha.set_xticks(range(len(labels)));ha.set_xticklabels(labels,fontsize=5.2,rotation=52,ha='right');clean(ha)
 ba=fig.add_axes([x+tw+nw+hw+3*gap,y,bw,h]);ev=effect[order];ba.barh(range(len(genes)),ev,height=.64,color=[BLUE if q<0 else RED for q in ev]);ba.set_ylim(len(genes)-.5,-.5);ba.set_xlim(-1.6,1.6);ba.axvline(0,color='#9B9FA5',lw=.4);ba.set_yticks([]);ba.set_xticks([-1,0,1]);ba.set_xlabel('log2FC',fontsize=5.5,labelpad=2);ba.spines['left'].set_visible(False);clean(ba)
 for i,q in enumerate(ev):
  if not np.isfinite(q):ba.plot(0,i,'x',ms=2,color='#999')
 ha.text(.5,1.012,'Relative expression',ha='center',transform=ha.transAxes,fontsize=5.3);ba.text(.5,1.012,contrast,ha='center',transform=ba.transAxes,fontsize=5.3);title(k,text);orders[k]=[genes[i] for i in order]
matched_gene_block('c',fz,fef,['WT female','HET WT-like','HET KO-like'],'Female neuronal candidates','KO-like / WT-like')
matched_gene_block('f',mz,mef,['WT males','KO males'],'Male neuronal candidates','KO / WT')
assert orders['c']==orders['f']
# d: restrained actual ligand/receptor components with a shared zero-centred scale.
ct=['C0002','C0001','C0003','C0011','C0006'];cn=['RG-like progenitors','Neurogenic progenitors','MGE IPCs','Cycling neuroblasts','Immature neurons'];de={(r['stratum'],r['gene_symbol']):r for r in read(SIG/'paired_DE_all_celltypes_with_signaling_annotations.tsv')};arr=np.array([[float(de[(c,g)]['logFC']) for g in ['Nrg1','Erbb4']] for c in ct]);x,y,w,h=rects['d'];na=fig.add_axes([x,y,w*.54,h]);na.set_xlim(0,1);na.set_ylim(4.5,-.5);na.axis('off')
for i,n in enumerate(cn):na.text(.98,i,n,ha='right',va='center',fontsize=5.5)
ax=fig.add_axes([x+w*.56,y,w*.36,h]);im=ax.imshow(arr,aspect='auto',cmap='RdBu_r',vmin=-1.2,vmax=1.2);ax.set_yticks([]);ax.set_xticks([0,1]);ax.set_xticklabels(['Nrg1\nligand','Erbb4\nreceptor'],fontsize=5.7);ax.tick_params(length=0,pad=3)
for i in range(5):
 for j in range(2):ax.text(j,i,'{:+.2f}'.format(arr[i,j]),ha='center',va='center',fontsize=5.7,color='white' if abs(arr[i,j])>.7 else INK)
ca=fig.add_axes([x+w*.95,y,w*.035,h]);cb=fig.colorbar(im,cax=ca,ticks=[-1,0,1]);cb.ax.tick_params(labelsize=5,length=2,pad=2);title('d','Female Nrg1–Erbb4 (log2FC)')
# One compact shared key; black labels remain readable independently of colour.
cm=plt.get_cmap('Purples');legend=fig.add_axes([.59,.030,.12,.010]);gradient=np.tile(np.linspace(-1.5,1.5,100),(2,1));legend.imshow(gradient,aspect='auto',cmap=cm,vmin=-1.5,vmax=1.5);legend.axis('off');fig.text(.59,.018,'Low',fontsize=5);fig.text(.71,.018,'High',fontsize=5,ha='right');fig.text(.75,.033,'Expression scaled within sex',fontsize=5.5)
fig.add_artist(plt.Rectangle((.06,.031),.012,.010,transform=fig.transFigure,fc=BLUE,ec='none'));fig.text(.078,.031,'Lower in KO-like / KO',fontsize=5.3);fig.add_artist(plt.Rectangle((.29,.031),.012,.010,transform=fig.transFigure,fc=RED,ec='none'));fig.text(.308,.031,'Higher in KO-like / KO',fontsize=5.3)
fig.text(.523,.031,'×',fontsize=6,color='#999');fig.text(.538,.031,'Not tested',fontsize=5)
for ext in ['png','pdf','svg']:fig.savefig(str(OUT/('matched_publication_figure.'+ext)),dpi=450,facecolor='white')
(OUT/'matched_publication_validation.json').write_text(json.dumps({'figure_size_inches':[7.2,6.5],'body_font_pt':[5.1,6.5],'B_E_same_measurement':'cell-level exon1 UMI / total cell UMI * 10000','B_E_same_y_range':[-.04,8],'B_E_same_KDE_bandwidth':.22,'B_E_above_range_fractions':display_fractions,'male_group_selection':'known genotype, all neurons; no probe gate','female_group_selection':'matched A-positive WT and HET WT-like; HET KO-like A-negative/downstream-positive','C_F_same_genes_and_order':orders['c'],'C_F_same_expression_scale':[-1.5,1.5],'C_F_same_effect_scale':[-1.6,1.6],'C_F_common_display_normalization':'joint TMM normalization of 15 sample/subset libraries; within-sex gene standardization across individual sample pseudobulks, then equal-sample group means','C_F_inferential_models':'unchanged paired HET and unpaired male edgeR fits','missing_effects':{'female':[g for g,q in zip(genes,fef) if not np.isfinite(q)],'male':[g for g,q in zip(genes,mef) if not np.isfinite(q)]},'exploratory_candidates_not_FDR_significant':True},indent=2))
with open(str(OUT/'matched_publication_gene_values.tsv'),'w') as f:
 wr=csv.writer(f,delimiter='\t');wr.writerow(['gene','WT_female','HET_WT_like','HET_KO_like','WT_male','KO_male','female_paired_logFC','male_logFC'])
 for i,g in enumerate(genes):wr.writerow([g]+fm[i].tolist()+mm[i].tolist()+[fef[i],mef[i]])
(OUT/'MATCHED_PUBLICATION_LEGEND.txt').write_text('A: saved female MGE UMAP with broad marker-supported annotations; unresolved identities remain gray. B/E: same exon-1 readout (raw probe A UMI/total cell UMI*10000), linear 0-8 display, identical KDE bandwidth 0.22, median-independent mean diamonds; full observations determine violin shapes and tails above 8 are documented in validation. Female WT and HET WT-like neurons are both A-positive; HET KO-like are A-negative/downstream-positive, so their zero exon1 signal is imposed by selection. Known male WT and KO neurons retain all cells, without probe gating. C/F: the same 18 genes, hierarchy, gene-row order, expression colour limits and effect-axis limits. Gene set is the union of nine previous female candidates and nine previous male candidates, not newly selected for matching directions. Joint TMM logCPM uses all 15 sample/subset libraries for display; group means weight samples equally and heatmaps standardize each gene across individual sample pseudobulks within sex before averaging by group. A single linkage clusters the concatenated five within-sex group profiles. Sidebars retain the existing paired HET and unpaired WT/KO male edgeR contrasts; they are not refitted by visual styling. D: fitted Nrg1 and Erbb4 paired HET gene log2FC by broad neural population, not measured interaction or signaling activity. Statistical findings remain exploratory and do not pass gene-wide FDR correction. The six-panel layout is prepared at 182.9 x 165.1 mm, with predominantly 5-7 pt editable text and native PDF/SVG exports.\n')
print('MATCHED_PUBLICATION_FIGURE_COMPLETE',flush=True)
