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
fig=plt.figure(figsize=(7.2,6.0))
rects={'a':[.055,.73,.62,.20],'b':[.075,.405,.17,.215],'c':[.42,.405,.555,.215],'d':[.75,.73,.225,.20],'e':[.075,.095,.17,.215],'f':[.42,.095,.555,.215]}
fig.text(.055,.978,'Pcdh19: female mosaic states and male controls',fontsize=7,weight='bold',color=INK)
def title(k,text):
 x,y,w,h=rects[k];fig.text(x-.027,y+h+.020,k,fontsize=8,weight='bold',color=INK);fig.text(x,y+h+.020,text,fontsize=6.5,color=INK)
def clean(ax):ax.tick_params(length=2,pad=2)
# a: contextual atlas with broad, provenance-supported labels.
broad=[('RG-like progenitors',['C0002'],'#7186A6'),('MGE IPCs',['C0003','C0008'],'#CE9667'),('Neurogenic progenitors',['C0001','C0007','C0012','C0024'],'#B8A363'),('Cycling neuroblasts',['C0011'],'#9B789D'),('Immature inhibitory neurons',['C0004','C0005','C0006','C0020'],'#4D928A'),('Mixed neural / vascular',['C0015'],'#B77D86'),('Unresolved populations',[],GRAY)]
lookup={c:i for i,(_,cc,_) in enumerate(broad) for c in cc};cache=read(OUT/'atlas_sex_genotype_display.tsv.gz');xy=np.array([[float(r['UMAP1']),float(r['UMAP2'])] for r in cache]);category=np.array([lookup.get(r['cluster'],6) for r in cache]);sex=np.array([r['sex'] for r in cache]);geno=np.array([r['genotype'] for r in cache]);lim=[xy[:,0].min()-.4,xy[:,0].max()+.4,xy[:,1].min()-.4,xy[:,1].max()+.4]
po=np.random.RandomState(20261031).permutation(len(cache))
for j,(view,lab) in enumerate([('cell','Cell type'),('sex','Sex'),('genotype','Genotype')]):
 ax=fig.add_axes([.055+j*.21,.745,.185,.17])
 if view=='cell':
  cols=np.array([broad[i][2] for i in category])
 else:
  entries=[('Female','#4D928A'),('Male','#7186A6')] if view=='sex' else [('WT','#969BA3'),('HET','#9F77AB'),('KO','#CE9667')]
  data=sex if view=='sex' else geno;colors=dict(entries);cols=np.array([colors[v] for v in data])
  hh=[Line2D([0],[0],marker='o',ls='',color=col,label=name,markersize=2.5) for name,col in entries];ax.legend(handles=hh,loc='upper center',bbox_to_anchor=(.5,-.02),ncol=len(entries),frameon=False,fontsize=5,handletextpad=.3,columnspacing=.5,handlelength=.5,borderaxespad=0)
 ax.scatter(xy[po,0],xy[po,1],s=.22,c=cols[po],alpha=.85,linewidths=0,rasterized=True)
 ax.set_xlim(lim[:2]);ax.set_ylim(lim[2:]);ax.set_aspect('equal');ax.axis('off');ax.set_title(lab,fontsize=6,pad=2)
lg=fig.add_axes([.055,.653,.62,.065]);lg.axis('off');hh=[Line2D([0],[0],marker='o',ls='',color=col,label=name,markersize=2.5) for name,cc,col in broad];lg.legend(handles=hh,ncol=3,frameon=False,fontsize=5,loc='upper left',handletextpad=.35,handlelength=.5,columnspacing=.9,labelspacing=.25,borderaxespad=0)
title('a','Atlas: cell identity, sex and genotype')
# b/e: identical cell-level metric, axis range, violin settings and mean markers.
fcache=np.load(str(FEMALE/'neuronal_exon1_signal.npz'));mcache=np.load(str(OUT/'male_exon1_display.npz'));display_fractions={}
def violin(k,data,group,labels,text):
 ax=fig.add_axes(rects[k]);ds=[data[group==i] for i in range(len(labels))];non=[i for i,d in enumerate(ds) if np.ptp(d)>0];vp=ax.violinplot([ds[i] for i in non],positions=non,widths=.74,showextrema=False,points=180,bw_method=.22)
 for body in vp['bodies']:body.set_facecolor('#B1B5BB');body.set_edgecolor('#656B74');body.set_alpha(1);body.set_linewidth(.55)
 for i in range(len(labels)):
  if i not in non:ax.plot([i-.27,i+.27],[0,0],color=INK,lw=1,clip_on=False)
 ax.plot(range(len(ds)),[d.mean() for d in ds],'-D',color=INK,ms=2,lw=.6);ax.set_ylim(-.04,8);ax.set_yticks([0,2,4,6,8]);ax.set_xlim(-.46,len(ds)-.54);ax.set_xticks(range(len(ds)));ax.set_xticklabels([t.replace(' ','\n',1) for t in labels],fontsize=5.2,rotation=0,ha='center');ax.set_ylabel('Exon-1 UMI per 10,000',fontsize=5.5,labelpad=2);clean(ax);title(k,text);display_fractions[k]=[float((d>8).mean()) for d in ds]
 ax.text(.5,1.012,'0–8 shown',transform=ax.transAxes,ha='center',fontsize=5,color='#6B7077')
violin('b',fcache['exon1_CP10K'],fcache['group'],['WT female','HET WT-like','HET KO-like'],'Female probe states')
violin('e',mcache['CP10K'],mcache['group'],['WT males','KO males'],'Male genotype groups')
# A common pre-existing candidate panel makes sex comparisons direct.
fgenes=['Ank3','Kcnc2','Gnao1','Ptprd','Epha5','Snap47','Pcdha5','Pcdha10','Pcdhga5'];mgenes=['Eps8l1','Serpine2','Syt9','Ecel1','Lrfn5','Grm7','Cntn1','Sphkap','Nrp2'];genes=fgenes+mgenes
common_female={r['gene_symbol'] for r in read(DEEP/'DE_C0006/within_HET_Aneg_BC1_vs_Apositive.tsv.gz')};common_male={r['gene_symbol'] for r in read(OUT/'male_C0006_DE.tsv')};dropped=[g for g in genes if g not in common_female or g not in common_male];genes=[g for g in genes if g in common_female and g in common_male]
rows=read(OUT/'matched_display_logCPM.tsv');values={r['gene_symbol']:r for r in rows};md=read(OUT/'matched_female_display_metadata.tsv');go=[('WT_F','A_positive'),('HET_F','A_positive'),('HET_F','Aneg_BC1')];md=sorted(md,key=lambda r:(go.index((r['design_group'],r['gate'])),r['sample']));flibs=[r['library'] for r in md];mlibs=['15662-JZ-'+str(i)+'__C0006' for i in [1,2,3,10,11,12]]
fm=np.array([[np.mean([float(values[g][l]) for l in flibs[k:k+3]]) for k in [0,3,6]] for g in genes]);mm=np.array([[np.mean([float(values[g][l]) for l in mlibs[k:k+3]]) for k in [0,3]] for g in genes])
def standardize(a):return (a-a.mean(1)[:,None])/np.maximum(a.std(1)[:,None],1e-9)
fx=np.array([[float(values[g][l]) for l in flibs] for g in genes]);mx=np.array([[float(values[g][l]) for l in mlibs] for g in genes])
fz=(fm-fx.mean(1)[:,None])/np.maximum(fx.std(1)[:,None],1e-9);mz=(mm-mx.mean(1)[:,None])/np.maximum(mx.std(1)[:,None],1e-9);link=linkage(np.column_stack([fz,mz]),method='average');fd={r['gene_symbol']:r for r in read(DEEP/'DE_C0006/within_HET_Aneg_BC1_vs_Apositive.tsv.gz')};mdict={r['gene_symbol']:r for r in read(OUT/'male_C0006_DE.tsv')};fef=np.array([float(fd[g]['logFC']) if g in fd else np.nan for g in genes]);mef=np.array([float(mdict[g]['logFC']) if g in mdict else np.nan for g in genes]);orders={}
def matched_gene_block(k,z,effect,labels,text,contrast):
 x,y,w,h=rects[k];tw=.14*w;nw=.12*w;hw=.42*w;bw=.245*w;gap=.025*w
 ta=fig.add_axes([x,y,tw,h]);order=dendrogram(link,no_labels=True,orientation='left',color_threshold=0,above_threshold_color='#555',ax=ta)['leaves'];ta.invert_yaxis();ta.axis('off')
 for col in ta.collections:col.set_linewidth(.6)
 na=fig.add_axes([x+tw+gap,y,nw,h]);na.set_xlim(0,1);na.set_ylim(len(genes)-.5,-.5);na.axis('off')
 for i,j in enumerate(order):na.text(.99,i,genes[j],ha='right',va='center',fontsize=5.6)
 ha=fig.add_axes([x+tw+nw+2*gap,y,hw,h]);ha.imshow(z[order],aspect='auto',cmap='Purples',vmin=-1.5,vmax=1.5);ha.set_yticks([]);ha.set_xticks(range(len(labels)));ha.set_xticklabels([t.replace(' ','\n',1) for t in labels],fontsize=5.2,rotation=0,ha='center');clean(ha)
 ba=fig.add_axes([x+tw+nw+hw+3*gap,y,bw,h]);ev=effect[order];ba.barh(range(len(genes)),ev,height=.64,color=[BLUE if q<0 else RED for q in ev]);ba.set_ylim(len(genes)-.5,-.5);ba.set_xlim(-1.6,1.6);ba.axvline(0,color='#9B9FA5',lw=.4);ba.set_yticks([]);ba.set_xticks([-1,0,1]);ba.set_xlabel('log2FC',fontsize=5.5,labelpad=2);ba.spines['left'].set_visible(False);clean(ba)
 for i,q in enumerate(ev):
  if not np.isfinite(q):ba.plot(0,i,'x',ms=2,color='#999')
 ha.text(.5,1.012,'Relative expression',ha='center',transform=ha.transAxes,fontsize=5.3);ba.text(.5,1.012,contrast,ha='center',transform=ba.transAxes,fontsize=5.3);title(k,text);orders[k]=[genes[i] for i in order]
matched_gene_block('c',fz,fef,['WT female','HET WT-like','HET KO-like'],'Female neuronal candidates','KO-like / WT-like')
matched_gene_block('f',mz,mef,['WT males','KO males'],'Male neuronal candidates','KO / WT')
assert orders['c']==orders['f']
# d: directional source-by-receiver expression potential within female HET tissue.
ct=['C0002','C0001','C0003','C0011','C0006'];short=['RG-like','Neurogenic','MGE IPCs','Neuroblasts','Neurons'];rows=['RG-like\nprogenitors','Neurogenic\nprogenitors','MGE IPCs','Cycling\nneuroblasts','Immature\nneurons'];ex={(r['stratum'],r['gene'],r['gate']):r for r in read(SIG/'signaling_gene_expression_support.tsv')};arr=np.zeros((5,5));network=[]
for i,source in enumerate(ct):
 for j,target in enumerate(ct):
  v=[float(ex[(source,'Nrg1',gate)]['mean_CPM']) for gate in ['A_positive','Aneg_BC1']]+[float(ex[(target,'Erbb4',gate)]['mean_CPM']) for gate in ['A_positive','Aneg_BC1']]
  assert all(int(ex[(st,g,gate)]['samples_CPM_ge1'])>=2 for st,g in [(source,'Nrg1'),(target,'Erbb4')] for gate in ['A_positive','Aneg_BC1'])
  wt=np.sqrt(v[0]*v[2]);ko=np.sqrt(v[1]*v[3]);arr[i,j]=np.log2(ko/wt);network.append(dict(Nrg1_source=rows[i].replace('\n',' '),Erbb4_receiver=rows[j].replace('\n',' '),WTlike_score=float(wt),KOlike_score=float(ko),log2_KOlike_vs_WTlike=float(arr[i,j])))
x,y,w,h=rects['d'];ax=fig.add_axes([x+w*.34,y+.015,w*.62,h-.035]);im=ax.imshow(arr,aspect='auto',cmap='RdBu_r',vmin=-.65,vmax=.65);ax.set_yticks(range(5));ax.set_yticklabels(rows,fontsize=5);ax.set_xticks(range(5));ax.set_xticklabels(short,fontsize=5,rotation=60,ha='right');ax.tick_params(length=0,pad=2);ax.set_xlabel('Erbb4 receiver cell type',fontsize=5,labelpad=1);ax.set_ylabel('Nrg1 source cell type',fontsize=5,labelpad=2)
for i in range(5):
 for j in range(5):ax.text(j,i,('0.0' if abs(arr[i,j])<.05 else '{:+.1f}'.format(arr[i,j])),ha='center',va='center',fontsize=5,color='white' if abs(arr[i,j])>.4 else INK)
ax.text(.5,1.025,'KO-like / WT-like potential (log2)',transform=ax.transAxes,ha='center',fontsize=5);title('d','Female Nrg1 → Erbb4')
with open(str(OUT/'directional_source_receiver_values.tsv'),'w') as f:
 wr=csv.DictWriter(f,fieldnames=list(network[0]),delimiter='\t');wr.writeheader();wr.writerows(network)
# One compact shared key; black labels remain readable independently of colour.
cm=plt.get_cmap('Purples');legend=fig.add_axes([.59,.030,.12,.010]);gradient=np.tile(np.linspace(-1.5,1.5,100),(2,1));legend.imshow(gradient,aspect='auto',cmap=cm,vmin=-1.5,vmax=1.5);legend.axis('off');fig.text(.59,.018,'Low',fontsize=5);fig.text(.71,.018,'High',fontsize=5,ha='right');fig.text(.75,.033,'Expression scaled within sex',fontsize=5.5)
fig.add_artist(plt.Rectangle((.06,.031),.012,.010,transform=fig.transFigure,fc=BLUE,ec='none'));fig.text(.078,.031,'Lower in KO-like / KO',fontsize=5.3);fig.add_artist(plt.Rectangle((.29,.031),.012,.010,transform=fig.transFigure,fc=RED,ec='none'));fig.text(.308,.031,'Higher in KO-like / KO',fontsize=5.3)

for ext in ['png','pdf','svg']:fig.savefig(str(OUT/('compact_directional_figure.'+ext)),dpi=450,facecolor='white')
(OUT/'compact_directional_validation.json').write_text(json.dumps({'figure_size_inches':[7.2,6.0],'atlas_views':['broad cell identity','sex','genotype'],'atlas_same_embedding':True,'atlas_render_order':'one fixed random interleaving, shared by all views; no sex/genotype drawn last as a block','atlas_cells':len(cache),'short_violin_height_inches':rects['b'][3]*6,'directional_panel':'Nrg1 source rows to Erbb4 receiver columns; KO-like versus WT-like sqrt(ligand CPM * receptor CPM) ratio','genes_dropped_for_missing_female_contrast':dropped,'body_font_pt':[5.1,6.5],'B_E_same_measurement':'cell-level exon1 UMI / total cell UMI * 10000','B_E_same_y_range':[-.04,8],'B_E_same_KDE_bandwidth':.22,'B_E_above_range_fractions':display_fractions,'male_group_selection':'known genotype, all neurons; no probe gate','female_group_selection':'matched A-positive WT and HET WT-like; HET KO-like A-negative/downstream-positive','C_F_same_genes_and_order':orders['c'],'C_F_same_expression_scale':[-1.5,1.5],'C_F_same_effect_scale':[-1.6,1.6],'C_F_common_display_normalization':'joint TMM normalization of 15 sample/subset libraries; within-sex gene standardization across individual sample pseudobulks, then equal-sample group means','C_F_inferential_models':'unchanged paired HET and unpaired male edgeR fits','missing_effects':{'female':[g for g,q in zip(genes,fef) if not np.isfinite(q)],'male':[g for g,q in zip(genes,mef) if not np.isfinite(q)]},'exploratory_candidates_not_FDR_significant':True},indent=2))
with open(str(OUT/'compact_directional_gene_values.tsv'),'w') as f:
 wr=csv.writer(f,delimiter='\t');wr.writerow(['gene','WT_female','HET_WT_like','HET_KO_like','WT_male','KO_male','female_paired_logFC','male_logFC'])
 for i,g in enumerate(genes):wr.writerow([g]+fm[i].tolist()+mm[i].tolist()+[fef[i],mef[i]])
(OUT/'COMPACT_DIRECTIONAL_LEGEND.txt').write_text('A: three views of the identical saved joint embedding from all 12 samples, showing broad marker-supported annotations, sex and registered genotype. Display is balanced to 6000 QC cells per sample; all views use the same points and limits. Unresolved identities remain gray. B/E: same exon-1 readout (raw probe A UMI/total cell UMI*10000), linear 0-8 display, identical KDE bandwidth 0.22, median-independent mean diamonds; full observations determine violin shapes and tails above 8 are documented in validation. Female WT and HET WT-like neurons are both A-positive; HET KO-like are A-negative/downstream-positive, so their zero exon1 signal is imposed by selection. Known male WT and KO neurons retain all cells, without probe gating. C/F: the same 14 genes, hierarchy, gene-row order, expression colour limits and effect-axis limits. Gene set is the intersection of genes tested in both contrasts from the union of nine previous female and nine previous male candidates; four female-untested genes are removed mechanically, without selecting effect directions. Joint TMM logCPM uses all 15 sample/subset libraries for display; group means weight samples equally and heatmaps standardize each gene across individual sample pseudobulks within sex before averaging by group. A single linkage clusters the concatenated five within-sex group profiles. Sidebars retain the existing paired HET and unpaired WT/KO male edgeR contrasts; they are not refitted by visual styling. D: directional Nrg1-source to Erbb4-receiver expression potential among five neural populations in female HET tissue. Each group score is sqrt(mean source ligand CPM * mean receiver receptor CPM); color and numbers show log2(KO-like score/WT-like score). Every component is expressed at >=1 CPM in at least two of three samples in each group. Matrix axes specify ligand source and receptor receiver. This is an expression-based candidate, not measured signaling, physical contact or a CellChat probability. Statistical findings remain exploratory and do not pass gene-wide FDR correction. The six-panel layout is prepared at 182.9 x 152.4 mm, with predominantly 5-7 pt editable text and native PDF/SVG exports.\n')
print('MATCHED_PUBLICATION_FIGURE_COMPLETE',flush=True)
