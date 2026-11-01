from pathlib import Path
from pcdh19_exploratory_paths import RUNTIME_ROOT
import csv,gzip,json,math,collections
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch,Patch
from scipy.cluster.hierarchy import linkage,dendrogram
ROOT=Path(__file__).resolve().parents[1];OUT=RUNTIME_ROOT/'results/female_male_signaling_20261031';FEMALE=RUNTIME_ROOT/'results/neuronal_matched_selection_20261031';DEEP=RUNTIME_ROOT/'results/mosaic_deep_dive_20261031';SIG=RUNTIME_ROOT/'results/paired_HET_signaling_20261031'
def read(p):
 with (gzip.open(str(p),'rt') if str(p).endswith('.gz') else open(str(p))) as f:return list(csv.DictReader(f,delimiter='\t'))
def write(p,rows):
 with open(str(p),'w') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':7,'pdf.fonttype':42,'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False})
blue='#2166AC';red='#B2182B';cts=['C0001','C0002','C0003','C0006','C0011'];ctnames={'C0001':'Neurogenic\nMGE progenitors','C0002':'RG-like\nprogenitors','C0003':'MGE intermediate\nprogenitors','C0006':'Immature inhibitory\nneurons','C0011':'Cycling\nneuroblasts'}
fig=plt.figure(figsize=(16,7.2));grid=fig.add_gridspec(2,5,left=.04,right=.98,top=.87,bottom=.16,hspace=.50,wspace=.62,height_ratios=[1.25,1],width_ratios=[1,1.25,1.1,.95,1.35])
fig.text(.04,.97,'Female mosaic tissue',fontsize=15,weight='bold');fig.text(.04,.935,'Matched WT females · WT-like and KO-like cells within HET females',fontsize=9)
# Row division conveys the independent comparisons.
bottom=grid[1,0].get_position(fig);fig.text(.04,bottom.y1+.045,'Known male controls · WT versus KO',fontsize=13,weight='bold')
def panel(ax,letter,title):
 ax.set_title(title,loc='left',fontsize=9,pad=9);ax.text(-.16,1.05,letter,transform=ax.transAxes,fontsize=12,weight='bold')
def block_title(slot,letter,title):
 b=slot.get_position(fig);fig.text(b.x0-.015,b.y1+.018,letter,fontsize=12,weight='bold');fig.text(b.x0+.02,b.y1+.018,title,fontsize=9)
# a: context, same saved coordinates.
r=[x for x in read(RUNTIME_ROOT/'results/mosaic_nature_figure_20261031/displayed_umap_cells.tsv.gz') if x['panel']=='background'];xy=np.array([[float(x['UMAP1']),float(x['UMAP2'])] for x in r]);mask=np.array([x['consensus_cluster']=='C0006' for x in r])
ax=fig.add_subplot(grid[0,0]);ax.scatter(xy[:,0],xy[:,1],s=.4,c='#ddd',linewidths=0,rasterized=True);ax.scatter(xy[mask,0],xy[mask,1],s=.8,c='#444',linewidths=0,rasterized=True);ax.set_aspect('equal');ax.axis('off');block_title(grid[0,0],'a','Inhibitory neuron population');ax.text(.02,-.02,'Dark: neurons analyzed',transform=ax.transAxes,fontsize=7)
# b: restored mixed-direction female gene list.
md=read(FEMALE/'nine_sample_metadata.tsv');order=[('WT_F','A_positive'),('HET_F','A_positive'),('HET_F','Aneg_BC1')];md=sorted(md,key=lambda x:(order.index((x['design_group'],x['gate'])),x['sample']));libs=[x['library'] for x in md];vals={x['gene_symbol']:np.array([float(x[l]) for l in libs]) for x in read(FEMALE/'nine_sample_log2CPM.tsv')}
fgenes=['Ank3','Kcnc2','Gnao1','Ptprd','Epha5','Snap47','Pcdha5','Pcdha10','Pcdhga5'];fm=np.array([[vals[g][k:k+3].mean() for k in [0,3,6]] for g in fgenes]);z=(fm-fm.mean(1)[:,None])/fm.std(1)[:,None]
def clustered(slot,mat,genes,columns,cmap,limits,letter,title,fs=6):
 inner=slot.subgridspec(1,4,width_ratios=[.40,.57,1,.07],wspace=.04);tree=fig.add_subplot(inner[0,0]);clean=np.nan_to_num(mat,nan=0);ro=dendrogram(linkage(clean,method='average'),orientation='left',no_labels=True,color_threshold=0,above_threshold_color='#333',ax=tree)['leaves'];tree.invert_yaxis();tree.axis('off')
 names=fig.add_subplot(inner[0,1]);names.set_ylim(len(genes)-.5,-.5);names.set_xlim(0,1);names.axis('off')
 for i,j in enumerate(ro):names.text(.97,i,genes[j],fontsize=fs,ha='right',va='center')
 ax=fig.add_subplot(inner[0,2]);cm=plt.get_cmap(cmap)
 import copy
 cm=copy.copy(cm);cm.set_bad('#ccc');im=ax.imshow(np.ma.masked_invalid(mat[ro]),aspect='auto',cmap=cm,vmin=limits[0],vmax=limits[1]);ax.set_yticks([]);ax.set_xticks(range(len(columns)));ax.set_xticklabels(columns,fontsize=6,rotation=40,ha='right')
 ca=fig.add_subplot(inner[0,3]);cb=fig.colorbar(im,cax=ca);cb.ax.tick_params(labelsize=6)
 block_title(slot,letter,title);return ro,cb
_,cb=clustered(grid[0,1],z,fgenes,['WT females','HET WT-like','HET KO-like'],'Purples',[-1.5,1.5],'b','Neuronal expression');cb.set_ticks([-1.5,1.5]);cb.set_ticklabels(['Low','High'])
# c: restored female effects.
fd={r['gene_symbol']:r for r in read(DEEP/'DE_C0006/within_HET_Aneg_BC1_vs_Apositive.tsv.gz')};fg=fgenes+['Lhx8','Sptbn2','Mctp1','Dlx1','Robo2','Erbb4'];fg=sorted(fg,key=lambda g:float(fd[g]['logFC']));fc=np.array([float(fd[g]['logFC']) for g in fg]);ax=fig.add_subplot(grid[0,2]);ax.barh(range(len(fg)),fc,color=[blue if v<0 else red for v in fc]);ax.set_yticks(range(len(fg)));ax.set_yticklabels(fg,fontsize=6);ax.invert_yaxis();ax.axvline(0,color='#888',lw=.6);ax.set_xlabel('KO-like / WT-like (log2FC)',fontsize=7);panel(ax,'c','Within HET females')
# d: unchanged matched exon-1 selection readout.
e=np.load(str(FEMALE/'neuronal_exon1_signal.npz'));ds=[e['exon1_CP10K'][e['group']==k] for k in range(3)];ax=fig.add_subplot(grid[0,3]);non=[k for k,d in enumerate(ds) if np.ptp(d)>0];v=ax.violinplot([ds[k] for k in non],positions=non,widths=.75,showextrema=False,points=150)
for body in v['bodies']:body.set_facecolor('#888');body.set_edgecolor('#555');body.set_alpha(.8)
for k in range(3):
 if k not in non:ax.plot([k-.28,k+.28],[0,0],color='#444',lw=2,clip_on=False)
upper=max(np.percentile(d,95) for d in ds)*1.08;ax.set_ylim(0,upper);ax.plot(range(3),[d.mean() for d in ds],'-D',color='#222',ms=3,lw=.8);ax.set_xticks(range(3));ax.set_xticklabels(['WT females','HET WT-like','HET KO-like'],rotation=45,ha='right',fontsize=6);ax.set_ylabel('Exon-1 UMI per 10,000',fontsize=7);panel(ax,'d','Pcdh19 exon 1');ax.text(.03,.90,'Matched selection',transform=ax.transAxes,fontsize=6,color='#666')
for k,d in enumerate(ds):ax.text(k,upper*1.02,'{:.0f}% above'.format(100*(d>upper).mean()),fontsize=4.5,ha='center')
# e: all supported non-gate protocadherin homophilic contacts, not only nominal hits.
expr={(r['stratum'],r['gene'],r['gate']):r for r in read(SIG/'signaling_gene_expression_support.tsv')}
def support(ct,genes,gate):
 rr=[expr.get((ct,g,gate)) for g in genes]
 if any(r is None or int(r['samples_CPM_ge1'])<2 for r in rr):return None
 return np.exp(np.mean([np.log(max(float(r['mean_CPM']),1e-8)) for r in rr]))
inter=[r for r in read(SIG/'expanded_signaling_pairs.tsv') if r['pathway']=='PCDH' and 'Pcdh19' not in r['ligand_genes']+';'+r['receptor_genes']]
net=[]
for i,sender in enumerate(cts):
 for receiver in cts[i+1:]:
  ws=[];ks=[];used=[]
  for r in inter:
   lg=r['ligand_genes'].split(';');rg=r['receptor_genes'].split(';');v=[support(sender,lg,g) for g in ['A_positive','Aneg_BC1']]+[support(receiver,rg,g) for g in ['A_positive','Aneg_BC1']]
   if any(x is None for x in v):continue
   ws.append(math.sqrt(v[0]*v[2]));ks.append(math.sqrt(v[1]*v[3]));used.append(r['interaction'])
  if not ws:continue
  w=sum(ws);k=sum(ks);net.append(dict(sender=sender,receiver=receiver,WTlike_score=w,KOlike_score=k,log2_score_change=math.log(k/w,2),interactions=';'.join(used)))
write(OUT/'displayed_PCDH_contact_network.tsv',net)
# A specific ligand/receptor pair is clearer than the nearly cancelling family total.
net_nrg=[]
for sender in cts:
 for receiver in cts:
  if sender==receiver:continue
  v=[support(sender,['Nrg1'],g) for g in ['A_positive','Aneg_BC1']]+[support(receiver,['Erbb4'],g) for g in ['A_positive','Aneg_BC1']]
  if any(x is None for x in v):continue
  w=math.sqrt(v[0]*v[2]);k=math.sqrt(v[1]*v[3]);net_nrg.append(dict(sender=sender,receiver=receiver,ligand='Nrg1',receptor='Erbb4',WTlike_score=w,KOlike_score=k,log2_score_change=math.log(k/w,2)))
write(OUT/'NRG1_ERBB4_all_supported_edges.tsv',net_nrg)
shown=sorted(net_nrg,key=lambda r:abs(r['log2_score_change']),reverse=True)[:8];write(OUT/'displayed_NRG1_ERBB4_edges.tsv',shown)
ax=fig.add_subplot(grid[0,4]);theta=np.linspace(np.pi/2,np.pi/2+2*np.pi,6)[:-1];pos={ct:np.array([np.cos(t),np.sin(t)]) for ct,t in zip(cts,theta)}
for r in sorted(shown,key=lambda r:abs(r['log2_score_change'])):
 a=pos[r['sender']];b=pos[r['receiver']];v=r['log2_score_change'];edge=FancyArrowPatch(a,b,connectionstyle='arc3,rad=.15',arrowstyle='-|>',mutation_scale=7,color=blue if v<0 else red,linewidth=.6+3*abs(v),alpha=.75,shrinkA=6,shrinkB=6,zorder=1);ax.add_patch(edge)
for ct in cts:
 q=pos[ct];ax.scatter(*q,s=85,color='#555',edgecolors='white',lw=.8,zorder=3);t=q*1.38;ax.text(*t,ctnames[ct],ha='center',va='center',fontsize=5.5)
ax.set_xlim(-1.9,1.9);ax.set_ylim(-1.85,1.85);ax.set_aspect('equal');ax.axis('off');block_title(grid[0,4],'e','NRG1 → ERBB4 potential');ax.text(.5,-.03,'KO-like vs WT-like · 8 largest changes\nExpression-based estimate',transform=ax.transAxes,ha='center',fontsize=6)
# f: independent known male genotypes, all neurons, no probe gate.
ma=read(OUT/'male_exon1_detection.tsv');ax=fig.add_subplot(grid[1,0]);means=[]
for k,group in enumerate(['WT males','KO males']):
 vv=np.array([100*float(r['fraction']) for r in ma if r['group']==group]);means.append(float(vv.mean()));ax.bar(k,vv.mean(),width=.55,color='#999' if k==0 else blue);ax.errorbar(k,vv.mean(),yerr=vv.std(ddof=1)/np.sqrt(len(vv)),fmt='none',color='#333',capsize=3)
ax.set_xticks([0,1]);ax.set_xticklabels(['WT males','KO males']);ax.set_ylabel('Exon-1-positive neurons (%)',fontsize=7);panel(ax,'f','Pcdh19 exon-1 detection')
# g,h: nominal neuronal leads in the same neuronal population; full test retained.
male={st:{r['gene_symbol']:r for r in read(OUT/('male_'+st+'_DE.tsv'))} for st in cts};m6=male['C0006'];neural=set('Eps8l1 Serpine2 Syt9 Ecel1 Lrfn5 Grm7 Cntn1 Magi3 Sphkap Nxph1 Ncam1 Spock2 Htr1f Hcn1 Sstr1 Slc6a1 Slc32a1 Gad1 Gad2 Lhx6 Dlx1 Dlx2 Nrxn1 Nlgn1 Nlgn2 Nlgn3 Kcnc2 Ank3 Gria2 Gria3 Grin1 Grin2b Pcdha5 Ptprd Sema3a Robo1 Robo2 Dcc Reln Dcx Tubb3 Map2 Stmn2 Cntnap2'.split());pg=read(DEEP/'marker_program_gene_coverage.tsv');neural.update(r['gene'] for r in pg if r['program'] in ['neural_core','adhesion_synapse','neurite_guidance','early_neuronal_maturation','later_neuronal_maturation','MGE_interneuron'])
selected=sorted([g for g in m6 if g in neural and abs(float(m6[g]['logFC']))>=.25],key=lambda g:float(m6[g]['PValue']))[:18];mat=np.array([[float(male[ct][g]['logFC']) if g in male[ct] else np.nan for ct in cts] for g in selected]);_,cb=clustered(grid[1,1:3],mat,selected,['Neurogenic\nprogenitors','RG-like\nprogenitors','MGE IPCs','Immature\nneurons','Cycling\nneuroblasts'],'RdBu_r',[-1.5,1.5],'g','Male neuronal candidates across cell types',fs=6);cb.set_label('KO / WT log2FC',fontsize=6)
ax=fig.add_subplot(grid[1,3:]);leaders=[g for g in selected if abs(float(m6[g]['logFC']))>=.5][:9];leaders=sorted(leaders,key=lambda g:float(m6[g]['logFC']));change=np.array([100*(2**float(m6[g]['logFC'])-1) for g in leaders]);ax.barh(range(len(leaders)),change,color=[blue if v<0 else red for v in change]);ax.set_yticks(range(len(leaders)));ax.set_yticklabels(leaders,fontsize=7);ax.invert_yaxis();ax.axvline(0,color='#888',lw=.6);ax.set_xlabel('KO change relative to WT (%)',fontsize=8);panel(ax,'h','Exploratory neuronal changes')
# A compact key belongs to the figure, not a methodological footer.
fig.text(.70,.96,'BLUE: lower in KO-like / KO',fontsize=8,color=blue,ha='right');fig.text(.98,.96,'RED: higher in KO-like / KO',fontsize=8,color=red,ha='right')
for ext in ['pdf','png','svg']:fig.savefig(str(OUT/('female_mosaic_and_male_controls.'+ext)),dpi=350,facecolor='white')
write(OUT/'displayed_male_neuronal_candidates.tsv',[dict(stratum=ct,**male[ct][g]) for ct in cts for g in selected if g in male[ct]])
(OUT/'figure_validation.json').write_text(json.dumps({'female_WT_selection':'same exon1-positive selection as HET WT-like','female_neuronal_panels':'restored mixed-direction list; same provisional neuronal cluster','male_samples_per_group':3,'male_genotype_selection':'known registered groups; no probe selection','male_exon1_means_percent':means,'male_neuronal_FDR05':sum(float(r['FDR'])<.05 for r in m6.values()),'male_leads_are_nominal':True,'male_selected_genes':selected,'network':'NRG1-ERBB4 expression potential; top 8 absolute differences of all supported directed routes; PCDH family calculation retained separately','network_edges':len(shown),'network_full_edges':len(net_nrg),'network_all_component_FDR05':False,'network_not_measured_communication':True,'female_E_F_expression_panels_removed':True,'individual_samples_displayed':False},indent=2))
(OUT/'FIGURE_LEGEND.txt').write_text('Female top row: matched exon-1-positive WT-female and HET WT-like neurons, versus A-negative downstream-positive HET KO-like neurons. A is saved UMAP context. B is gene-wise standardized equal-sample mean logCPM. C is paired HET edgeR log2FC with restored mixed-direction genes. D is defining exon-1 probe signal, not independent genotype validation; violins use actual single-cell UMI/10000 with marked above-range percentages. E is an exploratory directed NRG1-ERBB4 expression-potential network among provisional neural populations. Every component has CPM>=1 in at least two of three samples in each gate. Edge score is sqrt(ligand CPM * receptor CPM); color is log2 KO-like/WT-like score, width is 0.6+3*absolute log2 score change. The eight largest absolute changes of all 20 supported directed non-self routes are displayed. The initial PCDH family network is retained separately; its aggregate changes were modest. This is not a CellChat probability, a measured interaction, or evidence of physical contact. The initial nominal-component pathway screen is retained separately. Male bottom row: three known WT-male and three known KO-male samples, all QC-passing cells without probe selection. F gives equal-sample mean exon-1 detection fractions +/- SEM in the same neuronal cluster. G shows exploratory neuronal gene log2FC across five broad consensus populations; gray means not tested. H shows fitted percent changes for the nine best nominal-P neuronal candidates with absolute log2FC>=0.5 in the neuronal cluster. Gene lists are exploratory; no male gene passes per-stratum gene-wide FDR 0.05. Female LR component changes also lack corrected neural support. Samples, not cells, are differential-expression units. Original embryo/litter dependence and covariates remain unverified. Full source tables and selection methods are retained.\n')
print('COMBINED_FIGURE_COMPLETE',flush=True)
