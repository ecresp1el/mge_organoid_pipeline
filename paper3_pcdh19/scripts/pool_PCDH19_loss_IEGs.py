"""Descriptive pooled IEG hypothesis check, preserving sex and embryo provenance."""
from pathlib import Path
from pcdh19_exploratory_paths import RUNTIME_ROOT
import csv,json,itertools,shutil
import h5py,numpy as np
from scipy import sparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1];OUT=RUNTIME_ROOT/'results/pooled_PCDH19_IEGs_20261101';OUT.mkdir(exist_ok=True,parents=True)
BASE=RUNTIME_ROOT/'results';DIAG=BASE/'mosaic_deep_dive_20261031/per_cell_marker_diagnostics.h5';RAW=BASE/'primary_processing/06_technical_sample_batch_diagnostics/06_technical_sample_batch_diagnostics_20260903_114110_6087eeb/objects/pcdh19_step06_unintegrated_diagnostics.h5ad'
def strings(d):return np.array([x.decode() if isinstance(x,bytes) else str(x) for x in d[:]])
with h5py.File(str(DIAG),'r') as f:
 ids=strings(f['cell_id']);sample=strings(f['sample']);cluster=strings(f['cluster']);total=f['total_counts'][:];gn=list(strings(f['gate_names']));gates=f['gate_membership'][:]
key=list(csv.DictReader((ROOT/'config/sample_key.csv').open()));design={r['technical_sample_id']:r['design_group'] for r in key};dg=np.array([design[s] for s in sample]);neuron=cluster=='C0006';positive=gates[:,gn.index('informative')].astype(bool);apos=gates[:,gn.index('A_positive')].astype(bool);khet=gates[:,gn.index('Aneg_BC1')].astype(bool)
source_masks={'WT_F_positive':(dg=='WT_F')&positive,'WT_M_positive':(dg=='WT_M')&positive,'HET_WT':(dg=='HET_F')&apos,'HET_KO':(dg=='HET_F')&khet,'KO_M':dg=='KO_M'}
variants={'WT_all_positive':((np.isin(dg,['WT_F','WT_M'])&positive)|((dg=='HET_F')&apos)), 'pure_WT_positive_only':np.isin(dg,['WT_F','WT_M'])&positive,'exact_CE_WT':((dg=='WT_F')&apos)|((dg=='WT_M')&positive)|((dg=='HET_F')&apos),'all_WT_exon1_positive':np.isin(dg,['WT_F','WT_M','HET_F'])&apos}
ko=((dg=='KO_M')|((dg=='HET_F')&khet));genes=['Egr1','Npas4','Nr4a1'];selected=np.flatnonzero(neuron);cache=OUT/'C0006_IEG_raw_counts.npz'
(OUT/'PLAN.json').write_text(json.dumps({'genes':genes,'neuronal_cluster':'C0006, same provisional immature inhibitory neuronal cluster','primary_WT':'pure WT male/female with any Pcdh19 probe detected, plus A-positive HET WT cells','KO':'all KO male neurons plus HET A-negative/B+C>=1 cells','sensitivities':list(variants),'primary_expression':'mean per-cell raw IEG UMI / total cell UMI *10000; cell-weighted pooled mean','secondary_expression':'equal-embryo library means; raw sum counts / summed total UMI *1e6; fraction detected','no_pseudoreplication':'pooled cell averages are descriptive; do not treat individual cells as independent embryos','literature':'PMID 35613587, https://pmc.ncbi.nlm.nih.gov/articles/PMC9152703/'},indent=2))
if cache.exists():
 c=np.load(str(cache));assert np.array_equal(c['cell_id'],ids[selected]);counts=c['counts']
else:
 counts=np.zeros((len(selected),len(genes)),dtype=np.int32)
 with h5py.File(str(RAW),'r') as f:
  assert np.array_equal(strings(f['obs/cell_id']),ids);v=f['var/gene_symbol']
  if isinstance(v,h5py.Group):cats=strings(v['categories']);symbols=cats[v['codes'][:]]
  else:symbols=strings(v)
  gi=[list(symbols).index(g) for g in genes];x=f['X'];ptr=x['indptr'][:]
  for start in range(0,len(ids),4000):
   end=min(start+4000,len(ids));ix=np.flatnonzero((selected>=start)&(selected<end))
   if not len(ix):continue
   lo,hi=int(ptr[start]),int(ptr[end]);mat=sparse.csr_matrix((x['data'][lo:hi],x['indices'][lo:hi],ptr[start:end+1]-lo),shape=(end-start,len(symbols)));counts[ix]=mat[selected[ix]-start][:,gi].toarray()
   if start%100000==0:print('IEG_EXTRACTION',start,flush=True)
 np.savez_compressed(str(cache),cell_id=ids[selected],genes=np.array(genes),counts=counts)
ss=sample[selected];dd=total[selected];cp=counts/dd[:,None]*10000.;rows=[]
for source,mask in source_masks.items():
 for sn in sorted(set(ss[mask[selected]])):
  m=mask[selected]&(ss==sn)
  for j,gene in enumerate(genes):rows.append(dict(source=source,sample=sn,gene=gene,cells=int(m.sum()),raw_counts=int(counts[m,j].sum()),total_UMI=float(dd[m].sum()),mean_CP10K=float(cp[m,j].mean()),mean_log1p_CP10K=float(np.log1p(cp[m,j]).mean()),detection_percent=float((counts[m,j]>0).mean()*100),pseudobulk_CPM=float(counts[m,j].sum()/dd[m].sum()*1e6)))
def write(name,rr):
 with (OUT/name).open('w') as f:
  w=csv.DictWriter(f,fieldnames=list(rr[0]),delimiter='\t');w.writeheader();w.writerows(rr)
write('per_sample_source_IEGs.tsv',rows)
pooled=[]
for variant,wt in variants.items():
 for j,gene in enumerate(genes):
  for label,mask in [('WT',wt),('KO',ko)]:
   m=mask[selected];values=[]
   for sn in sorted(set(ss[m])):values.append(float(cp[m&(ss==sn),j].mean()))
   pooled.append(dict(variant=variant,gene=gene,group=label,cells=int(m.sum()),embryo_libraries=len(values),pooled_mean_CP10K=float(cp[m,j].mean()),pooled_mean_log1p_CP10K=float(np.log1p(cp[m,j]).mean()),equal_embryo_mean_CP10K=float(np.mean(values)),pooled_pseudobulk_CPM=float(counts[m,j].sum()/dd[m].sum()*1e6),detection_percent=float((counts[m,j]>0).mean()*100)))
write('pooled_IEG_summary.tsv',pooled)
contrasts=[]
for variant in variants:
 for gene in genes:
  rr=[next(r for r in pooled if r['variant']==variant and r['gene']==gene and r['group']==g) for g in ['WT','KO']]
  contrasts.append(dict(variant=variant,gene=gene,WT_cells=rr[0]['cells'],KO_cells=rr[1]['cells'],cell_weighted_KO_over_WT=rr[1]['pooled_mean_CP10K']/rr[0]['pooled_mean_CP10K'] if rr[0]['pooled_mean_CP10K'] else None,equal_embryo_KO_over_WT=rr[1]['equal_embryo_mean_CP10K']/rr[0]['equal_embryo_mean_CP10K'] if rr[0]['equal_embryo_mean_CP10K'] else None,pseudobulk_KO_over_WT=rr[1]['pooled_pseudobulk_CPM']/rr[0]['pooled_pseudobulk_CPM'] if rr[0]['pooled_pseudobulk_CPM'] else None))
write('all_pooling_sensitivities.tsv',contrasts)
source_summary=[]
for source in source_masks:
 for gene in genes:
  rr=[r for r in rows if r['source']==source and r['gene']==gene];source_summary.append(dict(source=source,gene=gene,cells=sum(r['cells'] for r in rr),equal_embryo_mean_CP10K=float(np.mean([r['mean_CP10K'] for r in rr])),minimum=float(min(r['mean_CP10K'] for r in rr)),maximum=float(max(r['mean_CP10K'] for r in rr))))
write('source_group_summary.tsv',source_summary)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'pdf.fonttype':42,'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False})
fig,axs=plt.subplots(2,3,figsize=(9,5.8));source_order=list(source_masks);source_labels=['WT female\nPcdh19+','WT male\nPcdh19+','HET\nWT cells','HET\nKO cells','KO\nmale'];colors=['#A4ABB4','#437C90','#5F997E','#B78366','#D18A43']
for j,gene in enumerate(genes):
 ax=axs[0,j];rr=[next(r for r in pooled if r['variant']=='WT_all_positive' and r['gene']==gene and r['group']==g) for g in ['WT','KO']];val=[100*r['pooled_mean_CP10K'] for r in rr]
 ax.plot([0,1],val,color='#BEC3C9',lw=1.3)
 for i,r in enumerate(rr):ax.scatter(i,val[i],s=55,color=['#5F997E','#B78366'][i]);ax.text(i,val[i],'%0.3g'%val[i],ha='center',va='bottom',fontsize=8)
 ax.set_title(gene,style='italic',weight='bold');ax.set_xticks([0,1]);ax.set_xticklabels(['Pooled WT\nPcdh19+','Pooled KO']);ax.set_xlim(-.3,1.3);ax.set_ylim(0,max(val)*1.3 if max(val)>0 else 1)
 if j==0:ax.set_ylabel('Mean per-cell expression (CPM)')
 ax=axs[1,j]
 for i,source in enumerate(source_order):
  r=next(r for r in source_summary if r['source']==source and r['gene']==gene);m=100*r['equal_embryo_mean_CP10K'];ax.errorbar(i,m,yerr=[[m-100*r['minimum']],[100*r['maximum']-m]],fmt='o',color=colors[i],capsize=2,markersize=5,lw=.8)
 ax.set_xticks(range(5));ax.set_xticklabels(source_labels,fontsize=6.5);ax.set_xlim(-.5,4.5);ax.set_ylim(bottom=0)
 if j==0:ax.set_ylabel('Equal-embryo mean expression (CPM)')
fig.suptitle('Same immature inhibitory neurons · pooled PCDH19-loss hypothesis',weight='bold',fontsize=11)
fig.tight_layout(rect=[0,.035,1,.95]);fig.text(.08,.012,'Top: cells pooled. Bottom: three-embryo means and sample ranges by source; shared HET embryos remain identified.',fontsize=7)
for ext in ['png','pdf','svg']:fig.savefig(str(OUT/('pooled_IEG_check.'+ext)),dpi=350,facecolor='white')
shutil.copy2(__file__,str(OUT/Path(__file__).name));print('POOLING_RESULTS',json.dumps(contrasts,indent=2),flush=True)
