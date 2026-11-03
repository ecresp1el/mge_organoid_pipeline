"""Focused sample-level Eomes/Tbr2, RG and IPC composition comparisons."""
from pathlib import Path
from pcdh19_exploratory_paths import RUNTIME_ROOT
import csv,json,itertools,hashlib,gzip,shutil
import h5py,numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1]
OUT=RUNTIME_ROOT/'results/nested_proliferation_20261101';OUT.mkdir(exist_ok=True,parents=True)
SRC=RUNTIME_ROOT/'results/mosaic_deep_dive_20261031/per_cell_marker_diagnostics.h5'
def strings(x):return np.array([v.decode() if isinstance(v,bytes) else str(v) for v in x[:]])
with h5py.File(str(SRC),'r') as f:
 sample=strings(f['sample']);cluster=strings(f['cluster']);genes=list(strings(f['marker_gene']));e=f['marker_raw_counts'][:,genes.index('Eomes')];gates=list(strings(f['gate_names']));membership=f['gate_membership'][:];ids=strings(f['cell_id'])
key=list(csv.DictReader((ROOT/'config/sample_key.csv').open()));groups=['WT_F','HET_WT','HET_KO','WT_M_positive','KO_M'];names=['WT female','HET WT cells','HET KO cells','WT males Pcdh19+','KO males']
lookup={v:i for i,v in enumerate(ids)};wt_positive=np.zeros(len(sample),bool)
for n in [1,2,3]:
 sn='15662-JZ-'+str(n);source=RUNTIME_ROOT/'results/pcdh19_probe_audit/per_sample'/sn/'pcdh19_probe_patterns.tsv'
 with source.open() as f:
  for r in csv.DictReader(f,delimiter='\t'):
   i=lookup.get(sn+'_'+r['barcode'])
   if i is not None:wt_positive[i]=int(r['Pcdh19_total_UMI'])>0
Apositive=membership[:,gates.index('A_positive')].astype(bool);AnegBC=membership[:,gates.index('Aneg_BC1')].astype(bool)
definitions={'WT_F':('WT_F',Apositive),'HET_WT':('HET_F',Apositive),'HET_KO':('HET_F',AnegBC),'WT_M_positive':('WT_M',wt_positive),'KO_M':('KO_M',np.ones(len(sample),bool))}
SCORE=RUNTIME_ROOT/'results/hicat/02_hierarchy_validation/hicat_validation_20260908_221625_d5d06a5d/outputs/full_data_projection/plot_program_scores.tsv.gz'
parents={'All progenitors / cycling neuroblasts':['C0001','C0002','C0003','C0007','C0008','C0011','C0012','C0024'],'RG':['C0002'],'IPCs':['C0003','C0008'],'IPC compartment 1':['C0003'],'IPC compartment 2':['C0008'],'Neurogenic progenitors':['C0001','C0007','C0012','C0024'],'Cycling neuroblasts':['C0011']}
plan={'parents':parents,'groups_CE':{k:v[0] for k,v in definitions.items()},'groups_genotype_sex':['WT_M','KO_M','WT_F','HET_F'],'cycle_definitions':['S-phase score >0.1','G2/M score >0.1','Either score >0.1','Either score >0.8','S-dominant with max score >0.1','G2/M-dominant with max score >0.1','Top quartile of max(S,G2/M) within each sample and parent, threshold fitted on all cells independent of probe gates'],'other_endpoints':['Eomes >=1 UMI','Ascl1 >=1 UMI','Mean S-phase score','Mean G2/M score','Child cluster shares within broad progenitor compartment'],'selection':'Report all endpoints; rank proportions with >=30 denominator cells/sample, >=5 numerator cells/sample/group; primary within-HET and within-sex contrasts; no threshold optimized against effect','interpretation':'saved descriptive marker-program scores, not experimentally measured proliferation rate or validated cell-cycle phase calls'}
(OUT/'EXPLORATORY_PLAN.json').write_text(json.dumps(plan,indent=2))
sc=np.full((len(sample),2),np.nan)
with gzip.open(str(SCORE),'rt') as f:
 for r in csv.DictReader(f,delimiter='\t'):
  i=lookup.get(r['cell_id'])
  if i is not None:sc[i]=[float(r['S_phase']),float(r['G2M'])]
assert np.isfinite(sc).all()
np.savez_compressed(str(OUT/'saved_cycle_scores.npz'),cell_id=ids,S_phase=sc[:,0],G2M=sc[:,1])
# Independent raw-marker sensitivity analysis; no phase calling.
RAW=RUNTIME_ROOT/'results/primary_processing/06_technical_sample_batch_diagnostics/06_technical_sample_batch_diagnostics_20260903_114110_6087eeb/objects/pcdh19_step06_unintegrated_diagnostics.h5ad'
cycle_genes=['Mki67','Top2a','Pcna','Mcm2','Ube2c','Birc5'];cache=OUT/'raw_cycle_markers.npz'
if cache.exists():
 raw_cache=np.load(str(cache));assert np.array_equal(raw_cache['cell_id'],ids);cycle_counts=raw_cache['counts']
else:
 from scipy import sparse
 cycle_counts=np.zeros((len(ids),len(cycle_genes)),dtype=np.int32)
 with h5py.File(str(RAW),'r') as f:
  assert np.array_equal(strings(f['obs/cell_id']),ids)
  v=f['var/gene_symbol']
  if isinstance(v,h5py.Group):
   cats=strings(v['categories']);symbols=cats[v['codes'][:]]
  else:symbols=strings(v)
  gi=[list(symbols).index(g) for g in cycle_genes];x=f['X'];ptr=x['indptr'][:]
  for start in range(0,len(ids),4000):
   end=min(start+4000,len(ids));lo,hi=int(ptr[start]),int(ptr[end]);mat=sparse.csr_matrix((x['data'][lo:hi],x['indices'][lo:hi],ptr[start:end+1]-lo),shape=(end-start,len(symbols)));cycle_counts[start:end]=mat[:,gi].toarray()
   if start%100000==0:print('RAW_MARKERS',start,flush=True)
 np.savez_compressed(str(cache),cell_id=ids,genes=np.array(cycle_genes),counts=cycle_counts)
a=sc[:,0];b=sc[:,1];mx=np.maximum(a,b)
ascl=None
with h5py.File(str(SRC),'r') as f:ascl=f['marker_raw_counts'][:,genes.index('Ascl1')]
allcells=np.ones(len(sample),bool)
modes={'CE_cell_groups':definitions,'genotype_sex':{g:(g,allcells) for g in ['WT_M','KO_M','WT_F','HET_F']}}
labels={'WT_F':'WT female','HET_WT':'HET WT cells','HET_KO':'HET KO cells','WT_M_positive':'WT males Pcdh19+','KO_M':'KO males','WT_M':'WT males','HET_F':'HET females'}
rows=[];metric_order=[]
for parent,cc in parents.items():
 pm=np.isin(cluster,cc);quartile=np.zeros(len(sample),bool)
 for r in key:
  ix=(sample==r['technical_sample_id'])&pm
  if ix.any():quartile[ix]=mx[ix]>np.percentile(mx[ix],75)
 endpoints={'S-phase score >0.1':a>.1,'G2/M score >0.1':b>.1,'Either cycle score >0.1':mx>.1,'Either cycle score >0.8':mx>.8,'S-dominant':(mx>.1)&(a>=b),'G2/M-dominant':(mx>.1)&(b>a),'High-cycle top quartile':quartile,'Eomes-positive':e>=1,'Ascl1-positive':ascl>=1,'Mean S-phase score':a,'Mean G2/M score':b}
 for i,gene in enumerate(cycle_genes):
  endpoints[gene+' detected']=cycle_counts[:,i]>=1
  if gene in ['Mki67','Top2a']:endpoints[gene+' >=2 UMI']=cycle_counts[:,i]>=2
 endpoints['At least 2 cycle markers']=(cycle_counts>0).sum(1)>=2
 if parent=='All progenitors / cycling neuroblasts':
  for child in ['RG','IPCs','Neurogenic progenitors','Cycling neuroblasts']:endpoints['Share: '+child]=np.isin(cluster,parents[child])
 for endpoint,v in endpoints.items():
  kind='mean_score' if endpoint.startswith('Mean ') else 'percent';metric_order.append((parent,endpoint))
  for mode,defs in modes.items():
   for group,(design,gate) in defs.items():
    for r in key:
     if r['design_group']!=design:continue
     m=(sample==r['technical_sample_id'])&pm&gate;den=int(m.sum());num=int(v[m].sum()) if kind=='percent' else -1
     val=(100*num/den if kind=='percent' else float(v[m].mean())) if den else float('nan')
     rows.append(dict(mode=mode,parent=parent,endpoint=endpoint,kind=kind,group=group,sample=r['technical_sample_id'],numerator=num,denominator=den,value=val))
def write(name,rr):
 with (OUT/name).open('w') as f:
  w=csv.DictWriter(f,fieldnames=list(rr[0]),delimiter='\t');w.writeheader();w.writerows(rr)
write('all_nested_sample_values.tsv',rows)
summary=[];tests=[]
for mode,defs in modes.items():
 for parent,endpoint in metric_order:
  rr=[r for r in rows if r['mode']==mode and r['parent']==parent and r['endpoint']==endpoint];kind=rr[0]['kind']
  vals={g:np.array([r['value'] for r in rr if r['group']==g]) for g in defs}
  for g,v in vals.items():summary.append(dict(mode=mode,parent=parent,endpoint=endpoint,kind=kind,group=g,mean=float(v.mean()),minimum=float(v.min()),maximum=float(v.max())))
  for ga,gb in itertools.combinations(defs,2):
   aa,bb=vals[ga],vals[gb]
   if not (np.isfinite(aa).all() and np.isfinite(bb).all()):continue
   effect=float(bb.mean()-aa.mean());paired=(ga,gb)==('HET_WT','HET_KO');d=bb-aa
   if paired:
    null=[float((d*np.array(signs)).mean()) for signs in itertools.product([-1,1],repeat=3)];loo=[float(np.delete(d,i).mean()) for i in range(3)];agree=int((d*effect>0).sum())
   else:
    joined=np.r_[aa,bb];null=[]
    for ix in itertools.combinations(range(6),3):null.append(float(joined[list(ix)].mean()-joined[[j for j in range(6) if j not in ix]].mean()))
    loo=[float(np.delete(bb,j).mean()-np.delete(aa,i).mean()) for i in range(3) for j in range(3)];agree=int(sum((y-x)*effect>0 for x in aa for y in bb))
   rrr=[r for r in rr if r['group'] in [ga,gb]];mind=min(r['denominator'] for r in rrr);minn=min(r['numerator'] for r in rrr);stable=effect!=0 and all(x*effect>0 for x in loo)
   primary=paired or set([ga,gb]) in [set(['WT_M','KO_M']),set(['WT_M_positive','KO_M']),set(['WT_F','HET_F'])]
   tests.append(dict(mode=mode,parent=parent,endpoint=endpoint,kind=kind,reference=ga,comparison=gb,reference_mean=float(aa.mean()),comparison_mean=float(bb.mean()),difference=effect,exact_p=sum(abs(x)>=abs(effect)-1e-10 for x in null)/len(null),paired=paired,primary_contrast=primary,agreeing_sample_differences=agree,leave_out_direction_stable=stable,minimum_abs_leave_out_effect=min(abs(x) for x in loo),minimum_denominator=mind,minimum_numerator=minn,eligible=bool(kind=='percent' and mind>=30 and minn>=5),all_pairs_separated=bool(aa.max()<bb.min() or bb.max()<aa.min())))
# BH across all tested contrasts, including overlapping endpoints, separately by analysis mode.
for mode in modes:
 ix=sorted([i for i,r in enumerate(tests) if r['mode']==mode],key=lambda i:tests[i]['exact_p']);q=1.
 for rank in range(len(ix)-1,-1,-1):
  i=ix[rank];q=min(q,tests[i]['exact_p']*len(ix)/(rank+1));tests[i]['BH_q']=q
write('all_nested_contrasts.tsv',tests);write('group_summary.tsv',summary)
ranked=sorted(tests,key=lambda r:(not r['primary_contrast'],not r['eligible'],not r['leave_out_direction_stable'],-r['minimum_abs_leave_out_effect'] if r['kind']=='percent' else 0))
write('ranked_exploratory_contrasts.tsv',ranked)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False})
selected={}
for mode,defs in modes.items():
 leads=[r for r in ranked if r['mode']==mode and r['primary_contrast'] and r['eligible'] and r['leave_out_direction_stable'] and r['endpoint']!='High-cycle top quartile' and (mode!='CE_cell_groups' or r['paired'])]
 chosen=[];used=set()
 for r in leads:
  identity=(r['parent'],r['endpoint'])
  if identity in used:continue
  chosen.append(r);used.add(identity)
  if len(chosen)==6:break
 selected[mode]=chosen
 if not chosen:continue
 fig,axs=plt.subplots(2,3,figsize=(10,6.3));group_list=list(defs);colors=['#AEB5BE','#5F997E','#B78366','#437C90','#D18A43'] if mode=='CE_cell_groups' else ['#437C90','#D18A43','#AEB5BE','#9A779C']
 for ax,r in zip(axs.flat,chosen):
  ss=[next(x for x in summary if x['mode']==mode and x['parent']==r['parent'] and x['endpoint']==r['endpoint'] and x['group']==g) for g in group_list];means=[x['mean'] for x in ss]
  ax.bar(range(len(defs)),means,color=colors,width=.65);ax.errorbar(range(len(defs)),means,yerr=[[x['mean']-x['minimum'] for x in ss],[x['maximum']-x['mean'] for x in ss]],fmt='none',color='#333333',capsize=2,lw=.7)
  parent=r['parent'].replace('All progenitors / cycling neuroblasts','Progenitors / cycling neuroblasts');ax.set_title(parent+'\n'+r['endpoint'],fontsize=8,weight='bold');ax.set_ylabel('% within parent population');ax.set_xticks(range(len(defs)));ax.set_xticklabels([labels[g].replace(' ','\n',1) for g in group_list],fontsize=6);ax.set_ylim(bottom=0)
 for ax in list(axs.flat)[len(chosen):]:ax.axis('off')
 fig.suptitle('Nested proliferation · '+('C/E cell groups' if mode=='CE_cell_groups' else 'genotype and sex'),weight='bold');fig.tight_layout(rect=[0,.025,1,.95]);fig.text(.04,.01,'Exploratory selected effects; bars = sample means, whiskers = sample range; n = 3 per group.',fontsize=7)
 for ext in ['png','pdf','svg']:fig.savefig(str(OUT/('nested_'+mode+'.'+ext)),dpi=300)
(OUT/'SELECTED_LEADS.json').write_text(json.dumps(selected,indent=2));shutil.copy2(__file__,str(OUT/Path(__file__).name))
print('NESTED_COMPLETE',len(tests),flush=True)
for mode in modes:
 print(mode,json.dumps([r for r in ranked if r['mode']==mode and r['primary_contrast'] and r['eligible'] and r['leave_out_direction_stable'] and r['endpoint']!='High-cycle top quartile' and (mode!='CE_cell_groups' or r['paired'])][:12],indent=2),flush=True)
