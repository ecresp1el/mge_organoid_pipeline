"""Focused sample-level Eomes/Tbr2, RG and IPC composition comparisons."""
from pathlib import Path
from pcdh19_exploratory_paths import RUNTIME_ROOT
import csv,json,itertools,hashlib
import h5py,numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1]
OUT=RUNTIME_ROOT/'results/Eomes_RG_proportions_20261101';OUT.mkdir(exist_ok=True,parents=True)
SRC=RUNTIME_ROOT/'results/mosaic_deep_dive_20261031/per_cell_marker_diagnostics.h5'
def strings(x):return np.array([v.decode() if isinstance(v,bytes) else str(v) for v in x[:]])
with h5py.File(str(SRC),'r') as f:
 sample=strings(f['sample']);cluster=strings(f['cluster']);genes=list(strings(f['marker_gene']));e=f['marker_raw_counts'][:,genes.index('Eomes')]
key=list(csv.DictReader((ROOT/'config/sample_key.csv').open()));groups=['WT_M','KO_M','WT_F','HET_F'];names=['WT males','KO males','WT females','HET females']
progen=np.isin(cluster,['C0001','C0002','C0003','C0008']);ipc=np.isin(cluster,['C0003','C0008']);rg=cluster=='C0002'
masks={'Eomes+ (>=1 UMI)':(e>=1,np.ones(len(e),bool)),'Eomes+ (>=2 UMI)':(e>=2,np.ones(len(e),bool)),'Radial glia':(rg,np.ones(len(e),bool)),'MGE IPCs':(ipc,np.ones(len(e),bool)),'Eomes+ within progenitors':(e>=1,progen),'Eomes+ within RG':(e>=1,rg),'Eomes+ within IPCs':(e>=1,ipc)}
rows=[]
for label,(hit,denom) in masks.items():
 for r in key:
  m=(sample==r['technical_sample_id'])&denom;n=int((hit&m).sum());total=int(m.sum());rows.append(dict(population=label,sample=r['technical_sample_id'],group=r['design_group'],count=n,total=total,percent=100*n/total if total else float('nan')))
def write(name,rr):
 with (OUT/name).open('w') as f:
  w=csv.DictWriter(f,fieldnames=list(rr[0]),delimiter='\t');w.writeheader();w.writerows(rr)
write('per_sample_proportions.tsv',rows)
summary=[];contrasts=[]
for label in masks:
 vals={g:np.array([r['percent'] for r in rows if r['population']==label and r['group']==g]) for g in groups}
 for g in groups:summary.append(dict(population=label,group=g,mean_percent=float(vals[g].mean()),min_percent=float(vals[g].min()),max_percent=float(vals[g].max()),n_samples=3))
 for ga,gb in itertools.combinations(groups,2):
  a,b=vals[ga],vals[gb];diff=b.mean()-a.mean();joined=np.r_[a,b];null=[]
  for ix in itertools.combinations(range(6),3):
   ii=np.array(ix);jj=np.array([i for i in range(6) if i not in ix]);null.append(joined[ii].mean()-joined[jj].mean())
  loo=[np.delete(b,j).mean()-np.delete(a,i).mean() for i in range(3) for j in range(3)]
  contrasts.append(dict(population=label,reference=ga,comparison=gb,difference_percentage_points=float(diff),exact_two_sided_p=float(np.mean(np.abs(null)>=abs(diff)-1e-10)),complete_sample_separation=bool(a.max()<b.min() or b.max()<a.min()),leave_one_per_group_direction_consistent=bool(all(v*diff>0 for v in loo)),sex_matched=bool(set([ga,gb]) in [set(['WT_M','KO_M']),set(['WT_F','HET_F'])])))
write('group_summary.tsv',summary);write('all_contrasts.tsv',contrasts)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False})
fig,axes=plt.subplots(1,3,figsize=(9.0,3.2));colors=['#437C90','#D18A43','#7F858C','#9A779C']
for ax,label,title in zip(axes,['Eomes+ (>=1 UMI)','Radial glia','MGE IPCs'],['Eomes / Tbr2-positive','Radial glia','MGE IPCs']):
 ss=[r for g in groups for r in summary if r['population']==label and r['group']==g];means=[r['mean_percent'] for r in ss];lo=[r['mean_percent']-r['min_percent'] for r in ss];hi=[r['max_percent']-r['mean_percent'] for r in ss]
 ax.bar(range(4),means,color=colors,width=.65);ax.errorbar(range(4),means,yerr=[lo,hi],fmt='none',ecolor='#333333',capsize=3,lw=1)
 ax.set_xticks(range(4));ax.set_xticklabels([n.replace(' ','\n') for n in names]);ax.set_title(title,fontweight='bold');ax.set_ylabel('% of all QC cells');ax.set_ylim(bottom=0)
fig.suptitle('Mouse MGE · E14.5 · 3 samples per group',fontweight='bold');fig.tight_layout(rect=[0,0,1,.94])
for ext in ['png','pdf','svg']:fig.savefig(str(OUT/('Eomes_RG_IPC_proportions.'+ext)),dpi=300)
(OUT/'METHODS.json').write_text(json.dumps(dict(source=str(SRC),sample_key=str(ROOT/'config/sample_key.csv'),n_cells=len(sample),probe_filter='none; all QC cells, including Pcdh19-negative cells',Eomes='Tbr2 alias; raw UMI >=1, >=2 sensitivity analysis',RG='existing provisional C0002 annotation',IPC='existing C0003+C0008 annotation; not equated to Eomes positivity',averaging='equal-weight mean of 3 sample percentages per group',plot_errorbars='sample minimum to maximum',statistics='all 20 exact assignments, two-sided; unadjusted exploratory P values; minimum possible P=0.1 with 3 vs 3; no cell-level replication',cross_sex='KO male vs female HET/WT confounds genotype with sex'),indent=2))
print(json.dumps(summary,indent=2));print('KO_VS_WT_MALE',json.dumps([r for r in contrasts if r['reference']=='WT_M' and r['comparison']=='KO_M'],indent=2))
