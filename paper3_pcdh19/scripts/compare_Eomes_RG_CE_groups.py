"""Focused sample-level Eomes/Tbr2, RG and IPC composition comparisons."""
from pathlib import Path
from pcdh19_exploratory_paths import RUNTIME_ROOT
import csv,json,itertools,hashlib
import h5py,numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1]
OUT=RUNTIME_ROOT/'results/Eomes_RG_CE_groups_20261101';OUT.mkdir(exist_ok=True,parents=True)
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
progen=np.isin(cluster,['C0001','C0002','C0003','C0008']);ipc=np.isin(cluster,['C0003','C0008']);rg=cluster=='C0002'
masks={'Eomes+ (>=1 UMI)':(e>=1,np.ones(len(e),bool)),'Eomes+ (>=2 UMI)':(e>=2,np.ones(len(e),bool)),'Radial glia':(rg,np.ones(len(e),bool)),'MGE IPCs':(ipc,np.ones(len(e),bool)),'Eomes+ within progenitors':(e>=1,progen),'Eomes+ within RG':(e>=1,rg),'Eomes+ within IPCs':(e>=1,ipc)}
rows=[]
for label,(hit,denom) in masks.items():
 for display_group,(design,gate) in definitions.items():
  for r in key:
   if r['design_group']!=design:continue
   m=(sample==r['technical_sample_id'])&denom&gate;n=int((hit&m).sum());total=int(m.sum());rows.append(dict(population=label,sample=r['technical_sample_id'],group=display_group,count=n,total=total,percent=100*n/total if total else float('nan')))

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
  if ga=='HET_WT' and gb=='HET_KO':
   diffs=b-a;null=[np.mean(diffs*np.array(signs)) for signs in itertools.product([-1,1],repeat=3)]
  loo=[np.delete(diffs,i).mean() for i in range(3)] if ga=='HET_WT' and gb=='HET_KO' else [np.delete(b,j).mean()-np.delete(a,i).mean() for i in range(3) for j in range(3)]
  contrasts.append(dict(population=label,reference=ga,comparison=gb,difference_percentage_points=float(diff),exact_two_sided_p=float(np.mean(np.abs(null)>=abs(diff)-1e-10)),complete_sample_separation=bool(a.max()<b.min() or b.max()<a.min()),leave_one_per_group_direction_consistent=bool(all(v*diff>0 for v in loo)),paired=bool(ga=='HET_WT' and gb=='HET_KO'),sex_matched=bool((ga in ['WT_F','HET_WT','HET_KO'])==(gb in ['WT_F','HET_WT','HET_KO']))))
write('group_summary.tsv',summary);write('all_contrasts.tsv',contrasts)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False})
fig,axes=plt.subplots(1,3,figsize=(10.5,3.6));colors=['#AEB5BE','#5F997E','#B78366','#437C90','#D18A43']
for ax,label,title in zip(axes,['Eomes+ (>=1 UMI)','Radial glia','MGE IPCs'],['Eomes / Tbr2-positive','Radial glia','MGE IPCs']):
 ss=[r for g in groups for r in summary if r['population']==label and r['group']==g];means=[r['mean_percent'] for r in ss];lo=[r['mean_percent']-r['min_percent'] for r in ss];hi=[r['max_percent']-r['mean_percent'] for r in ss]
 ax.bar(range(5),means,color=colors,width=.65);ax.errorbar(range(5),means,yerr=[lo,hi],fmt='none',ecolor='#333333',capsize=3,lw=1)
 ax.set_xticks(range(5));ax.set_xticklabels(['WT\nfemale','HET\nWT cells','HET\nKO cells','WT males\nPcdh19+','KO\nmales'],fontsize=7.5);ax.set_title(title,fontweight='bold');ax.set_ylabel('Cells (% of defined group)');ax.set_ylim(bottom=0)
fig.suptitle('Mouse MGE · groups defined as in panels C and E',fontweight='bold');fig.tight_layout(rect=[0,0,1,.94])
for ext in ['png','pdf','svg']:fig.savefig(str(OUT/('Eomes_RG_IPC_proportions.'+ext)),dpi=300)
(OUT/'METHODS.json').write_text(json.dumps(dict(source=str(SRC),sample_key=str(ROOT/'config/sample_key.csv'),n_cells=len(sample),probe_filter='C/E group definitions extended to all QC cell types: female WT A-positive, HET WT A-positive, HET KO A=0 and B+C>=1; male WT total Pcdh19>0, male KO all cells',Eomes='Tbr2 alias; raw UMI >=1, >=2 sensitivity analysis',RG='existing provisional C0002 annotation',IPC='existing C0003+C0008 annotation; not equated to Eomes positivity',averaging='equal-weight mean of 3 sample percentages per group',plot_errorbars='sample minimum to maximum',statistics='unpaired: all 20 exact assignments; within-HET: 8 paired sign flips; two-sided unadjusted exploratory P values; no cell-level replication',cross_sex='KO male vs female HET/WT confounds genotype with sex'),indent=2))
print(json.dumps(summary,indent=2));print('KO_VS_WT_MALE',json.dumps([r for r in contrasts if r['reference']=='WT_M' and r['comparison']=='KO_M'],indent=2))
