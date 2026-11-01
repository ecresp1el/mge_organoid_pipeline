from pathlib import Path
from pcdh19_exploratory_paths import RUNTIME_ROOT
import csv,json
import h5py,numpy as np
from scipy import sparse
ROOT=Path(__file__).resolve().parents[1];OUT=RUNTIME_ROOT/'results/female_male_signaling_20261031';BASE=Path('/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/paper3_pcdh19/results')
RAW=BASE/'primary_processing/06_technical_sample_batch_diagnostics/06_technical_sample_batch_diagnostics_20260903_114110_6087eeb/objects/pcdh19_step06_unintegrated_diagnostics.h5ad'
def strings(ds):
 if isinstance(ds,h5py.Group):
  cats=strings(ds['categories']);return np.array([cats[i] for i in ds['codes'][:]])
 return np.array([v.decode() if isinstance(v,bytes) else str(v) for v in ds[:]])
if (OUT/'MALE_PREPARATION_COMPLETE.json').exists():print('REUSE_MALE_PSEUDOBULKS');raise SystemExit
with h5py.File(str(RUNTIME_ROOT/'results/mosaic_deep_dive_20261031/per_cell_marker_diagnostics.h5'),'r') as f:
 samples=strings(f['sample']);cluster=strings(f['cluster'])
ms=['15662-JZ-'+str(i) for i in [1,2,3,10,11,12]];strata=['all_cells','C0001','C0002','C0003','C0006','C0011'];libs=[(s,st) for st in strata for s in ms]
with h5py.File(str(RAW),'r') as f:
 genes=strings(f['var/gene_id']);symbols=strings(f['var/gene_symbol']);x=f['X'];ptr=x['indptr'][:];counts=np.zeros((len(genes),len(libs)),dtype=np.int64);nc=np.zeros(len(libs),dtype=int)
 for start in range(0,len(samples),2000):
  end=min(start+2000,len(samples));ss=samples[start:end];cc=cluster[start:end]
  if not np.isin(ss,ms).any():continue
  lo,hi=int(ptr[start]),int(ptr[end]);mat=sparse.csr_matrix((x['data'][lo:hi],x['indices'][lo:hi],ptr[start:end+1]-lo),shape=(end-start,len(genes)))
  for k,(s,st) in enumerate(libs):
   mask=(ss==s)&((cc==st) if st!='all_cells' else True)
   if mask.any():counts[:,k]+=np.asarray(mat[mask].sum(0)).ravel();nc[k]+=int(mask.sum())
meta=[]
for k,(s,st) in enumerate(libs):meta.append(dict(library=s+'__'+st,sample=s,stratum=st,group='WT male' if s in ms[:3] else 'KO male',cells=int(nc[k]),UMI=int(counts[:,k].sum())))
with open(str(OUT/'male_metadata.tsv'),'w') as f:
 w=csv.DictWriter(f,fieldnames=list(meta[0]),delimiter='\t');w.writeheader();w.writerows(meta)
for st in strata:
 ix=[k for k,z in enumerate(libs) if z[1]==st]
 with open(str(OUT/('male_'+st+'_counts.tsv')),'w') as f:
  w=csv.writer(f,delimiter='\t');w.writerow(['gene_id','gene_symbol']+[meta[k]['library'] for k in ix]);w.writerows([[genes[j],symbols[j]]+counts[j,ix].tolist() for j in range(len(genes))])
comp=[]
for s in ms:
 sm=samples==s;den=int(sm.sum())
 for c in sorted(set(cluster)):
  n=int((sm&(cluster==c)).sum());comp.append(dict(sample=s,group='WT male' if s in ms[:3] else 'KO male',cluster=c,cells=n,total=den,fraction=n/den))
with open(str(OUT/'male_cluster_composition.tsv'),'w') as f:
 w=csv.DictWriter(f,fieldnames=list(comp[0]),delimiter='\t');w.writeheader();w.writerows(comp)
(OUT/'MALE_PREPARATION_COMPLETE.json').write_text(json.dumps({'samples':ms,'strata':strata,'cells_all':{s:int(sum(samples==s)) for s in ms},'genotype_source':'registered sample assignments; no probe selection','unit':'sample pseudobulk; three samples per group','raw_source':str(RAW)},indent=2));print('MALE_PREPARED',flush=True)
