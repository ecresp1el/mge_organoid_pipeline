from pathlib import Path
from pcdh19_exploratory_paths import RUNTIME_ROOT
import csv,json,h5py,numpy as np
from scipy import sparse
ROOT=Path(__file__).resolve().parents[1];OUT=RUNTIME_ROOT/'results/female_male_signaling_20261031';BASE=Path('/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/paper3_pcdh19/results');RAW=BASE/'primary_processing/06_technical_sample_batch_diagnostics/06_technical_sample_batch_diagnostics_20260903_114110_6087eeb/objects/pcdh19_step06_unintegrated_diagnostics.h5ad'
def strings(ds):
 if isinstance(ds,h5py.Group):
  cats=strings(ds['categories']);return np.array([cats[i] for i in ds['codes'][:]])
 return np.array([v.decode() if isinstance(v,bytes) else str(v) for v in ds[:]])
if (OUT/'male_cell_expression.npz').exists():print('REUSE_MALE_CELL_PANELS');raise SystemExit
with h5py.File(str(RUNTIME_ROOT/'results/mosaic_deep_dive_20261031/per_cell_marker_diagnostics.h5'),'r') as f:
 sample=strings(f['sample']);cluster=strings(f['cluster']);total=f['total_counts'][:];ids=strings(f['cell_id'])
ms=['15662-JZ-'+str(i) for i in [1,2,3,10,11,12]];rows=np.flatnonzero(np.isin(sample,ms)&(cluster=='C0006'));genes=['Eps8l1','Grm7','Lrfn5'];counts=np.zeros((len(rows),len(genes)),dtype=np.int32)
with h5py.File(str(RAW),'r') as f:
 symbols=strings(f['var/gene_symbol']);gi=[list(symbols).index(g) for g in genes];x=f['X'];ptr=x['indptr'][:]
 for start in range(0,len(sample),2000):
  end=min(start+2000,len(sample));sel=np.flatnonzero((rows>=start)&(rows<end))
  if not len(sel):continue
  lo,hi=int(ptr[start]),int(ptr[end]);mat=sparse.csr_matrix((x['data'][lo:hi],x['indices'][lo:hi],ptr[start:end+1]-lo),shape=(end-start,len(symbols)));counts[sel]=mat[rows[sel]-start][:,gi].toarray()
lookup={ids[row]:i for i,row in enumerate(rows)};a=np.full(len(rows),-1,dtype=int)
for s in ms:
 with open(str(BASE/'pcdh19_probe_audit/per_sample'/s/'pcdh19_probe_patterns.tsv')) as f:
  for r in csv.DictReader(f,delimiter='\t'):
   i=lookup.get(s+'_'+r['barcode'])
   if i is not None:a[i]=int(r['A_UMI'])
assert (a>=0).all()
group=np.array([0 if s in ms[:3] else 1 for s in sample[rows]])
np.savez_compressed(str(OUT/'male_cell_expression.npz'),genes=np.array(genes),CP10K=counts*10000./total[rows,None],group=group,sample=sample[rows],exon1_raw=a)
with open(str(OUT/'male_exon1_detection.tsv'),'w') as f:
 w=csv.writer(f,delimiter='\t');w.writerow(['sample','group','cells','exon1_positive','fraction'])
 for s in ms:
  mask=sample[rows]==s;n=int(mask.sum());pos=int((a[mask]>0).sum());w.writerow([s,'WT males' if s in ms[:3] else 'KO males',n,pos,pos/n])
print('MALE_CELL_PANELS_READY',np.bincount(group).tolist(),flush=True)
