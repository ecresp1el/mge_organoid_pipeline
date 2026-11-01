from pathlib import Path
from pcdh19_exploratory_paths import RUNTIME_ROOT
import csv,json,h5py
import numpy as np
ROOT=Path(__file__).resolve().parents[1];OUT=RUNTIME_ROOT/'results/female_male_signaling_20261031';BASE=Path('/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/paper3_pcdh19/results')
def strings(ds):return np.array([v.decode() if isinstance(v,bytes) else str(v) for v in ds[:]])
with h5py.File(str(RUNTIME_ROOT/'results/mosaic_deep_dive_20261031/per_cell_marker_diagnostics.h5'),'r') as f:
 s=strings(f['sample']);cl=strings(f['cluster']);ids=strings(f['cell_id']);total=f['total_counts'][:]
ms=['15662-JZ-'+str(i) for i in [1,2,3,10,11,12]];rows=np.flatnonzero(np.isin(s,ms)&(cl=='C0006'));lookup={ids[row]:i for i,row in enumerate(rows)};counts=np.full(len(rows),-1,dtype=int)
for sample in ms:
 with open(str(BASE/'pcdh19_probe_audit/per_sample'/sample/'pcdh19_probe_patterns.tsv')) as f:
  for r in csv.DictReader(f,delimiter='\t'):
   i=lookup.get(sample+'_'+r['barcode'])
   if i is None:continue
   counts[i]=int(r['Pcdh19_total_UMI']);assert counts[i]==sum(int(r[k]) for k in ['A_UMI','B_UMI','C_UMI'])
assert np.all(counts>=0)
group=np.array([0 if sample in ms[:3] else 1 for sample in s[rows]]);v=counts*10000./total[rows]
np.savez_compressed(str(OUT/'male_total_Pcdh19_display.npz'),CP10K=v,raw_counts=counts,group=group,sample=s[rows])
summary={}
for i,label in enumerate(['WT males','KO males']):
 mask=group==i;summary[label]=dict(cells=int(mask.sum()),detected_fraction=float((counts[mask]>0).mean()),mean_CP10K=float(v[mask].mean()),percentile99_CP10K=float(np.percentile(v[mask],99)))
(OUT/'MALE_TOTAL_PCDH19_METHODS.json').write_text(json.dumps(dict(selection='all QC C0006 immature inhibitory neurons from known male WT/KO genotypes; no A/B/C expression gate or female-state filter',measurement='total Pcdh19 UMI=A+B+C, divided by total cell UMI * 10000',zero_cells_retained=True,summary=summary),indent=2))
print('MALE_TOTAL_PCDH19',summary,flush=True)
