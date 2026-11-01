from pathlib import Path
from pcdh19_exploratory_paths import RUNTIME_ROOT
import numpy as np,h5py,json
ROOT=Path(__file__).resolve().parents[1];OUT=RUNTIME_ROOT/'results/female_male_signaling_20261031'
c=np.load(str(OUT/'male_cell_expression.npz'))
def strings(ds):return np.array([v.decode() if isinstance(v,bytes) else str(v) for v in ds[:]])
with h5py.File(str(RUNTIME_ROOT/'results/mosaic_deep_dive_20261031/per_cell_marker_diagnostics.h5'),'r') as f:
 s=strings(f['sample']);cl=strings(f['cluster']);total=f['total_counts'][:]
ms=['15662-JZ-'+str(i) for i in [1,2,3,10,11,12]];rows=np.flatnonzero(np.isin(s,ms)&(cl=='C0006'));assert np.array_equal(s[rows],c['sample']);v=c['exon1_raw']*10000./total[rows]
assert np.isfinite(v).all() and (v>=0).all()
np.savez_compressed(str(OUT/'male_exon1_display.npz'),CP10K=v,group=c['group'])
print('MALE_EXON1_SAME_UNITS_READY',flush=True)
