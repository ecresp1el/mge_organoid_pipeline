from pathlib import Path
from pcdh19_exploratory_paths import RUNTIME_ROOT
import csv,json
import numpy as np,h5py
ROOT=Path(__file__).resolve().parents[1];OUT=RUNTIME_ROOT/'results/neuronal_mosaic_focus_20261031'
BASE=Path('/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/paper3_pcdh19/results')
def strings(ds):return np.array([v.decode() if isinstance(v,bytes) else str(v) for v in ds[:]])
with h5py.File(str(RUNTIME_ROOT/'results/mosaic_deep_dive_20261031/per_cell_marker_diagnostics.h5'),'r') as f:
 sample=strings(f['sample']);cluster=strings(f['cluster']);gates=list(strings(f['gate_names']));gm=f['gate_membership'][:]
 group=np.full(len(sample),-1,dtype=np.int8)
 group[np.isin(sample,['15662-JZ-'+str(i) for i in [4,5,6]])&(cluster=='C0006')]=0
 het=np.isin(sample,['15662-JZ-'+str(i) for i in [7,8,9]])&(cluster=='C0006')
 group[het&gm[:,gates.index('A_positive')]]=1;group[het&gm[:,gates.index('Aneg_BC1')]]=2
 rows=np.flatnonzero(group>=0);ids=strings(f['cell_id'])[rows]
cache=np.load(str(OUT/'neuronal_violin_cells.npz'));assert np.array_equal(cache['group'],group[rows])
lookup={c:i for i,c in enumerate(ids)};a=np.full(len(rows),-1,dtype=np.int32);tot=np.full(len(rows),-1,dtype=np.int32)
for n in [4,5,6,7,8,9]:
 s='15662-JZ-'+str(n)
 with open(str(BASE/'pcdh19_probe_audit/per_sample'/s/'pcdh19_probe_patterns.tsv')) as f:
  for r in csv.DictReader(f,delimiter='\t'):
   i=lookup.get(s+'_'+r['barcode'])
   if i is not None:a[i]=int(r['A_UMI']);tot[i]=sum(int(r[k]) for k in ['A_UMI','B_UMI','C_UMI'])
assert (a>=0).all();assert np.array_equal(tot,cache['raw_counts'][:,0]);assert (a[group[rows]==2]==0).all()
np.savez_compressed(str(OUT/'neuronal_exon1_signal.npz'),group=cache['group'],exon1_CP10K=a*10000./cache['total_counts'],raw_exon1=a)
(OUT/'exon1_panel_provenance.json').write_text(json.dumps({'measurement':'exon-1 probe A UMI normalized per 10000 total cell UMI','total_gene_equals_A_plus_B_plus_C':True,'all_cells_verified_against_gene_counts':int(len(rows)),'KO_like_zero_by_definition':True,'not_independent_KO_validation':True,'known_KO_downstream_signal':'documented in HANDOFF.md cross-sample probe-pattern snapshot'},indent=2))
print('EXON1_VERIFIED',len(rows))
