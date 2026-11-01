from pathlib import Path
from pcdh19_exploratory_paths import RUNTIME_ROOT
import json
import h5py
import numpy as np
from scipy import sparse
ROOT=Path(__file__).resolve().parents[1]
OUT=RUNTIME_ROOT/'results/neuronal_matched_selection_20261031'
CACHE=OUT/'D_pattern_violin_cells.npz'
if CACHE.exists():
    print('REUSING_CELL_EXPRESSION');raise SystemExit
SIDE=RUNTIME_ROOT/'results/mosaic_deep_dive_20261031/per_cell_marker_diagnostics.h5'
RAW=Path('/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/paper3_pcdh19/results/primary_processing/06_technical_sample_batch_diagnostics/06_technical_sample_batch_diagnostics_20260903_114110_6087eeb/objects/pcdh19_step06_unintegrated_diagnostics.h5ad')
def strings(ds):
    if isinstance(ds,h5py.Group):
        cats=strings(ds['categories']);return np.array([cats[i] if i>=0 else '' for i in ds['codes'][:]])
    return np.array([v.decode() if isinstance(v,bytes) else str(v) for v in ds[:]])
with h5py.File(str(SIDE),'r') as f:
    sample=strings(f['sample']);cluster=strings(f['cluster']);gates=strings(f['gate_names']);total=f['total_counts'][:];gm=f['gate_membership'][:]
    group=np.full(len(sample),-1,dtype=np.int8)
    wt=np.isin(sample,['15662-JZ-'+str(i) for i in [4,5,6]])
    het=np.isin(sample,['15662-JZ-'+str(i) for i in [7,8,9]])
    group[wt & (cluster=='C0006')]=0
    group[het & (cluster=='C0006') & gm[:,list(gates).index('A_positive')]]=1
    group[het & (cluster=='C0006') & gm[:,list(gates).index('Aneg_BC1')]]=2
fullrows=np.flatnonzero(group>=0)
exon=np.load(str(RUNTIME_ROOT/'results/neuronal_mosaic_focus_20261031/neuronal_exon1_signal.npz'))
keep=(group[fullrows]!=0)|(exon['raw_exon1']>0);rows=fullrows[keep]
genes=['Pcdh19']+json.loads((OUT/'D_pattern_gene_selection.json').read_text())['genes'];counts=np.zeros((len(rows),len(genes)),dtype=np.int32)
with h5py.File(str(RAW),'r') as f:
    symbols=strings(f['var/gene_symbol']);gi=[list(symbols).index(g) for g in genes];x=f['X'];ptr=x['indptr'][:]
    for start in range(0,len(sample),2000):
        end=min(start+2000,len(sample));selected=np.flatnonzero((rows>=start)&(rows<end))
        if not len(selected):continue
        lo,hi=int(ptr[start]),int(ptr[end]);mat=sparse.csr_matrix((x['data'][lo:hi],x['indices'][lo:hi],ptr[start:end+1]-lo),shape=(end-start,len(symbols)))
        counts[selected]=mat[rows[selected]-start][:,gi].toarray()
normalized=counts.astype(float)*10000/total[rows,None]
assert np.isfinite(normalized).all() and (normalized>=0).all()
np.savez_compressed(str(CACHE),group=group[rows],raw_counts=counts,CP10K=normalized,genes=np.array(genes),total_counts=total[rows])
(OUT/'D_pattern_violin_data_provenance.json').write_text(json.dumps({'genes':genes,'cells_by_group':np.bincount(group[rows]).tolist(),'normalization':'raw gene UMI / cell total UMI * 10000; linear, no log transformation','zeros_included':True,'source':str(RAW),'subset':'same broad neuronal consensus cluster; WT and HET WT-like both A-positive; HET KO-like A-negative B/C-positive','cell_sampling':'all eligible cells; no distribution trimming'},indent=2))
print('CELL_EXPRESSION_COMPLETE',len(rows),flush=True)

for j,g in enumerate(genes):print('CELL_MEANS',g,[round(float(normalized[group[rows]==k,j].mean()),4) for k in range(3)],flush=True)
