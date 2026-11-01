"""Aggregate raw Step 06 counts by sample and frozen HET probe class."""
import csv
import gzip
import json
import hashlib
from pathlib import Path
from pcdh19_exploratory_paths import RUNTIME_ROOT
from collections import Counter
import h5py
import numpy as np
from scipy import sparse

ROOT=Path(__file__).resolve().parents[1]
OUT=RUNTIME_ROOT/'results/probe_subset_DE_20261031'
OUT.mkdir(exist_ok=True,parents=True)
BASE=Path('/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/paper3_pcdh19/results')
H5=BASE/'primary_processing/06_technical_sample_batch_diagnostics/06_technical_sample_batch_diagnostics_20260903_114110_6087eeb/objects/pcdh19_step06_unintegrated_diagnostics.h5ad'
CALLS=BASE/'step_07_pcdh19_het_female_wt_ko_like_classification/step_07_het_female_cell_classifications.tsv'
with CALLS.open() as f:
    calls={r['biological_sample_id']+'_'+r['cell_barcode']:r['final_classification'] for r in csv.DictReader(f,delimiter='\t')}
assert len(calls)==101102
libraries=[('WT_F_all','15662-JZ-'+str(i)) for i in [4,5,6]]
libraries += [(g,'15662-JZ-'+str(i)) for g in ['HET_WT_like','HET_KO_like'] for i in [7,8,9]]
lookup={v:i for i,v in enumerate(libraries)}
def strings(node):
    if isinstance(node,h5py.Group):
        cats=strings(node['categories']);return np.array([cats[i] if i>=0 else '' for i in node['codes'][:]])
    return np.array([x.decode() if isinstance(x,bytes) else str(x) for x in node[:]])
with h5py.File(H5,'r') as f:
    ids=strings(f['obs/cell_id']);genes=strings(f['var/gene_id']);symbols=strings(f['var/gene_symbol'])
    assert len(ids)==446349 and len(genes)==19071
    selected=np.full(len(ids),-1,dtype=np.int16)
    for i,cell in enumerate(ids):
        s=cell.split('_')[0]
        if s in ['15662-JZ-4','15662-JZ-5','15662-JZ-6']:selected[i]=lookup['WT_F_all',s]
        elif cell in calls and calls[cell] in ['WT_like','KO_like']:
            selected[i]=lookup['HET_'+calls[cell],s]
    cells=Counter(int(x) for x in selected if x>=0)
    assert [cells[i] for i in range(3,9)]==[467,1069,2411,26,10,103]
    x=f['X'];assert x.attrs['encoding-type']=='csr_matrix'
    assert x['data'].dtype.kind in 'iu'
    ptr=x['indptr'][:];sums=np.zeros((len(genes),9),dtype=np.int64)
    for start in range(0,len(ids),2000):
        end=min(start+2000,len(ids));groups=selected[start:end]
        if not np.any(groups>=0):continue
        lo,hi=int(ptr[start]),int(ptr[end]);data=x['data'][lo:hi]
        assert np.all(data>=0)
        mat=sparse.csr_matrix((data,x['indices'][lo:hi],ptr[start:end+1]-lo),shape=(end-start,len(genes)))
        for g in set(int(v) for v in groups if v>=0):sums[:,g]+=np.asarray(mat[groups==g].sum(axis=0)).ravel()
        print('AGGREGATED_ROWS',end,flush=True)
assert np.all(sums.sum(axis=0)>0)
with gzip.open(OUT/'raw_pseudobulk_counts.tsv.gz','wt') as f:
    w=csv.writer(f,delimiter='\t');w.writerow(['gene_id','gene_symbol']+['%s__%s'%lib for lib in libraries])
    for i,(gene,symbol) in enumerate(zip(genes,symbols)):w.writerow([gene,symbol]+list(sums[i]))
with (OUT/'sample_metadata.tsv').open('w') as f:
    w=csv.writer(f,delimiter='\t');w.writerow(['library','group','sample','cells','total_UMI'])
    for i,(g,s) in enumerate(libraries):w.writerow(['%s__%s'%(g,s),g,s,cells[i],int(sums[:,i].sum())])
(OUT/'aggregation_validation.json').write_text(json.dumps(dict(raw_integer_counts=True,genes=len(genes),cells={('%s__%s'%lib):cells[i] for i,lib in enumerate(libraries)},input_h5=str(H5),input_h5_bytes=H5.stat().st_size,calls_sha256=hashlib.sha256(CALLS.read_bytes()).hexdigest()),indent=2))
print('AGGREGATION_COMPLETE',str(OUT),flush=True)
