from pathlib import Path
from pcdh19_exploratory_paths import RUNTIME_ROOT
import csv,gzip,json,h5py
import numpy as np
ROOT=Path(__file__).resolve().parents[1];OUT=RUNTIME_ROOT/'results/female_male_signaling_20261031';SRC=Path('/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/paper3_pcdh19/results/hicat/02_hierarchy_validation/hicat_validation_20260908_221625_d5d06a5d/outputs/full_data_projection/plot_program_scores.tsv.gz')
with gzip.open(str(OUT/'atlas_sex_genotype_display.tsv.gz'),'rt') as f:cache=list(csv.DictReader(f,delimiter='\t'))
rows=np.array([int(r['row']) for r in cache])
with h5py.File(str(RUNTIME_ROOT/'results/mosaic_deep_dive_20261031/per_cell_marker_diagnostics.h5'),'r') as f:ids=[v.decode() if isinstance(v,bytes) else str(v) for v in f['cell_id'][:]]
selected=[ids[i] for i in rows];lookup={v:i for i,v in enumerate(selected)};assert len(lookup)==len(selected)
values=np.full((len(rows),3),np.nan);names=['MGE','LGE','CGE']
with gzip.open(str(SRC),'rt') as f:
 for r in csv.DictReader(f,delimiter='\t'):
  i=lookup.get(r['cell_id'])
  if i is not None:values[i]=[float(r[n]) for n in names]
assert np.isfinite(values).all()
limits=np.percentile(values,[2,98]);np.savez_compressed(str(OUT/'saved_MGE_LGE_CGE_scores.npz'),row=rows,cell_id=np.array(selected),names=np.array(names),scores=values,color_limits=limits)
method=dict(source=str(SRC),cells=len(rows),join='exact cell ID, using current atlas display rows',new_scores_computed=False,new_UMAP_computed=False,color_limits=limits.tolist(),shared_color_scale=True,score_means=dict(zip(names,values.mean(0).tolist())),score_medians=dict(zip(names,np.median(values,axis=0).tolist())),interpretation='saved regional marker-program scores; descriptive identity support, not proof of anatomical origin')
(OUT/'SAVED_REGIONAL_SCORES_PROVENANCE.json').write_text(json.dumps(method,indent=2));print('REGIONAL_SCORES_READY',method,flush=True)
