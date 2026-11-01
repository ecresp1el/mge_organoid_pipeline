from pathlib import Path
from pcdh19_exploratory_paths import RUNTIME_ROOT
import csv,gzip,json,h5py,numpy as np
ROOT=Path(__file__).resolve().parents[1];OUT=RUNTIME_ROOT/'results/female_male_signaling_20261031';BASE=Path('/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/paper3_pcdh19/results')
p=OUT/'atlas_sex_genotype_display.tsv.gz'
if p.exists():print('REUSE_ATLAS_SEX_GENOTYPE');raise SystemExit
def strings(ds):return np.array([v.decode() if isinstance(v,bytes) else str(v) for v in ds[:]])
with h5py.File(str(RUNTIME_ROOT/'results/mosaic_deep_dive_20261031/per_cell_marker_diagnostics.h5'),'r') as f:s=strings(f['sample']);cl=strings(f['cluster'])
with h5py.File(str(BASE/'primary_processing/06_technical_sample_batch_diagnostics/06_technical_sample_batch_diagnostics_20260903_114110_6087eeb/objects/pcdh19_step06_unintegrated_diagnostics.h5ad'),'r') as f:xy=f['obsm/X_umap'][:]
assert len(xy)==len(s);rng=np.random.RandomState(20261031);ix=[]
for n in range(1,13):
 rows=np.flatnonzero(s=='15662-JZ-'+str(n));ix.extend(rng.choice(rows,min(6000,len(rows)),replace=False).tolist())
with gzip.open(str(p),'wt') as f:
 w=csv.writer(f,delimiter='\t');w.writerow(['row','cluster','sex','genotype','UMAP1','UMAP2'])
 for i in ix:
  n=int(s[i].split('-')[-1]);w.writerow([i,cl[i],'Female' if n in [4,5,6,7,8,9] else 'Male','WT' if n<=6 else ('HET' if n<=9 else 'KO'),float(xy[i,0]),float(xy[i,1])])
(OUT/'atlas_sex_genotype_provenance.json').write_text(json.dumps({'cells':len(ix),'sample_balanced_display':6000,'all_12_registered_samples':True,'same_saved_joint_embedding_for_all_views':True,'embedding_refit':False,'sample_seed':20261031},indent=2));print('ATLAS_SEX_GENOTYPE_READY',len(ix),flush=True)
