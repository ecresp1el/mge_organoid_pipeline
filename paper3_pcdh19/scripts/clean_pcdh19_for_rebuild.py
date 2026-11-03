#!/usr/bin/env python3
"""Execute the user-authorized PCDH19 output reset after preserving essential inputs."""
from pathlib import Path
import datetime,hashlib,json,os,shutil
import anndata as ad
import numpy as np

ROOT=Path('/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/paper3_pcdh19')
OUT=ROOT/'results/rebuild_20261103'
BASELINE=OUT/'baseline';BASELINE.mkdir(parents=True,exist_ok=True)
STATUS=OUT/'cleanup_status.json'
def publish(value):
 p=STATUS.with_suffix('.tmp');p.write_text(json.dumps(value,indent=2));p.replace(STATUS)
def sha(path):
 h=hashlib.sha256()
 with path.open('rb') as f:
  for block in iter(lambda:f.read(8*1024*1024),b''):h.update(block)
 return h.hexdigest()
def copy(source,target):
 target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(str(source),str(target))
 assert sha(source)==sha(target)
def bounded_files(folder,depth=0):
 if depth>12:raise ValueError('Depth guard: '+str(folder))
 for p in folder.iterdir():
  if p.is_symlink():continue
  if p.stat().st_dev!=ROOT.stat().st_dev:raise ValueError('Cross-filesystem guard: '+str(p))
  if p.is_dir():yield from bounded_files(p,depth+1)
  elif p.is_file():yield p

if STATUS.exists():
 raise RuntimeError('Cleanup already started; inspect its exact manifest/status rather than running twice')
publish({'stage':'preserving','authorization':'User explicitly chose preserve final archive and essential inputs; remove replaceable outputs on 2026-11-03'})
archive=ROOT/'final_figures/mge_activity_IEG_FINAL_20261101'
for row in json.loads((archive/'frozen_manifest.json').read_text()):
 assert sha(archive/row['path'])==row['sha256'], row['path']
raw=ROOT/'results/primary_processing/00_input_validation_and_canonical_anndata/00_input_validation_and_canonical_anndata_20260830_113749_d8b6bf7/objects/pcdh19_step00_canonical_raw_counts.h5ad'
qc=ROOT/'results/primary_processing/02_qc_filtering/02_qc_filtering_20260830_124611_97e1bb5/objects/pcdh19_step02_qc_filtered.h5ad'
assert raw.is_file() and qc.is_file()
assert qc.stat().st_size==6554152927
assert sha(qc)=='fadba4a25a7b6b7320219b21c189b6325687493519ba0fe1bd27efd79606b103'
old=ROOT/'results/primary_processing/06_technical_sample_batch_diagnostics/06_technical_sample_batch_diagnostics_20260903_114110_6087eeb'
a=ad.read_h5ad(old/'objects/pcdh19_step06_unintegrated_diagnostics.h5ad',backed='r')
assert a.n_obs==446349 and a.n_vars==19071
np.savez_compressed(BASELINE/'original_step06_embeddings.npz',cell_id=np.asarray(a.obs_names,dtype=str),X_umap=a.obsm['X_umap'],X_pca=a.obsm['X_pca'])
a.obs.to_csv(BASELINE/'original_step06_obs.tsv.gz',sep='\t',compression='gzip')
a.var.to_csv(BASELINE/'original_step06_var.tsv.gz',sep='\t',compression='gzip')
a.file.close()
review=ROOT/'results/hicat/03_full_data_consensus_benchmark/hicat_consensus_review_20260909_142348_efd9e13d'
copy(review/'outputs/tables/all_cell_assignments.tsv.gz',BASELINE/'original_hicat_assignments.tsv.gz')
hicat=ROOT/'results/hicat/03_full_data_consensus_benchmark/consensus_restart_20260909_000256'
copy(hicat/'config/run.json',BASELINE/'previous_hicat_run.json')
# Save a bounded standalone upstream Python source snapshot before deleting old runs.
for relative in ['adapter','reference/allen_python']:
 for source in bounded_files(hicat/relative):
  rel=source.relative_to(hicat)
  if any(x in ['.git','__pycache__','logs'] for x in rel.parts):continue
  if source.suffix in ['.py','.json','.toml','.txt','.md','.rst','.cfg','.ini'] or source.name in ['LICENSE','setup.py','setup.cfg']:
   copy(source,BASELINE/'previous_hicat_source'/rel)
# Preserve lightweight processing contracts/parameters/ledgers without duplicate matrices.
for source in bounded_files(ROOT/'results/primary_processing'):
 if source.suffix in ['.json','.env'] or source.name=='APPROVAL_LEDGER.tsv':
  if source.stat().st_size<5*1024*1024:
   copy(source,BASELINE/'previous_processing_contracts'/source.relative_to(ROOT/'results/primary_processing'))
for name in ['APPROVAL_LEDGER.tsv']:
 copy(ROOT/'results/hicat'/name,BASELINE/'previous_hicat_approval_ledger.tsv')
keep_results={'primary_processing','rebuild_20261103','pcdh19_probe_audit','egfp_probe_audit','xgfp_probe_audit','00_developing_mouse_mge_reference_curation'}
keep_primary={'00_input_validation_and_canonical_anndata','02_qc_filtering'}
targets=[p for p in (ROOT/'results').iterdir() if p.is_dir() and p.name not in keep_results]
targets += [p for p in (ROOT/'results/primary_processing').iterdir() if p.is_dir() and p.name not in keep_primary]
targets += [p for p in (ROOT/'final_figures').iterdir() if p.is_dir() and p!=archive]
protected=[raw,qc,archive,OUT,ROOT/'inputs',ROOT/'software']
rows=[]
for target in sorted(targets):
 assert target.is_relative_to(ROOT) and not target.is_symlink()
 assert not any(target==p or target in p.parents or p in target.parents for p in protected), str(target)
 entries=list(bounded_files(target))
 rows.append({'path':str(target),'files':len(entries),'bytes':sum(p.stat().st_size for p in entries),
              'allocated_bytes':sum(p.stat().st_blocks*512 for p in entries)})
(OUT/'cleanup_manifest.json').write_text(json.dumps({'targets':rows,'preserved_raw':str(raw),'preserved_qc':str(qc),
 'preserved_archive':str(archive),'preserved_baseline':str(BASELINE),'preserve_inputs_and_software':True,
 'note':'Sizes are path totals; hard links may reduce actual quota space reclaimed.'},indent=2))
copy(Path(__file__),OUT/Path(__file__).name)
publish({'stage':'deleting','directories':len(rows),'path_bytes':sum(r['bytes'] for r in rows)})
for row in rows:
 target=Path(row['path'])
 # shutil removes this exact project subtree and never follows contained symlinks.
 shutil.rmtree(str(target))
 with (OUT/'cleanup_deleted.jsonl').open('a') as f:f.write(json.dumps(row)+'\n')
assert raw.exists() and qc.exists() and archive.exists()
for row in json.loads((archive/'frozen_manifest.json').read_text()):assert sha(archive/row['path'])==row['sha256']
publish({'stage':'complete','directories_removed':len(rows),'files_removed':sum(r['files'] for r in rows),
 'removed_path_bytes':sum(r['bytes'] for r in rows),'removed_path_GiB':sum(r['bytes'] for r in rows)/2**30,
 'removed_allocated_path_GiB':sum(r['allocated_bytes'] for r in rows)/2**30,
 'archive_hashes_verified':True,'qc_input_hash_verified':True,'raw_checkpoint_preserved':True,
 'completed_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat()})
print(STATUS.read_text(),flush=True)
