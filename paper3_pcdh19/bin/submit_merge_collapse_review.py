#!/usr/bin/env python3
"""Snapshot completed Step 07 evidence and submit a read-only partition audit.

The source report job may still be plotting. Only already completed numeric
assets are copied; source and copy hashes must agree. No original run is changed.
The resulting package records source paths and hashes and never imports HiCAT.
"""
import argparse
from datetime import datetime,timezone
import hashlib,json,shutil,subprocess
from pathlib import Path
import h5py
import pandas as pd


def sha(p):
    """Hash a finite source asset without loading an expression matrix."""
    h=hashlib.sha256()
    with p.open('rb') as f:
        for block in iter(lambda:f.read(8*1024*1024),b''):h.update(block)
    return h.hexdigest()


def decode(node):
    """Decode only AnnData categorical/string metadata used by this audit."""
    if isinstance(node,h5py.Group):
        cats=decode(node['categories']);return [cats[i] if i>=0 else None for i in node['codes'][:]]
    return [v.decode() if isinstance(v,bytes) else v for v in node[:]]


def main():
    """Freeze exact input bytes, metadata extraction and audit code before SLURM."""
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source-run',type=Path,required=True);args=p.parse_args()
    source=args.source_run;bundle=Path(__file__).resolve().parents[1]
    stamp=datetime.now(timezone.utc).strftime('%Y%m%d_%H%M%S')
    package=source/'review'/('merge_collapse_'+stamp);package.mkdir(parents=True)
    for sub in ['inputs','code','logs']:(package/sub).mkdir()
    mapping={'program_means.tsv':'marker_programs/fine_cluster_program_means.tsv','counts.tsv':'parameter_sensitivity/fine_counts.tsv','assignments.tsv.gz':'parameter_sensitivity/cell_assignments.tsv.gz','phase_pairs.tsv':'cell_cycle/matched_identity_pairs.tsv','phase_effects.tsv':'cell_cycle/matched_identity_phase_effects.tsv','seed.tsv':'stability/fine_clusters.tsv','sample_fractions.tsv':'sample_composition/fine_cluster_fractions.tsv','cluster_review.tsv':'annotation_review/cluster_validation_summary.tsv','gene_coverage.tsv':'marker_programs/gene_coverage.tsv','program_config.json':'marker_programs/canonical_program_config.json','canonical_expression.tsv.gz':'marker_programs/canonical_cell_expression.tsv.gz'}
    cfg=json.loads((source/'config/hicat_validation.json').read_text());baseline=Path(cfg['baseline_run']);records=[]
    for name,relative in mapping.items():
        file=source/'outputs'/relative
        if not file.exists():file=source/'staging/outputs'/relative
        before=sha(file);shutil.copy2(file,package/'inputs'/name)
        assert before==sha(file)==sha(package/'inputs'/name)
        records.append(dict(file=name,source=str(file),sha256=before))
    file=baseline/'baseline/global_fine_pairwise_evidence.tsv';shutil.copy2(file,package/'inputs/baseline_de.tsv');records.append(dict(file='baseline_de.tsv',source=str(file),sha256=sha(file)))
    for parent in ['C0001','C0002','C0003','C0004']:
        file=source/'outputs/parameter_sensitivity'/('fine_'+parent)/'candidate_summary.json'
        if not file.exists():file=source/'staging/outputs/parameter_sensitivity'/('fine_'+parent)/'candidate_summary.json'
        name=parent+'_candidate_summary.json';shutil.copy2(file,package/'inputs'/name);records.append(dict(file=name,source=str(file),sha256=sha(file)))
    pilot=baseline/'pcdh19_hicat_coarse_fine.h5ad'
    with h5py.File(pilot,'r') as h:
        o=h['obs'];key=o.attrs['_index'];index=decode(o[key]);names=['hicat_fine_baseline','hicat_coarse_baseline','technical_sample_id']
        frame=pd.DataFrame({name:decode(o[name]) for name in names},index=index)
    file=package/'inputs/baseline_cells.tsv';frame.to_csv(file,sep='\t',index_label='cell_id');records.append(dict(file=file.name,source=str(pilot)+'::obs (metadata extraction)',sha256=sha(file)))
    pd.DataFrame(records).to_csv(package/'input_manifest.tsv',sep='\t',index=False)
    thresholds=dict(partial_source_fraction=.1,majority_source_fraction=.5,near_complete_source_fraction=.9,effect_sd=.5,marker_mean_ln1pcpm_difference=.5,marker_detection_difference=.1,cycle_tv=.4,sample_tv=.4,new_fits=False,source_run=str(source),original_pilot_sha256=cfg['pilot_sha256'])
    (package/'review_config.json').write_text(json.dumps(thresholds,indent=2))
    shutil.copy2(bundle/'scripts/hicat/merge_collapse_review.py',package/'code/merge_collapse_review.py');shutil.copy2(__file__,package/'code/submit_merge_collapse_review.py')
    interpreter='/home/elcrespo/miniconda3/envs/pcdh19-hicat-py38/bin/python'
    (package/'environment.txt').write_text(subprocess.check_output([interpreter,'-m','pip','freeze'],text=True))
    jobfile=package/'code/run.sbatch'
    jobfile.write_text('#!/bin/bash\nset -Eeuo pipefail\nexport OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MPLBACKEND=Agg PYTHONDONTWRITEBYTECODE=1\n'+interpreter+' '+str(package/'code/merge_collapse_review.py')+' --package '+str(package)+'\n')
    (package/'frozen_code_manifest.json').write_text(json.dumps({str(p.relative_to(package)):sha(p) for p in (package/'code').iterdir()},indent=2))
    job=subprocess.check_output(['sbatch','--parsable','--account=parent0','--partition=standard','--cpus-per-task=1','--mem=8G','--time=00:30:00','--job-name=pcdh19-merge-review','--output='+str(package/'logs/scheduler.out'),'--error='+str(package/'logs/scheduler.err'),str(jobfile)],text=True).strip()
    (package/'job_id.txt').write_text(job+'\n');print(json.dumps(dict(package=str(package),job_id=job)))


if __name__=='__main__':main()
