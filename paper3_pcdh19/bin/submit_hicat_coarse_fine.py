#!/usr/bin/env python3
"""Freeze an isolated versioned HiCAT pilot and submit its bounded SLURM job.

This submitter performs no scientific fitting. It snapshots local source and
the pinned upstream checkout, records pending repository changes by content,
and submits only the explicitly authorized 12,000-cell coarse/fine pilot.
"""
import argparse
from datetime import datetime,timezone
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys


def command(*args):
    """Return checked subprocess stdout using argument arrays, without a shell."""
    return subprocess.check_output(list(args),text=True).strip()


def main():
    """Parse source location, freeze reproducibility inputs, then submit one job."""
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--upstream',type=Path,required=True)
    parser.add_argument('--prepare-only',action='store_true')
    args=parser.parse_args()
    bundle=Path(__file__).resolve().parents[1]
    repo=bundle.parent
    project=Path('/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/paper3_pcdh19')
    cfg=json.loads((bundle/'config/hicat_coarse_fine.json').read_text())
    if command('git','-C',str(args.upstream),'rev-parse','HEAD')!=cfg['upstream_commit']:
        raise ValueError('Upstream commit mismatch')
    if command('git','-C',str(args.upstream),'status','--porcelain'):
        raise ValueError('Upstream reference source must be unmodified')
    interpreter='/home/elcrespo/miniconda3/envs/pcdh19-hicat-py38/bin/python'
    command(interpreter,'-m','pip','check')
    stamp=datetime.now(timezone.utc).strftime('%Y%m%d_%H%M%S')
    config_hash=hashlib.sha256((bundle/'config/hicat_coarse_fine.json').read_bytes()).hexdigest()[:8]
    run=project/'results/hicat/01_coarse_fine_expanded_pilot'/('hicat_coarse_fine_'+stamp+'_'+config_hash)
    run.mkdir(parents=True)
    for name in ['code','config','inputs','logs','provenance','scratch']:(run/name).mkdir()
    shutil.copytree(bundle/'scripts/hicat',run/'code/hicat',ignore=shutil.ignore_patterns('__pycache__'))
    shutil.copy2(bundle/'slurm/hicat_coarse_fine.sbatch',run/'code/hicat_coarse_fine.sbatch')
    shutil.copy2(__file__,run/'code/submit_hicat_coarse_fine.py')
    # Preserve test-relative fixture paths so the frozen suite remains runnable.
    shutil.copytree(bundle/'tests',run/'code/tests',ignore=shutil.ignore_patterns('__pycache__'))
    (run/'code/config').mkdir()
    for fixture in ['hicat_pilot.json','hicat_coarse_fine.json']:
        shutil.copy2(bundle/'config'/fixture,run/'code/config'/fixture)
    shutil.copy2(bundle/'config/hicat_coarse_fine.json',run/'config/hicat_pilot.json')
    shutil.copy2(bundle/'config/hicat_pilot.requirements.txt',run/'config/hicat_pilot.requirements.txt')
    shutil.copy2(bundle/'HICAT_COARSE_FINE_PROTOCOL.md',run/'config/HICAT_COARSE_FINE_PROTOCOL.md')
    # Copy the checked commit's content; no upstream execution uses its mutable temp checkout.
    shutil.copytree(args.upstream,run/'upstream',ignore=shutil.ignore_patterns('.git','__pycache__'))
    source_state=dict(commit=command('git','-C',str(repo),'rev-parse','HEAD'),
                      status=command('git','-C',str(repo),'status','--short'),upstream_commit=cfg['upstream_commit'])
    (run/'provenance/repository_state.json').write_text(json.dumps(source_state,indent=2)+'\n')
    (run/'provenance/repository_tracked_diff.patch').write_text(command('git','-C',str(repo),'diff','HEAD','--binary')+'\n')
    (run/'config/python_pip_freeze.txt').write_text(command(interpreter,'-m','pip','freeze')+'\n')
    (run/'config/authorization.txt').write_text('2026-09-08 user: Increase cells per sample and run coarse then fine with dendrograms, DEGs and top-20 heatmaps.\nPooled 12 samples; cluster discovery/count review only; no annotations; no batch correction; 1000 cells/sample (12000); no full-data run.\n')
    rows=[]
    for group in ['code','config','upstream']:
        for path in sorted((run/group).rglob('*')):
            if path.is_file():rows.append(dict(path=str(path.relative_to(run)),bytes=path.stat().st_size,sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
    (run/'provenance/frozen_source_manifest.json').write_text(json.dumps(rows,indent=2)+'\n')
    if args.prepare_only:
        print(run);return
    job_file=project/'jobs'/(run.name+'.sbatch')
    if job_file.exists():raise FileExistsError(str(job_file))
    shutil.copy2(run/'code/hicat_coarse_fine.sbatch',job_file)
    job=command('sbatch','--parsable','--job-name=pcdh19-hicat-hierarchy','--account=parent0',
                '--partition=standard','--cpus-per-task=4','--mem=64G','--time=04:00:00',
                '--output='+str(run/'logs/scheduler.out'),'--error='+str(run/'logs/scheduler.err'),
                '--export=ALL,HICAT_RUN_DIR='+str(run)+',HICAT_PROJECT_ROOT='+str(project),
                str(job_file))
    (run/'provenance/job_id.txt').write_text(job+'\n')
    print(json.dumps(dict(run_dir=str(run),job_id=job,scope=cfg['scope'])),flush=True)


if __name__=='__main__':
    main()
