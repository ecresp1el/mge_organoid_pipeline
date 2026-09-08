#!/usr/bin/env python3
"""Freeze and submit the explicitly authorized Step 07 pilot-validation job.

No fitting occurs on the login node. Source/config/test/reference bytes and
environment are frozen before submission; scientific input hashes are checked
again on the compute node. Every invocation creates a new versioned run.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import shutil
import subprocess


def command(*args):
    """Run a checked argument-array command and return stripped stdout."""
    return subprocess.check_output(list(args), text=True).strip()


def main():
    """Freeze the complete executable package, then submit one bounded job."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--upstream', type=Path, required=True)
    parser.add_argument('--prepare-only', action='store_true')
    args = parser.parse_args()
    bundle = Path(__file__).resolve().parents[1]
    cfg = json.loads((bundle/'config/hicat_validation.json').read_text())
    if command('git', '-C', str(args.upstream), 'rev-parse', 'HEAD') != cfg['upstream_commit']:
        raise ValueError('Upstream commit mismatch')
    if command('git', '-C', str(args.upstream), 'status', '--porcelain'):
        raise ValueError('Pinned upstream checkout must be clean')
    interpreter = '/home/elcrespo/miniconda3/envs/pcdh19-hicat-py38/bin/python'
    command(interpreter, '-m', 'pip', 'check')
    project = Path('/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/paper3_pcdh19')
    stamp = datetime.now(timezone.utc).strftime('%Y%m%d_%H%M%S')
    digest = hashlib.sha256((bundle/'config/hicat_validation.json').read_bytes()).hexdigest()[:8]
    run = project/'results/hicat/02_hierarchy_validation'/('hicat_validation_'+stamp+'_'+digest)
    run.mkdir(parents=True)
    for name in ['code', 'config', 'inputs', 'logs', 'provenance', 'scratch']:
        (run/name).mkdir()
    shutil.copytree(bundle/'scripts/hicat', run/'code/hicat', ignore=shutil.ignore_patterns('__pycache__'))
    shutil.copytree(bundle/'tests', run/'code/tests', ignore=shutil.ignore_patterns('__pycache__'))
    shutil.copy2(__file__, run/'code/submit_hicat_validation.py')
    shutil.copy2(bundle/'slurm/hicat_validation.sbatch', run/'code/hicat_validation.sbatch')
    shutil.copytree(bundle/'config/hicat_validation_references', run/'config/hicat_validation_references')
    for name in ['hicat_validation.json', 'hicat_validation_programs.json',
                 'hicat_validation_hypotheses.json', 'hicat_coarse_fine.json',
                 'hicat_pilot.json', 'hicat_pilot.requirements.txt']:
        shutil.copy2(bundle/'config'/name, run/'config'/name)
    # Existing tests resolve ../config; preserve that exact fixture relationship.
    shutil.copytree(run/'config', run/'code/config')
    for name in ['HICAT_VALIDATION_PROTOCOL.md', 'EXPERIMENT_CONTEXT.md']:
        shutil.copy2(bundle/name, run/'config'/name)
    shutil.copytree(args.upstream, run/'upstream', ignore=shutil.ignore_patterns('.git', '__pycache__'))
    state = dict(commit=command('git', '-C', str(bundle.parent), 'rev-parse', 'HEAD'),
                 status=command('git', '-C', str(bundle.parent), 'status', '--short'),
                 upstream_commit=cfg['upstream_commit'])
    (run/'provenance/repository_state.json').write_text(json.dumps(state, indent=2)+'\n')
    (run/'provenance/repository_tracked_diff.patch').write_text(
        command('git', '-C', str(bundle.parent), 'diff', 'HEAD', '--binary')+'\n')
    (run/'config/python_pip_freeze.txt').write_text(command(interpreter, '-m', 'pip', 'freeze')+'\n')
    (run/'config/authorization.txt').write_text(
        'User explicitly requested Step 07 validation and then asked to pick it up and run it.\n'
        'Dissected E14.5 mouse MGE; same 12000 pilot cells; existing two seeds; one fixed-parent fine-DE sensitivity.\n'
        'Full 446349 cells: canonical scores and existing UMAP display only. No full-data clustering, label transfer, regression, batch correction or locked annotations.\n')
    records = []
    for folder in ['code', 'config', 'upstream']:
        for path in sorted((run/folder).rglob('*')):
            if path.is_file():
                records.append(dict(path=str(path.relative_to(run)), bytes=path.stat().st_size,
                                    sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
    (run/'provenance/frozen_source_manifest.json').write_text(json.dumps(records, indent=2)+'\n')
    if args.prepare_only:
        print(run)
        return
    job_file = project/'jobs'/(run.name+'.sbatch')
    if job_file.exists():
        raise FileExistsError(str(job_file))
    shutil.copy2(run/'code/hicat_validation.sbatch', job_file)
    job = command('sbatch', '--parsable', '--job-name=pcdh19-hicat-validation',
                  '--account=parent0', '--partition=standard', '--cpus-per-task=4',
                  '--mem=64G', '--time=04:00:00',
                  '--output='+str(run/'logs/scheduler.out'), '--error='+str(run/'logs/scheduler.err'),
                  '--export=ALL,HICAT_RUN_DIR='+str(run), str(job_file))
    (run/'provenance/job_id.txt').write_text(job+'\n')
    print(json.dumps(dict(run_dir=str(run), job_id=job, scope=cfg['scope'])), flush=True)


if __name__ == '__main__':
    main()
