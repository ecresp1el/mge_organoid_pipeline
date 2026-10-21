#!/usr/bin/env python3
"""Freeze the reviewed DIV90 downstream implementation and submit its NFS job.

--dry-run prints locations without creating files or submitting. Runtime artifacts
are always on Turbo; existing source runs and result packages are never replaced.
"""
import argparse
import csv
import datetime as dt
import hashlib
import io
import json
import os
from pathlib import Path
import shutil
import subprocess
import urllib.request
import time

os.environ['TZ'] = 'America/Detroit'
time.tzset()

ROOT = Path(__file__).resolve().parents[1]
PROJECT = Path('/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder')
SOURCE = PROJECT / 'results/div90_jia_lineage_urd/div90_allcells_jia_root10_neuron_s9_7tips_urd_resumable_v1/lineage_tree_jia_endpoint_tips_v1/div30_urd_lineage_tree_object.rds'
REFERENCE_URL = 'https://humantfs.ccbr.utoronto.ca/download/v_1.01/DatabaseExtract_v_1.01.csv'


def command(*args):
    return subprocess.check_output(args, cwd=str(ROOT), universal_newlines=True).strip()


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(8 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--dry-run', action='store_true')
    args = parser.parse_args()
    now = dt.datetime.now().astimezone()
    commit = command('git', 'rev-parse', 'HEAD')
    with (ROOT / 'config/div90_urd_branchpoint_genes.tsv').open() as f:
        settings = {r['parameter']: r['value'] for r in csv.DictReader(f, delimiter='\t')}
    prefix = ('div90_neuronal_BP11_BP9_single_window'
              if settings.get('neuronal_comparisons') == '1'
              else 'div90_BP12_OPC_vs_neurons_single_window')
    run_id = f'{prefix}_{now:%Y%m%d_%H%M%S_%Z}_{commit[:7]}'
    run = PROJECT / 'results/div90_urd_branchpoint_genes' / run_id
    job = PROJECT / 'jobs' / f'{run_id}.sbatch'
    sources = {
        'scripts/40_div90_urd_branchpoint_genes.R': 'code/40_div90_urd_branchpoint_genes.R',
        'scripts/40_submit_div90_urd_branchpoint_genes.py': 'code/40_submit_div90_urd_branchpoint_genes.py',
        'config/div90_urd_branchpoint_genes.tsv': 'config/div90_urd_branchpoint_genes.tsv',
        'slurm_templates/64_div90_urd_branchpoint_genes.sbatch.template': 'provenance/job_template.sbatch',
        'python_notebooks/HANDOFF_div90_jia_lineage_urd_plan.md': 'provenance/HANDOFF_div90_jia_lineage_urd_plan.md',
    }
    for path in sources:
        if not (ROOT / path).is_file():
            raise FileNotFoundError(path)
    if not SOURCE.is_file():
        raise FileNotFoundError(SOURCE)
    if run.exists() or job.exists():
        raise FileExistsError('Refusing timestamp/run collision')
    print(f'Run: {run}\nJob: {job}\nProtected source: {SOURCE}', flush=True)
    if args.dry_run:
        print('Dry run: no files created, no reference downloaded, no job submitted.')
        return
    for folder in ('code', 'config', 'provenance', 'tables', 'plots', 'checkpoints'):
        (run / folder).mkdir(parents=True, exist_ok=False)
    for src, dest in sources.items():
        shutil.copy2(ROOT / src, run / dest)
    ref = run / 'config/DatabaseExtract_v_1.01.csv'
    data = urllib.request.urlopen(REFERENCE_URL, timeout=60).read()
    records = list(csv.DictReader(io.StringIO(data.decode('utf-8-sig'))))
    if len(records) < 1000 or not {'HGNC symbol', 'Is TF?', 'DBD'}.issubset(records[0]):
        raise ValueError('Unexpected Lambert TF reference schema or row count')
    reference_hash = hashlib.sha256(data).hexdigest()
    reference_dir = PROJECT / 'reference' / 'human_tfs_lambert_v1.01' / reference_hash[:12]
    reference_dir.mkdir(parents=True, exist_ok=True)
    reference_file = reference_dir / ref.name
    if reference_file.exists():
        if sha(reference_file) != reference_hash:
            raise ValueError('Existing reference hash mismatch')
    else:
        with reference_file.open('xb') as f:
            f.write(data)
    shutil.copy2(reference_file, ref)
    source_hash = sha(SOURCE)
    mapping = SOURCE.parent / 'tables/tree_tip_mapping.tsv'
    # Analysis only reads this source-side mapping; hash it alongside the RDS.
    hashes = [(SOURCE, source_hash), (mapping, sha(mapping))]
    (run / 'provenance/source_sha256.txt').write_text(''.join(f'{h}  {p}\n' for p, h in hashes))
    text = (run / 'provenance/job_template.sbatch').read_text()
    text = text.replace('__RUN_DIR__', str(run)).replace('__SOURCE_RDS__', str(SOURCE))
    with job.open('x') as f:
        f.write(text)
    job.chmod(0o750)
    shutil.copy2(job, run / 'provenance/submitted.sbatch')
    subprocess.run(['bash', '-n', str(job)], check=True)
    (run / 'provenance/git_commit.txt').write_text(commit + '\n')
    (run / 'provenance/git_status.txt').write_text(command('git', 'status', '--short') + '\n')
    # git diff omits untracked implementation files, so the frozen files/hashes
    # remain authoritative; the diff captures changes to tracked documentation.
    (run / 'provenance/git_diff.patch').write_text(command('git', 'diff', '--', *sources) + '\n')
    freeze = [run / dest for dest in sources.values()] + [ref, run / 'provenance/submitted.sbatch']
    (run / 'provenance/frozen_sha256.txt').write_text(''.join(f'{sha(p)}  {p.relative_to(run)}\n' for p in freeze))
    manifest = {
        'run_id': run_id, 'run_dir': str(run), 'prepared_at': now.isoformat(),
        'prepared_at_utc': now.astimezone(dt.timezone.utc).isoformat(),
        'timezone': 'America/Detroit', 'git_commit': commit, 'source_rds': str(SOURCE),
        'source_sha256': source_hash, 'submitted_job': str(job),
        'reference_url': REFERENCE_URL, 'reference_sha256': reference_hash,
        'reference_retrieved_at': dt.datetime.now(dt.timezone.utc).isoformat(),
        'reference_rows': len(records), 'reference_tf_rows': sum(x['Is TF?'] == 'Yes' for x in records),
        'reference_data_path': str(reference_file), 'command': ['sbatch', '--parsable', str(job)],
        'scheduler_log_template': str(PROJECT / 'logs/64_div90_urd_branchpoint_genes_%j.log'),
    }
    manifest_path = run / 'provenance/submission.json'
    manifest_path.write_text(json.dumps(manifest, indent=2) + '\n')
    submitted = subprocess.run(manifest['command'], stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    manifest.update(submitted_at=dt.datetime.now().astimezone().isoformat(),
                    sbatch_stdout=submitted.stdout, sbatch_stderr=submitted.stderr,
                    submission_exit_code=submitted.returncode)
    if submitted.returncode == 0:
        job_id = submitted.stdout.strip().split(';')[0]
        if not job_id.isdigit():
            raise ValueError(f'Unexpected sbatch job ID: {submitted.stdout}')
        manifest['job_id'] = job_id
        (run / 'provenance/job_id.txt').write_text(job_id + '\n')
        print(f'SUBMITTED job_id={job_id}\nMASTER_CSV={run}/DIV90_URD_gene_results.csv', flush=True)
    else:
        (run / 'SUBMISSION_FAILED.txt').write_text(submitted.stderr)
    manifest_path.write_text(json.dumps(manifest, indent=2) + '\n')
    submitted.check_returncode()


if __name__ == '__main__':
    main()
