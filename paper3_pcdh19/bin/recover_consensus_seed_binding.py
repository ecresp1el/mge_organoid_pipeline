#!/usr/bin/env python3
"""Recover the initial array's Annoy binding error without changing frozen seeds.

This is a one-time operational recovery, not a new scientific ensemble. Preserve
all active fits and failed attempts; cancel only pending superseded tasks. The
replacement array verifies/reuses every completed checkpoint. It waits for the
original iteration writers, avoiding concurrent writes and keeping the total
number of active fitting iterations at or below ten. Final stages stay separate.
"""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import shlex
import shutil
import subprocess

from submit_consensus_production import PYTHON, batch_text, digest, submit_script, write_json


def main():
    """Require qualification, preserve prior records, and submit exactly once."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run-dir', type=Path, required=True)
    root = parser.parse_args().run_dir.resolve()
    old = json.loads((root / 'PRODUCTION_SUBMISSION.json').read_text())
    if old.get('seed_binding_recovery'):
        raise RuntimeError('Recovery already submitted; use individual task retries')
    for name, value in json.loads((root / 'FROZEN.json').read_text()).items():
        if digest(root / name) != value:
            raise ValueError('Frozen input changed: ' + name)
    if subprocess.check_output([PYTHON, '-m', 'pip', 'freeze'], text=True) != (root / 'config/environment.txt').read_text():
        raise ValueError('Scientific environment changed')
    compatibility = root / 'production/annoy_seed64_v1'
    qualification = json.loads((root / 'production/annoy_seed64_qualification.json').read_text())
    checkpoint_test = json.loads((root / 'production/annoy_seed64_checkpoint_qualification.json').read_text())
    if qualification['status'] != 'PASS' or qualification['frozen_seeds_tested'] != 100 or checkpoint_test['status'] != 'PASS':
        raise ValueError('Exact seed / checkpoint qualification missing')
    if qualification['compatibility_manifest_sha256'] != digest(compatibility / 'MANIFEST.json'):
        raise ValueError('Qualified compatibility manifest changed')
    for name, value in json.loads((compatibility / 'MANIFEST.json').read_text())['files'].items():
        if digest(compatibility / name) != value:
            raise ValueError('Qualified runtime changed: ' + name)
    ops = root / 'production' / ('seed_recovery_' + datetime.now(timezone.utc).strftime('%Y%m%d_%H%M%S'))
    (ops / 'logs').mkdir(parents=True)
    shutil.copy2(__file__, ops / Path(__file__).name)
    shutil.copy2(root / 'PRODUCTION_SUBMISSION.json', ops / 'previous_submission.json')
    for name in ['annoy_seed64_qualification.json', 'annoy_seed64_checkpoint_qualification.json']:
        shutil.copy2(root / 'production' / name, ops / name)
    live = subprocess.check_output(['squeue', '-r', '-h', '-j', old['array_job_id'], '-o', '%i|%T'], text=True)
    (ops / 'previous_array_before_recovery.tsv').write_text(live)
    account = subprocess.check_output(['sacct', '-X', '-n', '-P', '-j', old['array_job_id'], '--format=JobID%40,State%30,Elapsed,ExitCode'], text=True)
    (ops / 'previous_array_accounting.tsv').write_text(account)
    active = [line.split('|')[0] for line in live.splitlines() if line.split('|')[1] in ['RUNNING', 'COMPLETING']]
    jobs = json.loads((root / 'job_ids.json').read_text())
    legacy = Path(json.loads((root / 'config/run.json').read_text())['legacy_benchmark'])
    legacy_jobs = json.loads((legacy / 'job_ids.json').read_text())
    # SLURM applies this state filter atomically: a task that has begun running
    # is never cancelled by this command. Existing resources remain untouched.
    cancelled = [old['array_job_id'], old['release_pending_job_id'], old['aggregation_job_id'], old['final_DE_job_id'], '60615631', '60615750']
    for job in cancelled:
        subprocess.run(['scancel', '--state=PENDING', job], check=True)
    worker = shlex.join([PYTHON, str(compatibility / 'worker.py'), '--run-dir', str(root), '--production'])
    script = ops / 'iterations.sbatch'
    script.write_text(batch_text(root, '/usr/bin/time -v ' + worker + ' --stage iteration --iteration "${SLURM_ARRAY_TASK_ID:?}"'))
    dependency = 'afterany:' + ':'.join([old['array_job_id'], jobs['map0'], jobs['map1_retry']])
    array = submit_script(ops, script, 'iteration', ['--array=0-99%10', '--dependency=' + dependency])
    record = dict(old)
    record.update(status='SUBMITTED_SEED_BINDING_RECOVERY', array_job_id=array,
        operations_dir=str(ops), previous_array_job_id=old['array_job_id'],
        previous_array_tasks_preserved_running=active,
        seed_binding_recovery=True, compatibility_manifest=str(compatibility / 'MANIFEST.json'),
        compatibility_manifest_sha256=digest(compatibility / 'MANIFEST.json'),
        frozen_seeds_and_sample_lists_unchanged=True,
        completed_checkpoints_reused=True, replacement_dependency=dependency,
        initial_failed_iteration_ids=[2, 3, 4, 5, 6, 8, 9],
        recovery_reason='Annoy signed-int32 Python binding rejected exact frozen uint32 seeds; qualified bridge calls same native uint64 setter',
        submitted_utc=datetime.now(timezone.utc).isoformat(), release_pending_job_id=None)
    # Publish immediately so a subsequent failure never causes blind duplicate
    # submission; scheduler_actions.jsonl identifies every accepted job.
    write_json(root / 'PRODUCTION_SUBMISSION.json', record)
    scientific = shlex.join([PYTHON, '-m', 'hicat.consensus_restart', '--run-dir', str(root), '--production'])
    for stage in ['aggregate', 'merge']:
        (ops / (stage + '.sbatch')).write_text(batch_text(root, '/usr/bin/time -v ' + scientific + ' --stage ' + stage))
    aggregate = submit_script(ops, ops / 'aggregate.sbatch', 'aggregate',
        ['--dependency=afterany:' + array + ',afterok:' + legacy_jobs['fullfit'] + ':' + jobs['report']])
    merge = submit_script(ops, ops / 'merge.sbatch', 'merge', ['--dependency=afterok:' + aggregate])
    status = ops / 'consensus_production_status.py'
    shutil.copy2(Path(__file__).with_name(status.name), status)
    monitoring = {}
    for role, dep, extra in [('array-status', 'afterany:' + array, ''),
                             ('report-status', 'afterok:' + jobs['report'], ' --reconcile-benchmark')]:
        target = ops / (role + '.sbatch')
        target.write_text(batch_text(root, shlex.join([PYTHON, str(status), '--run-dir', str(root)]) + extra, cpus=1))
        monitoring[role] = submit_script(ops, target, role, ['--dependency=' + dep], cpus=1, memory='1G', wall='00:30:00')
    record.update(aggregation_job_id=aggregate, final_DE_job_id=merge, monitoring_jobs=monitoring)
    write_json(root / 'PRODUCTION_SUBMISSION.json', record)
    write_json(ops / 'submission.json', record)
    write_json(root / 'RUN_STATE.json', dict(status=record['status'], production_allowed=True,
        array_job_id=array, preserved_running_tasks=active, bigcat_jobs_submitted=0,
        production_submission='PRODUCTION_SUBMISSION.json'))
    print(json.dumps(record, indent=2))


if __name__ == '__main__':
    main()
