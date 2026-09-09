#!/usr/bin/env python3
"""Submit the authorized frozen 100-iteration HiCAT array without changing science.

Only orchestration and the separate authorization record are written. The
frozen worker, seeds, sampled IDs, expression and parameters remain unchanged.
Iterations 0/1 reuse the two existing benchmark pipelines. Eight new fits can
start immediately; remaining pending tasks are released when those existing
pipelines terminate, keeping at most ten active 80% iterations during overlap.
The array itself always has throttle 10. Final aggregation and DE are separate
jobs and can be submitted again from their saved inputs after individual retry.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import sys

PYTHON='/home/elcrespo/miniconda3/envs/pcdh19-hicat-py38/bin/python'


def digest(path):
    """Hash a closed artifact without loading it into memory."""
    value=hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda:handle.read(8*1024*1024),b''):value.update(block)
    return value.hexdigest()


def write_json(path,value):
    """Publish operational metadata atomically, outside frozen inputs."""
    path=Path(path);temporary=path.with_name(path.name+'.tmp')
    temporary.write_text(json.dumps(value,indent=2)+'\n');temporary.replace(path)


def batch_text(root,command,cpus=8):
    """Use the exact production interpreter, import paths and thread controls."""
    return ('#!/bin/bash\nset -Eeuo pipefail\n'+
        'export PYTHONPATH='+shlex.quote(str(root/'adapter')+':'+str(root/'reference/allen_python'))+'\n'+
        'export OPENBLAS_NUM_THREADS='+str(cpus)+' OMP_NUM_THREADS='+str(cpus)+' NUMBA_NUM_THREADS=1\n'+
        'export MPLBACKEND=Agg PYTHONDONTWRITEBYTECODE=1\n'+
        'export NUMBA_CACHE_DIR='+shlex.quote(str(root/'logs/numba_cache'))+'\n'+command+'\n')


def submit_script(ops,script,role,extra=(),cpus=8,memory='150G',wall='2-00:00:00'):
    """Submit one recorded scheduler action; never batch all fits in a process."""
    args=['sbatch','--parsable','--account=parent0','--partition=standard',
          '--cpus-per-task='+str(cpus),'--mem='+memory,'--time='+wall,
          '--job-name=pcdh19-production-'+role,
          '--output='+str(ops/'logs'/(role+'_%A_%a.out')),
          '--error='+str(ops/'logs'/(role+'_%A_%a.err'))]+list(extra)+[str(script)]
    job=subprocess.check_output(args,text=True).strip().split(';')[0]
    with (ops/'scheduler_actions.jsonl').open('a') as handle:
        handle.write(json.dumps(dict(action='submit',job_id=job,role=role,args=args,
                                    script_sha256=digest(script)))+'\n')
    return job


def main():
    """Validate the frozen contract, record authorization and submit once."""
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run-dir',type=Path,required=True)
    args=parser.parse_args();root=args.run_dir.resolve()
    pointer=root/'PRODUCTION_SUBMISSION.json'
    if pointer.exists():
        raise RuntimeError('Production already has a submission record; use individual retries, not another 100-task array')
    cfg=json.loads((root/'config/run.json').read_text())
    frozen=json.loads((root/'FROZEN.json').read_text())
    for name,expected in frozen.items():
        if digest(root/name)!=expected:raise ValueError('Frozen artifact changed: '+name)
    if cfg['production_iterations']!=100:raise ValueError('Expected 100 frozen iterations')
    threshold=cfg['fitting']['thresholds']
    if (threshold['q1_thresh'],threshold['qdiff_thresh'],threshold['score_thresh'])!=(.4,.7,150):
        raise ValueError('Frozen scientific thresholds differ from authorization')
    # The already-qualified worker independently verifies all seals/environment
    # again on each task. This preflight records the interpreter itself too.
    env=subprocess.check_output([PYTHON,'-m','pip','freeze'],text=True)
    if env!=(root/'config/environment.txt').read_text():raise ValueError('Scientific environment changed')
    jobs=json.loads((root/'job_ids.json').read_text())
    legacy=Path(cfg['legacy_benchmark'])
    old=json.loads((legacy/'job_ids.json').read_text())
    stamp=datetime.now(timezone.utc).strftime('%Y%m%d_%H%M%S')
    ops=root/'production'/('submission_'+stamp);(ops/'logs').mkdir(parents=True)
    shutil.copy2(__file__,ops/'submit_consensus_production.py')
    shutil.copy2(root/'AUTHORIZATION.json',ops/'authorization_before.json')
    write_json(ops/'preflight.json',dict(frozen_files_verified=len(frozen),
        frozen_manifest_sha256=digest(root/'FROZEN.json'),adapter_commit=cfg['adapter_commit'],
        all_settings_unchanged=True,all_100_frozen_sample_sets_reused=True,
        python_version=subprocess.check_output([PYTHON,'--version'],text=True).strip(),
        python_executable=PYTHON,python_executable_sha256=digest(Path(PYTHON).resolve()),
        package_freeze_sha256=digest(root/'config/environment.txt')))
    write_json(root/'AUTHORIZATION.json',dict(production_allowed=True,
        approved_utc=datetime.now(timezone.utc).isoformat(),approved_iterations=list(range(100)),
        concurrency_cap=10,cpus_per_task=8,memory_GiB_per_task=150,
        scientific_settings_must_remain_frozen=True,bigcat_production_allowed=False,
        authority='User explicitly authorizes 100 x 80% production now and supersedes the earlier benchmark-report approval gate',
        previously_running_jobs_must_not_be_changed_for_resource_optimization=True))
    worker=shlex.join([PYTHON,'-m','hicat.consensus_restart','--run-dir',str(root),'--production'])
    array_script=ops/'iterations.sbatch'
    array_script.write_text(batch_text(root,'/usr/bin/time -v '+worker+
        ' --stage iteration --iteration "${SLURM_ARRAY_TASK_ID:?}"'))
    # Hold before dependency edits so an array task cannot race an existing
    # writer. No already-running benchmark job is modified or interrupted.
    array=submit_script(ops,array_script,'iteration',['--array=0-99%10','--hold'])
    record=dict(status='CONFIGURING_ARRAY',run_dir=str(root),operations_dir=str(ops),
        array_job_id=array,array_spec='0-99%10',concurrency_cap=10,iteration_tasks=100,
        cpus_per_task=8,memory_GiB_per_task=150,time_limit='2-00:00:00',
        adapter_commit=cfg['adapter_commit'],bigcat_jobs_submitted=0,
        preexisting_iteration_jobs={0:old['fit80'],1:jobs['fit1']},
        preexisting_membership_jobs={0:jobs['map0'],1:jobs['map1_retry']})
    write_json(pointer,record)
    for index,dependency in [(0,jobs['map0']),(1,jobs['map1_retry'])]:
        subprocess.run(['scontrol','update','JobId='+array+'_'+str(index),
                        'Dependency=afterany:'+dependency],check=True)
    # Release the tail after preexisting iteration pipelines are terminal,
    # including failure, so one failure never blocks all unrelated iterations.
    # This changes only pending holds, never resources or numerical settings.
    gate=ops/'release_remaining.sbatch'
    gate.write_text('#!/bin/bash\nset -Eeuo pipefail\n'+
        'for task in $(seq 10 99); do scontrol release "'+array+'_${task}"; done\n')
    release=submit_script(ops,gate,'release-pending',
        ['--dependency=afterany:'+jobs['map0']+':'+jobs['map1_retry']],cpus=1,memory='1G',wall='00:30:00')
    for stage in ['aggregate','merge']:
        (ops/(stage+'.sbatch')).write_text(batch_text(root,'/usr/bin/time -v '+worker+' --stage '+stage))
    # Actual 100-block aggregation additionally waits for the real two-block
    # algorithm test. An array failure does not recompute successful iterations:
    # aggregate validates all 100 seals and fails clearly if any remain absent.
    aggregation=submit_script(ops,ops/'aggregate.sbatch','aggregate',
        ['--dependency=afterany:'+array+',afterok:'+old['fullfit']+':'+jobs['report']])
    merging=submit_script(ops,ops/'merge.sbatch','merge',['--dependency=afterok:'+aggregation])
    record.update(status='SUBMITTED',release_pending_job_id=release,
                  aggregation_job_id=aggregation,final_DE_job_id=merging,
                  submitted_utc=datetime.now(timezone.utc).isoformat())
    write_json(pointer,record);write_json(ops/'submission.json',record)
    # All array elements exist now. Eight can perform new fits alongside the
    # two existing fits; 0/1 wait to verify/reuse their own saved checkpoints.
    for index in range(10):
        subprocess.run(['scontrol','release',array+'_'+str(index)],check=True)
    write_json(root/'RUN_STATE.json',dict(status='PRODUCTION_ARRAY_SUBMITTED',
        production_submission='PRODUCTION_SUBMISSION.json',array_job_id=array,
        production_allowed=True,benchmark_validation_continues=True))
    print(json.dumps(record,indent=2))


if __name__=='__main__':main()
