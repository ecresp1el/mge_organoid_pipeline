#!/usr/bin/env python3
"""Read all 100 task/checkpoint states and save an operational status snapshot.

This command never fits, merges, submits, cancels or requeues anything. Task
state comes from SLURM; checkpoint state checks the published manifest seal.
The frozen scientific worker separately reopens and hashes every artifact on
reuse. Failed array elements are listed with individual fresh-job retry commands.

--reconcile-benchmark updates only the benchmark report's obsolete approval
wording after it runs, preserving its original report and measured results.
The user authorized production before benchmark completion; scientific code
and completed iteration checkpoints must not be edited to reflect that change.
"""
import argparse
from collections import Counter
from datetime import datetime,timezone
import hashlib
import json
from pathlib import Path
import re
import shlex
import subprocess


def read_json(path):
    """Read an existing JSON record without changing it."""
    return json.loads(Path(path).read_text())


def publish(path,data):
    """Atomically replace an explicitly mutable operational snapshot."""
    path=Path(path);tmp=path.with_name(path.name+'.tmp')
    tmp.write_text(json.dumps(data,indent=2)+'\n');tmp.replace(path)


def sealed(stage):
    """Check the completion manifest seal, leaving full asset checks to worker."""
    try:
        current=read_json(stage/'CURRENT.json');attempt=stage/current['attempt']
        complete=read_json(attempt/'COMPLETE.json')
        actual=hashlib.sha256((attempt/'ARTIFACTS.json').read_bytes()).hexdigest()
        return actual==current['manifest_sha256']==complete['manifest_sha256']
    except (OSError,ValueError,KeyError):return False


def snapshot(root):
    """Combine scheduler states with each independently published checkpoint."""
    submission=read_json(root/'PRODUCTION_SUBMISSION.json');array=submission['array_job_id']
    states={};pattern=re.compile(r'^'+re.escape(array)+r'_(\d+)$')
    account=subprocess.check_output(['sacct','-X','-n','-P','-j',array,
                                    '--format=JobID%40,State%30'],text=True)
    for line in account.splitlines():
        fields=line.split('|');match=pattern.match(fields[0].strip())
        if match:states[int(match[1])]=fields[1].strip().split()[0].rstrip('+')
    live=subprocess.check_output(['squeue','-r','-h','-j',array,'-o','%i|%T'],text=True)
    for line in live.splitlines():
        job,state=line.split('|');match=pattern.match(job.strip())
        if match:states[int(match[1])]=state.strip()
    failure={'FAILED','CANCELLED','TIMEOUT','OUT_OF_MEMORY','NODE_FAIL','BOOT_FAIL','DEADLINE','PREEMPTED'}
    def retry_command(indices):
        """Fresh job IDs let the frozen worker identify stale killed-job locks.

        Queue one array of failed indices behind the current array, preserving
        the global concurrency cap. Successful iteration indices are omitted.
        Do not use scontrol requeue: a SIGKILL can leave a lock with that same
        job ID, which would appear active when the job is requeued.
        """
        ops=Path(submission['operations_dir'])
        return shlex.join(['sbatch','--parsable','--account=parent0','--partition=standard',
            '--cpus-per-task=8','--mem=150G','--time=2-00:00:00',
            '--job-name=pcdh19-production-iteration-retry',
            '--array='+','.join(map(str,indices))+'%10','--dependency=afterany:'+array,
            '--output='+str(ops/'logs/retry_%A_%a.out'),
            '--error='+str(ops/'logs/retry_%A_%a.err'),str(ops/'iterations.sbatch')])
    # A compatibility recovery can preserve active writers from the previous
    # array. Report these explicitly instead of mistaking 100 replacement tasks
    # waiting on them for 100 idle scientific iterations.
    previous_states={};previous=submission.get('previous_array_job_id')
    if previous:
        old_pattern=re.compile(r'^'+re.escape(previous)+r'_(\d+)$')
        old_account=subprocess.check_output(['sacct','-X','-n','-P','-j',previous,
            '--format=JobID%40,State%30'],text=True)
        for line in old_account.splitlines():
            fields=line.split('|');match=old_pattern.match(fields[0].strip())
            if match:previous_states[int(match[1])]=fields[1].strip().split()[0].rstrip('+')
    preexisting_live={}
    for key in ['preexisting_iteration_jobs','preexisting_membership_jobs']:
        for index,job in submission.get(key,{}).items():
            query=subprocess.run(['squeue','-h','-j',job,'-o','%T'],
                                 text=True,capture_output=True)
            # Finished jobs eventually leave squeue; sacct retains accounting.
            # Do not suppress unrelated scheduler/connection errors.
            if query.returncode and 'Invalid job id specified' not in query.stderr:
                query.check_returncode()
            value=query.stdout.strip() if query.returncode==0 else ''
            if value in ['RUNNING','COMPLETING']:preexisting_live[int(index)]=job
    rows=[]
    for index in range(100):
        directory=root/'iterations'/('%03d'%index)
        state=states.get(index,'UNKNOWN_ACCOUNTING_NOT_AVAILABLE')
        fit=sealed(directory/'fit');mapping=sealed(directory/'mapping')
        logical=('COMPLETED' if fit and mapping else 'RUNNING' if
            state in ['RUNNING','COMPLETING'] or previous_states.get(index) in ['RUNNING','COMPLETING']
            or index in preexisting_live else 'FAILED' if state in failure else 'QUEUED')
        rows.append(dict(iteration=index,array_task=array+'_'+str(index),scheduler_state=state,
            logical_iteration_state=logical,previous_array_state=previous_states.get(index,''),
            preexisting_active_job=preexisting_live.get(index,''),
            fit_manifest_sealed=fit,mapping_manifest_sealed=mapping,
            complete_iteration=fit and mapping,
            retry_command=retry_command([index]) if state in failure and not (fit and mapping) else '',
            missing_fit=not fit,missing_membership=not mapping))
    counts=Counter(r['scheduler_state'] for r in rows)
    result=dict(utc=datetime.now(timezone.utc).isoformat(),array_job_id=array,
        concurrency_cap=10,states=dict(counts),
        queued=sum(counts[s] for s in ['PENDING','CONFIGURING','REQUEUED','REQUEUE_HOLD']),
        running=sum(counts[s] for s in ['RUNNING','COMPLETING']),
        completed=counts['COMPLETED'],failed=sum(counts[s] for s in failure),
        sealed_iterations=sum(r['complete_iteration'] for r in rows),
        final_aggregation_sealed=sealed(root/'final/production/aggregation'),
        final_DE_sealed=sealed(root/'final/production/DE'),
        full_artifact_verification='Each frozen worker reopens/hashes every input artifact before reuse; this status scan checks manifest seals only',
        preexisting_iteration_jobs=submission['preexisting_iteration_jobs'],
        production_authorized=read_json(root/'AUTHORIZATION.json')['production_allowed'],
        bigcat_production_jobs_submitted=0,iterations=rows)
    result.update(logical_iteration_counts=dict(Counter(r['logical_iteration_state'] for r in rows)),
        previous_array_job_id=previous,initial_failed_iteration_ids=submission.get('initial_failed_iteration_ids',[]),
        failed_attempts_preserved_and_queued_for_retry=bool(submission.get('seed_binding_recovery')))
    failed_indices=[r['iteration'] for r in rows if r['retry_command']]
    result.update(recommended_single_retry_array_command=retry_command(failed_indices) if failed_indices else '',
        retry_note='Use one fresh array of failed indices, capped at 10, after the current array; then resubmit final aggregation/DE with new dependencies. No same-job-ID requeue after a killed stage.')
    folder=root/'production';folder.mkdir(exist_ok=True)
    publish(folder/'STATUS.json',result)
    columns=list(rows[0])
    (folder/'iteration_status.tsv').write_text('\t'.join(columns)+'\n'+
        ''.join('\t'.join(str(row[c]) for c in columns)+'\n' for row in rows))
    return result


def reconcile(root,status):
    """Preserve benchmark measurements while reflecting the newer approval."""
    report=root/'report';source=report/'RESOURCE_REPORT.json'
    if not source.exists():raise ValueError('Benchmark resource report has not completed')
    archived=report/'original_benchmark_approval_metadata';archived.mkdir(exist_ok=True)
    for name in ['RESOURCE_REPORT.json','RESOURCE_REPORT.md']:
        target=archived/name
        if not target.exists():target.write_bytes((report/name).read_bytes())
    data=read_json(source)
    data.update(status='BENCHMARK_COMPLETE_PRODUCTION_AUTHORIZED',
        production_array_job_id=status['array_job_id'],
        production_iterations_executed=status['sealed_iterations'],
        production_progress_snapshot_utc=status['utc'],
        production_iteration_count_is_a_snapshot=True,
        approval_note='User superseded the earlier stop-before-production condition; scientific benchmark measurements unchanged')
    publish(source,data)
    note=('**Production was explicitly authorized before this benchmark finished. '
          'Array '+status['array_job_id']+' is the current production submission. '
          'This report measures the two-iteration benchmark; consult production/STATUS.json '
          'for production progress. Original approval metadata is archived.**\n\n')
    text=(archived/'RESOURCE_REPORT.md').read_text()
    text=text.replace('**Stop here for approval. No 100-run production submission.**',
                      '**Benchmark validation completed; production is already authorized.**')
    (report/'RESOURCE_REPORT.md').write_text(note+text)
    publish(root/'RUN_STATE.json',dict(status='BENCHMARK_COMPLETE_PRODUCTION_AUTHORIZED',
        array_job_id=status['array_job_id'],production_submission='PRODUCTION_SUBMISSION.json',
        production_status='production/STATUS.json',production_allowed=True))


def main():
    """Print a compact summary and save the complete per-iteration status table."""
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run-dir',type=Path,required=True)
    parser.add_argument('--reconcile-benchmark',action='store_true')
    args=parser.parse_args();status=snapshot(args.run_dir)
    if args.reconcile_benchmark:reconcile(args.run_dir,status)
    print(json.dumps({k:v for k,v in status.items() if k!='iterations'},indent=2))


if __name__=='__main__':main()
