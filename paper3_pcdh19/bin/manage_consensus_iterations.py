#!/usr/bin/env python3
"""Freeze 100 independent seeds/inputs; submit only the authorized two-fit test.

Actions: create-benchmark, status, submit-one. No monolithic fitting job exists.
Production permission is a separate operational file; changing permission later
does not invalidate scientific checkpoint contracts or recompute completed fits.
Do not set production_allowed before the user approves the benchmark report.
"""
import argparse
from datetime import datetime,timezone
import hashlib
import json
from pathlib import Path
import shlex
import shutil
import subprocess
import sys
import numpy as np
import pandas as pd
import h5py

PYTHON='/home/elcrespo/miniconda3/envs/pcdh19-hicat-py38/bin/python'


def sha(path):
    """Stream the exact artifact hash used in the frozen input contract."""
    digest=hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda:handle.read(8*1024*1024),b''):digest.update(block)
    return digest.hexdigest()


def submit(root,stage,index=None,dependency=None,inject=False,afterany=False):
    """Submit one deterministic job; additional production IDs require approval."""
    root=Path(root);cfg=json.loads((root/'config/run.json').read_text())
    authorization=json.loads((root/'AUTHORIZATION.json').read_text())
    if index is not None and index not in cfg['benchmark_iteration_ids'] and not authorization['production_allowed']:
        raise ValueError('Production iterations require explicit user approval')
    if stage=='report':cpus,memory,wall=1,'16G','04:00:00'
    elif stage=='mapping':cpus,memory,wall=8,'64G','08:00:00'
    else:cpus,memory,wall=8,'150G','2-00:00:00'
    label=stage+('-%03d'%index if index is not None else '')+('-injected' if inject else '')
    serial=datetime.now(timezone.utc).strftime('%Y%m%d_%H%M%S_%f')
    script=root/'jobs'/(label+'_'+serial+'.sbatch')
    if stage=='report':command=[PYTHON,'-m','hicat.restart_report','--run-dir',str(root)]
    else:
        command=[PYTHON,'-m','hicat.consensus_restart','--run-dir',str(root),'--stage',stage]
        if index is not None:command+=['--iteration',str(index)]
        if inject:command+=['--inject-failure']
    script.write_text('#!/bin/bash\nset -Eeuo pipefail\n'+
        'export PYTHONPATH='+shlex.quote(str(root/'adapter')+':'+str(root/'reference/allen_python'))+'\n'+
        'export OPENBLAS_NUM_THREADS='+str(cpus)+' OMP_NUM_THREADS='+str(cpus)+' NUMBA_NUM_THREADS=1\n'+
        'export MPLBACKEND=Agg PYTHONDONTWRITEBYTECODE=1\n'+
        'export NUMBA_CACHE_DIR='+shlex.quote(str(root/'logs/numba_cache'))+'\n'+
        '/usr/bin/time -v '+shlex.join(command)+'\n')
    args=['sbatch','--parsable','--account=parent0','--partition=standard','--cpus-per-task='+str(cpus),
          '--mem='+memory,'--time='+wall,'--job-name=pcdh19-resume-'+label,
          '--output='+str(root/'logs'/(label+'_'+serial+'.out')),
          '--error='+str(root/'logs'/(label+'_'+serial+'.err'))]
    if dependency:args+=['--dependency='+('afterany:' if afterany else 'afterok:')+dependency,'--kill-on-invalid-dep=yes']
    args.append(str(script));job=subprocess.check_output(args,text=True).strip().split(';')[0]
    with (root/'submitted_jobs.jsonl').open('a') as handle:
        handle.write(json.dumps(dict(stage=stage,iteration=index,injected_failure=inject,job_id=job,
            dependency=dependency,script=str(script),script_sha256=sha(script),cpus=cpus,memory=memory))+'\n')
    return job


def create(args):
    """Freeze all 100 seeds and ID sets, qualify snapshot, then launch two only."""
    bundle=Path(__file__).resolve().parents[1];legacy=args.legacy_benchmark
    for log in [args.qualification_log,args.restart_test_log]:
        if not log.read_text().rstrip().endswith('OK'):raise ValueError('Completed passing qualification logs required')
    identity=json.loads((legacy/'inputs/input_identity.json').read_text())
    stamp=datetime.now(timezone.utc).strftime('%Y%m%d_%H%M%S')
    root=legacy.parent/('consensus_restart_'+stamp);root.mkdir()
    for folder in ['adapter','config','reference','iterations','jobs','logs','qualification','job_metrics','final']:(root/folder).mkdir()
    fitting=json.loads((legacy/'config/benchmark.json').read_text())
    shutil.copytree(bundle/'scripts/hicat',root/'adapter/hicat',ignore=shutil.ignore_patterns('__pycache__'))
    shutil.copy2(__file__,root/'adapter/manage_consensus_iterations.py')
    for name in ['test_hicat_consensus.py','test_hicat_restart.py']:
        (root/'adapter/tests').mkdir(exist_ok=True);shutil.copy2(bundle/'tests'/name,root/'adapter/tests'/name)
    (root/'adapter/config').mkdir();shutil.copy2(bundle/'config/hicat_coarse_fine.json',root/'adapter/config/hicat_coarse_fine.json')
    shutil.copytree(bundle/'references/allen_consensus_audit_20260908',root/'adapter/references/allen_consensus_audit_20260908')
    shutil.copytree(legacy/'reference/allen_python',root/'reference/allen_python')
    for name in ['HICAT_ALLEN_CONSENSUS_ALGORITHM_AUDIT.md','HICAT_ALLEN_CONSENSUS_TECHNICAL_NOTES.md',
                 'FULL_DATA_CONSENSUS_BENCHMARK.md']:
        shutil.copy2(bundle/name,root/'config'/name)
    # A local snapshot commit binds the exact adapter independently of unrelated
    # edits/autocommits in the shared project worktree. Nothing is pushed.
    subprocess.run(['git','init','-q',str(root/'adapter')],check=True)
    subprocess.run(['git','-C',str(root/'adapter'),'add','.'],check=True)
    subprocess.run(['git','-C',str(root/'adapter'),'-c','user.name=Codex','-c','user.email=codex@local',
                    'commit','-q','-m','Freeze restart-safe HiCAT benchmark adapter and oracle tests'],check=True)
    commit=subprocess.check_output(['git','-C',str(root/'adapter'),'rev-parse','HEAD'],text=True).strip()
    universe=pd.read_csv(legacy/'inputs/full_cell_index.tsv.gz',sep='\t')
    universe.to_csv(root/'config/cell_universe.tsv.gz',sep='\t',index=False)
    with h5py.File(legacy/'inputs/full_log1p_cpm.h5ad','r') as f:
        node=f['var'];genes=node[node.attrs['_index']].asstr()[:]
    np.save(root/'config/gene_ids.npy',genes.astype(str))
    universe_hash=sha(root/'config/cell_universe.tsv.gz')
    # Iteration 0 preserves the already-running authorized benchmark. Remaining
    # streams use separate SeedSequence children; sampling and fit seeds differ.
    children=np.random.SeedSequence(20260908).spawn(99)
    seeds=[dict(iteration=0,sampling_seed=20260908,fit_seed=20260908)]
    for index,child in enumerate(children,1):
        state=child.generate_state(2);seeds.append(dict(iteration=index,sampling_seed=int(state[0]),fit_seed=int(state[1])))
    if len({row['sampling_seed'] for row in seeds})!=100:raise ValueError('Sampling seed collision')
    pd.DataFrame(seeds).to_csv(root/'config/seeds.tsv',sep='\t',index=False)
    cfg=dict(schema='restart_safe_consensus_v1',legacy_benchmark=str(legacy),production_iterations=100,
             benchmark_iteration_ids=[0,1],adapter_commit=commit,fitting=fitting,universe_sha256=universe_hash,
             normalized_sha256=identity['normalized_sha256'],consensus_format='per-iteration CSC B and integer labels in frozen universe order',
             same_seed_and_inputs_required_for_retry=True,independent_iteration_jobs=True,
             source_worktree_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=bundle,text=True).strip())
    (root/'config/run.json').write_text(json.dumps(cfg,indent=2)+'\n')
    (root/'config/environment.txt').write_text(subprocess.check_output([PYTHON,'-m','pip','freeze'],text=True))
    shutil.copy2(args.qualification_log,root/'qualification/numerical_checks.log')
    shutil.copy2(args.restart_test_log,root/'qualification/restart_checks.log')
    for record in seeds:
        directory=root/'iterations'/('%03d'%record['iteration']);directory.mkdir()
        rows=np.sort(np.random.default_rng(record['sampling_seed']).choice(len(universe),357079,replace=False)).astype(np.int32)
        np.save(directory/'sampled_source_rows.npy',rows)
        universe.iloc[rows][['cell_id','source_row']].to_csv(directory/'sampled_cell_ids.tsv.gz',sep='\t',index=False)
        settings=dict(record,universe_sha256=universe_hash,adapter_commit=commit,settings=fitting,
                      sampled_source_rows_sha256=sha(directory/'sampled_source_rows.npy'),
                      sampled_cell_ids_sha256=sha(directory/'sampled_cell_ids.tsv.gz'))
        (directory/'config.json').write_text(json.dumps(settings,indent=2)+'\n')
    frozen={}
    for base in ['adapter','config','reference','qualification','iterations']:
        for path in (root/base).rglob('*'):
            if path.is_file() and '.git' not in path.parts:frozen[str(path.relative_to(root))]=sha(path)
    (root/'FROZEN.json').write_text(json.dumps(frozen,indent=2)+'\n')
    (root/'AUTHORIZATION.json').write_text(json.dumps(dict(production_allowed=False,
        approved_benchmark_iterations=[0,1],must_stop_after='completed resource/restart report for user approval'),indent=2)+'\n')
    # Freeze the new adapter before any new execution. The first fit is adopted
    # after its legacy job completes; completed work is not recomputed.
    old=json.loads((legacy/'job_ids.json').read_text())
    jobs={}
    jobs['adopt_fit0']=submit(root,'fit',0,old['fit80'])
    jobs['fit1']=submit(root,'fit',1,old['prepare'])
    jobs['map0']=submit(root,'mapping',0,jobs['adopt_fit0'])
    jobs['map1_injected_failure']=submit(root,'mapping',1,jobs['fit1'],inject=True)
    jobs['map1_retry']=submit(root,'mapping',1,jobs['map1_injected_failure'],afterany=True)
    jobs['aggregate']=submit(root,'aggregate',dependency=':'.join([jobs['map0'],jobs['map1_retry'],old['fullfit']]))
    jobs['final_DE']=submit(root,'merge',dependency=jobs['aggregate'])
    jobs['report']=submit(root,'report',dependency=jobs['final_DE'])
    (root/'job_ids.json').write_text(json.dumps(jobs,indent=2)+'\n')
    (root/'RUN_STATE.json').write_text(json.dumps(dict(status='TWO_REAL_ITERATION_BENCHMARK_SUBMITTED',jobs=jobs,
        production_jobs=0,adapter_commit=commit,seed_sets_frozen=100),indent=2)+'\n')
    print(json.dumps(dict(run_dir=str(root),jobs=jobs,adapter_commit=commit,production_jobs=0),indent=2))


def status(root):
    """Detect each missing, failed, active or sealed stage independently."""
    rows=[]
    for index in range(100):
        for stage in ['fit','mapping']:
            directory=root/'iterations'/('%03d'%index)/stage
            if (directory/'CURRENT.json').exists():state='COMPLETE_SEALED_REVERIFY_ON_REUSE'
            elif (directory/'RUNNING.json').exists():state='RUNNING_OR_STALE_LOCK_CHECK_SCHEDULER'
            elif list(directory.glob('attempts/*/FAILURE.json')):state='FAILED_INDIVIDUALLY_RESUBMITTABLE'
            else:state='MISSING_NOT_SUBMITTED_OR_WAITING'
            rows.append(dict(iteration=index,stage=stage,status=state))
    print(pd.DataFrame(rows).to_csv(sep='\t',index=False))


def main():
    """Create the bounded benchmark or inspect/resubmit one independent stage."""
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action',choices=['create-benchmark','status','submit-one'])
    parser.add_argument('--run-dir',type=Path)
    parser.add_argument('--legacy-benchmark',type=Path)
    parser.add_argument('--qualification-log',type=Path)
    parser.add_argument('--restart-test-log',type=Path)
    parser.add_argument('--stage',choices=['iteration','fit','mapping','aggregate','merge','report'])
    parser.add_argument('--iteration',type=int)
    args=parser.parse_args()
    if args.action=='create-benchmark':create(args)
    elif args.action=='status':status(args.run_dir)
    else:print(submit(args.run_dir,args.stage,args.iteration))


if __name__=='__main__':main()
