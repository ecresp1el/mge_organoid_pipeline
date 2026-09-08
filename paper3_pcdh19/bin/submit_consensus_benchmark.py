#!/usr/bin/env python3
"""Freeze and submit one real full-data HiCAT benchmark, never production.

Only a fixed five-stage dependency chain is supported. There is deliberately
no array-size or production flag. Scientific inputs come from approved Step02;
all code/config/reference hashes are verified before every stage executes.
The final stage writes a measured resource report and stops IN_REVIEW.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import io
import json
from pathlib import Path
import shlex
import shutil
import subprocess
import tarfile


def sha(path):
    """Stream a file checksum used to freeze code and qualification evidence."""
    digest=hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda:handle.read(8*1024*1024),b''):digest.update(block)
    return digest.hexdigest()


def main():
    """Require passing checks, snapshot assets, then submit only benchmark jobs."""
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--qualification-log',type=Path,required=True)
    args=parser.parse_args()
    if not args.qualification_log.read_text().rstrip().endswith('OK'):
        raise ValueError('Passing completed qualification log required')
    bundle=Path(__file__).resolve().parents[1]
    base=Path('/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/paper3_pcdh19')
    stamp=datetime.now(timezone.utc).strftime('%Y%m%d_%H%M%S')
    run=base/'results/hicat/03_full_data_consensus_benchmark'/('consensus_benchmark_'+stamp)
    run.mkdir(parents=True)
    for key in ['config','inputs','code','reference','logs','progress','fits','qualification']:(run/key).mkdir()
    cfg=json.loads((bundle/'config/hicat_coarse_fine.json').read_text())
    cfg['thresholds']=dict(cfg.pop('fine_thresholds'),q1_thresh=.4,qdiff_thresh=.7,score_thresh=150)
    for key in ['coarse_thresholds','pilot_cells_per_sample','runs']:cfg.pop(key,None)
    cfg.update(schema_version='full_data_benchmark_v1',scope='one_real_80percent_fit_then_measured_consensus_branch',
               source_path=str(base/'results/primary_processing/02_qc_filtering'/cfg['step02_run_id']/'objects/pcdh19_step02_qc_filtered.h5ad'),
               sampling_seed=20260908,fit_seed=20260908,block_cells=4096,
               pca_fit_cells='all_cells_at_each_node',production_iterations=100,executed_subsample_iterations=1,
               production_launch_authorized=False,requires_user_approval_after_report=True,
               consensus_reference_commit='9af2f04cb837a7b53b005dee88a92eae337afead',
               heldout_assignment='mean_prototype_Pearson',denominator='all_successful_iterations_after_mapping',
               dense_full_pair_matrix=False,fixed_pilot_parents=False,target_K=None,
               qualified_consensus_branch='large_graph; stop if actual first-iteration G selects unqualified smaller branch')
    (run/'config/benchmark.json').write_text(json.dumps(cfg,indent=2)+'\n')
    (run/'config/LATEST_USER_AUTHORIZATION.txt').write_text(
        'Definitive target: 100 x 80% of full 446349 cells, q1=.4 qdiff=.7 DEscore150. '
        'User explicitly instructed: do not launch 100 production iterations until completed audit and real full-data '
        'benchmark are reported; stop for approval after measured runtime/RSS/CPU/disk/resources/concurrency report.\n')
    shutil.copytree(bundle/'scripts/hicat',run/'code/hicat',ignore=shutil.ignore_patterns('__pycache__'))
    shutil.copy2(__file__,run/'code/submit_consensus_benchmark.py')
    shutil.copy2(bundle/'tests/test_hicat_consensus.py',run/'qualification/test_hicat_consensus.py')
    shutil.copy2(args.qualification_log,run/'qualification/passed_checks.log')
    for name in ['HICAT_ALLEN_CONSENSUS_ALGORITHM_AUDIT.md','HICAT_ALLEN_CONSENSUS_TECHNICAL_NOTES.md',
                 'FULL_DATA_CONSENSUS_BENCHMARK.md']:
        shutil.copy2(bundle/name,run/'config'/name)
    shutil.copytree(bundle/'references/allen_consensus_audit_20260908',run/'reference/audit')
    repository=Path('/tmp/pcdh19_hicat_validation_reference')
    archive=subprocess.check_output(['git','archive',cfg['upstream_commit']],cwd=repository)
    destination=run/'reference/allen_python';destination.mkdir()
    with tarfile.open(fileobj=io.BytesIO(archive)) as tar:
        for member in tar.getmembers():
            if member.name.startswith('/') or '..' in Path(member.name).parts or member.issym() or member.islnk():
                raise ValueError('Unsafe reference archive member')
        tar.extractall(destination)
    python='/home/elcrespo/miniconda3/envs/pcdh19-hicat-py38/bin/python'
    (run/'config/environment.txt').write_text(subprocess.check_output([python,'-m','pip','freeze'],text=True))
    resources=dict(prepare=(1,'32G','02:00:00'),fit80=(8,'150G','2-00:00:00'),
                   fullfit=(8,'150G','2-00:00:00'),consensus=(8,'150G','08:00:00'),report=(1,'8G','01:00:00'))
    for stage,(cpus,memory,wall) in resources.items():
        command=[python,'-m','hicat.consensus_benchmark','--run-dir',str(run),'--stage',stage]
        jobfile=run/'code'/('run_'+stage+'.sbatch')
        jobfile.write_text('#!/bin/bash\nset -Eeuo pipefail\n'+
            'export OPENBLAS_NUM_THREADS='+str(cpus)+' OMP_NUM_THREADS='+str(cpus)+' NUMBA_NUM_THREADS=1\n'+
            'export MPLBACKEND=Agg PYTHONDONTWRITEBYTECODE=1\n'+
            'export PYTHONPATH='+shlex.quote(str(run/'code')+':'+str(destination))+'\n'+
            'export NUMBA_CACHE_DIR='+shlex.quote(str(run/'logs/numba_cache'))+'\n'+
            '/usr/bin/time -v '+shlex.join(command)+'\n')
    frozen={str(p.relative_to(run)):sha(p) for directory in ['code','config','reference','qualification']
            for p in (run/directory).rglob('*') if p.is_file()}
    (run/'frozen_manifest.json').write_text(json.dumps(frozen,indent=2)+'\n')
    jobs={};previous=None
    # Sequential dependency chain: no simultaneous fit or array submission.
    for stage,(cpus,memory,wall) in resources.items():
        command=['sbatch','--parsable','--account=parent0','--partition=standard',
                 '--cpus-per-task='+str(cpus),'--mem='+memory,'--time='+wall,
                 '--job-name=pcdh19-bench-'+stage,'--output='+str(run/'logs'/(stage+'.out')),
                 '--error='+str(run/'logs'/(stage+'.err'))]
        if previous:command+=['--dependency=afterok:'+previous,'--kill-on-invalid-dep=yes']
        command+=[str(run/'code'/('run_'+stage+'.sbatch'))]
        previous=subprocess.check_output(command,text=True).strip().split(';')[0];jobs[stage]=previous
        (run/'job_ids.json').write_text(json.dumps(jobs,indent=2)+'\n')
    (run/'RUN_STATE.json').write_text(json.dumps(dict(status='BENCHMARK_SUBMITTED',jobs=jobs,production_jobs=0,
        stopped_for_approval_after='outputs/RESOURCE_REPORT.md'),indent=2)+'\n')
    print(json.dumps(dict(run_dir=str(run),jobs=jobs,production_jobs=0),indent=2))


if __name__=='__main__':main()
