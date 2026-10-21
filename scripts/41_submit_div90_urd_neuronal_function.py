#!/usr/bin/env python3
"""Prepare timestamped carbon copies and submit the independent functional extension."""
import argparse
import datetime as dt
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import time
from importlib.machinery import SourceFileLoader

os.environ['TZ']='America/Detroit';time.tzset()
ROOT=Path(__file__).resolve().parents[1]
PROJECT=Path('/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder')
PREVIOUS=PROJECT/'results/div90_urd_branchpoint_genes/div90_neuronal_BP11_BP9_single_window_20261021_160849_EDT_5baa285'
SOURCE=PROJECT/'results/div90_jia_lineage_urd/div90_allcells_jia_root10_neuron_s9_7tips_urd_resumable_v1/lineage_tree_jia_endpoint_tips_v1/div30_urd_lineage_tree_object.rds'
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
    return h.hexdigest()
def git(*args):return subprocess.check_output(['git']+list(args),cwd=str(ROOT)).decode().strip()
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--dry-run',action='store_true');ap.add_argument('--prepare-only',action='store_true');a=ap.parse_args()
    now=dt.datetime.now().astimezone();commit=git('rev-parse','HEAD')
    rid='div90_neuronal_function_{:%Y%m%d_%H%M%S_%Z}_{}'.format(now,commit[:7])
    run=PROJECT/'results/div90_urd_neuronal_function'/rid;job=PROJECT/'jobs'/(rid+'.sbatch')
    print('RUN',run,'\nJOB',job,flush=True)
    if a.dry_run:return
    if run.exists() or job.exists():raise FileExistsError('Run collision')
    for d in ['code','config','provenance','tables','plots']:(run/d).mkdir(parents=True,exist_ok=False)
    sources=['scripts/41_div90_urd_neuronal_function.R','scripts/41_div90_function_references.py',
      'scripts/41_submit_div90_urd_neuronal_function.py','slurm_templates/65_div90_urd_neuronal_function.sbatch.template',
      'python_notebooks/HANDOFF_div90_jia_lineage_urd_plan.md']
    copies=[]
    for s in sources:
        dest=run/('code' if s.startswith('scripts/') else 'provenance')/Path(s).name
        shutil.copy2(ROOT/s,dest);copies.append(dest)
    ref=SourceFileLoader('functional_refs',str(run/'code/41_div90_function_references.py')).load_module()
    # Avoid Python bytecode runtime artifacts in the checkout; frozen helper runs on NFS.
    ref.prepare(run/'config')
    copies+=list((run/'config').iterdir())
    protected=[SOURCE,SOURCE.parent/'tables/tree_tip_mapping.tsv',PREVIOUS/'DIV90_URD_gene_results.csv',
       PREVIOUS/'DIV90_URD_gene_summary.csv',PREVIOUS/'tables/comparison_cells.tsv.gz',PREVIOUS/'tables/window_support.tsv']
    (run/'provenance/source_sha256.txt').write_text(''.join('{}  {}\n'.format(sha(p),p) for p in protected))
    rendered=(ROOT/'slurm_templates/65_div90_urd_neuronal_function.sbatch.template').read_text()
    for key,value in [('__RUN_DIR__',run),('__SOURCE_RDS__',SOURCE),('__PREVIOUS__',PREVIOUS)]:rendered=rendered.replace(key,str(value))
    with job.open('x') as f:f.write(rendered)
    shutil.copy2(job,run/'provenance/submitted.sbatch');copies.append(run/'provenance/submitted.sbatch')
    subprocess.run(['bash','-n',str(job)],check=True)
    (run/'provenance/git_commit.txt').write_text(commit+'\n')
    (run/'provenance/git_status.txt').write_text(git('status','--short')+'\n')
    (run/'provenance/git_diff.patch').write_text(git('diff','--')+'\n')
    (run/'provenance/frozen_sha256.txt').write_text(''.join('{}  {}\n'.format(sha(p),p.relative_to(run)) for p in copies))
    manifest={'run_id':rid,'prepared_at':now.isoformat(),'prepared_at_utc':now.astimezone(dt.timezone.utc).isoformat(),
      'git_commit':commit,'previous_run':str(PREVIOUS),'tree':str(SOURCE),'command':['sbatch','--parsable',str(job)],
      'analysis_scope':'reuse six existing DE contrasts; no DE rerun; functional expression and enrichment downstream',
      'log_template':str(PROJECT/'logs/65_div90_urd_neuronal_function_%j.log')}
    target=run/'provenance/submission.json';target.write_text(json.dumps(manifest,indent=2)+'\n')
    if a.prepare_only:print('PREPARED_ONLY',run,flush=True);return
    result=subprocess.run(manifest['command'],stdout=subprocess.PIPE,stderr=subprocess.PIPE,universal_newlines=True)
    manifest.update(submitted_at=dt.datetime.now().astimezone().isoformat(),submission_exit_code=result.returncode,stdout=result.stdout,stderr=result.stderr)
    if result.returncode==0:
        jid=result.stdout.strip().split(';')[0];assert jid.isdigit();manifest['job_id']=jid
        (run/'provenance/job_id.txt').write_text(jid+'\n');print('SUBMITTED',jid,flush=True)
    else:(run/'SUBMISSION_FAILED.txt').write_text(result.stderr)
    target.write_text(json.dumps(manifest,indent=2)+'\n');result.check_returncode()
if __name__=='__main__':main()
