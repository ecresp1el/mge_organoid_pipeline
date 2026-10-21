#!/usr/bin/env python3
"""Submit the bounded activity-program extension using NFS carbon copies."""
import datetime as dt
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import time

os.environ['TZ']='America/Detroit';time.tzset()
ROOT=Path(__file__).resolve().parents[1]
PROJECT=Path('/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder')
EARLY=PROJECT/'results/div90_urd_branchpoint_genes/div90_neuronal_BP11_BP9_single_window_20261021_160849_EDT_5baa285'
FUNCTION=PROJECT/'results/div90_urd_neuronal_function/div90_neuronal_function_20261021_163824_EDT_5baa285'
TREE=PROJECT/'results/div90_jia_lineage_urd/div90_allcells_jia_root10_neuron_s9_7tips_urd_resumable_v1/lineage_tree_jia_endpoint_tips_v1/div30_urd_lineage_tree_object.rds'
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(8388608),b''):h.update(b)
 return h.hexdigest()
def main():
 now=dt.datetime.now().astimezone();commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=str(ROOT)).decode().strip()
 rid='div90_activity_programs_{:%Y%m%d_%H%M%S_%Z}_{}'.format(now,commit[:7]);run=PROJECT/'results/div90_urd_neuronal_function'/rid
 job=PROJECT/'jobs'/(rid+'.sbatch')
 if run.exists() or job.exists():raise FileExistsError('Run collision')
 for d in ['code','config','provenance','tables','plots']:(run/d).mkdir(parents=True,exist_ok=False)
 copies=[]
 for name in ['scripts/42_div90_activity_programs.R','scripts/42_submit_div90_activity_programs.py','python_notebooks/HANDOFF_div90_jia_lineage_urd_plan.md','slurm_templates/66_div90_activity_programs.sbatch.template']:
  dest=run/('code' if name.startswith('scripts/')else 'provenance')/Path(name).name;shutil.copy2(ROOT/name,dest);copies.append(dest)
 # Execute the exact native helpers used by the completed neuronal DE run.
 dest=run/'code/40_div90_urd_branchpoint_genes.R';shutil.copy2(EARLY/'code/40_div90_urd_branchpoint_genes.R',dest);copies.append(dest)
 config={'early_run':str(EARLY),'function_run':str(FUNCTION),'tree':str(TREE),'near_split_FDR':'unchanged',
  'later_window_fractions':[.05,.10],'minimum_cells':30,'minimum_bin_cells':5,'minimum_joint_bins':3,'matching_bins':5,
  'minimum_detection':.10,'minimum_log2fc':.25,'maximum_FDR':.05,'AUCPR_gate':False,'detection_difference_gate':False,
  'maximum_later_comparisons':6,'preflight_supported_pairs':['BP11 5/3','BP9 1/2'],'IEG_coexpression_primary':False}
 dest=run/'config/resolved_settings.json';dest.write_text(json.dumps(config,indent=2)+'\n');copies.append(dest)
 refs=[FUNCTION/'config'/n for n in ['go-basic.obo','goa_human.gaf.gz','functional_gene_sets.tsv','reference_manifest.json']]
 protected=[TREE,EARLY/'DIV90_URD_gene_results.csv',EARLY/'code/40_div90_urd_branchpoint_genes.R']+refs
 (run/'provenance/source_sha256.txt').write_text(''.join('{}  {}\n'.format(sha(p),p)for p in protected))
 text=(ROOT/'slurm_templates/66_div90_activity_programs.sbatch.template').read_text()
 for key,val in [('__RUN__',run),('__EARLY__',EARLY),('__FUNCTION__',FUNCTION),('__TREE__',TREE)]:text=text.replace(key,str(val))
 with job.open('x')as f:f.write(text)
 dest=run/'provenance/submitted.sbatch';shutil.copy2(job,dest);copies.append(dest)
 subprocess.run(['bash','-n',str(job)],check=True)
 for n,cmd in [('git_status.txt',['git','status','--short']),('git_diff.patch',['git','diff','--'])]:
  (run/'provenance'/n).write_bytes(subprocess.check_output(cmd,cwd=str(ROOT)))
 (run/'provenance/git_commit.txt').write_text(commit+'\n')
 (run/'provenance/frozen_sha256.txt').write_text(''.join('{}  {}\n'.format(sha(p),p.relative_to(run))for p in copies))
 manifest={'run_id':rid,'run_dir':str(run),'prepared_at':now.isoformat(),'prepared_at_utc':now.astimezone(dt.timezone.utc).isoformat(),
  'git_commit':commit,'command':['sbatch','--parsable',str(job)],'settings':config,'reference_reuse':'Frozen references protected by hashes; no downloads or source writes'}
 p=run/'provenance/submission.json';p.write_text(json.dumps(manifest,indent=2)+'\n')
 r=subprocess.run(manifest['command'],stdout=subprocess.PIPE,stderr=subprocess.PIPE,universal_newlines=True)
 manifest.update(submitted_at=dt.datetime.now().astimezone().isoformat(),submission_exit_code=r.returncode,stdout=r.stdout,stderr=r.stderr)
 if r.returncode==0:
  jid=r.stdout.strip().split(';')[0];assert jid.isdigit();manifest['job_id']=jid;(run/'provenance/job_id.txt').write_text(jid+'\n')
 else:(run/'SUBMISSION_FAILED.txt').write_text(r.stderr)
 p.write_text(json.dumps(manifest,indent=2)+'\n');r.check_returncode();print('SUBMITTED',manifest['job_id'],'\nRUN',run,flush=True)
if __name__=='__main__':main()
