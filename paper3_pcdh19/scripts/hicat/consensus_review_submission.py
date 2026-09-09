"""Freeze and submit a same-stage consensus review with established Turbo layout."""
import argparse
from datetime import datetime,timezone
import hashlib
import json
from pathlib import Path
import shutil
import subprocess


class ConsensusReviewSubmission:
    """Prepare versioned HiCAT-03 code/config, job and logs before submission."""

    def __init__(self, bundle):
        """Resolve canonical source, existing project root and review configuration."""
        self.bundle=Path(bundle)
        self.project=Path('/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/paper3_pcdh19')
        self.cfg=json.loads((self.bundle/'config/hicat_consensus_review.json').read_text())

    @staticmethod
    def digest(path):
        """Hash input bytes without retaining large files in memory."""
        h=hashlib.sha256()
        with Path(path).open('rb') as handle:
            for block in iter(lambda:handle.read(8*1024*1024),b''):h.update(block)
        return h.hexdigest()

    def prepare(self):
        """Freeze code and input identity into a new numbered-stage package."""
        stamp=datetime.now(timezone.utc).strftime('%Y%m%d_%H%M%S')
        config=self.bundle/'config/hicat_consensus_review.json'
        self.root=self.project/'results/hicat'/self.cfg['stage']/('hicat_consensus_review_'+stamp+'_'+self.digest(config)[:8])
        self.root.mkdir(parents=True)
        for folder in ['code','config','inputs','provenance','scratch']:(self.root/folder).mkdir()
        shutil.copytree(self.bundle/'scripts/hicat',self.root/'code/hicat',ignore=shutil.ignore_patterns('__pycache__'))
        shutil.copy2(self.bundle/'bin/submit_hicat_consensus_review.py',self.root/'code/submit_hicat_consensus_review.py')
        shutil.copy2(self.bundle/'slurm/hicat_consensus_review.sbatch',self.root/'code/hicat_consensus_review.sbatch')
        shutil.copy2(config,self.root/'config'/config.name)
        for name in ['HICAT_CONSENSUS_REVIEW_PROTOCOL.md','HICAT_CONSENSUS_RESTART_PROTOCOL.md','EXPERIMENT_CONTEXT.md']:
            shutil.copy2(self.bundle/name,self.root/'config'/name)
        source=Path(self.cfg['source_run'])
        shutil.copytree(source/'reference/allen_python',self.root/'upstream',ignore=shutil.ignore_patterns('.git','__pycache__'))
        shutil.copy2(source/'config/environment.txt',self.root/'config/python_pip_freeze.txt')
        relatives=['config/run.json','config/cell_universe.tsv.gz','config/gene_ids.npy',
                   self.cfg['source_consensus']+'/SUMMARY.json',self.cfg['source_consensus']+'/aggregation/CURRENT.json',
                   self.cfg['source_consensus']+'/DE/CURRENT.json']
        identity=dict(source_run=str(source),files={name:self.digest(source/name) for name in relatives})
        (self.root/'inputs/source_identity.json').write_text(json.dumps(identity,indent=2)+'\n')
        (self.root/'config/authorization.txt').write_text(
            'User requested the 98-iteration consensus and then directed following the existing Markdown workflow and numbering. '
            'This package validates and republishes saved results within HiCAT stage 03; no refit, new major step or approval promotion.\n')
        state=subprocess.check_output(['git','-C',str(self.bundle),'status','--short'],text=True)
        commit=subprocess.check_output(['git','-C',str(self.bundle),'rev-parse','HEAD'],text=True).strip()
        (self.root/'provenance/repository_state.json').write_text(json.dumps(dict(commit=commit,status=state),indent=2)+'\n')
        diff=subprocess.check_output(['git','-C',str(self.bundle),'diff','HEAD','--binary'],text=True)
        (self.root/'provenance/repository_tracked_diff.patch').write_text(diff)
        rows=[]
        for folder in ['code','config','upstream','inputs']:
            for path in sorted((self.root/folder).rglob('*')):
                if path.is_file():rows.append(dict(path=str(path.relative_to(self.root)),bytes=path.stat().st_size,sha256=self.digest(path)))
        (self.root/'provenance/frozen_source_manifest.json').write_text(json.dumps(rows,indent=2)+'\n')
        log=self.project/'logs/hicat'/self.root.name;log.mkdir(parents=True)
        (self.root/'logs').symlink_to(log,target_is_directory=True)
        return self.root

    def submit(self):
        """Submit the frozen SLURM wrapper via the canonical jobs directory."""
        jobfile=self.project/'jobs'/(self.root.name+'.sbatch')
        shutil.copy2(self.root/'code/hicat_consensus_review.sbatch',jobfile)
        job=subprocess.check_output(['sbatch','--parsable','--account=parent0','--partition=standard',
             '--cpus-per-task=2','--mem=16G','--time=01:00:00','--job-name=pcdh19-hicat03-review',
             '--output='+str(self.root/'logs/scheduler.out'),'--error='+str(self.root/'logs/scheduler.err'),
             '--export=ALL,HICAT_RUN_DIR='+str(self.root),str(jobfile)],text=True).strip()
        (self.root/'provenance/job_id.txt').write_text(job+'\n')
        return dict(run_dir=str(self.root),job_id=job,status='SUBMITTED',stage=self.cfg['stage'])


def main(bundle=None):
    """Prepare or submit the single configured existing-consensus review."""
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--prepare-only',action='store_true')
    args=parser.parse_args();submission=ConsensusReviewSubmission(bundle or Path(__file__).resolve().parents[2])
    root=submission.prepare()
    print(json.dumps(dict(run_dir=str(root)) if args.prepare_only else submission.submit()),flush=True)
