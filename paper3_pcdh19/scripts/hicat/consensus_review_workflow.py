"""Publish an existing consensus as a versioned IN_REVIEW HiCAT stage-03 package."""
import argparse
import json
from pathlib import Path
import os
import pandas as pd
from .consensus_review_io import ConsensusReviewSettings,ConsensusReviewInputs
from .consensus_review_report import ConsensusReviewReport
from .provenance import Progress,write_json,manifest,sha256


class ConsensusReviewWorkflow:
    """Own verification, numbered reporting and atomic publication; no fitting."""

    def __init__(self, root):
        """Initialize frozen settings and a single-writer staging directory."""
        self.settings=ConsensusReviewSettings.read(root)
        self.root=self.settings.root
        self.output=self.root/'staging/outputs'
        self.output.mkdir(parents=True)
        (self.output/'tables').mkdir();(self.output/'validation').mkdir()
        self.progress=Progress(self.root/'provenance')

    def run(self):
        """Validate source checkpoints, render figures, verify copies and publish."""
        with self.progress.track('consensus_review.verify_saved_checkpoints'):
            data=ConsensusReviewInputs(self.settings).read(self.output)
        with self.progress.track('consensus_review.numbered_report'):
            ConsensusReviewReport(self.output,self.settings.config).render(data)
        with self.progress.track('consensus_review.validate_and_publish'):
            for name in ['all_cell_assignments.tsv.gz','pairwise_DE.tsv','cluster_mean_log1p_cpm.tsv.gz','cluster_detection.tsv.gz']:
                data.check('published_copy:'+name,sha256(data.final/name)==sha256(self.output/'tables'/name))
            figures=pd.read_csv(self.output/'figures/figure_index.tsv',sep='\t')
            data.check('contiguous_figure_numbering',figures.page.to_list()==list(range(1,len(figures)+1)))
            for row in figures.itertuples():
                data.check('figure:'+row.png,(self.output/'figures'/row.png).is_file())
                for source in row.source_tables.split(';'):data.check('figure_source:'+source,(self.output/source).is_file())
            status=dict(stage=self.settings.config['stage'],run_id=self.root.name,status='IN_REVIEW',
                        iterations=self.settings.config['iterations'],cells=len(data.assignments),genes=data.means.shape[1],
                        clusters=len(data.means),residual_DE_pairs=data.scope['residual_DE_pairs'],
                        source_consensus=str(data.endpoint),refitted=False,annotations_locked=False,
                        source_compute_job='60707938',publication_job=os.environ.get('SLURM_JOB_ID'))
            write_json(self.output/'STEP_STATUS.json',status)
            pd.DataFrame(data.checks).to_csv(self.output/'validation/validation_checks.tsv',sep='\t',index=False)
            readme=('# HiCAT 03: 98-iteration consensus review\n\n'
                    '**IN_REVIEW — 34 clusters across 446,349 dissected E14.5 mouse MGE cells.**\n\n'
                    'This versioned package reviews the saved consensus of iterations 0–97, denominator 98. '
                    'No iteration, aggregation, refinement or DE fit was rerun during publication. '
                    'The source computation is preserved at '+str(data.endpoint)+'.\n\n'
                    'Open [the review PDF](figures/hicat_consensus_review.pdf). '
                    'Numbered PNGs 001–003 and their precise source tables are indexed in '
                    '`figures/figure_index.tsv`. Cell assignments are in `tables/all_cell_assignments.tsv.gz`; '
                    '`diagnostic_cluster` is the inherited worker column containing final-DE consensus IDs.\n\n'
                    'Aggregation/refinement produced 35 groups, final DE produced 34. Three of 561 audited '
                    'pairs do not meet the configured DE separation criteria. Cluster IDs remain unannotated.\n\n'
                    'The report has cluster sizes/sample fractions, a display-only canonical marker panel '
                    'and an all-gene cluster-mean expression dendrogram. It does not contain a new UMAP, '
                    'a new AnnData, ranked top-20 markers or an accepted taxonomy. The dendrogram is '
                    'not the saved recursive fitting tree.\n\n'
                    'All source aggregation/DE artifacts were hashed and reopened before publication; '
                    'cell order, gene order, sample count, source hashes, denominator and parent seals '
                    'were checked. The source computation used focused membership-input checks; this '
                    'publication does not retroactively claim unused fitting models were reopened. '
                    'See `tables/source_checkpoints.json` and `validation/validation_checks.tsv`.\n\n'
                    'Frozen code/config are in the parent run. Canonical results live on Turbo under '
                    'HiCAT stage 03; this is a same-stage review, not primary-processing Step 03 or a new Step 08.\n')
            (self.output/'REVIEW_README.md').write_text(readme)
            manifest(self.output).to_csv(self.output/'output_manifest.tsv',sep='\t',index=False)
            for row in pd.read_csv(self.output/'output_manifest.tsv',sep='\t').itertuples():
                path=self.output/row.relative_path
                if path.stat().st_size!=row.bytes or sha256(path)!=row.sha256:raise ValueError('Published manifest mismatch')
            self.output.rename(self.root/'outputs')
            (self.root/'COMPUTATION_SUCCESS.txt').write_text('COMPLETED\nstatus=IN_REVIEW\n')
        print(json.dumps(status),flush=True)


def main():
    """Run only the selected frozen review package."""
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--run-dir',type=Path,required=True)
    args=parser.parse_args();ConsensusReviewWorkflow(args.run_dir).run()

if __name__=='__main__':main()
