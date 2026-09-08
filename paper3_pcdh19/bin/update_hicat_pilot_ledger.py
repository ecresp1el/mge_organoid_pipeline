#!/usr/bin/env python3
"""Reconcile explicit HiCAT pilot success/failure records into a review ledger.

Only terminal records are included. Existing human APPROVED/REJECTED decisions
are preserved. A successful computation is inserted as IN_REVIEW, never approved.
This post-run operation reads scientific assets without rewriting run packages.
"""
import argparse
import csv
import json
import os
from pathlib import Path


class PilotStatusLedger:
    """Maintain one human-readable table of versioned pilot dispositions."""

    fields=['stage','run_id','status','job_id','scope','source_step02_run_id',
            'cells','genes','output_checkpoint','approval_decision','notes']

    def __init__(self,workflow_root):
        """Resolve the HiCAT workflow root, not a scientific output directory."""
        self.root=Path(workflow_root)
        self.path=self.root/'APPROVAL_LEDGER.tsv'

    def update(self):
        """Merge terminal run statuses and atomically replace the ledger.

        Returns
        -------
        pathlib.Path
            Updated TSV. No active or pending run is inferred to have failed.
        """
        rows={}
        if self.path.exists():
            with self.path.open() as handle:
                rows={row['run_id']:row for row in csv.DictReader(handle,delimiter='\t')}
        for run in sorted((self.root/'00_method_and_technical_pilot').glob('hicat_pilot_*')):
            if rows.get(run.name,{}).get('status') in ('APPROVED','REJECTED'):
                continue
            statusfile=run/'outputs/STEP_STATUS.json'
            failure=run/'COMPUTATION_FAILED.txt'
            if not statusfile.exists() and not failure.exists():continue
            cfg=json.loads((run/'config/hicat_pilot.json').read_text())
            success=statusfile.exists()
            identityfile=run/'inputs/input_identity.json'
            shape=json.loads(identityfile.read_text())['shape'] if identityfile.exists() else ['','']
            row={key:'' for key in self.fields}
            row.update(stage='00_method_and_technical_pilot',run_id=run.name,
                       status='IN_REVIEW' if success else 'FAILED',
                       job_id=(run/'provenance/job_id.txt').read_text().strip(),
                       scope=cfg['scope'],source_step02_run_id=cfg['step02_run_id'],
                       cells=shape[0],genes=shape[1],
                       output_checkpoint=str(run/'outputs/pcdh19_hicat_pilot.h5ad') if success else '',
                       notes='Technical pilot only; no full-data K or annotation.' if success else 'See immutable failure marker and scheduler logs; no output checkpoint published.')
            rows[run.name]=row
        temporary=self.path.with_suffix('.tsv.tmp')
        with temporary.open('w',newline='') as handle:
            writer=csv.DictWriter(handle,fieldnames=self.fields,delimiter='\t')
            writer.writeheader()
            writer.writerows(rows[key] for key in sorted(rows))
        os.replace(temporary,self.path)
        return self.path


def main():
    """Parse the explicit workflow root and reconcile only recorded terminal runs."""
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--workflow-root',type=Path,required=True)
    args=parser.parse_args()
    print(PilotStatusLedger(args.workflow_root).update())


if __name__=='__main__':
    main()
