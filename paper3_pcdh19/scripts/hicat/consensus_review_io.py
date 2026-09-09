"""Settings and verified readers for a saved consensus review within HiCAT 03.

This module reads immutable completed aggregation/DE checkpoints. It never fits,
merges, edits the source run, or changes the denominator of an existing result.
"""
from dataclasses import dataclass
import json
from pathlib import Path
import shutil
import h5py
import numpy as np
import pandas as pd
from .checkpoints import complete_attempt
from .provenance import sha256, write_json


@dataclass(frozen=True)
class ConsensusReviewSettings:
    """Bind one review to its frozen stage, source and input membership count."""
    root: Path
    config: dict

    @classmethod
    def read(cls, root):
        """Read copied configuration; refuse scientific execution or promotion."""
        root=Path(root)
        cfg=json.loads((root/'config/hicat_consensus_review.json').read_text())
        if cfg['stage']!='03_full_data_consensus_benchmark' or cfg['status']!='IN_REVIEW':
            raise ValueError('Review stage/status mismatch')
        if cfg['refit'] or cfg['lock_annotations']:
            raise ValueError('Review cannot refit or lock annotations')
        return cls(root,cfg)


class ConsensusReviewInputs:
    """Fully validate saved consensus outputs, then expose aligned review tables."""

    def __init__(self, settings):
        """Resolve frozen source paths independently of the publication directory."""
        self.settings=settings
        self.cfg=settings.config
        self.source=Path(self.cfg['source_run'])
        self.endpoint=self.source/self.cfg['source_consensus']
        self.checks=[]

    def check(self, name, condition):
        """Record and enforce an input or saved-output invariant."""
        self.checks.append(dict(check=name,status='PASS' if condition else 'FAIL'))
        if not condition:raise ValueError(name)

    def read(self, output):
        """Verify all aggregate/DE artifacts and exact cells, genes and R.

        Historical fitting models are not inputs to this reporting stage. The
        source consensus's original validation scope is carried into provenance
        without claiming it has been retroactively changed.
        """
        output=Path(output)
        pins=json.loads((self.settings.root/'inputs/source_identity.json').read_text())
        for relative,digest in pins['files'].items():
            self.check('source_hash:'+relative,sha256(self.source/relative)==digest)
        self.aggregation,aggseal=complete_attempt(self.endpoint/'aggregation')
        self.final,deseal=complete_attempt(self.endpoint/'DE')
        self.check('DE_aggregation_parent',deseal['contract']['parents']['aggregation']==aggseal['manifest_sha256'])
        membership=json.loads((self.aggregation/'membership_index.json').read_text())
        self.check('denominator',membership['denominator']==self.cfg['iterations'])
        self.check('exact_iteration_ids',set(membership['parent_seals'])==
                   set(map(str,range(self.cfg['iterations'])))|{'fullfit_labels'})
        self.universe=pd.read_csv(self.source/'config/cell_universe.tsv.gz',sep='\t')
        self.assignments=pd.read_csv(self.final/'all_cell_assignments.tsv.gz',sep='\t')
        self.check('unique_cells',self.assignments.cell_id.is_unique and self.universe.cell_id.is_unique)
        self.check('cell_order',np.array_equal(self.assignments.cell_id,self.universe.cell_id))
        self.check('population',len(self.assignments)==self.cfg['cells'])
        self.check('all_cells_assigned',not self.assignments.diagnostic_cluster.isna().any())
        self.check('samples',self.universe.technical_sample_id.nunique()==self.cfg['samples'])
        self.means=pd.read_csv(self.final/'cluster_mean_log1p_cpm.tsv.gz',sep='\t',index_col=0)
        genes=np.load(self.source/'config/gene_ids.npy',allow_pickle=False).astype(str)
        self.check('gene_order',np.array_equal(self.means.columns.astype(str),genes))
        self.check('genes',len(genes)==self.cfg['genes'])
        self.check('cluster_labels',set(self.means.index)==set(self.assignments.diagnostic_cluster))
        self.scope=json.loads((self.final/'scope.json').read_text())
        self.check('source_scope',self.scope['real_iterations']==self.cfg['iterations'] and not self.scope['annotations_locked'])
        cfg=json.loads((self.source/'config/run.json').read_text())
        normalized=Path(cfg['legacy_benchmark'])/'inputs/full_log1p_cpm.h5ad'
        with h5py.File(normalized,'r') as handle:
            values=handle['var/gene_symbol']
            if isinstance(values,h5py.Group):
                symbols=values['categories'].asstr()[:][values['codes'][:]]
            else:symbols=values.asstr()[:]
        self.check('symbol_count',len(symbols)==len(genes))
        self.symbols=dict(zip(genes,symbols))
        for name in ['all_cell_assignments.tsv.gz','pairwise_DE.tsv','cluster_mean_log1p_cpm.tsv.gz','cluster_detection.tsv.gz']:
            shutil.copy2(self.final/name,output/'tables'/name)
        for name in ['cluster_consensus_diagnostics.tsv','cell_consensus_diagnostics.tsv.gz','refinement_trace.json','merge_actions.json']:
            shutil.copy2(self.aggregation/name,output/'tables'/name)
        write_json(output/'tables/source_checkpoints.json',dict(
            aggregation=str(self.aggregation),DE=str(self.final),aggregation_manifest=aggseal['manifest_sha256'],
            DE_manifest=deseal['manifest_sha256'],source_summary=json.loads((self.endpoint/'SUMMARY.json').read_text()),
            review_validation='All saved aggregation and DE artifacts hashed and reopened; exact source identities and alignment checked.'))
        return self
