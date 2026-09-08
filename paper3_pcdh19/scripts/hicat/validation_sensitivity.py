"""One controlled fine-clustering comparison inside fixed pilot coarse parents.

Only q1 and qdiff change from the saved baseline configuration. The baseline
seed, exact cells, normalized expression, graph/PCA settings and DE score are
retained. Candidate merges affect this separate sensitivity partition only.
This module never fits coarse clusters, chooses a desired K, or reads labels
from the provisional biological review reference.
"""
from copy import deepcopy
from pathlib import Path
import json
import numpy as np
import pandas as pd
import anndata as ad
from .engine import PilotEngine
from .provenance import write_json


def resolve_sensitivity_config(baseline_config, validation_config):
    """Return the baseline engine controls with exactly two fine-DE changes.

    Reject a changed seed, DE score or unlocked-parent request. The explicit
    returned diff is saved for review so unrelated parameter drift is visible.
    """
    cfg = deepcopy(baseline_config)
    if validation_config['seed'] != cfg['baseline_seed']:
        raise ValueError('Sensitivity must retain the baseline fitting seed')
    if not validation_config['fixed_baseline_parents']:
        raise ValueError('Refitting coarse parents is outside Step 07')
    before = dict(cfg['fine_thresholds'])
    if any(before.get(key) != value for key, value in
           {'q1_thresh': 0.3, 'qdiff_thresh': 0.5, 'score_thresh': 150}.items()):
        raise ValueError('Baseline must be the specified q1=.3, qdiff=.5, score=150 pilot')
    changes = validation_config['fine_changes']
    if changes != {'q1_thresh': 0.4, 'qdiff_thresh': 0.7}:
        raise ValueError('Only the specified Allen-reference comparison is authorized')
    if before['score_thresh'] != validation_config['fixed_fine_score']:
        raise ValueError('Fine DE score must be unchanged')
    cfg['thresholds'] = dict(before, **changes)
    diff = [{'parameter': key, 'baseline': before[key], 'sensitivity': value}
            for key, value in changes.items()]
    return cfg, diff


class FixedParentSensitivity:
    """Fit only the alternative fine partition and save every Allen node asset."""

    def __init__(self, directory, baseline_config, validation_config, progress):
        """Resolve immutable controls without allocating expression or fitting."""
        self.directory = Path(directory)
        self.directory.mkdir(parents=True, exist_ok=True)
        self.cfg, self.diff = resolve_sensitivity_config(baseline_config, validation_config)
        self.seed = validation_config['seed']
        self.progress = progress

    def run(self, pilot):
        """Return candidate labels and parent counts using exactly the pilot rows.

        Input ``pilot.X`` remains raw counts; only the existing ``log1p_cpm``
        layer is copied into bounded, dense parent matrices for Allen's API.
        Expression is not renormalized, regressed or batch corrected. Candidate
        IDs have a separate ``.A`` namespace to avoid implying label identity.
        """
        import transcriptomic_clustering as tc
        tc.memory.set_memory_limit(GB=self.cfg['pilot_scratch_memory_gb'])
        tc.memory.allow_chunking = False
        write_json(self.directory/'resolved_engine_config.json', self.cfg)
        pd.DataFrame(self.diff).to_csv(self.directory/'parameter_changes.tsv', sep='\t', index=False)
        parent = pilot.obs['hicat_coarse_baseline'].astype(str)
        result = pd.Series(index=pilot.obs_names, dtype=object, name='hicat_fine_allen_reference')
        summaries = []
        for key in sorted(parent.unique()):
            rows = np.flatnonzero(parent.to_numpy() == key)
            with self.progress.track('validation.fixed_parent_sensitivity', parent=key,
                                     cells=len(rows), seed=self.seed,
                                     thresholds=self.cfg['thresholds']) as event:
                # Do not pass metadata/hypotheses to the clustering engine.
                obj = ad.AnnData(pilot.layers['log1p_cpm'][rows].toarray(),
                                 obs=pd.DataFrame(index=pilot.obs_names[rows]),
                                 var=pilot.var.copy(), dtype=np.float64)
                candidate = dict(name='allen_reference/'+key, seed=self.seed,
                                 score_thresh=self.cfg['thresholds']['score_thresh'])
                engine = PilotEngine(self.cfg, candidate, self.directory/('fine_'+key), self.progress)
                labels, summary, _, _ = engine.run(obj)
                labels = [key+'.A'+label[1:] for label in labels]
                result.iloc[rows] = labels
                repeat = pilot.obs.iloc[rows]['hicat_fine_seed_repeat'].astype(str)
                summaries.append(dict(parent=key, n_cells=len(rows),
                    n_fine_current=pilot.obs.iloc[rows]['hicat_fine_baseline'].nunique(),
                    n_fine_Allen_reference=len(set(labels)),
                    n_repeat_labels_intersecting_fixed_parent=repeat.nunique(),
                    unresolved_candidate_pairs=summary['unresolved_final_pairs']))
                event.update(summary)
                del obj, engine
        if result.isna().any():
            raise ValueError('Candidate does not cover every pilot cell')
        if not all(result.str.split('.').str[0] == parent):
            raise ValueError('Candidate fine partition crossed a fixed coarse parent')
        table = pd.DataFrame(summaries)
        result.to_csv(self.directory/'cell_assignments.tsv.gz', sep='\t', index_label='cell_id')
        table.to_csv(self.directory/'parent_counts.tsv', sep='\t', index=False)
        return result, table
