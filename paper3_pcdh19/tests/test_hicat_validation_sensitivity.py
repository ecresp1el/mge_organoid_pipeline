"""Controlled-contrast and saved fixed-parent execution checks for Step 07."""
from copy import deepcopy
import json
from pathlib import Path
import tempfile
import unittest

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

from hicat.provenance import Progress
from hicat.validation_sensitivity import FixedParentSensitivity, resolve_sensitivity_config


class SensitivityContracts(unittest.TestCase):
    """Exercise real Allen terminal-node outputs without a long graph benchmark."""

    def setUp(self):
        """Load the actual baseline controls and the explicitly authorized contrast."""
        self.baseline = json.loads((Path(__file__).parents[1] / 'config/hicat_coarse_fine.json').read_text())
        self.validation = dict(seed=self.baseline['baseline_seed'], fixed_baseline_parents=True,
                               fine_changes=dict(q1_thresh=.4, qdiff_thresh=.7), fixed_fine_score=150)

    def test_exact_contrast_preserves_all_unrelated_controls_and_inputs(self):
        """Resolved controls differ only at q1/qdiff, with no caller-dict mutation."""
        before = deepcopy(self.baseline)
        validation_before = deepcopy(self.validation)
        resolved, changes = resolve_sensitivity_config(self.baseline, self.validation)
        self.assertEqual(self.baseline, before)
        self.assertEqual(self.validation, validation_before)
        for key in before:
            self.assertEqual(resolved[key], before[key], key)
        desired = dict(before['fine_thresholds'], q1_thresh=.4, qdiff_thresh=.7)
        self.assertEqual(resolved['thresholds'], desired)
        self.assertEqual({x['parameter'] for x in changes}, {'q1_thresh', 'qdiff_thresh'})
        self.assertEqual(resolved['thresholds']['score_thresh'], 150)
        resolved['thresholds']['min_genes'] = 99
        self.assertEqual(self.baseline['fine_thresholds']['min_genes'], before['fine_thresholds']['min_genes'])

    def test_seed_parent_score_and_additional_parameter_drift_are_rejected(self):
        """An unrelated seed/score/graph change cannot be disguised as sensitivity."""
        cases = [dict(seed=self.validation['seed'] + 1), dict(fixed_baseline_parents=False),
                 dict(fixed_fine_score=151), dict(fine_changes=dict(q1_thresh=.4, qdiff_thresh=.7, graph_k=30)),
                 dict(fine_changes=dict(q1_thresh=.5, qdiff_thresh=.7))]
        for changed in cases:
            with self.subTest(changed=changed):
                config = dict(self.validation, **changed)
                with self.assertRaises(ValueError):
                    resolve_sensitivity_config(self.baseline, config)

    def test_changed_baseline_contrast_is_rejected_even_with_matching_score(self):
        """The reference contrast must start from exactly q1=.3/qdiff=.5/score150."""
        for key, value in [('q1_thresh', .2), ('qdiff_thresh', .6), ('score_thresh', 151)]:
            with self.subTest(parameter=key):
                baseline = deepcopy(self.baseline)
                baseline['fine_thresholds'][key] = value
                validation = dict(self.validation)
                if key == 'score_thresh':
                    validation['fixed_fine_score'] = value
                with self.assertRaises(ValueError):
                    resolve_sensitivity_config(baseline, validation)

    def test_fixed_parent_run_saves_separate_labels_and_preserves_expression(self):
        """The actual adapter retains 80 cells and both fixed parents without mutation.

        Each parent contains 40 cells while the declared synthetic minimum for
        recursive fitting is 41. Allen's real terminal-node means/evidence and
        the real adapter's manifests/labels are exercised; graph splitting is
        intentionally not benchmarked by this contract test.
        """
        baseline = deepcopy(self.baseline)
        baseline.update(min_recursive_cells=41, pilot_scratch_memory_gb=1,
                        max_hvgs=60, pca_components=10, max_retained_pcs=6, graph_k=5)
        rng = np.random.default_rng(23)
        raw = sparse.csr_matrix(rng.poisson(3, size=(80, 80)).astype(np.int32))
        normalized = raw.astype(np.float64).multiply(1e6 / np.asarray(raw.sum(axis=1))).tocsr()
        normalized.data = np.log1p(normalized.data)
        parent = np.repeat(['C0001', 'C0002'], 40)
        fine = np.repeat(['C0001.F0001', 'C0001.F0002', 'C0002.F0001', 'C0002.F0002'], 20)
        obs = pd.DataFrame(dict(hicat_coarse_baseline=parent, hicat_fine_baseline=fine,
                                hicat_fine_seed_repeat=fine, unrelated_metadata=np.arange(80)),
                           index=['cell%03d' % i for i in range(80)])
        var = pd.DataFrame(dict(gene_symbol=['Gene%d' % i for i in range(80)]),
                           index=['g%d' % i for i in range(80)])
        pilot = ad.AnnData(raw.copy(), obs=obs.copy(), var=var.copy(), dtype=np.int32)
        pilot.layers['log1p_cpm'] = normalized.copy()
        with tempfile.TemporaryDirectory(prefix='hicat_fixed_parent_test_', dir='/tmp') as temporary:
            directory = Path(temporary)
            runner = FixedParentSensitivity(directory / 'sensitivity', baseline,
                                            self.validation, Progress(directory / 'progress'))
            result, counts = runner.run(pilot)
            self.assertEqual(result.index.tolist(), pilot.obs_names.tolist())
            self.assertFalse(result.isna().any())
            self.assertEqual(sorted(result.unique()), ['C0001.A0001', 'C0002.A0001'])
            self.assertTrue((result.str.split('.').str[0].to_numpy() == parent).all())
            self.assertFalse(set(result).intersection(fine))
            self.assertEqual(counts.n_cells.sum(), 80)
            self.assertEqual(counts.n_fine_Allen_reference.tolist(), [1, 1])
            self.assertEqual(counts.n_fine_current.tolist(), [2, 2])
            saved = pd.read_csv(directory / 'sensitivity/cell_assignments.tsv.gz', sep='\t', index_col=0)
            self.assertEqual(saved.index.tolist(), pilot.obs_names.tolist())
            self.assertEqual(saved.iloc[:, 0].tolist(), result.tolist())
            for key in sorted(set(parent)):
                base = directory / ('sensitivity/fine_' + key)
                tree = pd.read_csv(base / 'recursion_tree.tsv', sep='\t')
                self.assertTrue(tree.terminal_reason.eq('below_min_recursive_cells').all())
                self.assertTrue((base / 'cluster_mean_log1p_cpm.tsv.gz').is_file())
                self.assertTrue((base / 'final_pairwise_evidence.tsv').is_file())
            actual_cfg = json.loads((directory / 'sensitivity/resolved_engine_config.json').read_text())
            self.assertEqual(actual_cfg['thresholds'], dict(baseline['fine_thresholds'], q1_thresh=.4, qdiff_thresh=.7))
        self.assertEqual((pilot.X != raw).nnz, 0)
        self.assertEqual((pilot.layers['log1p_cpm'] != normalized).nnz, 0)
        pd.testing.assert_frame_equal(pilot.obs, obs)
        pd.testing.assert_frame_equal(pilot.var, var)


if __name__ == '__main__':
    unittest.main()
