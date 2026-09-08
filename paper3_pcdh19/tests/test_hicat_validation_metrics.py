"""Synthetic contracts for same-cell overlap and descriptive Step 07 diagnostics."""
import json
import unittest

import numpy as np
import pandas as pd

from hicat.validation_metrics import (partition_overlap, sample_composition,
                                      phase_matched_identity)


def series(values):
    """Construct explicit cell-indexed labels for synthetic cases."""
    return pd.Series(values, index=['cell%d' % i for i in range(len(values))])


class ValidationMetricsContracts(unittest.TestCase):
    """Verify known split/merge and confounding cases, not implementation details."""

    def test_split_and_merge_are_distinguished_without_renaming(self):
        """A splits into X/Y while B and C merge into Z; matched metrics are exact."""
        baseline = series(['A'] * 4 + ['B'] * 4 + ['C'] * 4)
        comparison = series(['X'] * 2 + ['Y'] * 2 + ['Z'] * 8)
        result = partition_overlap(baseline, comparison.iloc[::-1])
        table = result['clusters']
        self.assertEqual(table.loc['A', 'best_match'], 'X')
        self.assertEqual(table.loc['A', 'jaccard'], .5)
        self.assertEqual(table.loc['A', 'precision'], 1)
        self.assertEqual(table.loc['A', 'recall'], .5)
        self.assertTrue(table.loc['A', 'splits'])
        self.assertFalse(table.loc['A', 'merges'])
        self.assertEqual(json.loads(table.loc['A', 'splits_into']), ['X', 'Y'])
        self.assertEqual(json.loads(table.loc['B', 'merges_with']), ['C'])
        self.assertTrue(table.loc['B', 'merges'])
        self.assertFalse(table.loc['B', 'splits'])
        self.assertEqual(table.loc['B', 'precision'], .5)
        self.assertEqual(table.loc['B', 'recall'], 1)
        self.assertEqual(table.loc['B', 'boundary_survives'], 'partial')
        np.testing.assert_allclose(result['row_fractions'].sum(axis=1), 1)
        np.testing.assert_allclose(result['column_fractions'].sum(axis=0), 1)
        self.assertEqual(result['counts'].to_numpy().sum(), 12)

    def test_strict_ten_percent_and_identical_partition(self):
        """A target receiving exactly ten percent does not count as a split."""
        labels = series(['A'] * 10)
        result = partition_overlap(labels, series(['X'] * 9 + ['Y']))
        self.assertFalse(result['clusters'].loc['A', 'splits'])
        self.assertEqual(result['clusters'].loc['A', 'n_targets_over_split_fraction'], 1)
        identical = partition_overlap(series(['A', 'A', 'B', 'B']), series(['X', 'X', 'Y', 'Y']))
        self.assertEqual(identical['summary']['ARI'], 1)
        self.assertEqual(identical['summary']['NMI'], 1)
        self.assertTrue(identical['clusters'].stability.eq('HIGH_STABILITY').all())

    def test_candidate_display_order_preserves_actual_parent_blocks(self):
        """Separate .A sensitivity namespaces preserve parents as .F labels do."""
        baseline = series(['C0001.F0001'] * 4 + ['C0002.F0001'] * 4)
        comparison = series(['C0002.A0001'] * 4 + ['C0001.A0001'] * 4)
        result = partition_overlap(baseline, comparison)
        self.assertEqual(result['summary']['comparison_order'], ['C0001.A0001', 'C0002.A0001'])
        self.assertEqual(result['clusters'].loc['C0001.F0001', 'best_match'], 'C0002.A0001')

    def test_mismatched_or_duplicate_ids_raise(self):
        """An overlap must not silently intersect two different cell selections."""
        a = series(['A', 'B'])
        b = pd.Series(['X', 'Y'], index=['cell0', 'other_cell'])
        with self.assertRaises(ValueError):
            partition_overlap(a, b)
        b.index = ['cell0', 'cell0']
        with self.assertRaises(ValueError):
            partition_overlap(a, b)

    def test_balanced_twelve_sample_cluster(self):
        """A balanced design has unit enrichment, entropy log12, effective12."""
        samples = series(['JZ-%02d' % (i + 1) for i in range(12)] * 10)
        labels = series(['A'] * 60 + ['B'] * 60)
        result = sample_composition(labels, samples)
        np.testing.assert_allclose(result['enrichment'], 1)
        np.testing.assert_allclose(result['cluster_fractions'].sum(axis=1), 1)
        np.testing.assert_allclose(result['sample_fractions'].sum(axis=0), 1)
        np.testing.assert_allclose(result['clusters'].sample_entropy, np.log(12))
        np.testing.assert_allclose(result['clusters'].effective_samples, 12)
        np.testing.assert_allclose(result['clusters'].chi_square, 0)
        self.assertFalse(result['clusters'].sample_concern.any())

    def test_sample_concentration_is_descriptive_and_small_clusters_remain(self):
        """Single-sample concentration is flagged only after the size rule."""
        samples = series(['JZ-%02d' % (i + 1) for i in range(12) for _ in range(100)])
        labels = series(['A'] * 50 + ['tiny'] * 5 + ['other'] * 1145)
        table = sample_composition(labels, samples)['clusters']
        self.assertTrue(table.loc['A', 'sample_concern'])
        self.assertEqual(table.loc['A', 'effective_samples'], 1)
        self.assertFalse(table.loc['tiny', 'sample_concern'])
        self.assertTrue(table.loc['tiny', 'limited_by_small_cluster'])
        self.assertEqual(table.loc['tiny', 'n_cells'], 5)

    def test_phase_composition_can_explain_unstratified_score_difference(self):
        """Identical identity scores within phases yield no retained difference."""
        labels = series(['A'] * 100 + ['B'] * 100)
        parents = series(['C0001'] * 200)
        phases = series(['S'] * 80 + ['G1-like'] * 20 + ['S'] * 20 + ['G1-like'] * 80)
        scores = pd.DataFrame({'identity': phases.eq('S').astype(float) * 2}, index=labels.index)
        result = phase_matched_identity(labels, parents, phases, scores)
        pair = result['pairs'].iloc[0]
        self.assertEqual(pair.identity_difference_persists_after_cell_cycle_stratification, 'no')
        self.assertAlmostEqual(pair.cell_cycle_total_variation, .6)
        self.assertTrue(pair.cell_cycle_boundary_review)
        self.assertEqual(pair.n_eligible_phases, 2)
        self.assertGreater(pair.max_abs_unstratified_identity_effect, .5)
        stratified = result['phase_effects'].query('phase != "all"')
        np.testing.assert_allclose(stratified.standardized_effect, 0)

    def test_identity_difference_persists_in_both_matched_phases(self):
        """An independent group difference survives adjustment by stratification."""
        labels = series(['A'] * 100 + ['B'] * 100)
        parents = series(['C'] * 200)
        phases = series((['S'] * 50 + ['G1-like'] * 50) * 2)
        scores = pd.DataFrame({'identity': labels.eq('A').astype(float) * 2
                              + phases.eq('S').astype(float) * .1}, index=labels.index)
        result = phase_matched_identity(labels, parents, phases, scores)
        pair = result['pairs'].iloc[0]
        self.assertEqual(pair.identity_difference_persists_after_cell_cycle_stratification, 'yes')
        self.assertEqual(pair.n_retained_phases, 2)
        self.assertFalse(pair.cell_cycle_boundary_review)

    def test_nonoverlapping_phases_are_insufficient_not_negative(self):
        """No matched phase can exclude an independent identity difference."""
        labels = series(['A'] * 40 + ['B'] * 40)
        parents = series(['C'] * 80)
        phases = series(['S'] * 40 + ['G1-like'] * 40)
        scores = pd.DataFrame({'identity': labels.eq('A').astype(float)}, index=labels.index)
        result = phase_matched_identity(labels, parents, phases, scores)
        pair = result['pairs'].iloc[0]
        self.assertEqual(pair.identity_difference_persists_after_cell_cycle_stratification,
                         'insufficient_cells')
        self.assertFalse(pair.cell_cycle_boundary_review)
        self.assertEqual(pair.n_eligible_phases, 0)

    def test_single_eligible_phase_is_partial_and_cross_parent_labels_rejected(self):
        """One matched phase supplies limited evidence, not multi-phase support."""
        labels = series(['A'] * 40 + ['B'] * 40)
        phases = series(['S'] * 80)
        parents = series(['C'] * 80)
        scores = pd.DataFrame({'identity': labels.eq('A').astype(float)}, index=labels.index)
        result = phase_matched_identity(labels, parents, phases, scores)
        self.assertEqual(result['pairs'].iloc[0].identity_difference_persists_after_cell_cycle_stratification,
                         'partial')
        parents.iloc[0] = 'other'
        with self.assertRaises(ValueError):
            phase_matched_identity(labels, parents, phases, scores)


if __name__ == '__main__':
    unittest.main()
