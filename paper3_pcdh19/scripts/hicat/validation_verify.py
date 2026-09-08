"""Independently verify a finalized Step 07 review package without modifying it.

Run ``python -m hicat.validation_verify --run-dir /path/to/versioned/run``.
The CLI prints JSON checks to stdout and returns a nonzero exit status for a
failed check. Redirect stdout to a provenance log if desired; this module never
writes to the package. No normalization, gene scoring, clustering, statistical
refitting, or source-raw-count reads are performed.

Verification streams each output file once for the declared SHA-256 manifest,
then reads small HDF5 metadata/observation arrays. It checks original/sensitivity
membership bookkeeping, cell-cycle source-table alignment, independently
recomputed seed agreement, exact AnnData inventories, full-data sidecar scope,
and all individual figure/source and fine-review-page references.
"""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import time

import h5py
import numpy as np
import pandas as pd
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score


EXPECTED_PILOT_CELLS = 12000
EXPECTED_PILOT_GENES = 19071
EXPECTED_COARSE = 4
EXPECTED_FINE = 38
EXPECTED_FULL_CELLS = 446349
EXPECTED_CANONICAL_GENES = 68


def _native(value):
    """Decode HDF5/NumPy values into deterministic JSON-compatible primitives."""
    if isinstance(value, bytes):
        return value.decode('utf-8')
    if isinstance(value, np.ndarray):
        return [_native(v) for v in value.tolist()] if value.ndim else _native(value.item())
    if isinstance(value, np.generic):
        return _native(value.item())
    if isinstance(value, dict):
        return {str(k): _native(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_native(v) for v in value]
    return value


def _element(node):
    """Read the simple arrays/dicts/categoricals used in these AnnData assets.

    Unknown group encodings fail rather than being silently flattened into an
    inaccurate inventory. Expression CSR groups are never passed to this reader.
    """
    if isinstance(node, h5py.Dataset):
        return _native(node[()])
    encoding = _native(node.attrs.get('encoding-type', 'dict'))
    if encoding == 'categorical':
        categories = _element(node['categories'])
        return [categories[int(code)] if code >= 0 else None for code in node['codes'][:]]
    if encoding == 'dict':
        return {key: _element(node[key]) for key in node}
    if encoding == 'nullable-boolean' or encoding == 'nullable-integer':
        values = np.asarray(node['values'][:], dtype=object)
        values[node['mask'][:]] = None
        return _native(values)
    raise ValueError('Unsupported inventory encoding %r at %s' % (encoding, node.name))


def _frame_metadata(group):
    """Return dataframe index and ordered columns from HDF5 encoding metadata."""
    index_key = _native(group.attrs['_index'])
    return _element(group[index_key]), _native(group.attrs['column-order'])


def _shape(node):
    """Return sparse or dense matrix shape without reading matrix contents."""
    return list(node.shape) if isinstance(node, h5py.Dataset) else _native(node.attrs['shape'])


def _digest(path):
    """Stream one saved output's SHA-256 with bounded memory."""
    result = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b''):
            result.update(block)
    return result.hexdigest()


def _json_equal(left, right):
    """Compare inventories including explicit NaN values and nested arrays."""
    return json.dumps(_native(left), sort_keys=True, allow_nan=True) == json.dumps(
        _native(right), sort_keys=True, allow_nan=True)


class PackageVerifier:
    """Accumulate independent, read-only checks for one finalized review run."""

    def __init__(self, run_dir):
        """Record paths and an empty check collection; create no files."""
        self.root = Path(run_dir).resolve()
        self.out = self.root / 'outputs'
        self.checks = []
        self.started = time.monotonic()

    def check(self, name, condition, detail=None):
        """Record a Boolean contract and return its truth value for dependencies."""
        passed = bool(condition)
        self.checks.append(dict(check=name, status='PASS' if passed else 'FAIL',
                                detail=_native(detail)))
        return passed

    def _inside_output(self, relative):
        """Resolve only nonempty relative paths contained in this output package."""
        relative = str(relative)
        if not relative or Path(relative).is_absolute():
            raise ValueError('Expected a package-relative asset path: %r' % relative)
        path = (self.out / relative).resolve()
        if self.out.resolve() not in path.parents:
            raise ValueError('Asset escapes the finalized output directory: %r' % relative)
        return path

    def verify_manifest(self):
        """Hash all declared outputs and detect unmanifested/missing package files."""
        manifest_path = self.out / 'output_manifest.tsv'
        table = pd.read_csv(manifest_path, sep='\t')
        self.check('manifest_required_columns', {'relative_path', 'bytes', 'sha256'}.issubset(table))
        self.check('manifest_unique_paths', table.relative_path.is_unique)
        declared = set(table.relative_path.astype(str))
        actual = {str(p.relative_to(self.out)) for p in self.out.rglob('*')
                  if p.is_file() and p.name != 'output_manifest.tsv'}
        self.check('manifest_complete_file_set', declared == actual,
                   dict(unmanifested=sorted(actual - declared), missing=sorted(declared - actual)))
        for row in table.itertuples():
            try:
                path = self._inside_output(row.relative_path)
                present = path.is_file()
                correct_size = present and path.stat().st_size == int(row.bytes)
                correct_hash = present and _digest(path) == str(row.sha256)
                self.check('manifest_file:' + str(row.relative_path), correct_size and correct_hash,
                           dict(exists=present, size_matches=correct_size, sha256_matches=correct_hash))
            except (OSError, ValueError) as error:
                self.check('manifest_file:' + str(row.relative_path), False, str(error))
        return len(table)

    def verify_pilot(self):
        """Check saved memberships and complete slots/uns using HDF5 metadata."""
        path = self.out / 'pcdh19_hicat_hierarchy_validation.h5ad'
        slots = json.loads((self.out / 'anndata_slot_inventory.json').read_text())
        inventory = json.loads((self.out / 'anndata_uns_inventory.json').read_text())
        with h5py.File(path, 'r') as handle:
            shape = _shape(handle['X'])
            cells, obs_columns = _frame_metadata(handle['obs'])
            genes, var_columns = _frame_metadata(handle['var'])
            obs = pd.DataFrame({column: _element(handle['obs'][column]) for column in obs_columns}, index=cells)
            actual_uns = _element(handle['uns'])
            self.check('pilot_shape', shape == [EXPECTED_PILOT_CELLS, EXPECTED_PILOT_GENES], shape)
            self.check('pilot_unique_cell_gene_ids', len(set(cells)) == len(cells) and len(set(genes)) == len(genes))
            self.check('pilot_integer_raw_CSR', _native(handle['X'].attrs.get('encoding-type')) == 'csr_matrix'
                       and handle['X']['data'].dtype.kind in 'iu')
            self.check('pilot_normalized_layer_shape', _shape(handle['layers']['log1p_cpm']) == shape)
            self.check('slot_inventory_shape', slots['shape'] == shape)
            for name, actual in [('obs', obs_columns), ('var', var_columns),
                                 ('layers', list(handle['layers'])), ('obsm', list(handle['obsm'])),
                                 ('obsp', list(handle['obsp'])), ('uns', list(handle['uns']))]:
                self.check('slot_inventory_' + name, set(slots[name]) == set(actual),
                           dict(declared=slots[name], actual=actual))
            self.check('slot_inventory_raw_present', bool(slots['raw_present']) == ('raw' in handle))
            self.check('complete_uns_inventory_exact', _json_equal(inventory, actual_uns))
            self.check('pilot_no_locked_annotations', not actual_uns['hicat_validation']['annotation_locked'])
            self.check('pilot_no_full_clustering', not actual_uns['hicat_validation']['full_data_clustering'])
            self.check('pilot_saved_program_score_shape', handle['obsm']['X_validation_program_scores'].shape[0] == EXPECTED_PILOT_CELLS)
        self.check('pilot_baseline_coarse_count', obs.hicat_coarse_baseline.nunique() == EXPECTED_COARSE)
        self.check('pilot_baseline_fine_count', obs.hicat_fine_baseline.nunique() == EXPECTED_FINE)
        self.check('pilot_repeat_counts', obs.hicat_coarse_seed_repeat.nunique() == 5 and obs.hicat_fine_seed_repeat.nunique() == 39)
        counts = obs.technical_sample_id.value_counts()
        self.check('pilot_exactly_1000_per_12_samples', len(counts) == 12 and counts.eq(1000).all())
        fixed = obs.hicat_fine_allen_reference.astype(str).str.split('.').str[0]
        self.check('candidate_retains_fixed_parents', fixed.equals(obs.hicat_coarse_baseline.astype(str)))
        self.check('candidate_separate_A_namespace', obs.hicat_fine_allen_reference.astype(str).str.match(r'^C[0-9]+\.A[0-9]+$').all())
        self.check('candidate_complete_membership', obs.hicat_fine_allen_reference.notna().all())
        original = pd.read_csv(self.root / 'inputs/original_pilot_assignments.tsv.gz', sep='\t', index_col=0)
        self.check('original_cell_order_retained', original.index.tolist() == obs.index.tolist())
        for name in ['hicat_coarse_baseline', 'hicat_fine_baseline', 'hicat_coarse_seed_repeat', 'hicat_fine_seed_repeat']:
            self.check('original_partition_retained_' + name, original[name].astype(str).equals(obs[name].astype(str)))
        return obs

    def verify_phase(self, obs):
        """Recompute phase rules and compare every cell to the saved source table."""
        table = pd.read_csv(self.out / 'marker_programs/cell_cycle_scores.tsv.gz', sep='\t', index_col=0)
        ids_match = table.index.is_unique and set(table.index) == set(obs.index)
        self.check('phase_source_exact_pilot_cell_set', ids_match, dict(rows=len(table)))
        if not ids_match:
            return
        table = table.loc[obs.index]
        finite = np.isfinite(table[['S_phase', 'G2M']].to_numpy()).all()
        self.check('phase_scores_finite', finite)
        expected = pd.Series('S', index=table.index)
        expected.loc[table.G2M > table.S_phase] = 'G2/M'
        expected.loc[(table.G2M < 0) & (table.S_phase < 0)] = 'G1-like'
        self.check('phase_rule_recomputed_from_scores', expected.equals(obs.validation_phase.astype(str)))
        self.check('phase_table_matches_saved_object', table.validation_phase.astype(str).equals(obs.validation_phase.astype(str)))
        self.check('phase_S_values_match_object', np.allclose(table.S_phase, obs.validation_S_score, atol=1e-12, rtol=1e-12))
        self.check('phase_G2M_values_match_object', np.allclose(table.G2M, obs.validation_G2M_score, atol=1e-12, rtol=1e-12))
        for level in ['coarse', 'fine']:
            phase_counts = pd.read_csv(self.out / ('cell_cycle/' + level + '_phase_counts.tsv'), sep='\t', index_col=0)
            actual = pd.crosstab(obs['hicat_' + level + '_baseline'], obs.validation_phase.replace({'G2/M': 'G2M'}))
            actual = actual.reindex(index=phase_counts.index, columns=phase_counts.columns, fill_value=0)
            self.check('phase_count_table_' + level, np.array_equal(actual.to_numpy(), phase_counts.to_numpy())
                       and phase_counts.to_numpy().sum() == EXPECTED_PILOT_CELLS)

    def verify_partitions(self, obs, report_summary):
        """Independently recompute seed agreement and every raw overlap count."""
        for level in ['coarse', 'fine']:
            left = obs['hicat_' + level + '_baseline'].astype(str)
            right = obs['hicat_' + level + '_seed_repeat'].astype(str)
            metrics = json.loads((self.out / ('stability/' + level + '_summary.json')).read_text())
            ari = float(adjusted_rand_score(left, right))
            nmi = float(normalized_mutual_info_score(left, right, average_method='arithmetic'))
            self.check('seed_' + level + '_ARI', np.isclose(ari, metrics['ARI'], atol=1e-12, rtol=0), dict(recomputed=ari, saved=metrics['ARI']))
            self.check('seed_' + level + '_NMI', np.isclose(nmi, metrics['NMI'], atol=1e-12, rtol=0), dict(recomputed=nmi, saved=metrics['NMI']))
            raw = pd.read_csv(self.out / ('stability/' + level + '_counts.tsv'), sep='\t', index_col=0)
            expected = pd.crosstab(left, right).reindex(index=raw.index, columns=raw.columns, fill_value=0)
            self.check('seed_' + level + '_raw_overlap_exact', np.array_equal(expected.to_numpy(), raw.to_numpy())
                       and raw.to_numpy().sum() == EXPECTED_PILOT_CELLS)
            fractions = pd.read_csv(self.out / ('stability/' + level + '_row_fractions.tsv'), sep='\t', index_col=0)
            expected_fraction = expected.div(expected.sum(axis=1), axis=0)
            self.check('seed_' + level + '_row_fraction_exact', fractions.index.equals(expected_fraction.index)
                       and fractions.columns.equals(expected_fraction.columns)
                       and np.allclose(expected_fraction, fractions, atol=1e-12, rtol=0))
            per_cluster = pd.read_csv(self.out / ('stability/' + level + '_clusters.tsv'), sep='\t', index_col=0)
            row_sizes, column_sizes = expected.sum(axis=1), expected.sum(axis=0)
            for cluster in expected.index:
                overlap = expected.loc[cluster]
                jaccard = overlap / (row_sizes[cluster] + column_sizes - overlap)
                best = sorted(expected.columns, key=lambda target: (-jaccard[target], -overlap[target], target))[0]
                saved = per_cluster.loc[cluster]
                precision, recall = overlap[best] / column_sizes[best], overlap[best] / row_sizes[cluster]
                f1 = 2 * overlap[best] / (row_sizes[cluster] + column_sizes[best])
                valid = saved.best_match == best and np.allclose(
                    [saved.jaccard, saved.precision, saved.recall, saved.F1, saved.fraction_retained],
                    [jaccard[best], precision, recall, f1, recall], atol=1e-12, rtol=0)
                self.check('seed_' + level + '_best_match_metrics:' + cluster, valid)
            if level == 'fine':
                summary_seed = report_summary.get('readiness_assessment_axes', {}).get('seed_stability', {})
                self.check('main_summary_seed_agreement', np.isclose(summary_seed.get('ARI', np.nan), ari, atol=1e-12, rtol=0)
                           and np.isclose(summary_seed.get('NMI', np.nan), nmi, atol=1e-12, rtol=0))
        for name, column in [('cluster_validation_summary.tsv', 'hicat_fine_baseline'),
                             ('coarse_cluster_validation_summary.tsv', 'hicat_coarse_baseline')]:
            table = pd.read_csv(self.out / ('annotation_review/' + name), sep='\t').set_index('cluster')
            counts = obs[column].value_counts().sort_index()
            self.check('review_membership_counts_' + name, set(table.index) == set(counts.index)
                       and table.n_cells.reindex(counts.index).equals(counts.rename('n_cells')))

    def verify_full_sidecar(self, pilot_obs):
        """Verify full display scope using small metadata, never full expression."""
        path = self.out / 'full_data_projection/canonical_full_data.h5ad'
        summary = json.loads((self.out / 'full_data_projection/full_data_projection_summary.json').read_text())
        with h5py.File(path, 'r') as handle:
            shape = _shape(handle['X'])
            cells, columns = _frame_metadata(handle['obs'])
            genes, _ = _frame_metadata(handle['var'])
            selected = np.asarray(_element(handle['obs']['is_pilot']), dtype=bool)
            uns = _element(handle['uns'])
            self.check('full_sidecar_exact_shape', shape == [EXPECTED_FULL_CELLS, EXPECTED_CANONICAL_GENES], shape)
            self.check('full_sidecar_unique_ids', len(set(cells)) == len(cells) and len(set(genes)) == len(genes))
            self.check('full_sidecar_exact_12000_pilot_cells', int(selected.sum()) == EXPECTED_PILOT_CELLS)
            self.check('full_sidecar_pilot_id_set', set(np.asarray(cells)[selected]) == set(pilot_obs.index))
            self.check('full_sidecar_no_clustering_labels', not any('hicat' in name.lower() or 'cluster' in name.lower() for name in columns), columns)
            self.check('full_sidecar_no_new_graphs', not len(handle.get('obsp', {})))
            self.check('full_sidecar_no_layers_or_raw', not len(handle.get('layers', {})) and 'raw' not in handle)
            allowed = {'X_umap_step06_display', 'X_validation_program_scores', 'X_validation_signal_detected_genes'}
            self.check('full_sidecar_only_declared_embeddings', set(handle['obsm']) == allowed, list(handle['obsm']))
            self.check('full_sidecar_UMAP_shape', list(handle['obsm']['X_umap_step06_display'].shape) == [EXPECTED_FULL_CELLS, 2])
            names = uns['validation_program_names']
            for field in ['X_validation_program_scores', 'X_validation_signal_detected_genes']:
                self.check('full_sidecar_' + field + '_shape', list(handle['obsm'][field].shape) == [EXPECTED_FULL_CELLS, len(names)])
            declared = json.loads(uns['asset_inventory_json'])
            self.check('full_sidecar_inventory_exact', _json_equal(declared, summary['asset_inventory']))
            self.check('full_sidecar_summary_dimensions', summary['n_cells'] == shape[0] and summary['n_selected_genes'] == shape[1]
                       and summary['n_pilot_cells'] == EXPECTED_PILOT_CELLS)
            self.check('full_sidecar_source_scope', not summary['full_data_clustering'] and not summary['label_transfer']
                       and not summary['umap_refitted'] and not summary['batch_correction'] and not summary['regression'])

    def verify_figures(self, pilot_obs, summary):
        """Check every figure/source reference and one page for each fine cluster."""
        index = pd.read_csv(self.out / 'figures/figure_index.tsv', sep='\t')
        self.check('figure_names_unique', index.name.is_unique)
        self.check('figure_summary_count', len(index) == int(summary['individual_figures']))
        self.check('detailed_page_sequence', index.detailed_page.tolist() == list(range(1, len(index) + 1)))
        self.check('main_page_sequence', index.main_page.dropna().astype(int).tolist() == list(range(1, int(summary['main_pages']) + 1)))
        references = set()
        for row in index.itertuples():
            sources = str(row.source_tables).split(';')
            self.check('figure_has_numeric_sources:' + row.name, bool(sources) and all(s not in ['', 'nan'] for s in sources))
            references.update([row.png, row.pdf] + sources)
        for relative in sorted(references):
            try:
                path = self._inside_output(relative)
                self.check('figure_or_source:' + relative, path.is_file() and path.stat().st_size > 0)
            except (OSError, ValueError) as error:
                self.check('figure_or_source:' + relative, False, str(error))
        for filename in ['hicat_hierarchy_validation_main.pdf', 'hicat_hierarchy_validation_detailed.pdf']:
            path = self.out / 'figures' / filename
            with path.open('rb') as handle:
                signature = handle.read(5)
            self.check('PDF_signature:' + filename, signature == b'%PDF-')
        pages = pd.read_csv(self.out / 'annotation_review/fine_review_page_index.tsv', sep='\t')
        fine_counts = pilot_obs.hicat_fine_baseline.value_counts()
        self.check('fine_review_38_unique_pages', len(pages) == EXPECTED_FINE and pages.cluster.is_unique)
        self.check('fine_review_exact_cluster_set', set(pages.cluster) == set(fine_counts.index))
        self.check('fine_review_header_counts', np.array_equal(pages.n_cells, fine_counts.reindex(pages.cluster).to_numpy()))
        self.check('fine_review_parent_headers', pages.cluster.str.split('.').str[0].equals(pages.parent))
        indexed_fine = set(index.loc[index.name.str.startswith('fine_review_'), 'name'])
        self.check('fine_review_figure_index_complete', indexed_fine == {'fine_review_' + name for name in fine_counts.index})
        for relative in pages.figure:
            self.check('fine_review_page:' + relative, self._inside_output(relative).is_file())

    def run(self):
        """Run all independent checks and return a serializable verification record."""
        manifested = 0
        try:
            self.check('finalized_outputs_directory', self.out.is_dir())
            manifested = self.verify_manifest()
            summary = json.loads((self.out / 'summary.json').read_text())
            status = json.loads((self.out / 'STEP_STATUS.json').read_text())
            self.check('checkpoint_IN_REVIEW', status['status'] == 'IN_REVIEW'
                       and not status['full_data_clustering_started'] and not status['annotations_locked'])
            obs = self.verify_pilot()
            self.verify_phase(obs)
            self.verify_partitions(obs, summary)
            self.verify_full_sidecar(obs)
            self.verify_figures(obs, summary)
            original_checks = pd.read_csv(self.out / 'validation/checks.tsv', sep='\t')
            self.check('workflow_validation_checks_all_pass', len(original_checks) > 0 and original_checks.status.eq('PASS').all())
        except Exception as error:
            self.check('verification_exception', False, dict(type=type(error).__name__, message=str(error)))
        failed = [row for row in self.checks if row['status'] != 'PASS']
        return dict(status='PASS' if not failed else 'FAIL', run_dir=str(self.root),
                    verifier='independent_read_only_step07_package_verifier_v1',
                    elapsed_seconds=time.monotonic() - self.started, manifest_files=manifested,
                    checks_passed=len(self.checks) - len(failed), checks_failed=len(failed),
                    source_raw_counts_read=False, scientific_refitting_performed=False,
                    writes_performed=False, failed_checks=failed, checks=self.checks)


def main():
    """Print verification JSON only; report failure via exit status without writes."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run-dir', required=True, type=Path)
    args = parser.parse_args()
    result = PackageVerifier(args.run_dir).run()
    print(json.dumps(result, indent=2, default=_native), flush=True)
    return 0 if result['status'] == 'PASS' else 1


if __name__ == '__main__':
    sys.exit(main())
