"""Independent canonical program scoring for dissected E14.5 mouse MGE.

The scorer receives measured expression and a literature-based configuration,
never cluster labels or a DEG table. It fits expression-bin matched control
sets once on the 12,000-cell pilot, records every coefficient, and reuses the
same linear score model for full-data display. A score is a mean difference
on natural-log(1 + counts per million), not a cell identity probability.

Tuning is centralized in ``config/hicat_validation_programs.json``: gene lists,
source references, control-bin count, random seed, minimum assay coverage and
an explicitly operational positive-score rule. Changing these settings creates
a new review run; no scores are regressed from expression or used to fit HiCAT.
"""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy import sparse
from .provenance import write_json


def normalize_log1p_cpm(raw, target=1000000.0):
    """Return a new float64 CSR matrix normalized over every measured gene.

    Parameters
    ----------
    raw : scipy sparse matrix
        Unmodified integer counts, cells by all measured genes.
    target : float
        Counts-per-million denominator target; default 1,000,000.

    Returns
    -------
    scipy.sparse.csr_matrix
        Natural-log(1 + target * count / total measured counts). Empty cells
        retain zeros. The caller's count matrix is never changed.
    """
    x = sparse.csr_matrix(raw, dtype=np.float64).copy()
    totals = np.asarray(x.sum(axis=1)).ravel()
    multiplier = np.divide(target, totals, out=np.zeros_like(totals), where=totals > 0)
    x = sparse.diags(multiplier).dot(x).tocsr()
    np.log1p(x.data, out=x.data)
    return x


class CanonicalProgramScorer:
    """Fit auditable matched-control signatures independently of cluster labels.

    Parameters
    ----------
    config_path : path or dict
        Frozen program configuration. Dictionary input is useful in tests.

    Notes
    -----
    Symbols are matched to ``adata.var['gene_symbol']`` case-insensitively;
    three known renamed cell-cycle symbols have explicit aliases. Ambiguous
    symbols matching multiple features are excluded and recorded rather than
    silently choosing a feature. This is conserved-symbol matching, not a new
    orthology mapping. Unavailable programs retain NaN scores, not zeros.
    """

    def __init__(self, config_path):
        """Load settings without reading expression, annotations, or DEG files."""
        self.config = dict(config_path) if isinstance(config_path, dict) else json.loads(Path(config_path).read_text())
        self.settings = self.config['scoring']
        self.names = [p['name'] for p in self.config['programs']]
        self.weights = None
        self.identity_weights = None
        self.result = None

    def _resolve(self, requested, symbols):
        """Resolve one configured signature and report every missing/alias match."""
        aliases = self.config.get('cycle_reference', {}).get('renamed_aliases', {})
        lower_to_indices = {}
        for index, symbol in enumerate(symbols):
            lower_to_indices.setdefault(str(symbol).lower(), []).append(index)
        indices, audit = [], []
        for symbol in requested:
            resolved = aliases.get(symbol, symbol)
            hits = lower_to_indices.get(resolved.lower(), [])
            status = 'available' if len(hits) == 1 else ('absent_from_assay' if not hits else 'ambiguous_multiple_features')
            if len(hits) == 1:
                indices.append(hits[0])
            audit.append(dict(requested_symbol=symbol, lookup_symbol=resolved,
                              measured_symbol=str(symbols[hits[0]]) if len(hits) == 1 else '',
                              feature_index=int(hits[0]) if len(hits) == 1 else -1,
                              status=status, mapping='documented_renamed_alias' if symbol in aliases else 'case_insensitive_symbol'))
        return sorted(set(indices)), audit

    def _coefficients(self, signal, averages, ctrl_size):
        """Build Scanpy-compatible bin-matched mean-signal minus mean-control weights.

        Bins use average measured pilot expression and minimum-rank ties, as in
        Scanpy 1.9.8. The same random seed starts each signature independently.
        Controls can include other programs' genes; own signal genes are removed
        after the per-bin draws, matching the reference implementation.
        """
        weights = np.zeros(len(averages), dtype=np.float64)
        if len(signal) < self.settings['min_available_genes']:
            return weights, [], 'insufficient_assay_coverage'
        n_items = max(1, int(np.round(len(averages) / (self.settings['n_bins'] - 1))))
        cuts = np.asarray(pd.Series(averages).rank(method='min') // n_items)
        rng = np.random.RandomState(self.settings['random_seed'])
        controls = set()
        for cut in np.unique(cuts[signal]):
            candidates = np.flatnonzero(cuts == cut)
            rng.shuffle(candidates)
            controls.update(candidates[:ctrl_size].tolist())
        controls.difference_update(signal)
        controls = sorted(controls)
        if not controls:
            return weights, [], 'no_control_genes'
        weights[signal] = 1.0 / len(signal)
        weights[controls] = -1.0 / len(controls)
        return weights, controls, 'available'

    def transform(self, log_expression, identity=False):
        """Apply frozen pilot weights to identically ordered normalized genes.

        This operation estimates no new model, bin, threshold, or cluster. It
        supports small matrices or streamed full-data chunks. Scores with
        insufficient measured genes are overwritten with NaN.
        """
        if self.weights is None:
            raise RuntimeError('Call fit_transform on the pilot before transform')
        weights = self.identity_weights if identity else self.weights
        if log_expression.shape[1] != weights.shape[0]:
            raise ValueError('Expression feature count differs from fitted model')
        values = np.asarray(log_expression.dot(weights))
        valid = self.identity_valid if identity else self.valid
        values[:, ~valid] = np.nan
        return values

    def fit_transform(self, adata, output_dir):
        """Score one fixed pilot and save coefficients, coverage and cell values.

        Parameters
        ----------
        adata : anndata.AnnData
            Raw counts in X, normalized values in layers['log1p_cpm'], gene
            IDs in var_names and mouse symbols in var['gene_symbol'].
        output_dir : path
            New run's marker_programs directory, created if absent.

        Returns
        -------
        dict
            Cell-indexed DataFrames: scores, identity_scores, signal_means,
            detection_counts, positive, canonical_expression; phase Series;
            gene-level coverage table and JSON-safe model. Identity scores use
            separate signatures after removing every standard S/G2M gene and
            Ccnd1/2/3, with separately fitted controls and no explicit cycle
            genes in either signal or controls.

        Saved data
        ----------
        Coefficients and gene coverage TSV; model/config JSON; score/phase/
        expression TSV.gz; raw signal summaries; cluster-level summaries are
        intentionally delegated to the review workflow. The object receives
        validation score/phase fields, but no biological annotations.
        """
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        symbols = adata.var['gene_symbol'].astype(str).to_numpy()
        if not adata.var_names.is_unique:
            raise ValueError('Feature IDs must be unique')
        self.gene_ids = adata.var_names.astype(str).to_numpy()
        self.symbols = symbols
        x = sparse.csr_matrix(adata.layers['log1p_cpm'], dtype=np.float64)
        if np.any(x.data < 0) or not np.all(np.isfinite(x.data)):
            raise ValueError('Expected finite nonnegative log1p CPM expression')
        averages = np.asarray(x.mean(axis=0)).ravel()
        self.signal_indices, audit_rows = [], []
        for program in self.config['programs']:
            indices, rows = self._resolve(program['genes'], symbols)
            self.signal_indices.append(indices)
            for row in rows:
                row.update(program=program['name'], feature_id=self.gene_ids[row['feature_index']] if row['feature_index'] >= 0 else '')
            audit_rows.extend(rows)
        cycle_counts = [len(self.signal_indices[self.names.index(name)]) for name in ['S_phase', 'G2M']]
        cycle_control_size = min(cycle_counts)
        exclusion_indices, _ = self._resolve(self.config['identity_cycle_exclusions'], symbols)
        excluded = set(exclusion_indices)
        columns, identity_columns, models, identity_models = [], [], [], []
        valid, identity_valid = [], []
        for program, signal in zip(self.config['programs'], self.signal_indices):
            ctrl_size = cycle_control_size if program['group'] == 'cycle' else self.settings['control_size']
            weights, controls, status = self._coefficients(signal, averages, ctrl_size)
            columns.append(weights)
            valid.append(status == 'available')
            models.append(dict(program=program['name'], status=status, signal_feature_ids=self.gene_ids[signal].tolist(),
                               control_feature_ids=self.gene_ids[controls].tolist(), control_size_requested=ctrl_size,
                               n_signal=len(signal), n_controls=len(controls)))
            # Identity comparisons exclude the full cycle list from signals AND
            # controls. Merely dropping Ccnd2 from an IPC label would not remove
            # the effect if it remained in its control pool.
            identity_signal = [i for i in signal if i not in excluded]
            iw, icontrols, istatus = self._coefficients(identity_signal, averages, self.settings['control_size'])
            icontrols = [i for i in icontrols if i not in excluded]
            iw[:] = 0
            if len(identity_signal) >= self.settings['min_available_genes'] and icontrols and program['group'] not in ['cycle', 'quality']:
                iw[identity_signal] = 1.0 / len(identity_signal)
                iw[icontrols] = -1.0 / len(icontrols)
                istatus = 'available'
            else:
                istatus = 'not_an_identity_signature' if program['group'] in ['cycle', 'quality'] else 'insufficient_assay_coverage'
            identity_columns.append(iw)
            identity_valid.append(istatus == 'available')
            identity_models.append(dict(program=program['name'], status=istatus, signal_feature_ids=self.gene_ids[identity_signal].tolist(), control_feature_ids=self.gene_ids[icontrols].tolist()))
        self.weights = np.column_stack(columns)
        self.identity_weights = np.column_stack(identity_columns)
        self.valid = np.asarray(valid)
        self.identity_valid = np.asarray(identity_valid)
        scores = pd.DataFrame(self.transform(x), index=adata.obs_names, columns=self.names)
        identity_scores = pd.DataFrame(self.transform(x, identity=True), index=adata.obs_names, columns=self.names)
        means, detected = {}, {}
        for name, indices in zip(self.names, self.signal_indices):
            means[name] = np.asarray(x[:, indices].mean(axis=1)).ravel() if indices else np.full(adata.n_obs, np.nan)
            detected[name] = np.asarray((x[:, indices] > 0).sum(axis=1)).ravel() if indices else np.full(adata.n_obs, np.nan)
        means = pd.DataFrame(means, index=adata.obs_names)
        detected = pd.DataFrame(detected, index=adata.obs_names)
        positive = ((scores > self.settings['positive_score_threshold']) &
                    (detected >= self.settings['min_detected_genes_positive'])).astype(float)
        positive.loc[:, ~self.valid] = np.nan
        phase = pd.Series('S', index=adata.obs_names, name='validation_phase')
        phase[scores['G2M'] > scores['S_phase']] = 'G2/M'
        phase[(scores['G2M'] < 0) & (scores['S_phase'] < 0)] = 'G1-like'
        if scores[['S_phase','G2M']].isna().any().any():
            raise ValueError('Cell-cycle reference does not have sufficient measured coverage')
        adata.obs['validation_S_score'] = scores['S_phase']
        adata.obs['validation_G2M_score'] = scores['G2M']
        adata.obs['validation_phase'] = pd.Categorical(phase, categories=['G1-like','S','G2/M'])
        adata.obsm['X_validation_program_scores'] = scores.to_numpy(dtype=np.float32)
        adata.uns['validation_program_names'] = np.asarray(self.names, dtype=str)
        requested = list(dict.fromkeys(g for genes in self.config['dotplot_groups'].values() for g in genes))
        self.canonical_indices, canonical_audit = self._resolve(requested, symbols)
        # Keep biological display order, not genomic order.
        self.canonical_indices = [r['feature_index'] for r in canonical_audit if r['status'] == 'available']
        canonical = pd.DataFrame(x[:, self.canonical_indices].toarray(), index=adata.obs_names, columns=symbols[self.canonical_indices])
        coverage = pd.DataFrame(audit_rows)
        self.model = dict(schema_version='pcdh19_fixed_program_scores_v1', feature_ids=self.gene_ids.tolist(),
                          fitted_cells=adata.n_obs, normalization='natural log(1 + counts per million), denominator all measured genes',
                          annotation_or_cluster_fields_used=False, settings=self.settings, programs=models,
                          identity_programs=identity_models, identity_cycle_exclusion_feature_ids=self.gene_ids[exclusion_indices].tolist(),
                          phase_rule='S by default; G2/M when G2M > S; G1-like when both < 0. This is a transcriptional proxy.',
                          score_interpretation='Mean signature expression minus expression-bin matched control mean; not identity probabilities.',
                          positive_rule='score > %.3f and at least %d measured signature genes detected' % (self.settings['positive_score_threshold'], self.settings['min_detected_genes_positive']))
        write_json(output_dir/'score_model.json', self.model)
        write_json(output_dir/'canonical_program_config.json', self.config)
        coverage.to_csv(output_dir/'gene_coverage.tsv', sep='\t', index=False)
        pd.DataFrame(canonical_audit).to_csv(output_dir/'canonical_dotplot_coverage.tsv', sep='\t', index=False)
        pd.DataFrame(models).to_csv(output_dir/'program_coverage.tsv', sep='\t', index=False)
        coefficients = []
        for identity, matrix in [(False,self.weights),(True,self.identity_weights)]:
            for column, name in enumerate(self.names):
                for row in np.flatnonzero(matrix[:,column]):
                    coefficients.append(dict(model='cycle_excluded_identity' if identity else 'standard',program=name,feature_id=self.gene_ids[row],gene_symbol=symbols[row],weight=matrix[row,column]))
        pd.DataFrame(coefficients).to_csv(output_dir/'score_coefficients.tsv', sep='\t', index=False)
        for name, frame in [('cell_program_scores', scores), ('cell_identity_scores', identity_scores),
                            ('cell_signal_means',means), ('cell_signal_detected_genes',detected),
                            ('cell_program_positive',positive), ('canonical_cell_expression', canonical)]:
            frame.rename_axis('cell_id').to_csv(output_dir/(name+'.tsv.gz'), sep='\t', compression='gzip')
        pd.concat([phase, scores[['S_phase','G2M']]],axis=1).rename_axis('cell_id').to_csv(output_dir/'cell_cycle_scores.tsv.gz',sep='\t',compression='gzip')
        self.result = dict(scores=scores,identity_scores=identity_scores,signal_means=means,detection_counts=detected,
                           positive=positive,phase=phase,coverage=coverage,model=self.model,canonical_expression=canonical,
                           program_metadata=self.config['programs'])
        return self.result
