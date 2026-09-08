"""Explain sample structure using descriptive summaries of full-cell results.

``Step06MetricBuilder.build`` combines preserved raw-count QC metadata with
PCA/UMAP coordinates, provisional states, the unintegrated distance graph,
and sample pseudobulks. It returns tables; it does not change expression.
``NeighborhoodDiversity.calculate`` measures sample-label entropy on graph
edges. ``_sample_decisions`` supplies per-sample flags; ``_outcome`` assigns
the report's automatic headline using hard-coded heuristic thresholds.

These flags do not establish technical causation or absence of batch effects.
No inferential tests, null-model calibration, integration, or filtering occur.
The grouping field is the registered technical sample ID by default; calling
it a biological sample does not establish embryo/litter independence.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.spatial.distance import pdist, squareform
from sklearn.metrics import pairwise_distances

from .step06_models import Step06Artifacts, Step06Results, Step06Settings


def _long_square(matrix: np.ndarray, labels: list[str], value_name: str) -> pd.DataFrame:
    """Convert a labeled square matrix to a stable long-form table.

    Parameters
    ----------
    matrix : numpy.ndarray
        Square matrix aligned to labels on both axes.
    labels : list[str]
        Sample order corresponding to both matrix axes.
    value_name : str
        Name of the numeric output column.

    Returns
    -------
    pandas.DataFrame
        sample_1, sample_2 and value_name columns; diagonal and both pair directions
        are retained.

    Notes
    -----
    Reshape only; no distance or correlation is calculated here.
    """

    return pd.DataFrame(matrix, index=labels, columns=labels).rename_axis(
        "sample_1"
    ).reset_index().melt(id_vars="sample_1", var_name="sample_2", value_name=value_name)


class NeighborhoodDiversity:
    """Measure local biological-sample diversity on the unintegrated graph.

    Notes
    -----
    Own graph-neighborhood sample diversity metrics. The distance graph supplies
    adjacency only; registered sample labels supply the categories. Output metrics
    are descriptive and are not corrected for unequal sample sizes or biological
    composition.
    """

    def calculate(
        self, distances: sparse.csr_matrix, samples: pd.Series
    ) -> tuple[pd.DataFrame, pd.DataFrame]:
        """Return per-cell entropy/effective-sample metrics and summaries.

        Parameters
        ----------
        distances : scipy.sparse.csr_matrix
            Full-cell distance graph with rows aligned to samples.
        samples : pandas.Series
            Registered sample labels, indexed by cell ID, in graph row order.

        Returns
        -------
        tuple[pandas.DataFrame, pandas.DataFrame]
            Cell-indexed diversity metrics and an all-cell describe() summary with
            5/25/50/75/95 percentiles.

        Notes
        -----
        Treat each STORED distance-graph entry as an equally weighted neighbor;
        distance magnitude and fuzzy connectivity weight are ignored. Do not assume
        exactly k neighbors: the stored graph can omit self or zero-distance entries.
        For sample proportions p_s, H=-sum(p_s*ln(p_s)); normalized H=H/ln(S), effective
        samples=exp(H), inverse Simpson=1/sum(p_s**2), where S is the global number
        of registered samples. Normalized H=1 denotes equal representation, NOT
        mixing in proportion to unequal input sample sizes. No null model or sample
        size adjustment is applied. Requires multiple samples and nonempty graph
        rows for scientific interpretation; empty rows are numerically protected but
        the resulting effective-sample statistics would not be meaningful.
        """

        categories = sorted(samples.astype(str).unique())
        codes = pd.Categorical(samples.astype(str), categories=categories).codes
        adjacency = distances.copy().tocsr()
        # Count stored neighbors equally; remove distance magnitudes without modifying the input graph.
        adjacency.data = np.ones_like(adjacency.data, dtype=np.float32)
        one_hot = sparse.csr_matrix(
            (
                np.ones(len(codes), dtype=np.float32),
                (np.arange(len(codes)), codes),
            ),
            shape=(len(codes), len(categories)),
        )
        counts = (adjacency @ one_hot).toarray()
        totals = counts.sum(axis=1)
        proportions = counts / np.maximum(totals[:, None], 1)
        log_proportions = np.zeros_like(proportions)
        np.log(proportions, out=log_proportions, where=proportions > 0)
        # Shannon entropy with natural logs; 0*log(0) is represented by zero.
        entropy = -(proportions * log_proportions).sum(axis=1)
        inverse_simpson = 1.0 / np.maximum((proportions**2).sum(axis=1), 1e-12)
        per_cell = pd.DataFrame(
            {
                "local_neighbor_count": totals.astype(int),
                "local_sample_entropy": entropy,
                "local_sample_entropy_normalized": entropy / np.log(len(categories)),
                "local_effective_samples": np.exp(entropy),
                "local_inverse_simpson": inverse_simpson,
            },
            index=samples.index,
        )
        summary = per_cell.describe(percentiles=[0.05, 0.25, 0.5, 0.75, 0.95]).T
        summary.insert(0, "scope", "all_cells")
        summary.insert(1, "metric", summary.index)
        return per_cell, summary.reset_index(drop=True)


class Step06MetricBuilder:
    """Build every numerical companion table required by the review contract.

    Notes
    -----
    Own full-cell companion tables and automatic reporting heuristics. build()
    combines preserved QC metadata with analysis artifacts; _sample_decisions()
    flags individual samples and _outcome() separately assigns the global headline.
    None of these operations filters cells.
    """

    def __init__(self, settings: Step06Settings):
        """Store immutable settings.

        Parameters
        ----------
        settings : Step06Settings
            Resolved scientific/rendering controls; see the settings class and tuning
            guide.

        Notes
        -----
        Store grouping names and centroid-component count; no summaries are computed.
        """

        self.settings = settings

    def build(self, obs: pd.DataFrame, artifacts: Step06Artifacts) -> Step06Results:
        """Calculate composition, QC, sample similarity, mixing, and decisions.

        Parameters
        ----------
        obs : pandas.DataFrame
            All retained cells in artifact row order, with sample/design and raw QC
            columns.
        artifacts : Step06Artifacts
            Full-cell coordinates, graphs, scores and HVG/pseudobulk results in input
            order.

        Returns
        -------
        Step06Results
            Sixteen filename-keyed tables and the heuristic outcome code/label/rationale.

        Notes
        -----
        Use ALL cells. Preserve the original raw-count QC fields; do not recompute QC
        on log expression. Compute sample means in every PC, but Euclidean sample
        center distances in centroid_components PCs (20 by default). UMAP centers
        are descriptive display-space summaries. Within-sample state percentages
        sum across states; within-state percentages sum across samples.

        Pearson pseudobulk correlations use HVG-restricted log1p CPM. Mixing uses
        distance-graph edges, not UMAP geometry. _sample_decisions and _outcome apply
        separate hard-coded rules; the global headline can differ from individual
        sample flags. No inferential tests occur (recorded explicitly in a table).
        Grouping uses registered technical sample IDs; independence is not inferred.
        """

        sample = self.settings.sample_field
        order = sorted(obs[sample].astype(str).unique())
        metadata = obs.groupby(sample, observed=True).agg(
            submitted_sample_name=("submitted_sample_name", "first"),
            genotype=(self.settings.genotype_field, "first"),
            sex=(self.settings.sex_field, "first"),
            design_group=(self.settings.design_field, "first"),
            cells=(sample, "size"),
        ).reindex(order).reset_index()

        qc = obs.groupby(sample, observed=True).agg(
            cells=("total_counts", "size"),
            total_counts_mean=("total_counts", "mean"),
            total_counts_median=("total_counts", "median"),
            total_counts_q25=("total_counts", lambda x: x.quantile(0.25)),
            total_counts_q75=("total_counts", lambda x: x.quantile(0.75)),
            n_genes_mean=("n_genes_by_counts", "mean"),
            n_genes_median=("n_genes_by_counts", "median"),
            n_genes_q25=("n_genes_by_counts", lambda x: x.quantile(0.25)),
            n_genes_q75=("n_genes_by_counts", lambda x: x.quantile(0.75)),
            pct_mt_mean=("pct_counts_mt", "mean"),
            pct_mt_median=("pct_counts_mt", "median"),
            pct_mt_q25=("pct_counts_mt", lambda x: x.quantile(0.25)),
            pct_mt_q75=("pct_counts_mt", lambda x: x.quantile(0.75)),
        ).reindex(order).reset_index()

        pca_columns = [f"PC{i + 1}" for i in range(artifacts.pca.shape[1])]
        pca_frame = pd.DataFrame(artifacts.pca, index=obs.index, columns=pca_columns)
        pca_frame[sample] = obs[sample].astype(str).to_numpy()
        centroids = pca_frame.groupby(sample, observed=True)[pca_columns].mean().reindex(order)
        centroid_table = centroids.reset_index()
        # Numeric distances use the first centroid_components PCs, not just the two displayed PCs.
        centroid_distance = pairwise_distances(
            centroids.iloc[:, : self.settings.centroid_components], metric="euclidean"
        )
        centroid_long = _long_square(centroid_distance, order, "pca_centroid_distance")

        umap_columns = ["UMAP1", "UMAP2"]
        umap_frame = pd.DataFrame(artifacts.umap, index=obs.index, columns=umap_columns)
        umap_frame[sample] = obs[sample].astype(str).to_numpy()
        umap_centroids = umap_frame.groupby(sample, observed=True)[umap_columns].mean().reindex(order)

        states = artifacts.annotations["step06_provisional_state"]
        state_order = sorted(states.unique())
        state_counts = pd.crosstab(obs[sample].astype(str), states).reindex(
            index=order, columns=state_order, fill_value=0
        )
        # Denominators differ: composition of a sample vs contribution of samples to one state.
        within_sample = state_counts.div(state_counts.sum(axis=1), axis=0) * 100
        within_state = state_counts.div(state_counts.sum(axis=0), axis=1) * 100
        state_count_long = state_counts.rename_axis(sample).reset_index().melt(
            id_vars=sample, var_name="provisional_state", value_name="cells"
        )
        sample_pct_long = within_sample.rename_axis(sample).reset_index().melt(
            id_vars=sample, var_name="provisional_state", value_name="within_sample_pct"
        )
        state_pct_long = within_state.rename_axis(sample).reset_index().melt(
            id_vars=sample, var_name="provisional_state", value_name="within_state_pct"
        )

        diversity, diversity_summary = NeighborhoodDiversity().calculate(
            artifacts.distances, obs[sample]
        )
        diversity[sample] = obs[sample].astype(str).to_numpy()
        diversity["provisional_state"] = states.to_numpy()
        by_state = diversity.groupby("provisional_state", observed=True).agg(
            cells=("local_effective_samples", "size"),
            effective_samples_mean=("local_effective_samples", "mean"),
            effective_samples_median=("local_effective_samples", "median"),
            effective_samples_q05=("local_effective_samples", lambda x: x.quantile(0.05)),
            effective_samples_q95=("local_effective_samples", lambda x: x.quantile(0.95)),
            normalized_entropy_mean=("local_sample_entropy_normalized", "mean"),
            inverse_simpson_mean=("local_inverse_simpson", "mean"),
        ).reset_index()

        correlations = artifacts.pseudobulk_logcpm.T.corr(method="pearson").reindex(
            index=order, columns=order
        )
        correlation_long = _long_square(correlations.to_numpy(), order, "pseudobulk_pearson_r")
        decision = self._sample_decisions(
            qc, centroids, umap_centroids, state_counts, diversity, order
        )
        outcome_code, outcome_label, rationale = self._outcome(
            centroid_distance, state_counts, diversity
        )
        parameters = pd.DataFrame(
            [{"parameter": key, "value": value} for key, value in vars(self.settings).items()]
        )
        tables = {
            "sample_metadata_and_cells.tsv": metadata,
            "qc_summary_by_sample.tsv": qc,
            "pca_centroids.tsv": centroid_table,
            "pca_centroid_distances.tsv": centroid_long,
            "umap_centroids.tsv": umap_centroids.reset_index(),
            "pseudobulk_logcpm_hvgs.tsv.gz": artifacts.pseudobulk_logcpm.reset_index(names=sample),
            "pseudobulk_correlations.tsv": correlation_long,
            "provisional_state_counts.tsv": state_count_long,
            "provisional_state_within_sample_percentages.tsv": sample_pct_long,
            "provisional_state_within_state_percentages.tsv": state_pct_long,
            "local_sample_diversity_per_cell.tsv.gz": diversity.reset_index(names="cell_id"),
            "local_sample_diversity_summary.tsv": diversity_summary,
            "local_sample_diversity_by_state.tsv": by_state,
            "preintegration_sample_decisions.tsv": decision,
            "analysis_parameters.tsv": parameters,
            "statistical_comparisons.tsv": pd.DataFrame([{
                "comparison": "none",
                "scope": "descriptive_diagnostics_only",
                "result": "No inferential statistical tests were used in the report.",
            }]),
        }
        return Step06Results(tables, outcome_code, outcome_label, rationale)

    def _sample_decisions(
        self,
        qc: pd.DataFrame,
        centroids: pd.DataFrame,
        umap_centroids: pd.DataFrame,
        state_counts: pd.DataFrame,
        diversity: pd.DataFrame,
        order: list[str],
    ) -> pd.DataFrame:
        """Create a transparent per-sample review summary without integrating.

        Parameters
        ----------
        qc : pandas.DataFrame
            Sample-level original QC medians and sample identifier.
        centroids : pandas.DataFrame
            Sample-indexed PC means; first centroid_components columns define
            displacement.
        umap_centroids : pandas.DataFrame
            Sample-indexed UMAP1/UMAP2 means.
        state_counts : pandas.DataFrame
            Samples by provisional states, raw cell counts.
        diversity : pandas.DataFrame
            Full-cell sample labels and local_effective_samples.
        order : list[str]
            Stable sample order for output rows.

        Returns
        -------
        pandas.DataFrame
            One row per sample with QC, displacement, state representation, review flags
            and unresolved interpretation.

        Notes
        -----
        Code-only rules: flag PCA/UMAP displacement if distance from the equally
        weighted mean of sample centers exceeds median + 3*UNSCALED MAD across
        samples. This differs from Step 02's 1.4826-scaled MAD QC filtering.
        A represented state needs >=100 cells AND >=1% of that sample. A missing-state
        flag requires >=1% globally AND (<50 cells OR <0.1% of that sample).
        These literal thresholds trigger review only; none filter cells or prove
        technical causation. UMAP displacement is sensitive to embedding geometry.
        """

        center = centroids.iloc[:, : self.settings.centroid_components].mean(axis=0)
        displacement = np.sqrt(
            ((centroids.iloc[:, : self.settings.centroid_components] - center) ** 2).sum(axis=1)
        )
        # CODE-ONLY rule: median + 3 unscaled MAD. This is NOT the Step 02 QC MAD formula.
        cutoff = float(np.median(displacement) + 3 * np.median(np.abs(displacement - np.median(displacement))))
        umap_center = umap_centroids.mean(axis=0)
        umap_displacement = np.sqrt(((umap_centroids - umap_center) ** 2).sum(axis=1))
        umap_cutoff = float(
            np.median(umap_displacement)
            + 3 * np.median(np.abs(umap_displacement - np.median(umap_displacement)))
        )
        global_fraction = state_counts.sum(axis=0) / state_counts.values.sum()
        rows = []
        for sample_id in order:
            counts = state_counts.loc[sample_id]
            # CODE-ONLY thresholds: >=100 AND >=1%; missing uses a different rule below.
            represented = counts[(counts >= 100) & (counts / counts.sum() >= 0.01)].index.tolist()
            missing = [
                state
                for state in state_counts.columns
                if global_fraction[state] >= 0.01 and (counts[state] < 50 or counts[state] / counts.sum() < 0.001)
            ]
            local = diversity.loc[diversity[self.settings.sample_field] == sample_id]
            rows.append(
                {
                    self.settings.sample_field: sample_id,
                    "cells": int(counts.sum()),
                    "median_total_counts": float(qc.set_index(self.settings.sample_field).loc[sample_id, "total_counts_median"]),
                    "median_n_genes": float(qc.set_index(self.settings.sample_field).loc[sample_id, "n_genes_median"]),
                    "median_pct_mt": float(qc.set_index(self.settings.sample_field).loc[sample_id, "pct_mt_median"]),
                    "pca_centroid_displacement": float(displacement.loc[sample_id]),
                    "major_pca_displacement": bool(displacement.loc[sample_id] > cutoff),
                    "umap_centroid_displacement": float(umap_displacement.loc[sample_id]),
                    "major_umap_displacement": bool(umap_displacement.loc[sample_id] > umap_cutoff),
                    "broad_states_represented": "; ".join(represented),
                    "evidence_of_missing_states": "; ".join(missing) if missing else "none",
                    "median_local_effective_samples": float(local["local_effective_samples"].median()),
                    "unusual_sample_specific_structure": "review" if displacement.loc[sample_id] > cutoff or umap_displacement.loc[sample_id] > umap_cutoff or missing else "not_flagged",
                    "interpretation": "unresolved; compare biology and technical covariates",
                }
            )
        return pd.DataFrame(rows)

    @staticmethod
    def _outcome(
        centroid_distance: np.ndarray,
        state_counts: pd.DataFrame,
        diversity: pd.DataFrame,
    ) -> tuple[int, str, str]:
        """Choose one transparent provisional outcome for human review.

        Parameters
        ----------
        centroid_distance : numpy.ndarray
            Symmetric sample-by-sample Euclidean distances across configured centroid PCs.
        state_counts : pandas.DataFrame
            Samples by provisional states in the same sample order as the distance matrix.
        diversity : pandas.DataFrame
            All-cell metrics containing local_sample_entropy_normalized.

        Returns
        -------
        tuple[int, str, str]
            Automatic outcome code, headline and rationale; all require human review.

        Notes
        -----
        Evaluate in this order (literal code thresholds, not CLI settings):
        1. Median normalized entropy >=0.65 AND largest sample-center distance <=2
           times median off-diagonal distance -> outcome 1, "Minimal...".
        2. Otherwise Pearson r between sample-center distances and Jensen-Shannon
           distances of state proportions >=0.5 -> outcome 3, composition/confounding.
        3. Otherwise -> outcome 2, whose current headline says "technical structure".

        These are UNVALIDATED DESCRIPTIVE HEURISTICS. Rule 1 measures relative spread,
        not absolute separation, and ignores sample-size-adjusted expected mixing.
        Pairwise distances are not independent replicates. No p-value is calculated.
        Rule 3 does not establish technical causation merely because other conditions
        failed; its wording is stronger than its evidence. The first passing branch
        wins, even if the composition correlation is high. Never treat any headline
        as approval, a calibrated classifier, or an integration decision.
        """

        composition = state_counts.div(state_counts.sum(axis=1), axis=0)
        composition_distance = squareform(pdist(composition, metric="jensenshannon"))
        triangle = np.triu_indices_from(centroid_distance, k=1)
        association = float(np.corrcoef(centroid_distance[triangle], composition_distance[triangle])[0, 1])
        median_mixing = float(diversity["local_sample_entropy_normalized"].median())
        # CODE-ONLY heuristic branch 1. Relative distance ratio does not measure absolute separation.
        if median_mixing >= 0.65 and np.nanmax(centroid_distance) <= 2 * np.nanmedian(centroid_distance[triangle]):
            return 1, "Minimal sample-associated structure", f"Median normalized local entropy={median_mixing:.3f}; centroid dispersion was limited."
        # Branch 2 only runs if branch 1 failed; correlation is descriptive, with no p-value.
        if np.isfinite(association) and association >= 0.5:
            return 3, "Strong sample differences potentially caused by biological composition or confounding", f"PCA-centroid and state-composition distances were correlated (r={association:.3f}); integration should wait for review."
        # Fallback headline overstates evidence: failing the other rules does NOT prove technical cause.
        return 2, "Clear sample-associated technical structure superimposed on shared biology", f"Median normalized local entropy={median_mixing:.3f}; centroid separation was not primarily explained by provisional-state composition."
