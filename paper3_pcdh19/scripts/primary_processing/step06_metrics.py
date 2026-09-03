"""Machine-readable sample, state, neighborhood, and decision diagnostics."""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.spatial.distance import pdist, squareform
from sklearn.metrics import pairwise_distances

from .step06_models import Step06Artifacts, Step06Results, Step06Settings


def _long_square(matrix: np.ndarray, labels: list[str], value_name: str) -> pd.DataFrame:
    """Convert a labeled square matrix to a stable long-form table."""

    return pd.DataFrame(matrix, index=labels, columns=labels).rename_axis(
        "sample_1"
    ).reset_index().melt(id_vars="sample_1", var_name="sample_2", value_name=value_name)


class NeighborhoodDiversity:
    """Measure local biological-sample diversity on the unintegrated graph."""

    def calculate(
        self, distances: sparse.csr_matrix, samples: pd.Series
    ) -> tuple[pd.DataFrame, pd.DataFrame]:
        """Return per-cell entropy/effective-sample metrics and summaries."""

        categories = sorted(samples.astype(str).unique())
        codes = pd.Categorical(samples.astype(str), categories=categories).codes
        adjacency = distances.copy().tocsr()
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
    """Build every numerical companion table required by the review contract."""

    def __init__(self, settings: Step06Settings):
        """Store immutable settings."""

        self.settings = settings

    def build(self, obs: pd.DataFrame, artifacts: Step06Artifacts) -> Step06Results:
        """Calculate composition, QC, sample similarity, mixing, and decisions."""

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
        """Create a transparent per-sample review summary without integrating."""

        center = centroids.iloc[:, : self.settings.centroid_components].mean(axis=0)
        displacement = np.sqrt(
            ((centroids.iloc[:, : self.settings.centroid_components] - center) ** 2).sum(axis=1)
        )
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
        """Choose one transparent provisional outcome for human review."""

        composition = state_counts.div(state_counts.sum(axis=1), axis=0)
        composition_distance = squareform(pdist(composition, metric="jensenshannon"))
        triangle = np.triu_indices_from(centroid_distance, k=1)
        association = float(np.corrcoef(centroid_distance[triangle], composition_distance[triangle])[0, 1])
        median_mixing = float(diversity["local_sample_entropy_normalized"].median())
        if median_mixing >= 0.65 and np.nanmax(centroid_distance) <= 2 * np.nanmedian(centroid_distance[triangle]):
            return 1, "Minimal sample-associated structure", f"Median normalized local entropy={median_mixing:.3f}; centroid dispersion was limited."
        if np.isfinite(association) and association >= 0.5:
            return 3, "Strong sample differences potentially caused by biological composition or confounding", f"PCA-centroid and state-composition distances were correlated (r={association:.3f}); integration should wait for review."
        return 2, "Clear sample-associated technical structure superimposed on shared biology", f"Median normalized local entropy={median_mixing:.3f}; centroid separation was not primarily explained by provisional-state composition."
