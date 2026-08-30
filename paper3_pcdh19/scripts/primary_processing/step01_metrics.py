"""Scanpy QC calculation, feature-set definition, and descriptive summaries."""

from __future__ import annotations

from dataclasses import dataclass

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

from .step01_models import Step01Settings


@dataclass(frozen=True)
class QCFeatureSetResult:
    """Describe mitochondrial and ribosomal gene flags observed in the panel."""

    mitochondrial_genes: tuple[str, ...]
    ribosomal_genes: tuple[str, ...]


class QCFeatureSetAnnotator:
    """Derive QC gene sets from delivered Cell Ranger gene symbols."""

    def __init__(self, settings: Step01Settings):
        """Store configured symbol prefixes without hard-coding human syntax."""

        self.settings = settings

    def annotate(self, adata: ad.AnnData) -> QCFeatureSetResult:
        """Add boolean `mt` and `ribo` flags to `var` and return their genes."""

        if "gene_symbol" not in adata.var.columns:
            raise ValueError("Step 01 requires var['gene_symbol'] from the approved Step 00 checkpoint")
        symbols = adata.var["gene_symbol"].astype(str)
        mt = symbols.str.startswith(self.settings.mitochondrial_prefix)
        ribo = symbols.str.startswith(self.settings.ribosomal_prefixes)
        adata.var["mt"] = mt.to_numpy(dtype=bool)
        adata.var["ribo"] = ribo.to_numpy(dtype=bool)
        return QCFeatureSetResult(
            mitochondrial_genes=tuple(symbols.loc[mt]),
            ribosomal_genes=tuple(symbols.loc[ribo]),
        )


class ScanpyQCMetricsCalculator:
    """Apply Scanpy's standard QC calculation without changing expression."""

    def __init__(self, settings: Step01Settings):
        """Store the explicit percent-top and gene-set calculation contract."""

        self.settings = settings

    def calculate(self, adata: ad.AnnData, feature_sets: QCFeatureSetResult) -> tuple[str, ...]:
        """Run `sc.pp.calculate_qc_metrics` and return the QC variables used."""

        qc_vars = ["mt"]
        if feature_sets.ribosomal_genes:
            qc_vars.append("ribo")
        sc.pp.calculate_qc_metrics(
            adata,
            qc_vars=qc_vars,
            percent_top=self.settings.percent_top,
            layer=None,
            use_raw=False,
            inplace=True,
            log1p=True,
        )
        return tuple(qc_vars)


class QCSummaryBuilder:
    """Build pooled, sample-level, and design-group descriptive QC tables."""

    METRICS = (
        "total_counts",
        "n_genes_by_counts",
        "total_counts_mt",
        "pct_counts_mt",
    )
    QUANTILES = (0.01, 0.05, 0.25, 0.50, 0.75, 0.95, 0.99)

    def pooled(self, obs: pd.DataFrame) -> pd.DataFrame:
        """Summarize QC distributions across every retained cell."""

        return self._summarize_groups(obs, None)

    def by_sample(self, obs: pd.DataFrame) -> pd.DataFrame:
        """Summarize QC distributions separately for all 12 samples."""

        return self._summarize_groups(obs, "technical_sample_id")

    def by_design_group(self, obs: pd.DataFrame) -> pd.DataFrame:
        """Summarize QC distributions for the four registered design groups."""

        return self._summarize_groups(obs, "design_group")

    def _summarize_groups(self, obs: pd.DataFrame, group_column: str | None) -> pd.DataFrame:
        """Calculate non-thresholding distribution statistics for each group."""

        grouped = [("all", obs)] if group_column is None else list(obs.groupby(group_column, observed=True, sort=True))
        rows: list[dict[str, object]] = []
        for group, frame in grouped:
            row: dict[str, object] = {"grouping": group_column or "pooled", "group": str(group), "n_cells": len(frame)}
            for metric in self.METRICS:
                values = frame[metric].to_numpy(dtype=float)
                row[f"{metric}_min"] = float(np.min(values))
                for quantile, value in zip(self.QUANTILES, np.quantile(values, self.QUANTILES)):
                    row[f"{metric}_p{int(quantile * 100):02d}"] = float(value)
                row[f"{metric}_mean"] = float(np.mean(values))
                row[f"{metric}_max"] = float(np.max(values))
            rows.append(row)
        return pd.DataFrame(rows)


class QCMetadataDictionary:
    """Describe the observation and variable fields added by Scanpy Step 01."""

    CORE_DESCRIPTIONS = {
        "total_counts": "Total raw Cell Ranger counts across all retained genes for this cell.",
        "log1p_total_counts": "Natural log1p of total_counts; metadata only, not a transformed expression matrix.",
        "n_genes_by_counts": "Number of genes with a nonzero raw count in this cell.",
        "log1p_n_genes_by_counts": "Natural log1p of n_genes_by_counts; metadata only.",
        "total_counts_mt": "Raw counts assigned to the 13 panel genes whose symbols begin with mt-.",
        "log1p_total_counts_mt": "Natural log1p of total_counts_mt; metadata only.",
        "pct_counts_mt": "100 × total_counts_mt / total_counts.",
        "mt": "True for delivered gene symbols beginning with mt-.",
        "ribo": "True for delivered symbols beginning with Rpl or Rps; no such genes occur in this panel.",
    }

    def build(
        self,
        before_obs: tuple[str, ...],
        before_var: tuple[str, ...],
        adata: ad.AnnData,
    ) -> pd.DataFrame:
        """Return a dictionary for every field newly created by Step 01."""

        rows = []
        for axis, before, frame in (("obs", before_obs, adata.obs), ("var", before_var, adata.var)):
            for column in frame.columns:
                if column in before:
                    continue
                rows.append({
                    "axis": axis,
                    "field": column,
                    "dtype": str(frame[column].dtype),
                    "description": self.CORE_DESCRIPTIONS.get(column, "Standard field created by scanpy.pp.calculate_qc_metrics()."),
                })
        return pd.DataFrame(rows)
