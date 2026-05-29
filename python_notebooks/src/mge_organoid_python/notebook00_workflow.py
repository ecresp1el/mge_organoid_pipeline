"""Reusable workflow pieces for Notebook 00."""

from __future__ import annotations

import time
from dataclasses import dataclass
from typing import Sequence

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import median_abs_deviation as mad


DEFAULT_QC_RULES = (
    {
        "metric": "log1p_total_counts",
        "flag_column": "qc_mad_flag_log1p_total_counts",
        "n_mads": 5,
        "upper_only": False,
    },
    {
        "metric": "pct_counts_mt",
        "flag_column": "qc_mad_flag_pct_counts_mt",
        "n_mads": 3,
        "upper_only": True,
    },
)


@dataclass
class PreprocessSettings:
    """Basic Scanpy preprocessing settings for Notebook 00."""

    target_sum: float = 1e4
    n_top_genes: int = 2000
    hvg_batch_key: str = "run_sample_id"
    hvg_flavor: str = "seurat"
    pca_solver: str = "arpack"


def calculate_qc_metrics(adata_list: Sequence[ad.AnnData], mito_prefix: str = "MT-", overwrite: bool = False) -> None:
    """Add Scanpy QC metric columns to each AnnData object's `.obs`."""
    required_obs_columns = {"total_counts", "log1p_total_counts", "pct_counts_mt"}

    for one_sample_adata in adata_list:
        missing_qc_columns = required_obs_columns.difference(one_sample_adata.obs.columns)
        if overwrite or missing_qc_columns:
            one_sample_adata.var["mt"] = one_sample_adata.var_names.str.startswith(mito_prefix)
            sc.pp.calculate_qc_metrics(
                one_sample_adata,
                qc_vars=["mt"],
                percent_top=None,
                log1p=True,
                inplace=True,
            )


def _mad_bounds(values: pd.Series | np.ndarray, n_mads: float) -> tuple[float, float, float, float]:
    finite_values = pd.Series(values, dtype="float64").replace([np.inf, -np.inf], np.nan).dropna().to_numpy()
    median_value = float(np.median(finite_values))
    mad_value = float(mad(finite_values))
    lower_bound = median_value - n_mads * mad_value
    upper_bound = median_value + n_mads * mad_value
    return median_value, mad_value, lower_bound, upper_bound


def annotate_mad_qc(
    adata_names: Sequence[str],
    adata_list: Sequence[ad.AnnData],
    rules: Sequence[dict] = DEFAULT_QC_RULES,
) -> pd.DataFrame:
    """Write 0/1 MAD QC flags into per-sample `.obs` and return thresholds."""
    threshold_records = []

    for run_sample_id, one_sample_adata in zip(adata_names, adata_list):
        combined_review_mask = np.zeros(one_sample_adata.n_obs, dtype=bool)
        sample_thresholds = []

        for rule in rules:
            metric = rule["metric"]
            flag_column = rule["flag_column"]
            n_mads = float(rule["n_mads"])
            upper_only = bool(rule["upper_only"])
            if metric not in one_sample_adata.obs.columns:
                raise KeyError(f"{metric!r} is not in .obs for {run_sample_id}. Run calculate_qc_metrics() first.")

            values = one_sample_adata.obs[metric].astype(float)
            median_value, mad_value, lower_bound, upper_bound = _mad_bounds(values, n_mads=n_mads)
            if upper_only:
                mask = values > upper_bound
            else:
                mask = (values < lower_bound) | (values > upper_bound)

            one_sample_adata.obs[flag_column] = mask.astype("int8")
            combined_review_mask |= mask.to_numpy(dtype=bool)

            record = {
                "run_sample_id": str(run_sample_id),
                "metric": metric,
                "flag_column": flag_column,
                "n_mads": n_mads,
                "upper_only": upper_only,
                "n_cells": int(one_sample_adata.n_obs),
                "median": median_value,
                "mad": mad_value,
                "lower_bound": lower_bound,
                "upper_bound": upper_bound,
            }
            threshold_records.append(record)
            sample_thresholds.append(record)

        one_sample_adata.obs["qc_mad_any_flag"] = combined_review_mask.astype("int8")
        one_sample_adata.obs["qc_mad_pass"] = ~combined_review_mask
        one_sample_adata.obs["qc_mad_status"] = pd.Categorical(
            np.where(combined_review_mask, "flagged_for_review", "passes_mad_qc"),
            categories=["passes_mad_qc", "flagged_for_review"],
            ordered=True,
        )
        one_sample_adata.uns["qc_mad_annotation_rules"] = sample_thresholds
        one_sample_adata.uns["qc_mad_annotation_note"] = "Annotation only: no cells were filtered or removed."

    return pd.DataFrame(threshold_records)


def qc_annotation_summary(adata_names: Sequence[str], adata_list: Sequence[ad.AnnData]) -> pd.DataFrame:
    """Summarize QC flags per sample."""
    records = []
    for run_sample_id, one_sample_adata in zip(adata_names, adata_list):
        if "qc_mad_any_flag" not in one_sample_adata.obs.columns:
            raise KeyError(f"Missing qc_mad_any_flag for {run_sample_id}. Run annotate_mad_qc() first.")

        flagged = one_sample_adata.obs["qc_mad_any_flag"].astype(bool)
        records.append(
            {
                "run_sample_id": str(run_sample_id),
                "data_source": str(one_sample_adata.obs["data_source"].iloc[0])
                if "data_source" in one_sample_adata.obs.columns
                else "",
                "starting_n_cells": int(one_sample_adata.n_obs),
                "flagged_n_cells": int(flagged.sum()),
                "retained_n_cells": int((~flagged).sum()),
                "retained_pct": float(100.0 * (~flagged).sum() / one_sample_adata.n_obs)
                if one_sample_adata.n_obs
                else np.nan,
            }
        )
    return pd.DataFrame(records)


def filter_qc_pass_samples(
    adata_names: Sequence[str],
    adata_list: Sequence[ad.AnnData],
    pass_column: str = "qc_mad_pass",
) -> tuple[list[str], list[ad.AnnData], pd.DataFrame]:
    """Return retained-cell sample copies plus a per-sample filtering summary."""
    filtered_names: list[str] = []
    filtered_adatas: list[ad.AnnData] = []
    records = []

    for run_sample_id, one_sample_adata in zip(adata_names, adata_list):
        if pass_column not in one_sample_adata.obs.columns:
            raise KeyError(f"Missing {pass_column!r} for {run_sample_id}. Run annotate_mad_qc() first.")

        keep_mask = one_sample_adata.obs[pass_column].astype(bool).to_numpy()
        filtered_sample = one_sample_adata[keep_mask].copy()
        filtered_names.append(str(run_sample_id))
        filtered_adatas.append(filtered_sample)
        records.append(
            {
                "run_sample_id": str(run_sample_id),
                "starting_n_cells": int(one_sample_adata.n_obs),
                "retained_n_cells": int(keep_mask.sum()),
                "removed_n_cells": int((~keep_mask).sum()),
                "retained_pct": float(100.0 * keep_mask.sum() / one_sample_adata.n_obs)
                if one_sample_adata.n_obs
                else np.nan,
            }
        )

    return filtered_names, filtered_adatas, pd.DataFrame(records)


def concat_samples(
    adata_names: Sequence[str],
    adata_list: Sequence[ad.AnnData],
    label: str = "batch",
    join: str = "outer",
) -> ad.AnnData:
    """Concatenate per-sample AnnData objects with stable sample metadata."""
    combined_adata = ad.concat(
        {str(run_sample_id): one_sample_adata for run_sample_id, one_sample_adata in zip(adata_names, adata_list)},
        axis=0,
        join=join,
        merge="same",
        fill_value=0,
        index_unique=None,
        label=label,
    )

    if label in combined_adata.obs.columns:
        combined_adata.obs[label] = combined_adata.obs[label].astype("string")
        combined_adata.obs["run_sample_id"] = combined_adata.obs[label].astype("string")

    for column in ["DIV", "run_sample_id", "cell_line", "data_source"]:
        if column in combined_adata.obs.columns:
            combined_adata.obs[column] = combined_adata.obs[column].astype("string")

    if not combined_adata.obs_names.is_unique:
        duplicated = combined_adata.obs_names[combined_adata.obs_names.duplicated()].unique()[:10].tolist()
        raise ValueError(f"Combined AnnData has duplicated obs_names. Examples: {duplicated}")

    return combined_adata


def preprocess_basic(adata: ad.AnnData, settings: PreprocessSettings = PreprocessSettings()) -> dict:
    """Normalize, log-transform, select HVGs, and run PCA in place."""
    started_at = time.perf_counter()
    sc.pp.normalize_total(adata, target_sum=settings.target_sum)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(
        adata,
        n_top_genes=settings.n_top_genes,
        batch_key=settings.hvg_batch_key,
        flavor=settings.hvg_flavor,
    )
    sc.tl.pca(
        adata,
        mask_var="highly_variable",
        svd_solver=settings.pca_solver,
    )
    return {
        "target_sum": settings.target_sum,
        "n_top_genes": settings.n_top_genes,
        "hvg_batch_key": settings.hvg_batch_key,
        "elapsed_seconds": time.perf_counter() - started_at,
    }


def run_neighbors_umap(adata: ad.AnnData, neighbors_key: str | None = None) -> dict:
    """Run neighbors and UMAP in place."""
    started_at = time.perf_counter()
    if neighbors_key:
        sc.pp.neighbors(adata, key_added=neighbors_key)
        sc.tl.umap(adata, neighbors_key=neighbors_key)
    else:
        sc.pp.neighbors(adata)
        sc.tl.umap(adata)
    return {"neighbors_key": neighbors_key or "", "elapsed_seconds": time.perf_counter() - started_at}
