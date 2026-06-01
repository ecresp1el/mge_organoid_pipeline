"""Reusable manual_ec workflow pieces for Notebook 00."""

from __future__ import annotations

import time
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Sequence

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc


GROUP_COLUMNS = ("data_source", "run_sample_id", "cell_line")


@dataclass(frozen=True)
class ManualECFilterSettings:
    """Manual EC QC cutoffs requested for the Scanpy conversion path."""

    min_genes: int = 20
    min_cells: int = 3
    min_genes_by_counts: int = 1500
    max_genes_by_counts: int = 9000
    max_total_counts: int = 30000
    max_pct_counts_mt: float = 5.0
    mito_prefix: str = "MT-"


@dataclass(frozen=True)
class ManualECPreprocessSettings:
    """Manual EC normalization and HVG settings."""

    target_sum: float = 1e4
    hvg_min_mean: float = 0.0125
    hvg_max_mean: float = 3.0
    hvg_min_disp: float = 0.5


def settings_to_frame(settings: object, method: str) -> pd.DataFrame:
    """Return a stable two-column settings table."""
    records = [{"method": method, "parameter": key, "value": value} for key, value in asdict(settings).items()]
    return pd.DataFrame(records)


def calculate_qc_metrics(
    adata: ad.AnnData,
    mito_prefix: str = "MT-",
    log1p: bool = False,
    overwrite: bool = True,
) -> dict:
    """Add Scanpy QC metric columns to one AnnData object."""
    required_obs_columns = {"total_counts", "n_genes_by_counts", "pct_counts_mt"}
    missing_qc_columns = required_obs_columns.difference(adata.obs.columns)

    if overwrite or missing_qc_columns:
        adata.var["mt"] = adata.var_names.str.startswith(mito_prefix)
        sc.pp.calculate_qc_metrics(
            adata,
            qc_vars=["mt"],
            percent_top=None,
            log1p=log1p,
            inplace=True,
        )

    return {
        "mito_prefix": mito_prefix,
        "log1p": log1p,
        "n_cells": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
    }


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

    for column in ["DIV", "run_sample_id", "biological_label", "cell_line", "data_source"]:
        if column in combined_adata.obs.columns:
            combined_adata.obs[column] = combined_adata.obs[column].astype("string")

    if not combined_adata.obs_names.is_unique:
        duplicated = combined_adata.obs_names[combined_adata.obs_names.duplicated()].unique()[:10].tolist()
        raise ValueError(f"Combined AnnData has duplicated obs_names. Examples: {duplicated}")

    return combined_adata


def _available_group_columns(adata: ad.AnnData) -> list[str]:
    return [column for column in GROUP_COLUMNS if column in adata.obs.columns]


def _cell_counts_by_group(adata: ad.AnnData, count_column: str) -> pd.DataFrame:
    group_columns = _available_group_columns(adata)
    if not group_columns:
        return pd.DataFrame({count_column: [int(adata.n_obs)]})

    obs_df = adata.obs[group_columns].copy()
    for column in group_columns:
        obs_df[column] = obs_df[column].astype("string")

    return (
        obs_df.groupby(group_columns, observed=True)
        .size()
        .rename(count_column)
        .reset_index()
    )


def _merge_group_counts(frames: Sequence[pd.DataFrame]) -> pd.DataFrame:
    if not frames:
        return pd.DataFrame()

    group_columns = [column for column in GROUP_COLUMNS if any(column in frame.columns for frame in frames)]
    merged = frames[0].copy()
    for frame in frames[1:]:
        if group_columns:
            merged = merged.merge(frame, on=group_columns, how="outer")
        else:
            merged = pd.concat([merged, frame], axis=1)

    count_columns = [column for column in merged.columns if column.endswith("_n_cells")]
    for column in count_columns:
        merged[column] = merged[column].fillna(0).astype(int)
    return merged


def manual_ec_qc_summary(adata: ad.AnnData, stage: str) -> pd.DataFrame:
    """Summarize manual_ec QC metrics by source/sample/cell line."""
    required_metrics = ["n_genes_by_counts", "total_counts", "pct_counts_mt"]
    missing_metrics = [metric for metric in required_metrics if metric not in adata.obs.columns]
    if missing_metrics:
        raise KeyError(f"Missing QC metrics: {missing_metrics}. Run calculate_qc_metrics() first.")

    group_columns = _available_group_columns(adata)
    obs_df = adata.obs[[*group_columns, *required_metrics]].copy()
    for column in group_columns:
        obs_df[column] = obs_df[column].astype("string")

    if group_columns:
        summary_df = (
            obs_df.groupby(group_columns, observed=True)
            .agg(
                n_cells=("total_counts", "size"),
                median_n_genes_by_counts=("n_genes_by_counts", "median"),
                median_total_counts=("total_counts", "median"),
                median_pct_counts_mt=("pct_counts_mt", "median"),
                mean_n_genes_by_counts=("n_genes_by_counts", "mean"),
                mean_total_counts=("total_counts", "mean"),
                mean_pct_counts_mt=("pct_counts_mt", "mean"),
            )
            .reset_index()
        )
    else:
        summary_df = pd.DataFrame(
            {
                "n_cells": [int(adata.n_obs)],
                "median_n_genes_by_counts": [float(obs_df["n_genes_by_counts"].median())],
                "median_total_counts": [float(obs_df["total_counts"].median())],
                "median_pct_counts_mt": [float(obs_df["pct_counts_mt"].median())],
                "mean_n_genes_by_counts": [float(obs_df["n_genes_by_counts"].mean())],
                "mean_total_counts": [float(obs_df["total_counts"].mean())],
                "mean_pct_counts_mt": [float(obs_df["pct_counts_mt"].mean())],
            }
        )

    summary_df.insert(0, "stage", stage)
    return summary_df


def apply_manual_ec_filter(
    adata: ad.AnnData,
    settings: ManualECFilterSettings = ManualECFilterSettings(),
) -> tuple[ad.AnnData, pd.DataFrame, pd.DataFrame]:
    """Apply preliminary Scanpy filters and requested manual_ec QC cutoffs."""
    started_at = time.perf_counter()
    working = adata.copy()
    calculate_qc_metrics(working, mito_prefix=settings.mito_prefix, log1p=False, overwrite=True)

    starting_counts = _cell_counts_by_group(working, "starting_n_cells")
    starting_n_genes = int(working.n_vars)

    sc.pp.filter_cells(working, min_genes=settings.min_genes)
    sc.pp.filter_genes(working, min_cells=settings.min_cells)
    calculate_qc_metrics(working, mito_prefix=settings.mito_prefix, log1p=False, overwrite=True)

    post_prelim_counts = _cell_counts_by_group(working, "post_min_gene_cell_filter_n_cells")
    post_prelim_n_genes = int(working.n_vars)

    keep_mask = (
        (working.obs["n_genes_by_counts"] > settings.min_genes_by_counts)
        & (working.obs["n_genes_by_counts"] < settings.max_genes_by_counts)
        & (working.obs["total_counts"] < settings.max_total_counts)
        & (working.obs["pct_counts_mt"] < settings.max_pct_counts_mt)
    )
    filtered = working[keep_mask.to_numpy()].copy()
    filtered.obs["manual_ec_pass"] = True
    filtered.uns["manual_ec_filter_settings"] = asdict(settings)
    filtered.uns["manual_ec_filter_note"] = "manual_ec cutoffs applied after sc.pp.filter_cells/min_genes and sc.pp.filter_genes/min_cells."

    retained_counts = _cell_counts_by_group(filtered, "manual_ec_retained_n_cells")
    summary_df = _merge_group_counts([starting_counts, post_prelim_counts, retained_counts])
    summary_df["manual_ec_removed_n_cells"] = summary_df["starting_n_cells"] - summary_df["manual_ec_retained_n_cells"]
    summary_df["manual_ec_retained_pct"] = np.where(
        summary_df["starting_n_cells"] > 0,
        100.0 * summary_df["manual_ec_retained_n_cells"] / summary_df["starting_n_cells"],
        np.nan,
    )
    summary_df["starting_n_genes"] = starting_n_genes
    summary_df["post_min_cell_filter_n_genes"] = post_prelim_n_genes
    summary_df["manual_ec_final_n_genes"] = int(filtered.n_vars)
    summary_df["elapsed_seconds"] = time.perf_counter() - started_at

    parameters_df = settings_to_frame(settings, method="manual_ec_filter")
    return filtered, summary_df, parameters_df


def manual_ec_preprocess(
    adata: ad.AnnData,
    settings: ManualECPreprocessSettings = ManualECPreprocessSettings(),
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Normalize, log-transform, and select HVGs with manual_ec settings."""
    started_at = time.perf_counter()
    sc.pp.normalize_total(adata, target_sum=settings.target_sum)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(
        adata,
        min_mean=settings.hvg_min_mean,
        max_mean=settings.hvg_max_mean,
        min_disp=settings.hvg_min_disp,
    )
    adata.uns["manual_ec_preprocess_settings"] = asdict(settings)

    hvg_mask = adata.var["highly_variable"].fillna(False).astype(bool)
    hvg_genes_df = adata.var.loc[hvg_mask].copy()
    hvg_genes_df.insert(0, "gene", hvg_genes_df.index.astype(str))
    hvg_genes_df = hvg_genes_df.reset_index(drop=True)

    report_df = pd.DataFrame(
        [
            {
                "n_cells": int(adata.n_obs),
                "n_genes": int(adata.n_vars),
                "n_highly_variable_genes": int(hvg_mask.sum()),
                "target_sum": settings.target_sum,
                "hvg_min_mean": settings.hvg_min_mean,
                "hvg_max_mean": settings.hvg_max_mean,
                "hvg_min_disp": settings.hvg_min_disp,
                "elapsed_seconds": time.perf_counter() - started_at,
            }
        ]
    )
    parameters_df = settings_to_frame(settings, method="manual_ec_preprocess")
    return report_df, parameters_df, hvg_genes_df


def load_manual_ec_run_outputs(results_root: Path | str, run_labels: Sequence[str]) -> dict[str, pd.DataFrame]:
    """Load manual_ec output tables for completed run labels."""
    records = []
    qc_records = []
    hvg_by_run: dict[str, set[str]] = {}
    results_path = Path(results_root).expanduser().resolve()

    for run_label in run_labels:
        table_dir = results_path / run_label / "tables"
        filter_summary_path = table_dir / "manual_ec_filter_summary.tsv"
        qc_summary_path = table_dir / "manual_ec_qc_summary_after_filter.tsv"
        hvg_path = table_dir / "manual_ec_hvg_genes.tsv"

        if not filter_summary_path.exists():
            raise FileNotFoundError(f"Missing manual_ec filter summary for {run_label}: {filter_summary_path}")
        if not qc_summary_path.exists():
            raise FileNotFoundError(f"Missing manual_ec QC summary for {run_label}: {qc_summary_path}")
        if not hvg_path.exists():
            raise FileNotFoundError(f"Missing manual_ec HVG table for {run_label}: {hvg_path}")

        filter_df = pd.read_csv(filter_summary_path, sep="\t")
        filter_df.insert(0, "run_label", run_label)
        records.append(filter_df)

        qc_df = pd.read_csv(qc_summary_path, sep="\t")
        qc_df.insert(0, "run_label", run_label)
        qc_records.append(qc_df)

        hvg_df = pd.read_csv(hvg_path, sep="\t")
        hvg_by_run[run_label] = set(hvg_df["gene"].astype(str)) if "gene" in hvg_df.columns else set()

    filter_summary = pd.concat(records, ignore_index=True) if records else pd.DataFrame()
    qc_summary = pd.concat(qc_records, ignore_index=True) if qc_records else pd.DataFrame()

    overlap_records = []
    labels = list(run_labels)
    for i, left in enumerate(labels):
        for right in labels[i + 1 :]:
            left_genes = hvg_by_run[left]
            right_genes = hvg_by_run[right]
            intersection = left_genes & right_genes
            union = left_genes | right_genes
            overlap_records.append(
                {
                    "left_run_label": left,
                    "right_run_label": right,
                    "left_n_hvg": len(left_genes),
                    "right_n_hvg": len(right_genes),
                    "intersection_n_hvg": len(intersection),
                    "union_n_hvg": len(union),
                    "jaccard": (len(intersection) / len(union)) if union else np.nan,
                }
            )

    return {
        "manual_ec_source_comparison_summary": filter_summary,
        "manual_ec_qc_metric_comparison_by_sample": qc_summary,
        "manual_ec_hvg_overlap": pd.DataFrame(overlap_records),
    }
