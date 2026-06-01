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
import scipy.sparse as sp


GROUP_COLUMNS = ("data_source", "run_sample_id", "cell_line")


def _safe_filename_stem(value: object) -> str:
    """Return a filesystem-safe stem while preserving sample IDs like 9853-MW-1."""
    safe = str(value).strip().replace(" ", "_").replace("/", "_")
    return "".join(ch for ch in safe if ch.isalnum() or ch in {"_", "-", "."})


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
class ManualECCheckpointSettings:
    """Notebook 00 checkpoint output settings."""

    target_sum: float = 1e4
    counts_layer: str = "counts"
    h5ad_dirname: str = "h5ad"
    per_sample_dirname: str = "per_sample"
    combined_counts_filename: str = "manual_ec_filtered_counts.h5ad"
    combined_normalized_log1p_filename: str = "manual_ec_filtered_normalized_log1p.h5ad"


@dataclass(frozen=True)
class ManualECCheckpointPaths:
    """Run-scoped Notebook 00 checkpoint paths."""

    run_dir: Path
    h5ad_dir: Path
    per_sample_h5ad_dir: Path
    combined_counts: Path
    combined_normalized_log1p: Path

    @classmethod
    def from_run_dir(
        cls,
        run_dir: Path | str,
        settings: ManualECCheckpointSettings = ManualECCheckpointSettings(),
    ) -> "ManualECCheckpointPaths":
        """Build all checkpoint paths from the dynamic Notebook 00 run directory."""
        resolved_run_dir = Path(run_dir).expanduser().resolve()
        h5ad_dir = resolved_run_dir / settings.h5ad_dirname
        return cls(
            run_dir=resolved_run_dir,
            h5ad_dir=h5ad_dir,
            per_sample_h5ad_dir=h5ad_dir / settings.per_sample_dirname,
            combined_counts=h5ad_dir / settings.combined_counts_filename,
            combined_normalized_log1p=h5ad_dir / settings.combined_normalized_log1p_filename,
        )

    def per_sample_counts(self, run_sample_id: object) -> Path:
        """Return the counts checkpoint path for one sample."""
        sample_stem = _safe_filename_stem(run_sample_id)
        return self.per_sample_h5ad_dir / f"{sample_stem}.manual_ec_filtered_counts.h5ad"

    def per_sample_normalized_log1p(self, run_sample_id: object) -> Path:
        """Return the normalized/log1p checkpoint path for one sample."""
        sample_stem = _safe_filename_stem(run_sample_id)
        return self.per_sample_h5ad_dir / f"{sample_stem}.manual_ec_filtered_normalized_log1p.h5ad"


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


def _copy_matrix(matrix):
    """Return an independent matrix copy for AnnData X/layers values."""
    return matrix.copy() if hasattr(matrix, "copy") else np.array(matrix, copy=True)


def _matrix_equal(left, right) -> bool:
    """Return True when two dense/sparse matrices have identical values."""
    if left.shape != right.shape:
        return False
    if sp.issparse(left) or sp.issparse(right):
        left_sparse = left if sp.issparse(left) else sp.csr_matrix(left)
        right_sparse = right if sp.issparse(right) else sp.csr_matrix(right)
        difference = left_sparse != right_sparse
        return int(difference.nnz) == 0
    return bool(np.array_equal(np.asarray(left), np.asarray(right)))


def _matrix_nonnegative(matrix) -> bool:
    """Return True when a dense/sparse matrix has no negative values."""
    if sp.issparse(matrix):
        data = matrix.data
        return True if data.size == 0 else bool(np.nanmin(data) >= 0)
    values = np.asarray(matrix)
    return True if values.size == 0 else bool(np.nanmin(values) >= 0)


def _annotate_checkpoint(
    adata: ad.AnnData,
    *,
    checkpoint_name: str,
    x_state: str,
    scope: str,
    settings: ManualECCheckpointSettings,
    run_sample_id: object | None = None,
) -> None:
    """Store lightweight checkpoint metadata in AnnData.uns."""
    metadata = {
        "checkpoint_name": checkpoint_name,
        "x_state": x_state,
        "scope": scope,
        "target_sum": settings.target_sum,
        "counts_layer": settings.counts_layer,
    }
    if run_sample_id is not None:
        metadata["run_sample_id"] = str(run_sample_id)
    adata.uns["notebook00_checkpoint"] = metadata


def _prepare_for_h5ad_write(adata: ad.AnnData) -> None:
    """Normalize string columns to AnnData's stable H5AD representation."""
    adata.strings_to_categoricals()


def validate_manual_ec_checkpoint_pair(
    counts_adata: ad.AnnData,
    normalized_log1p_adata: ad.AnnData,
    *,
    counts_path: Path | str,
    normalized_log1p_path: Path | str,
    scope: str,
    settings: ManualECCheckpointSettings = ManualECCheckpointSettings(),
    run_sample_id: object | None = None,
) -> dict:
    """Validate the Notebook 00 counts and normalized/log1p checkpoint contract."""
    counts_path = Path(counts_path)
    normalized_log1p_path = Path(normalized_log1p_path)

    obs_names_match = counts_adata.obs_names.equals(normalized_log1p_adata.obs_names)
    var_names_match = counts_adata.var_names.equals(normalized_log1p_adata.var_names)
    counts_layer_present = settings.counts_layer in normalized_log1p_adata.layers
    counts_layer_matches_counts_x = (
        _matrix_equal(counts_adata.X, normalized_log1p_adata.layers[settings.counts_layer])
        if counts_layer_present
        else False
    )

    counts_metadata = counts_adata.uns.get("notebook00_checkpoint", {})
    normalized_metadata = normalized_log1p_adata.uns.get("notebook00_checkpoint", {})
    counts_x_state = counts_metadata.get("x_state", "")
    normalized_x_state = normalized_metadata.get("x_state", "")

    checks = {
        "obs_names_match": bool(obs_names_match),
        "var_names_match": bool(var_names_match),
        "counts_x_nonnegative": _matrix_nonnegative(counts_adata.X),
        "normalized_x_nonnegative": _matrix_nonnegative(normalized_log1p_adata.X),
        "counts_layer_present": bool(counts_layer_present),
        "counts_layer_matches_counts_x": bool(counts_layer_matches_counts_x),
        "counts_x_state_is_filtered_counts": counts_x_state == "filtered_counts",
        "normalized_x_state_is_normalized_log1p": normalized_x_state == "normalized_log1p",
        "counts_file_exists": counts_path.exists(),
        "normalized_log1p_file_exists": normalized_log1p_path.exists(),
    }

    record = {
        "scope": scope,
        "run_sample_id": "" if run_sample_id is None else str(run_sample_id),
        "counts_path": str(counts_path),
        "normalized_log1p_path": str(normalized_log1p_path),
        "n_cells": int(counts_adata.n_obs),
        "n_genes": int(counts_adata.n_vars),
        "counts_x_state": counts_x_state,
        "normalized_x_state": normalized_x_state,
        "counts_layer": settings.counts_layer,
        "counts_file_size_bytes": counts_path.stat().st_size if counts_path.exists() else 0,
        "normalized_log1p_file_size_bytes": (
            normalized_log1p_path.stat().st_size if normalized_log1p_path.exists() else 0
        ),
        **checks,
    }
    record["checkpoint_validation_passed"] = all(checks.values())
    return record


def _first_seen_sample_ids(adata: ad.AnnData, sample_column: str) -> list[str]:
    """Return sample IDs in first-seen obs order."""
    if sample_column not in adata.obs.columns:
        raise KeyError(f"{sample_column!r} not found in adata.obs")
    return [str(value) for value in adata.obs[sample_column].astype(str).drop_duplicates()]


def _checkpoint_summary_record(
    *,
    scope: str,
    run_sample_id: object | None,
    counts_path: Path,
    normalized_log1p_path: Path,
    counts_adata: ad.AnnData,
    settings: ManualECCheckpointSettings,
) -> dict:
    """Return one stable checkpoint summary row."""
    return {
        "scope": scope,
        "run_sample_id": "" if run_sample_id is None else str(run_sample_id),
        "counts_path": str(counts_path),
        "normalized_log1p_path": str(normalized_log1p_path),
        "n_cells": int(counts_adata.n_obs),
        "n_genes": int(counts_adata.n_vars),
        "counts_x": "filtered_counts",
        "normalized_log1p_x": "normalized_log1p",
        "normalized_log1p_counts_layer": settings.counts_layer,
        "target_sum": settings.target_sum,
    }


def save_manual_ec_checkpoints(
    filtered_counts_adata: ad.AnnData,
    run_dir: Path | str,
    *,
    sample_column: str = "run_sample_id",
    settings: ManualECCheckpointSettings = ManualECCheckpointSettings(),
) -> tuple[ad.AnnData, ad.AnnData, pd.DataFrame, pd.DataFrame, ManualECCheckpointPaths]:
    """Write combined and per-sample Notebook 00 counts/normalized checkpoints."""
    started_at = time.perf_counter()
    paths = ManualECCheckpointPaths.from_run_dir(run_dir, settings=settings)
    paths.h5ad_dir.mkdir(parents=True, exist_ok=True)
    paths.per_sample_h5ad_dir.mkdir(parents=True, exist_ok=True)

    counts_adata = filtered_counts_adata.copy()
    _annotate_checkpoint(
        counts_adata,
        checkpoint_name="manual_ec_filtered_counts",
        x_state="filtered_counts",
        scope="combined",
        settings=settings,
    )

    normalized_log1p_adata = counts_adata.copy()
    normalized_log1p_adata.layers[settings.counts_layer] = _copy_matrix(counts_adata.X)
    sc.pp.normalize_total(normalized_log1p_adata, target_sum=settings.target_sum)
    sc.pp.log1p(normalized_log1p_adata)
    _annotate_checkpoint(
        normalized_log1p_adata,
        checkpoint_name="manual_ec_filtered_normalized_log1p",
        x_state="normalized_log1p",
        scope="combined",
        settings=settings,
    )

    _prepare_for_h5ad_write(counts_adata)
    _prepare_for_h5ad_write(normalized_log1p_adata)
    counts_adata.write_h5ad(paths.combined_counts)
    normalized_log1p_adata.write_h5ad(paths.combined_normalized_log1p)

    summary_records = [
        _checkpoint_summary_record(
            scope="combined",
            run_sample_id=None,
            counts_path=paths.combined_counts,
            normalized_log1p_path=paths.combined_normalized_log1p,
            counts_adata=counts_adata,
            settings=settings,
        )
    ]
    validation_records = [
        validate_manual_ec_checkpoint_pair(
            counts_adata,
            normalized_log1p_adata,
            counts_path=paths.combined_counts,
            normalized_log1p_path=paths.combined_normalized_log1p,
            scope="combined",
            settings=settings,
        )
    ]

    sample_ids = _first_seen_sample_ids(counts_adata, sample_column=sample_column)
    for sample_id in sample_ids:
        sample_mask = counts_adata.obs[sample_column].astype(str).to_numpy() == sample_id
        sample_counts = counts_adata[sample_mask].copy()
        sample_normalized_log1p = normalized_log1p_adata[sample_mask].copy()
        _annotate_checkpoint(
            sample_counts,
            checkpoint_name="manual_ec_filtered_counts",
            x_state="filtered_counts",
            scope="per_sample",
            settings=settings,
            run_sample_id=sample_id,
        )
        _annotate_checkpoint(
            sample_normalized_log1p,
            checkpoint_name="manual_ec_filtered_normalized_log1p",
            x_state="normalized_log1p",
            scope="per_sample",
            settings=settings,
            run_sample_id=sample_id,
        )

        sample_counts_path = paths.per_sample_counts(sample_id)
        sample_normalized_log1p_path = paths.per_sample_normalized_log1p(sample_id)
        _prepare_for_h5ad_write(sample_counts)
        _prepare_for_h5ad_write(sample_normalized_log1p)
        sample_counts.write_h5ad(sample_counts_path)
        sample_normalized_log1p.write_h5ad(sample_normalized_log1p_path)

        summary_records.append(
            _checkpoint_summary_record(
                scope="per_sample",
                run_sample_id=sample_id,
                counts_path=sample_counts_path,
                normalized_log1p_path=sample_normalized_log1p_path,
                counts_adata=sample_counts,
                settings=settings,
            )
        )
        validation_records.append(
            validate_manual_ec_checkpoint_pair(
                sample_counts,
                sample_normalized_log1p,
                counts_path=sample_counts_path,
                normalized_log1p_path=sample_normalized_log1p_path,
                scope="per_sample",
                settings=settings,
                run_sample_id=sample_id,
            )
        )

    summary_df = pd.DataFrame(summary_records)
    summary_df["elapsed_seconds"] = time.perf_counter() - started_at
    validation_df = pd.DataFrame(validation_records)
    return counts_adata, normalized_log1p_adata, summary_df, validation_df, paths


def load_manual_ec_run_outputs(results_root: Path | str, run_labels: Sequence[str]) -> dict[str, pd.DataFrame]:
    """Load manual_ec output tables for completed run labels."""
    records = []
    qc_records = []
    results_path = Path(results_root).expanduser().resolve()

    for run_label in run_labels:
        table_dir = results_path / run_label / "tables"
        filter_summary_path = table_dir / "manual_ec_filter_summary.tsv"
        qc_summary_path = table_dir / "manual_ec_qc_summary_after_filter.tsv"

        if not filter_summary_path.exists():
            raise FileNotFoundError(f"Missing manual_ec filter summary for {run_label}: {filter_summary_path}")
        if not qc_summary_path.exists():
            raise FileNotFoundError(f"Missing manual_ec QC summary for {run_label}: {qc_summary_path}")

        filter_df = pd.read_csv(filter_summary_path, sep="\t")
        filter_df.insert(0, "run_label", run_label)
        records.append(filter_df)

        qc_df = pd.read_csv(qc_summary_path, sep="\t")
        qc_df.insert(0, "run_label", run_label)
        qc_records.append(qc_df)

    filter_summary = pd.concat(records, ignore_index=True) if records else pd.DataFrame()
    qc_summary = pd.concat(qc_records, ignore_index=True) if qc_records else pd.DataFrame()

    return {
        "manual_ec_source_comparison_summary": filter_summary,
        "manual_ec_qc_metric_comparison_by_sample": qc_summary,
    }
