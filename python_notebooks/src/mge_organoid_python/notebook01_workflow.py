"""Reusable Notebook 01 workflow scaffolding.

Notebook 01 starts from frozen Notebook 00 checkpoints. This module keeps path
resolution, input validation, and analysis-branch bookkeeping out of notebook
cells so combined and per-sample runs can be compared reproducibly.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable, Sequence

import anndata as ad
import pandas as pd


NOTEBOOK01_SUPPORTED_SCOPES = ("combined", "per_sample")
NOTEBOOK01_REGRESSION_BRANCHES = ("not_regressed", "regressed_qc")
NOTEBOOK01_DEFAULT_REGRESS_KEYS = ("total_counts", "pct_counts_mt")


@dataclass(frozen=True)
class Notebook01InputSettings:
    """Notebook 00 checkpoint selection for Notebook 01."""

    notebook00_run_label: str = "cellranger_filtered_manual_ec_div30_core_samples_freeze"
    counts_layer: str = "counts"
    normalized_log1p_filename: str = "manual_ec_filtered_normalized_log1p.h5ad"
    counts_filename: str = "manual_ec_filtered_counts.h5ad"
    notebook00_results_dirname: str = "notebook00"
    per_sample_dirname: str = "per_sample"


@dataclass(frozen=True)
class Notebook01RunSettings:
    """Notebook 01 run-scoped output settings."""

    run_label: str
    notebook01_results_dirname: str = "notebook01"
    scopes: tuple[str, ...] = NOTEBOOK01_SUPPORTED_SCOPES
    regression_branches: tuple[str, ...] = NOTEBOOK01_REGRESSION_BRANCHES


@dataclass(frozen=True)
class Notebook01RegressionVariant:
    """One branch of a Notebook 01 regression comparison."""

    branch: str
    regress_keys: tuple[str, ...] = ()
    regress_cell_cycle: bool = False
    description: str = ""

    @classmethod
    def not_regressed(cls) -> "Notebook01RegressionVariant":
        """Return the non-regressed comparison branch."""
        return cls(
            branch="not_regressed",
            regress_keys=(),
            regress_cell_cycle=False,
            description="No covariates are regressed from .X before scaling/PCA.",
        )

    @classmethod
    def qc_regressed(
        cls,
        regress_keys: Sequence[str] = NOTEBOOK01_DEFAULT_REGRESS_KEYS,
    ) -> "Notebook01RegressionVariant":
        """Return the initial QC-covariate regressed comparison branch."""
        return cls(
            branch="regressed_qc",
            regress_keys=tuple(regress_keys),
            regress_cell_cycle=False,
            description="QC covariates are regressed from .X before scaling/PCA.",
        )

    @classmethod
    def qc_and_cell_cycle_regressed(
        cls,
        regress_keys: Sequence[str] = (*NOTEBOOK01_DEFAULT_REGRESS_KEYS, "CC.Difference"),
    ) -> "Notebook01RegressionVariant":
        """Return the planned later branch after cell-cycle scoring is added."""
        return cls(
            branch="regressed_qc_cc",
            regress_keys=tuple(regress_keys),
            regress_cell_cycle=True,
            description="QC and cell-cycle covariates are regressed from .X before scaling/PCA.",
        )


@dataclass(frozen=True)
class Notebook01InputPaths:
    """Resolved Notebook 00 checkpoint inputs for Notebook 01."""

    data_root: Path
    notebook00_run_label: str
    notebook00_run_dir: Path
    h5ad_dir: Path
    per_sample_h5ad_dir: Path
    combined_counts: Path
    combined_normalized_log1p: Path

    @classmethod
    def from_data_root(
        cls,
        data_root: Path | str,
        settings: Notebook01InputSettings = Notebook01InputSettings(),
    ) -> "Notebook01InputPaths":
        """Resolve Notebook 00 checkpoint paths from the large data root."""
        resolved_data_root = Path(data_root).expanduser().resolve()
        run_dir = (
            resolved_data_root
            / "results"
            / settings.notebook00_results_dirname
            / settings.notebook00_run_label
        )
        h5ad_dir = run_dir / "h5ad"
        per_sample_h5ad_dir = h5ad_dir / settings.per_sample_dirname
        return cls(
            data_root=resolved_data_root,
            notebook00_run_label=settings.notebook00_run_label,
            notebook00_run_dir=run_dir,
            h5ad_dir=h5ad_dir,
            per_sample_h5ad_dir=per_sample_h5ad_dir,
            combined_counts=h5ad_dir / settings.counts_filename,
            combined_normalized_log1p=h5ad_dir / settings.normalized_log1p_filename,
        )

    def per_sample_normalized_log1p(self, run_sample_id: object) -> Path:
        """Return the per-sample normalized/log1p checkpoint path."""
        sample_stem = _safe_filename_stem(run_sample_id)
        return self.per_sample_h5ad_dir / f"{sample_stem}.manual_ec_filtered_normalized_log1p.h5ad"

    def per_sample_counts(self, run_sample_id: object) -> Path:
        """Return the per-sample counts checkpoint path."""
        sample_stem = _safe_filename_stem(run_sample_id)
        return self.per_sample_h5ad_dir / f"{sample_stem}.manual_ec_filtered_counts.h5ad"


@dataclass(frozen=True)
class Notebook01OutputPaths:
    """Run-scoped Notebook 01 output paths."""

    data_root: Path
    run_label: str
    run_dir: Path
    table_dir: Path
    plot_dir: Path
    h5ad_dir: Path
    executed_dir: Path

    @classmethod
    def from_data_root(
        cls,
        data_root: Path | str,
        settings: Notebook01RunSettings,
    ) -> "Notebook01OutputPaths":
        """Resolve Notebook 01 output paths from the large data root."""
        resolved_data_root = Path(data_root).expanduser().resolve()
        run_dir = resolved_data_root / "results" / settings.notebook01_results_dirname / settings.run_label
        return cls(
            data_root=resolved_data_root,
            run_label=settings.run_label,
            run_dir=run_dir,
            table_dir=run_dir / "tables",
            plot_dir=run_dir / "plots",
            h5ad_dir=run_dir / "h5ad",
            executed_dir=resolved_data_root / "results" / settings.notebook01_results_dirname / "executed",
        )

    def ensure_dirs(self) -> None:
        """Create Notebook 01 output directories."""
        for path in [self.table_dir, self.plot_dir, self.h5ad_dir, self.executed_dir]:
            path.mkdir(parents=True, exist_ok=True)

    def branch_h5ad_dir(self, scope: str, branch: str, run_sample_id: object | None = None) -> Path:
        """Return the h5ad directory for one scope/branch combination."""
        if scope == "combined":
            return self.h5ad_dir / scope / branch
        if scope == "per_sample":
            if run_sample_id is None:
                raise ValueError("run_sample_id is required for per_sample branch output.")
            return self.h5ad_dir / scope / _safe_filename_stem(run_sample_id) / branch
        raise ValueError(f"Unsupported scope: {scope!r}")


def _safe_filename_stem(value: object) -> str:
    """Return a filesystem-safe filename stem."""
    safe = str(value).strip().replace(" ", "_").replace("/", "_")
    return "".join(ch for ch in safe if ch.isalnum() or ch in {"_", "-", "."})


def settings_to_frame(settings: object, method: str) -> pd.DataFrame:
    """Return a stable settings table."""
    return pd.DataFrame(
        [{"method": method, "parameter": key, "value": value} for key, value in asdict(settings).items()]
    )


def default_regression_variants(
    regress_keys: Sequence[str] = NOTEBOOK01_DEFAULT_REGRESS_KEYS,
) -> tuple[Notebook01RegressionVariant, ...]:
    """Return the default non-regressed and QC-regressed branches."""
    return (
        Notebook01RegressionVariant.not_regressed(),
        Notebook01RegressionVariant.qc_regressed(regress_keys=regress_keys),
    )


def validate_notebook01_input(
    adata: ad.AnnData,
    *,
    input_path: Path | str,
    counts_layer: str = "counts",
    scope: str = "combined",
    run_sample_id: object | None = None,
) -> dict:
    """Validate one Notebook 00 normalized/log1p checkpoint for Notebook 01 use."""
    input_path = Path(input_path)
    counts_layer_present = counts_layer in adata.layers
    counts_layer_shape_matches_x = bool(counts_layer_present and adata.layers[counts_layer].shape == adata.X.shape)
    checkpoint_metadata = adata.uns.get("notebook00_checkpoint", {})

    record = {
        "scope": scope,
        "run_sample_id": "" if run_sample_id is None else str(run_sample_id),
        "input_path": str(input_path),
        "input_exists": input_path.exists(),
        "n_cells": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "x_state": checkpoint_metadata.get("x_state", ""),
        "counts_layer": counts_layer,
        "counts_layer_present": bool(counts_layer_present),
        "counts_layer_shape_matches_x": counts_layer_shape_matches_x,
        "obs_names_unique": bool(adata.obs_names.is_unique),
        "var_names_unique": bool(adata.var_names.is_unique),
    }
    checks = [
        record["input_exists"],
        record["x_state"] in {"normalized_log1p", ""},
        record["counts_layer_present"],
        record["counts_layer_shape_matches_x"],
        record["obs_names_unique"],
        record["var_names_unique"],
    ]
    record["notebook01_input_validation_passed"] = all(checks)
    return record


def infer_run_sample_ids(adata: ad.AnnData, sample_column: str = "run_sample_id") -> list[str]:
    """Return run_sample_id values in first-seen order."""
    if sample_column not in adata.obs.columns:
        raise KeyError(f"{sample_column!r} not found in adata.obs")
    return [str(value) for value in adata.obs[sample_column].astype(str).drop_duplicates()]


def planned_analysis_table(
    *,
    run_sample_ids: Iterable[str],
    variants: Sequence[Notebook01RegressionVariant] | None = None,
    scopes: Sequence[str] = NOTEBOOK01_SUPPORTED_SCOPES,
) -> pd.DataFrame:
    """Return the planned combined/per-sample branch table for Notebook 01."""
    selected_variants = tuple(variants or default_regression_variants())
    records = []
    for scope in scopes:
        if scope not in NOTEBOOK01_SUPPORTED_SCOPES:
            raise ValueError(f"Unsupported Notebook 01 scope: {scope!r}")
        sample_ids = [""] if scope == "combined" else list(run_sample_ids)
        for run_sample_id in sample_ids:
            for variant in selected_variants:
                records.append(
                    {
                        "scope": scope,
                        "run_sample_id": run_sample_id,
                        "branch": variant.branch,
                        "regress_cell_cycle": variant.regress_cell_cycle,
                        "regress_keys": ",".join(variant.regress_keys),
                        "description": variant.description,
                    }
                )
    return pd.DataFrame(records)
