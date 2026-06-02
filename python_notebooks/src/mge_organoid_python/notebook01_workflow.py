"""Reusable Notebook 01 workflow scaffolding.

Notebook 01 starts from frozen Notebook 00 checkpoints. This module keeps path
resolution, input validation, and analysis-branch bookkeeping out of notebook
cells so combined and per-sample runs can be compared reproducibly.
"""

from __future__ import annotations

import time
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable, Sequence

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc


NOTEBOOK01_SUPPORTED_SCOPES = ("combined", "per_sample")
NOTEBOOK01_REGRESSION_BRANCHES = ("not_regressed", "regressed_ccdifference")
NOTEBOOK01_DEFAULT_REGRESS_KEYS = ("CCDifference",)
NOTEBOOK01_DEFAULT_QC_REGRESS_KEYS = ("total_counts", "pct_counts_mt")
NOTEBOOK01_BASE_OPERATION_ORDER = (
    "optional_regress_out_obs_covariates_from_full_x",
    "select_hvgs_from_counts_layer_after_regression_step",
    "subset_to_hvgs",
    "scale_x",
    "pca_from_scaled_x",
    "neighbors_from_pca",
    "umap_from_neighbors",
    "leiden_from_neighbors",
)


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
class Notebook01HVGSettings:
    """Seurat-v3 HVG settings for Notebook 01."""

    n_top_genes: int = 4000
    flavor: str = "seurat_v3"
    layer: str = "counts"
    batch_key: str | None = None


@dataclass(frozen=True)
class Notebook01EmbeddingSettings:
    """PCA/neighbors/UMAP/clustering settings for Notebook 01 branches."""

    n_pcs: int = 50
    n_neighbors: int = 15
    leiden_resolution: float = 0.5
    random_state: int = 0
    scale_max_value: float = 10.0


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
        regress_keys: Sequence[str] = NOTEBOOK01_DEFAULT_QC_REGRESS_KEYS,
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
        regress_keys: Sequence[str] = (*NOTEBOOK01_DEFAULT_QC_REGRESS_KEYS, "CCDifference"),
    ) -> "Notebook01RegressionVariant":
        """Return the planned later branch after cell-cycle scoring is added."""
        return cls(
            branch="regressed_qc_cc",
            regress_keys=tuple(regress_keys),
            regress_cell_cycle=True,
            description="QC and cell-cycle covariates are regressed from .X before scaling/PCA.",
        )

    @classmethod
    def ccdifference_regressed(
        cls,
        regress_keys: Sequence[str] = ("CCDifference",),
    ) -> "Notebook01RegressionVariant":
        """Return the CCDifference-regressed comparison branch."""
        return cls(
            branch="regressed_ccdifference",
            regress_keys=tuple(regress_keys),
            regress_cell_cycle=True,
            description="CCDifference is regressed from .X before scaling/PCA.",
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
    """Return the default non-regressed and CCDifference-regressed branches."""
    return (
        Notebook01RegressionVariant.not_regressed(),
        Notebook01RegressionVariant.ccdifference_regressed(regress_keys=regress_keys),
    )


def branch_operation_order(variant: Notebook01RegressionVariant) -> tuple[str, ...]:
    """Return the concrete operation order for one regression branch."""
    regression_step = (
        "regress_out_obs_covariates_from_full_x"
        if variant.regress_keys
        else "skip_regress_out_keep_normalized_log1p_full_x"
    )
    return (
        regression_step,
        "select_hvgs_from_counts_layer_after_regression_step",
        "subset_to_hvgs",
        "scale_x",
        "pca_from_scaled_x",
        "neighbors_from_pca",
        "umap_from_neighbors",
        "leiden_from_neighbors",
    )


def parse_csv(value: str | Sequence[str] | None, default: Sequence[str] = ()) -> tuple[str, ...]:
    """Parse comma/semicolon/colon-delimited env-style strings."""
    if value is None:
        return tuple(default)
    if isinstance(value, str):
        if not value.strip():
            return tuple(default)
        normalized = value.replace(";", ",").replace(":", ",")
        return tuple(part.strip() for part in normalized.split(",") if part.strip())
    return tuple(str(part).strip() for part in value if str(part).strip())


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


def run_hvg_selection(
    adata: ad.AnnData,
    *,
    settings: Notebook01HVGSettings = Notebook01HVGSettings(),
    scope: str = "combined",
    run_sample_id: object | None = None,
    branch: str | None = None,
) -> tuple[np.ndarray, pd.DataFrame, pd.DataFrame]:
    """Select HVGs from the configured layer and return mask, gene table, and parameters."""
    started_at = time.perf_counter()
    if settings.layer not in adata.layers:
        raise KeyError(f"Missing counts layer {settings.layer!r} for HVG selection.")

    sc.pp.highly_variable_genes(
        adata,
        n_top_genes=settings.n_top_genes,
        flavor=settings.flavor,
        layer=settings.layer,
        batch_key=settings.batch_key,
        subset=False,
    )

    hvg_mask = adata.var["highly_variable"].fillna(False).to_numpy(dtype=bool)
    hvg_table = adata.var.loc[hvg_mask].copy()
    hvg_table.insert(0, "gene", hvg_table.index.astype(str))
    if branch is not None:
        hvg_table.insert(0, "branch", branch)
    hvg_table.insert(0, "run_sample_id", "" if run_sample_id is None else str(run_sample_id))
    hvg_table.insert(0, "scope", scope)
    hvg_table = hvg_table.reset_index(drop=True)

    parameter_df = settings_to_frame(settings, method="notebook01_hvg")
    if branch is not None:
        parameter_df.insert(0, "branch", branch)
    parameter_df.insert(0, "run_sample_id", "" if run_sample_id is None else str(run_sample_id))
    parameter_df.insert(0, "scope", scope)
    parameter_df["n_cells"] = int(adata.n_obs)
    parameter_df["n_genes"] = int(adata.n_vars)
    parameter_df["n_highly_variable_genes"] = int(hvg_mask.sum())
    parameter_df["elapsed_seconds"] = time.perf_counter() - started_at
    return hvg_mask, hvg_table, parameter_df


def _branch_log_prefix(scope: str, branch: str, run_sample_id: object | None = None) -> str:
    """Return a compact stdout prefix for branch-level execution logs."""
    sample = "" if run_sample_id is None else f" sample={run_sample_id}"
    return f"[Notebook01 matrix flow] scope={scope}{sample} branch={branch}"


def _print_branch_matrix_flow(
    *,
    scope: str,
    branch: str,
    run_sample_id: object | None = None,
    message: str,
) -> None:
    """Print branch matrix provenance messages into notebook and Slurm logs."""
    print(f"{_branch_log_prefix(scope, branch, run_sample_id)} {message}", flush=True)


def run_regression_embedding_branch(
    adata: ad.AnnData,
    *,
    variant: Notebook01RegressionVariant,
    hvg_settings: Notebook01HVGSettings = Notebook01HVGSettings(),
    settings: Notebook01EmbeddingSettings = Notebook01EmbeddingSettings(),
    scope: str = "combined",
    run_sample_id: object | None = None,
) -> tuple[ad.AnnData, dict, pd.DataFrame, pd.DataFrame]:
    """Run one regression/HVG/PCA/UMAP/clustering branch.

    Input `adata.X` is expected to be normalized/log1p expression from the
    Notebook 00 normalized checkpoint. The branch first works on the full gene
    matrix. If `variant.regress_keys` is nonempty, `sc.pp.regress_out` replaces
    full `.X` with residuals after modeling each gene against the named `.obs`
    covariates. HVG selection is then run on the branch object using the
    configured layer, followed by HVG subsetting, scaling, PCA, neighbors, UMAP,
    and Leiden.
    """
    started_at = time.perf_counter()
    operation_order = branch_operation_order(variant)
    regression_covariate_source = (
        ",".join(f".obs[{key}]" for key in variant.regress_keys) if variant.regress_keys else ""
    )
    branch_adata = adata.copy()
    _print_branch_matrix_flow(
        scope=scope,
        branch=variant.branch,
        run_sample_id=run_sample_id,
        message=(
            f"branch start: .X=normalized_log1p_full_matrix shape={branch_adata.X.shape}; "
            f".layers[{hvg_settings.layer!r}] shape={branch_adata.layers[hvg_settings.layer].shape}"
        ),
    )
    missing_regress_keys = [key for key in variant.regress_keys if key not in branch_adata.obs.columns]
    if missing_regress_keys:
        raise KeyError(f"Missing regression covariates in adata.obs: {missing_regress_keys}")

    if variant.regress_keys:
        _print_branch_matrix_flow(
            scope=scope,
            branch=variant.branch,
            run_sample_id=run_sample_id,
            message=f"regress_out input: full branch_adata.X; keys={list(variant.regress_keys)}; n_jobs=Scanpy default",
        )
        sc.pp.regress_out(branch_adata, keys=list(variant.regress_keys))
        _print_branch_matrix_flow(
            scope=scope,
            branch=variant.branch,
            run_sample_id=run_sample_id,
            message=f"regress_out saved: branch_adata.X=regressed_residuals_full_matrix shape={branch_adata.X.shape}",
        )
    else:
        _print_branch_matrix_flow(
            scope=scope,
            branch=variant.branch,
            run_sample_id=run_sample_id,
            message="regress_out skipped: branch_adata.X remains normalized_log1p_full_matrix",
        )

    _print_branch_matrix_flow(
        scope=scope,
        branch=variant.branch,
        run_sample_id=run_sample_id,
        message=(
            "HVG input: "
            f"branch_adata.layers[{hvg_settings.layer!r}] raw counts; "
            f"flavor={hvg_settings.flavor}; n_top_genes={hvg_settings.n_top_genes}; "
            f"batch_key={hvg_settings.batch_key}; subset=False"
        ),
    )
    hvg_mask, hvg_table, hvg_params = run_hvg_selection(
        branch_adata,
        settings=hvg_settings,
        scope=scope,
        run_sample_id=run_sample_id,
        branch=variant.branch,
    )
    _print_branch_matrix_flow(
        scope=scope,
        branch=variant.branch,
        run_sample_id=run_sample_id,
        message=f"HVG saved: branch_adata.var['highly_variable']; n_highly_variable_genes={int(hvg_mask.sum())}",
    )
    branch_adata = branch_adata[:, hvg_mask].copy()
    _print_branch_matrix_flow(
        scope=scope,
        branch=variant.branch,
        run_sample_id=run_sample_id,
        message=(
            f"HVG subset saved: branch_adata.X shape={branch_adata.X.shape}; "
            f"branch_adata.layers[{hvg_settings.layer!r}] shape={branch_adata.layers[hvg_settings.layer].shape}"
        ),
    )

    _print_branch_matrix_flow(
        scope=scope,
        branch=variant.branch,
        run_sample_id=run_sample_id,
        message=f"scale input: branch_adata.X shape={branch_adata.X.shape}; max_value={settings.scale_max_value}",
    )
    sc.pp.scale(branch_adata, max_value=settings.scale_max_value)
    _print_branch_matrix_flow(
        scope=scope,
        branch=variant.branch,
        run_sample_id=run_sample_id,
        message=f"scale saved: branch_adata.X=scaled matrix shape={branch_adata.X.shape}",
    )
    n_comps = min(settings.n_pcs, max(1, branch_adata.n_vars - 1), max(1, branch_adata.n_obs - 1))
    _print_branch_matrix_flow(
        scope=scope,
        branch=variant.branch,
        run_sample_id=run_sample_id,
        message=f"PCA input: scaled branch_adata.X shape={branch_adata.X.shape}; n_comps={n_comps}",
    )
    sc.pp.pca(branch_adata, n_comps=n_comps, svd_solver="arpack", random_state=settings.random_state)
    _print_branch_matrix_flow(
        scope=scope,
        branch=variant.branch,
        run_sample_id=run_sample_id,
        message=(
            "PCA saved: branch_adata.obsm['X_pca'] "
            f"shape={branch_adata.obsm['X_pca'].shape}; branch_adata.varm['PCs'] "
            f"shape={branch_adata.varm['PCs'].shape}; branch_adata.uns['pca']"
        ),
    )
    _print_branch_matrix_flow(
        scope=scope,
        branch=variant.branch,
        run_sample_id=run_sample_id,
        message=f"Neighbors input: branch_adata.obsm['X_pca'] shape={branch_adata.obsm['X_pca'].shape}",
    )
    sc.pp.neighbors(
        branch_adata,
        n_neighbors=settings.n_neighbors,
        n_pcs=n_comps,
        random_state=settings.random_state,
    )
    _print_branch_matrix_flow(
        scope=scope,
        branch=variant.branch,
        run_sample_id=run_sample_id,
        message=(
            "Neighbors saved: branch_adata.uns['neighbors']; "
            "branch_adata.obsp['distances']; branch_adata.obsp['connectivities']"
        ),
    )
    _print_branch_matrix_flow(
        scope=scope,
        branch=variant.branch,
        run_sample_id=run_sample_id,
        message="UMAP input: neighbors graph from branch_adata.uns['neighbors']/branch_adata.obsp",
    )
    sc.tl.umap(branch_adata, random_state=settings.random_state)
    _print_branch_matrix_flow(
        scope=scope,
        branch=variant.branch,
        run_sample_id=run_sample_id,
        message=f"UMAP saved: branch_adata.obsm['X_umap'] shape={branch_adata.obsm['X_umap'].shape}",
    )
    _print_branch_matrix_flow(
        scope=scope,
        branch=variant.branch,
        run_sample_id=run_sample_id,
        message="Leiden input: neighbors graph from branch_adata.uns['neighbors']/branch_adata.obsp",
    )
    sc.tl.leiden(
        branch_adata,
        resolution=settings.leiden_resolution,
        key_added="leiden",
        random_state=settings.random_state,
    )
    _print_branch_matrix_flow(
        scope=scope,
        branch=variant.branch,
        run_sample_id=run_sample_id,
        message=f"Leiden saved: branch_adata.obs['leiden']; n_clusters={int(branch_adata.obs['leiden'].nunique())}",
    )

    matrix_flow_validation = {
        "pca_input": "scaled branch_adata.X",
        "pca_saved": "branch_adata.obsm['X_pca']; branch_adata.varm['PCs']; branch_adata.uns['pca']",
        "neighbors_input": "branch_adata.obsm['X_pca']",
        "neighbors_saved": "branch_adata.uns['neighbors']; branch_adata.obsp['distances']; branch_adata.obsp['connectivities']",
        "umap_input": "neighbors graph",
        "umap_saved": "branch_adata.obsm['X_umap']",
        "leiden_input": "neighbors graph",
        "leiden_saved": "branch_adata.obs['leiden']",
    }

    branch_adata.uns["notebook01_branch"] = {
        "scope": scope,
        "run_sample_id": "" if run_sample_id is None else str(run_sample_id),
        "branch": variant.branch,
        "regress_keys": list(variant.regress_keys),
        "regression_covariate_source": regression_covariate_source,
        "regress_cell_cycle": variant.regress_cell_cycle,
        "operation_order": list(operation_order),
        "x_state_at_branch_start": "normalized_log1p_full_matrix",
        "x_state_after_regression_step": (
            "regressed_residuals_full_matrix" if variant.regress_keys else "normalized_log1p_full_matrix"
        ),
        "hvg_selection_layer": hvg_settings.layer,
        "hvg_selection_after_regression_step": True,
        "x_state_after_hvg_subset": (
            "regressed_residuals_hvg_subset" if variant.regress_keys else "normalized_log1p_hvg_subset"
        ),
        "x_state_after_scale": "scaled",
        "hvg_n_genes": int(hvg_mask.sum()),
        "hvg_settings": asdict(hvg_settings),
        "embedding_settings": asdict(settings),
        "matrix_flow_validation": matrix_flow_validation,
        **asdict(settings),
    }

    record = {
        "scope": scope,
        "run_sample_id": "" if run_sample_id is None else str(run_sample_id),
        "branch": variant.branch,
        "regress_keys": ",".join(variant.regress_keys),
        "regression_covariate_source": regression_covariate_source,
        "operation_order": " -> ".join(operation_order),
        "x_state_at_branch_start": "normalized_log1p_full_matrix",
        "x_state_after_regression_step": (
            "regressed_residuals_full_matrix" if variant.regress_keys else "normalized_log1p_full_matrix"
        ),
        "hvg_selection_layer": hvg_settings.layer,
        "hvg_selection_after_regression_step": True,
        "x_state_after_hvg_subset": (
            "regressed_residuals_hvg_subset" if variant.regress_keys else "normalized_log1p_hvg_subset"
        ),
        "x_state_after_scale": "scaled",
        "n_cells": int(branch_adata.n_obs),
        "n_hvg_genes": int(branch_adata.n_vars),
        "n_pcs": int(n_comps),
        "n_neighbors": settings.n_neighbors,
        "leiden_resolution": settings.leiden_resolution,
        "n_leiden_clusters": int(branch_adata.obs["leiden"].nunique()),
        **matrix_flow_validation,
        "elapsed_seconds": time.perf_counter() - started_at,
    }
    return branch_adata, record, hvg_table, hvg_params


def embedding_table(
    adata: ad.AnnData,
    *,
    scope: str,
    branch: str,
    run_sample_id: object | None = None,
) -> pd.DataFrame:
    """Return UMAP coordinates and selected metadata for one branch."""
    if "X_umap" not in adata.obsm:
        raise KeyError("Missing adata.obsm['X_umap']; run UMAP first.")
    coords = pd.DataFrame(adata.obsm["X_umap"], index=adata.obs_names, columns=["UMAP1", "UMAP2"])
    coords.insert(0, "cell_id", coords.index.astype(str))
    coords.insert(0, "branch", branch)
    coords.insert(0, "run_sample_id", "" if run_sample_id is None else str(run_sample_id))
    coords.insert(0, "scope", scope)
    for column in ["data_source", "run_sample_id", "cell_line", "total_counts", "pct_counts_mt", "n_genes_by_counts", "leiden"]:
        if column in adata.obs.columns and column not in coords.columns:
            coords[column] = adata.obs[column].to_numpy()
    return coords.reset_index(drop=True)


def save_umap_plot(
    adata: ad.AnnData,
    *,
    color: str,
    output_path: Path | str,
    title: str,
    show: bool = False,
    dpi: int = 160,
) -> Path:
    """Save one UMAP plot for a Notebook 01 branch."""
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    sc.pl.umap(adata, color=color, title=title, show=False)
    fig = plt.gcf()
    fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
    if show:
        plt.show()
    plt.close(fig)
    return output_path


def save_pca_variance_plot(
    adata: ad.AnnData,
    *,
    output_path: Path | str,
    title: str,
    show: bool = False,
    dpi: int = 160,
) -> Path:
    """Save one PCA variance-ratio plot for a Notebook 01 branch."""
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    sc.pl.pca_variance_ratio(adata, n_pcs=min(50, adata.uns["pca"]["variance_ratio"].shape[0]), show=False)
    fig = plt.gcf()
    fig.suptitle(title)
    fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
    if show:
        plt.show()
    plt.close(fig)
    return output_path
