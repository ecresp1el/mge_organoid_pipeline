#!/usr/bin/env python
"""Run Notebook 01 CCDifference before/after PCA diagnostics."""

from __future__ import annotations

import gc
import os
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from pandas.api.types import is_numeric_dtype


def _add_repo_src_to_path() -> Path:
    start = Path(__file__).resolve()
    for candidate in [start.parent, *start.parents]:
        src_dir = candidate / "python_notebooks" / "src"
        if src_dir.exists():
            sys.path.insert(0, str(src_dir))
            return candidate
    raise FileNotFoundError("Could not find python_notebooks/src from script path.")


REPO_ROOT = _add_repo_src_to_path()

from mge_organoid_python.cell_cycle import (  # noqa: E402
    CCDIFFERENCE_KEY,
    cell_cycle_score_summary,
    score_cell_cycle_and_ccdifference,
    select_cell_cycle_genes,
)
from mge_organoid_python.data_sources import resolve_data_root  # noqa: E402
from mge_organoid_python.notebook01_workflow import (  # noqa: E402
    Notebook01EmbeddingSettings,
    Notebook01InputPaths,
    Notebook01InputSettings,
    Notebook01OutputPaths,
    Notebook01RunSettings,
    infer_run_sample_ids,
    parse_csv,
    validate_notebook01_input,
)


PCA_PAIRS = ((1, 2), (1, 3), (2, 3))


def env_int(name: str, default: int) -> int:
    raw = os.environ.get(name)
    return default if raw is None or not raw.strip() else int(raw)


def safe_part(value: object) -> str:
    safe = str(value).strip().replace(" ", "_").replace("/", "_")
    return "".join(ch for ch in safe if ch.isalnum() or ch in {"_", "-", "."})


def pca_colors(scope: str) -> tuple[str, ...]:
    base = ("phase", "S_score", "G2M_score", CCDIFFERENCE_KEY)
    if scope == "combined":
        return (*base, "run_sample_id")
    return base


def prepare_pca_object(
    adata: sc.AnnData,
    *,
    cell_cycle_genes: tuple[str, ...],
    state: str,
    regress_key: str | None,
    n_pcs: int,
    random_state: int,
) -> sc.AnnData:
    """Create a cell-cycle-gene PCA object before or after CCDifference regression."""
    pca_adata = adata[:, list(cell_cycle_genes)].copy()
    pca_adata.uns["notebook01_ccdifference_pca_state"] = {
        "state": state,
        "regress_key": "" if regress_key is None else regress_key,
        "operation_order": [
            "subset_to_present_cell_cycle_genes",
            "optional_regress_out_CCDifference_from_x",
            "scale_x",
            "pca_from_scaled_x",
        ],
        "x_state_at_start": "normalized_log1p_cell_cycle_gene_subset",
        "x_state_after_regression": (
            "ccdifference_regressed_residuals" if regress_key else "normalized_log1p_cell_cycle_gene_subset"
        ),
        "x_state_after_scale": "scaled",
    }
    if regress_key is not None:
        sc.pp.regress_out(pca_adata, keys=[regress_key])
    sc.pp.scale(pca_adata, max_value=10.0)
    sc.tl.pca(
        pca_adata,
        n_comps=min(n_pcs, pca_adata.n_vars - 1, pca_adata.n_obs - 1),
        svd_solver="arpack",
        random_state=random_state,
    )
    return pca_adata


def save_pca_pair_plot(
    pca_adata: sc.AnnData,
    *,
    pc_x: int,
    pc_y: int,
    color: str,
    output_path: Path,
    title: str,
) -> None:
    """Save one PCA pair scatter plot."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    coords = pca_adata.obsm["X_pca"]
    x = coords[:, pc_x - 1]
    y = coords[:, pc_y - 1]
    values = pca_adata.obs[color]

    fig, ax = plt.subplots(figsize=(6.2, 5.2))
    point_size = 2.0 if pca_adata.n_obs > 50000 else 4.0
    if is_numeric_dtype(values):
        scatter = ax.scatter(x, y, c=values.to_numpy(), s=point_size, cmap="viridis", linewidths=0, alpha=0.75)
        fig.colorbar(scatter, ax=ax, label=color)
    else:
        categories = pd.Categorical(values.astype(str))
        palette = plt.get_cmap("tab20")
        codes = categories.codes
        for code, category in enumerate(categories.categories):
            mask = codes == code
            ax.scatter(
                x[mask],
                y[mask],
                s=point_size,
                color=palette(code % palette.N),
                linewidths=0,
                alpha=0.75,
                label=str(category),
            )
        ax.legend(markerscale=4, frameon=False, fontsize=8, loc="best")

    ax.set_xlabel(f"PC{pc_x}")
    ax.set_ylabel(f"PC{pc_y}")
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(output_path, dpi=160)
    plt.close(fig)


def pca_summary_record(
    pca_adata: sc.AnnData,
    *,
    scope: str,
    run_sample_id: object | None,
    state: str,
    regress_key: str | None,
) -> dict:
    """Return one PCA diagnostic summary row."""
    variance = pca_adata.uns["pca"]["variance_ratio"]
    record = {
        "scope": scope,
        "run_sample_id": "" if run_sample_id is None else str(run_sample_id),
        "state": state,
        "regress_key": "" if regress_key is None else regress_key,
        "n_cells": int(pca_adata.n_obs),
        "n_cell_cycle_genes": int(pca_adata.n_vars),
        "n_pcs": int(pca_adata.obsm["X_pca"].shape[1]),
    }
    for i in range(min(3, len(variance))):
        record[f"PC{i + 1}_variance_ratio"] = float(variance[i])
    return record


def run_scope(
    adata: sc.AnnData,
    *,
    input_path: Path,
    output_paths: Notebook01OutputPaths,
    scope: str,
    run_sample_id: object | None,
    embedding_settings: Notebook01EmbeddingSettings,
    score_records: list[pd.DataFrame],
    pca_records: list[dict],
    plot_records: list[dict],
    input_validation_records: list[dict],
) -> None:
    """Run scoring and before/after CCDifference PCA for one scope/sample."""
    validation = validate_notebook01_input(
        adata,
        input_path=input_path,
        counts_layer="counts",
        scope=scope,
        run_sample_id=run_sample_id,
    )
    input_validation_records.append(validation)
    if not validation["notebook01_input_validation_passed"]:
        raise AssertionError(f"Notebook 01 input validation failed: {validation}")

    selection = select_cell_cycle_genes(adata.var_names)
    score_cell_cycle_and_ccdifference(adata, selection=selection)
    score_records.append(cell_cycle_score_summary(adata, scope=scope, run_sample_id=run_sample_id))

    for state, regress_key in [
        ("before_ccdifference_regression", None),
        ("after_ccdifference_regression", CCDIFFERENCE_KEY),
    ]:
        pca_adata = prepare_pca_object(
            adata,
            cell_cycle_genes=selection.cell_cycle_genes,
            state=state,
            regress_key=regress_key,
            n_pcs=3,
            random_state=embedding_settings.random_state,
        )
        pca_records.append(
            pca_summary_record(
                pca_adata,
                scope=scope,
                run_sample_id=run_sample_id,
                state=state,
                regress_key=regress_key,
            )
        )

        plot_dir = output_paths.plot_dir / "cell_cycle_pca" / scope
        if scope == "per_sample":
            plot_dir = plot_dir / safe_part(run_sample_id)
        plot_dir = plot_dir / state

        for pc_x, pc_y in PCA_PAIRS:
            for color in pca_colors(scope):
                if color not in pca_adata.obs.columns:
                    continue
                plot_path = plot_dir / f"pca_pc{pc_x}_pc{pc_y}_{safe_part(color)}.png"
                save_pca_pair_plot(
                    pca_adata,
                    pc_x=pc_x,
                    pc_y=pc_y,
                    color=color,
                    output_path=plot_path,
                    title=f"{scope} {run_sample_id or ''} {state}: PC{pc_x} vs PC{pc_y} by {color}".strip(),
                )
                plot_records.append(
                    {
                        "scope": scope,
                        "run_sample_id": "" if run_sample_id is None else str(run_sample_id),
                        "state": state,
                        "pc_x": pc_x,
                        "pc_y": pc_y,
                        "color": color,
                        "path": str(plot_path),
                    }
                )
        del pca_adata
        gc.collect()


def main() -> None:
    data_root = resolve_data_root()
    notebook00_run_label = os.environ.get(
        "NOTEBOOK01_NOTEBOOK00_RUN_LABEL",
        "cellranger_filtered_manual_ec_div30_core_samples_freeze",
    )
    run_label = os.environ.get(
        "NOTEBOOK01_RUN_LABEL",
        f"{notebook00_run_label}_ccdifference_seurat_order_v1",
    )
    scopes = parse_csv(os.environ.get("NOTEBOOK01_SCOPES"), default=("combined", "per_sample"))
    embedding_settings = Notebook01EmbeddingSettings(
        n_pcs=3,
        random_state=env_int("NOTEBOOK01_RANDOM_STATE", 0),
    )

    input_settings = Notebook01InputSettings(notebook00_run_label=notebook00_run_label)
    run_settings = Notebook01RunSettings(run_label=run_label, scopes=scopes)
    input_paths = Notebook01InputPaths.from_data_root(data_root, settings=input_settings)
    output_paths = Notebook01OutputPaths.from_data_root(data_root, settings=run_settings)
    output_paths.ensure_dirs()

    print(f"Repo root: {REPO_ROOT}")
    print(f"Data root: {data_root}")
    print(f"Notebook 00 run label: {notebook00_run_label}")
    print(f"Notebook 01 run label: {run_label}")
    print(f"Scopes: {scopes}")
    print(f"Input checkpoint: {input_paths.combined_normalized_log1p}")
    print(f"Output root: {output_paths.run_dir}")

    input_validation_records: list[dict] = []
    score_records: list[pd.DataFrame] = []
    pca_records: list[dict] = []
    plot_records: list[dict] = []

    combined_adata = sc.read_h5ad(input_paths.combined_normalized_log1p)
    run_sample_ids = infer_run_sample_ids(combined_adata)

    if "combined" in scopes:
        run_scope(
            combined_adata,
            input_path=input_paths.combined_normalized_log1p,
            output_paths=output_paths,
            scope="combined",
            run_sample_id=None,
            embedding_settings=embedding_settings,
            score_records=score_records,
            pca_records=pca_records,
            plot_records=plot_records,
            input_validation_records=input_validation_records,
        )

    del combined_adata
    gc.collect()

    if "per_sample" in scopes:
        for sample_id in run_sample_ids:
            sample_path = input_paths.per_sample_normalized_log1p(sample_id)
            sample_adata = sc.read_h5ad(sample_path)
            run_scope(
                sample_adata,
                input_path=sample_path,
                output_paths=output_paths,
                scope="per_sample",
                run_sample_id=sample_id,
                embedding_settings=embedding_settings,
                score_records=score_records,
                pca_records=pca_records,
                plot_records=plot_records,
                input_validation_records=input_validation_records,
            )
            del sample_adata
            gc.collect()

    input_validation_df = pd.DataFrame(input_validation_records)
    score_summary_df = pd.concat(score_records, ignore_index=True) if score_records else pd.DataFrame()
    pca_summary_df = pd.DataFrame(pca_records)
    plot_manifest_df = pd.DataFrame(plot_records)

    input_validation_df.to_csv(output_paths.table_dir / "cell_cycle_input_validation.tsv", sep="\t", index=False)
    score_summary_df.to_csv(output_paths.table_dir / "cell_cycle_score_summary.tsv", sep="\t", index=False)
    pca_summary_df.to_csv(output_paths.table_dir / "cell_cycle_pca_diagnostic_summary.tsv", sep="\t", index=False)
    plot_manifest_df.to_csv(output_paths.table_dir / "cell_cycle_pca_plot_manifest.tsv", sep="\t", index=False)

    if not input_validation_df["notebook01_input_validation_passed"].all():
        raise AssertionError("One or more Notebook 01 input validations failed.")

    print(f"Input validation rows: {len(input_validation_df)}")
    print(f"Score summary rows: {len(score_summary_df)}")
    print(f"PCA summary rows: {len(pca_summary_df)}")
    print(f"Plot rows: {len(plot_manifest_df)}")
    print(pca_summary_df.to_string(index=False))


if __name__ == "__main__":
    main()
