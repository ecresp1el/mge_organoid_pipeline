#!/usr/bin/env python3
"""Run Seurat-based Shi label transfer smoke plots for DIV30."""

from __future__ import annotations

import os
import subprocess
from pathlib import Path

import anndata as ad
import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from mge_organoid_python.shi_label_transfer import (
    attach_reference_labels,
    broad_region_from_label,
    developmental_class_from_label,
    get_umap_coordinates,
    load_shi_table_s2_labels,
    natural_sort_key,
    plot_overlaid_density_by_group,
    plot_stacked_bar,
    plot_umap_categorical,
    plot_umap_continuous,
    safe_token,
    summarize_predictions_by_cluster,
)


MARKER = "[ShiSeuratTransfer]"


def log(message: str) -> None:
    print(f"{MARKER} {message}", flush=True)


def env_path(name: str, default: str) -> Path:
    return Path(os.environ.get(name, default)).expanduser()


def env_int(name: str, default: int) -> int:
    return int(os.environ.get(name, str(default)))


def env_bool(name: str, default: bool) -> bool:
    raw = os.environ.get(name)
    if raw is None:
        return default
    return str(raw).strip().lower() in {"1", "true", "yes", "y"}


def sanitize_score_label(label: str) -> str:
    return safe_token(label).replace("__", "_")


def prepare_reference_labels(
    reference_h5ad: Path,
    table_s2_xlsx: Path,
    labels_path: Path,
    table_dir: Path,
    min_match_fraction: float,
) -> list[str]:
    log(f"Reading Shi Table S2 labels: {table_s2_xlsx}")
    table_s2 = load_shi_table_s2_labels(table_s2_xlsx)

    log(f"Reading Shi reference AnnData metadata: {reference_h5ad}")
    reference = ad.read_h5ad(reference_h5ad, backed="r")
    try:
        label_obs, join_summary = attach_reference_labels(
            reference,
            table_s2,
            min_match_fraction=min_match_fraction,
        )
        for source_col, output_col in [
            ("cell_id", "reference_cell_id"),
            ("raw_cell_id", "reference_raw_cell_id"),
            ("barcode", "reference_barcode"),
        ]:
            if source_col in reference.obs.columns:
                label_obs[output_col] = reference.obs[source_col].astype(str).to_numpy()
            else:
                label_obs[output_col] = ""
    finally:
        reference.file.close()

    labels_path.parent.mkdir(parents=True, exist_ok=True)
    table_dir.mkdir(parents=True, exist_ok=True)
    label_obs.reset_index(drop=True).to_csv(labels_path, sep="\t", index=False)
    join_summary.to_csv(table_dir / "shi_reference_table_s2_label_join_summary.tsv", sep="\t", index=False)
    counts = (
        label_obs["shi_label"]
        .dropna()
        .astype(str)
        .value_counts()
        .rename_axis("shi_label")
        .reset_index(name="n_reference_cells")
    )
    counts.to_csv(table_dir / "shi_reference_table_s2_label_counts.tsv", sep="\t", index=False)
    log(f"Wrote Seurat label bridge: {labels_path}")
    return counts["shi_label"].astype(str).tolist()


def run_seurat_transfer(
    repo_root: Path,
    rscript_bin: Path,
    reference_seurat: Path,
    query_seurat: Path,
    labels_path: Path,
    seurat_dir: Path,
    dims: int,
    nfeatures: int,
    npcs: int,
    seed: int,
) -> None:
    script = repo_root / "python_notebooks" / "scripts" / "seurat_shi_label_transfer_export.R"
    cmd = [
        str(rscript_bin),
        str(script),
        "--reference",
        str(reference_seurat),
        "--query",
        str(query_seurat),
        "--labels",
        str(labels_path),
        "--outdir",
        str(seurat_dir),
        "--reference_assay",
        os.environ.get("SHI_SEURAT_REFERENCE_ASSAY", "RNA"),
        "--query_assay",
        os.environ.get("SHI_SEURAT_QUERY_ASSAY", "RNA"),
        "--transfer_name",
        "shi_label",
        "--normalization_method",
        os.environ.get("SHI_SEURAT_NORMALIZATION_METHOD", "LogNormalize"),
        "--dims",
        str(dims),
        "--nfeatures",
        str(nfeatures),
        "--npcs",
        str(npcs),
        "--seed",
        str(seed),
    ]
    log("Running Seurat transfer/export")
    log("Command: " + " ".join(cmd))
    subprocess.run(cmd, check=True)


def load_predictions(prediction_path: Path) -> pd.DataFrame:
    predictions = pd.read_csv(prediction_path, sep="\t")
    if "cell_id" not in predictions.columns:
        raise ValueError(f"Missing cell_id in Seurat predictions: {prediction_path}")
    if "predicted.id" not in predictions.columns:
        raise ValueError(f"Missing predicted.id in Seurat predictions: {prediction_path}")
    if "prediction.score.max" not in predictions.columns:
        raise ValueError(f"Missing prediction.score.max in Seurat predictions: {prediction_path}")
    predictions["cell_id"] = predictions["cell_id"].astype(str)
    return predictions


def align_predictions_to_adata(adata: ad.AnnData, predictions: pd.DataFrame) -> pd.DataFrame:
    pred_by_cell = predictions.set_index("cell_id", drop=False)
    obs_names = pd.Index(adata.obs_names.astype(str))
    direct_matches = obs_names.isin(pred_by_cell.index).sum()
    if direct_matches == adata.n_obs:
        log("Joining Seurat predictions by AnnData obs_names")
        return pred_by_cell.loc[obs_names].reset_index(drop=True)

    if "cell_id" in adata.obs.columns:
        cell_ids = pd.Index(adata.obs["cell_id"].astype(str))
        cell_matches = cell_ids.isin(pred_by_cell.index).sum()
        if cell_matches == adata.n_obs:
            log('Joining Seurat predictions by adata.obs["cell_id"]')
            return pred_by_cell.loc[cell_ids].reset_index(drop=True)

    diagnostics = pd.DataFrame(
        [
            {"join_key": "obs_names", "matched_cells": int(direct_matches), "query_cells": int(adata.n_obs)},
            {
                "join_key": "obs.cell_id",
                "matched_cells": int(adata.obs["cell_id"].astype(str).isin(pred_by_cell.index).sum())
                if "cell_id" in adata.obs.columns
                else 0,
                "query_cells": int(adata.n_obs),
            },
        ]
    )
    raise ValueError(
        "Could not fully join Seurat predictions to AnnData cells:\n"
        + diagnostics.to_string(index=False)
    )


def add_predictions_to_obs(adata: ad.AnnData, aligned: pd.DataFrame, label_order: list[str]) -> list[str]:
    prefix = "shi_seurat_full"
    adata.obs[f"{prefix}_predicted_shi_label"] = aligned["predicted.id"].astype(str).to_numpy()
    adata.obs[f"{prefix}_prediction_score"] = pd.to_numeric(
        aligned["prediction.score.max"],
        errors="coerce",
    ).to_numpy()
    adata.obs[f"{prefix}_uncertainty_score"] = 1.0 - adata.obs[f"{prefix}_prediction_score"].astype(float)
    adata.obs[f"{prefix}_broad_region_class"] = adata.obs[f"{prefix}_predicted_shi_label"].map(broad_region_from_label)
    adata.obs[f"{prefix}_developmental_class"] = adata.obs[f"{prefix}_predicted_shi_label"].map(
        developmental_class_from_label
    )

    score_cols = [col for col in aligned.columns if col.startswith("prediction.score.") and col != "prediction.score.max"]
    label_score_cols = []
    records = []
    for col in score_cols:
        label = col.removeprefix("prediction.score.")
        token = sanitize_score_label(label)
        out_col = f"{prefix}_prediction_score_{token}"
        adata.obs[out_col] = pd.to_numeric(aligned[col], errors="coerce").to_numpy()
        label_score_cols.append(out_col)
        records.append({"shi_label": label, "score_column": out_col})

    score_map = pd.DataFrame.from_records(records)
    adata.uns["shi_seurat_full_score_column_map"] = score_map.to_dict(orient="list")
    if label_order:
        adata.obs[f"{prefix}_predicted_shi_label"] = pd.Categorical(
            adata.obs[f"{prefix}_predicted_shi_label"].astype(str),
            categories=[str(label) for label in label_order],
        )
    return label_score_cols


def plot_score_panel(
    coords: np.ndarray,
    obs: pd.DataFrame,
    score_cols: list[str],
    labels: list[str],
    path: Path,
) -> None:
    fig, axes = plt.subplots(1, len(score_cols), figsize=(5.2 * len(score_cols), 4.8), squeeze=False)
    for ax, col, label in zip(axes.ravel(), score_cols, labels, strict=True):
        values = pd.to_numeric(obs[col], errors="coerce")
        artist = ax.scatter(
            coords[:, 0],
            coords[:, 1],
            c=values,
            s=1.0,
            cmap="Blues",
            vmin=0.0,
            vmax=1.0,
            linewidths=0,
            rasterized=True,
        )
        ax.set_title(label)
        ax.set_xlabel("UMAP 1")
        ax.set_ylabel("UMAP 2")
        ax.set_xticks([])
        ax.set_yticks([])
        fig.colorbar(artist, ax=ax, label="Seurat prediction score")
    fig.tight_layout()
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def write_long_score_table(obs: pd.DataFrame, score_cols: list[str], path: Path) -> None:
    id_cols = ["cell_id_for_join"]
    if "orig.ident" in obs.columns:
        id_cols.append("orig.ident")
    if "seurat_clusters" in obs.columns:
        id_cols.append("seurat_clusters")
    work = obs[id_cols + score_cols].copy()
    long = work.melt(
        id_vars=id_cols,
        value_vars=score_cols,
        var_name="score_column",
        value_name="score",
    )
    long["shi_label_token"] = long["score_column"].str.removeprefix("shi_seurat_full_prediction_score_")
    long.to_csv(path, sep="\t", index=False)


def main() -> None:
    repo_root = env_path("REPO_ROOT", "/home/elcrespo/Desktop/githubprojects/mge_organoid_pipeline")
    project_root = env_path("PROJECT_ROOT", "/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder")
    conda_env_bin = env_path("CONDA_ENV_BIN", "/home/elcrespo/miniconda3/envs/mge-organoid-python/bin")

    results_dirname = os.environ.get("SHI_SEURAT_RESULTS_DIRNAME", "shi_reference_div30_seurat_label_transfer")
    run_label = os.environ.get("SHI_SEURAT_RUN_LABEL", "shi_reference_div30_seurat_label_transfer_v1")
    run_dir = project_root / "results" / results_dirname / run_label
    seurat_dir = run_dir / "seurat"
    table_dir = run_dir / "tables"
    plot_dir = run_dir / "plots"
    h5ad_dir = run_dir / "h5ad"
    for directory in [seurat_dir, table_dir, plot_dir, h5ad_dir]:
        directory.mkdir(parents=True, exist_ok=True)

    reference_h5ad = env_path("SHI_SEURAT_REFERENCE_H5AD", str(project_root / "results/python_anndata/shi_2019_paper_qc.h5ad"))
    query_h5ad = env_path("SHI_SEURAT_QUERY_H5AD", str(project_root / "results/python_anndata/varela_div30.h5ad"))
    reference_seurat = env_path(
        "SHI_SEURAT_REFERENCE_RDS",
        str(project_root / "results/shi_2019_paper_qc/shi_2019_seurat.rds"),
    )
    query_seurat = env_path(
        "SHI_SEURAT_QUERY_RDS",
        str(project_root / "results/varela_this_paper/varela_this_paper_seurat.rds"),
    )
    table_s2_xlsx = env_path(
        "SHI_SEURAT_TABLE_S2_XLSX",
        str(project_root / "reference/shi_2021_tables_s2_to_s9/science.abj6641_table_s2.xlsx"),
    )
    rscript_bin = env_path("RSCRIPT", str(conda_env_bin / "Rscript"))
    dims = env_int("SHI_SEURAT_DIMS", 50)
    nfeatures = env_int("SHI_SEURAT_NFEATURES", 3000)
    npcs = env_int("SHI_SEURAT_NPCS", 50)
    seed = env_int("SHI_SEURAT_RANDOM_STATE", 0)
    min_match_fraction = float(os.environ.get("SHI_SEURAT_MIN_LABEL_MATCH_FRACTION", "0.95"))
    overwrite = env_bool("SHI_SEURAT_OVERWRITE", False)
    write_h5ad = env_bool("SHI_SEURAT_WRITE_H5AD", True)

    labels_path = seurat_dir / "shi_reference_labels_for_seurat.tsv"
    prediction_path = seurat_dir / "div30_shi_seurat_full_predictions.tsv.gz"
    score_path = seurat_dir / "div30_shi_seurat_full_prediction_scores.tsv.gz"

    log(f"Run label: {run_label}")
    log(f"Run dir: {run_dir}")
    label_order = prepare_reference_labels(
        reference_h5ad=reference_h5ad,
        table_s2_xlsx=table_s2_xlsx,
        labels_path=labels_path,
        table_dir=table_dir,
        min_match_fraction=min_match_fraction,
    )

    if overwrite or not prediction_path.exists() or not score_path.exists():
        run_seurat_transfer(
            repo_root=repo_root,
            rscript_bin=rscript_bin,
            reference_seurat=reference_seurat,
            query_seurat=query_seurat,
            labels_path=labels_path,
            seurat_dir=seurat_dir,
            dims=dims,
            nfeatures=nfeatures,
            npcs=npcs,
            seed=seed,
        )
    else:
        log(f"Using existing Seurat predictions: {prediction_path}")

    log(f"Reading DIV30 AnnData for plotting: {query_h5ad}")
    adata = ad.read_h5ad(query_h5ad)
    predictions = load_predictions(prediction_path)
    aligned = align_predictions_to_adata(adata, predictions)
    score_cols = add_predictions_to_obs(adata, aligned, label_order=label_order)
    if "cell_id" in adata.obs.columns:
        adata.obs["cell_id_for_join"] = adata.obs["cell_id"].astype(str)
    else:
        adata.obs["cell_id_for_join"] = adata.obs_names.astype(str)

    prefix = "shi_seurat_full"
    obs_cols = [
        "cell_id_for_join",
        "orig.ident",
        "seurat_clusters",
        f"{prefix}_predicted_shi_label",
        f"{prefix}_prediction_score",
        f"{prefix}_uncertainty_score",
        f"{prefix}_broad_region_class",
        f"{prefix}_developmental_class",
    ]
    obs_cols = [col for col in obs_cols if col in adata.obs.columns]
    adata.obs[obs_cols + score_cols].to_csv(
        table_dir / "div30_shi_seurat_label_transfer_obs.tsv.gz",
        sep="\t",
        index=True,
    )
    write_long_score_table(
        adata.obs,
        score_cols,
        table_dir / "div30_shi_seurat_label_transfer_label_scores_long.tsv.gz",
    )

    cluster_tables = summarize_predictions_by_cluster(
        adata.obs,
        label_col=f"{prefix}_predicted_shi_label",
        score_col=f"{prefix}_prediction_score",
        uncertainty_col=f"{prefix}_uncertainty_score",
        cluster_col="seurat_clusters",
        label_categories=label_order,
    )
    cluster_tables["cluster_summary"].to_csv(
        table_dir / "div30_shi_seurat_full_cluster_summaries.tsv",
        sep="\t",
        index=False,
    )
    cluster_tables["label_counts"].to_csv(
        table_dir / "div30_shi_seurat_full_label_fractions_by_cluster.tsv",
        sep="\t",
        index=False,
    )

    coords, umap_key = get_umap_coordinates(adata)
    log(f"Using UMAP key: {umap_key}")
    plot_umap_categorical(
        coords,
        adata.obs[f"{prefix}_predicted_shi_label"],
        "DIV30 Seurat-transfer Shi label",
        plot_dir / "div30_umap_shi_seurat_full_predicted_shi_label.png",
        point_size=1.0,
        category_order=label_order,
    )
    plot_umap_continuous(
        coords,
        adata.obs[f"{prefix}_prediction_score"],
        "DIV30 Seurat-transfer max prediction score",
        plot_dir / "div30_umap_shi_seurat_full_prediction_score.png",
        point_size=1.0,
        cmap="viridis",
    )

    panel_labels = ["MGE", "LGE", "CGE"]
    panel_cols = [f"{prefix}_prediction_score_{sanitize_score_label(label)}" for label in panel_labels]
    missing_panel = [col for col in panel_cols if col not in adata.obs.columns]
    if missing_panel:
        raise ValueError("Missing expected MGE/LGE/CGE score columns: " + ", ".join(missing_panel))
    plot_score_panel(
        coords,
        adata.obs,
        panel_cols,
        panel_labels,
        plot_dir / "div30_umap_shi_seurat_full_mge_lge_cge_score_panel.png",
    )

    if "orig.ident" in adata.obs.columns:
        plot_overlaid_density_by_group(
            adata.obs,
            value_col=f"{prefix}_prediction_score",
            group_col="orig.ident",
            title="Seurat-transfer max prediction score by DIV30 sample",
            path=plot_dir / "div30_shi_seurat_full_prediction_score_density_by_sample.png",
            bins=np.linspace(0, 1, 41),
            x_label="Seurat max prediction score",
        )

    plot_stacked_bar(
        cluster_tables["label_counts"],
        cluster_col="seurat_clusters",
        label_col=f"{prefix}_predicted_shi_label",
        path=plot_dir / "div30_shi_seurat_full_shi_label_stacked_bar_by_seurat_clusters.png",
        title="Seurat-transfer Shi label fractions by DIV30 seurat_clusters",
    )

    manifest = []
    for root in [seurat_dir, table_dir, plot_dir, h5ad_dir]:
        for path in sorted(root.glob("*")):
            if path.is_file():
                manifest.append(
                    {
                        "path": str(path),
                        "relative_path": str(path.relative_to(run_dir)),
                        "size_bytes": path.stat().st_size,
                    }
                )
    pd.DataFrame(manifest).to_csv(table_dir / "shi_seurat_label_transfer_output_manifest.tsv", sep="\t", index=False)

    if write_h5ad:
        out_h5ad = h5ad_dir / "div30_shi_seurat_label_transfer_predictions.h5ad"
        log(f"Writing annotated AnnData: {out_h5ad}")
        adata.write_h5ad(out_h5ad, compression="gzip")

    completion = pd.DataFrame(
        [
            {
                "run_label": run_label,
                "status": "complete",
                "n_query_cells": int(adata.n_obs),
                "n_score_columns": int(len(score_cols)),
                "n_plots": len(list(plot_dir.glob("*.png"))),
            }
        ]
    )
    completion.to_csv(table_dir / "shi_seurat_label_transfer_complete.tsv", sep="\t", index=False)
    log("complete")


if __name__ == "__main__":
    main()
