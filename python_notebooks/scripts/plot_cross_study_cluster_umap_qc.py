#!/usr/bin/env python3
"""Plot cross-study UMAP cluster QC grids from prepared marker-expression tables."""

from __future__ import annotations

import argparse
from pathlib import Path
import re
from typing import Iterable

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd


DEFAULT_PROJECT_ROOT = Path("/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder")
DEFAULT_RUN_LABEL = "cross_study_marker_expression_v12"
BASE_COLUMNS = ["cell_id", "study_id", "study_label", "sample", "cluster", "umap_1", "umap_2"]

DIV30_PAPER_CLUSTER_MAP = {
    "0": ("1", "Radial glia"),
    "3": ("1", "Radial glia"),
    "7": ("1", "Radial glia"),
    "6": ("2", "Inhibitory progenitors"),
    "1": ("3", "SST+ cIN"),
    "4": ("4", "PV neuron precursor"),
    "2": ("5", "MGE subpallial neurons"),
}


def parse_study_ids(raw: Iterable[str] | None) -> list[str]:
    values: list[str] = []
    for item in raw or []:
        values.extend(part for part in re.split(r"[,;\s]+", item) if part and part.lower() not in {"none", "null", "na"})
    return values


def natural_sort_key(value: object) -> list[object]:
    parts = re.split(r"(\d+)", str(value))
    return [int(part) if part.isdigit() else part.lower() for part in parts]


def run_dir(project_root: Path, run_label: str) -> Path:
    return project_root / "results" / "cross_study_marker_expression" / run_label


def table_dir(project_root: Path, run_label: str) -> Path:
    return run_dir(project_root, run_label) / "tables"


def per_study_table_path(project_root: Path, run_label: str, study_id: str) -> Path:
    return table_dir(project_root, run_label) / "per_study" / f"{study_id}_marker_expression.tsv.gz"


def div90_mapping_path(project_root: Path) -> Path:
    return (
        project_root
        / "results"
        / "div90_umap_cluster_label_audit"
        / "div90_umap_cluster_label_audit_v1"
        / "tables"
        / "div90_cluster_number_name_to_biology_mapping.tsv"
    )


def load_div90_mapping(project_root: Path) -> pd.DataFrame:
    path = div90_mapping_path(project_root)
    if not path.exists():
        raise FileNotFoundError(f"Missing DIV90 mapping table: {path}")
    mapping = pd.read_csv(path, sep="\t")
    mapping["cluster"] = mapping["cluster_id_numeric"].astype(str)
    mapping["mapped_cluster_id"] = mapping["cluster_id_numeric"].astype(str)
    mapping["mapped_cluster_name"] = mapping["biology_name"].astype(str)
    mapping["mapped_cluster_label"] = mapping["exact_metadata_name"].astype(str)
    return mapping[["cluster", "mapped_cluster_id", "mapped_cluster_name", "mapped_cluster_label"]]


def load_study_table(project_root: Path, run_label: str) -> pd.DataFrame:
    path = table_dir(project_root, run_label) / "cross_study_marker_expression_studies.tsv"
    studies = pd.read_csv(path, sep="\t")
    studies["study_id"] = studies["study_id"].astype(str)
    studies["study_label"] = studies["study_label"].astype(str)
    studies["study_order"] = pd.to_numeric(studies.get("study_order", np.arange(1, studies.shape[0] + 1)), errors="coerce")
    return studies


def apply_mapped_cluster_labels(data: pd.DataFrame, study_id: str, div90_mapping: pd.DataFrame) -> pd.DataFrame:
    out = data.copy()
    out["raw_cluster"] = out["cluster"].astype(str)
    out["mapped_cluster_id"] = out["raw_cluster"]
    out["mapped_cluster_name"] = out["raw_cluster"]
    out["mapped_cluster_label"] = out["raw_cluster"]
    out["cluster_mapping_source"] = "raw_cluster"

    if study_id == "varela_div30":
        mapped = out["raw_cluster"].map(DIV30_PAPER_CLUSTER_MAP)
        if mapped.isna().any():
            missing = sorted(out.loc[mapped.isna(), "raw_cluster"].unique(), key=natural_sort_key)
            raise ValueError(f"Unmapped DIV30 clusters in marker table: {missing}")
        out["mapped_cluster_id"] = mapped.map(lambda value: value[0])
        out["mapped_cluster_name"] = mapped.map(lambda value: value[1])
        out["mapped_cluster_label"] = out["mapped_cluster_id"] + " - " + out["mapped_cluster_name"]
        out["cluster_mapping_source"] = "div30_paper_cluster_annotation_mapping"
    elif study_id == "varela_div90":
        before = out.shape[0]
        out = out.merge(div90_mapping, on="cluster", how="left", suffixes=("", "_div90"))
        if out.shape[0] != before:
            raise ValueError("DIV90 mapping merge changed row count.")
        if out["mapped_cluster_label_div90"].isna().any():
            missing = sorted(out.loc[out["mapped_cluster_label_div90"].isna(), "raw_cluster"].unique(), key=natural_sort_key)
            raise ValueError(f"Unmapped DIV90 clusters in marker table: {missing}")
        out["mapped_cluster_id"] = out["mapped_cluster_id_div90"]
        out["mapped_cluster_name"] = out["mapped_cluster_name_div90"]
        out["mapped_cluster_label"] = out["mapped_cluster_label_div90"]
        out["cluster_mapping_source"] = str(div90_mapping_path(DEFAULT_PROJECT_ROOT))
        out = out.drop(columns=[col for col in out.columns if col.endswith("_div90")])
    return out


def load_cluster_tables(project_root: Path, run_label: str, studies: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    frames = []
    manifest_rows = []
    div90_mapping = load_div90_mapping(project_root)
    for _, row in studies.sort_values("study_order").iterrows():
        study_id = str(row["study_id"])
        path = per_study_table_path(project_root, run_label, study_id)
        exists = path.exists()
        manifest = {
            "study_id": study_id,
            "study_label": row["study_label"],
            "study_order": row["study_order"],
            "include_in_first_plot": row.get("include_in_first_plot", ""),
            "marker_table_path": str(path),
            "marker_table_exists": exists,
            "source_seurat_path": row.get("seurat_path", ""),
            "source_h5ad_path": row.get("h5ad_path", ""),
            "reduction": row.get("reduction", ""),
            "cluster_col_source": row.get("cluster_col", ""),
            "sample_col_source": row.get("sample_col", ""),
            "n_cells_loaded": 0,
        }
        if not exists:
            manifest_rows.append(manifest)
            continue
        data = pd.read_csv(path, sep="\t", usecols=BASE_COLUMNS, low_memory=False)
        data["study_id"] = data["study_id"].astype(str)
        data["study_label"] = data["study_label"].astype(str)
        data["sample"] = data["sample"].astype(str)
        data["cluster"] = data["cluster"].astype(str)
        data["umap_1"] = pd.to_numeric(data["umap_1"], errors="coerce")
        data["umap_2"] = pd.to_numeric(data["umap_2"], errors="coerce")
        data["study_order"] = row["study_order"]
        data = apply_mapped_cluster_labels(data, study_id, div90_mapping)
        manifest["n_cells_loaded"] = int(data.shape[0])
        manifest["mapping_source"] = str(data["cluster_mapping_source"].iloc[0])
        manifest_rows.append(manifest)
        frames.append(data)
    if not frames:
        raise FileNotFoundError("No per-study marker-expression tables could be loaded.")
    return pd.concat(frames, ignore_index=True), pd.DataFrame(manifest_rows)


def apply_figure_filters(data: pd.DataFrame, exclude_study_ids: list[str]) -> tuple[pd.DataFrame, pd.DataFrame]:
    summaries = []
    parts = []
    excluded = set(exclude_study_ids)
    for (study_id, study_label), group in data.groupby(["study_id", "study_label"], sort=False):
        rule = "none"
        filter_label = ""
        kept = group
        excluded_by_study = study_id in excluded
        if excluded_by_study:
            rule = "excluded_study"
            filter_label = "Excluded from v12 figure"
            kept = group.iloc[0:0].copy()
        elif study_id == "samarasinghe_2021":
            rule = "samarasinghe_2021_controls_only"
            filter_label = "Controls only"
            kept = group.loc[group["sample"].astype(str).str.contains("Ctrl", case=False, na=False)].copy()
        summaries.append(
            {
                "study_id": study_id,
                "study_label": study_label,
                "plot_filter": rule,
                "plot_filter_label": filter_label,
                "n_cells_before_filter": int(group.shape[0]),
                "n_cells_after_filter": int(kept.shape[0]),
                "n_cells_removed": int(group.shape[0] - kept.shape[0]),
                "sample_values_after_filter": ";".join(sorted(kept["sample"].astype(str).unique())) if not kept.empty else "",
            }
        )
        if not kept.empty:
            parts.append(kept)
    if not parts:
        raise ValueError("All cells were removed by figure filters.")
    return pd.concat(parts, ignore_index=True), pd.DataFrame(summaries)


def cluster_counts(data: pd.DataFrame, scope: str) -> pd.DataFrame:
    counts = (
        data.groupby(["study_id", "study_label", "raw_cluster", "mapped_cluster_id", "mapped_cluster_name", "mapped_cluster_label", "cluster_mapping_source"], sort=False)
        .size()
        .reset_index(name="n_cells")
    )
    counts["scope"] = scope
    totals = counts.groupby(["study_id", "study_label"], sort=False)["n_cells"].transform("sum")
    counts["fraction_cells"] = counts["n_cells"] / totals
    counts["cluster_sort_key"] = counts["mapped_cluster_id"].map(lambda x: repr(natural_sort_key(x)))
    counts = counts.sort_values(["study_id", "cluster_sort_key", "raw_cluster"]).drop(columns=["cluster_sort_key"])
    return counts


def sample_cluster_counts(data: pd.DataFrame, scope: str) -> pd.DataFrame:
    out = (
        data.groupby(["study_id", "study_label", "sample", "raw_cluster", "mapped_cluster_id", "mapped_cluster_name", "mapped_cluster_label"], sort=False)
        .size()
        .reset_index(name="n_cells")
    )
    out["scope"] = scope
    return out


def study_summary(data: pd.DataFrame, scope: str) -> pd.DataFrame:
    rows = []
    for (study_id, study_label), group in data.groupby(["study_id", "study_label"], sort=False):
        clusters = sorted(group["mapped_cluster_label"].astype(str).unique(), key=natural_sort_key)
        rows.append(
            {
                "scope": scope,
                "study_id": study_id,
                "study_label": study_label,
                "n_cells": int(group.shape[0]),
                "n_clusters": len(clusters),
                "cluster_labels": ";".join(clusters),
                "raw_cluster_labels": ";".join(sorted(group["raw_cluster"].astype(str).unique(), key=natural_sort_key)),
                "n_samples": int(group["sample"].nunique()),
                "sample_values": ";".join(sorted(group["sample"].astype(str).unique())),
            }
        )
    return pd.DataFrame(rows)


def cluster_palette(labels: list[str]) -> dict[str, tuple[float, float, float, float]]:
    cmap_names = ["tab20", "tab20b", "tab20c", "Set3", "Dark2", "Accent"]
    colors = []
    for name in cmap_names:
        cmap = plt.get_cmap(name)
        n = getattr(cmap, "N", 8)
        colors.extend(cmap(i) for i in range(n))
    return {label: colors[i % len(colors)] for i, label in enumerate(labels)}


def label_positions(group: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for cluster, subset in group.groupby("mapped_cluster_label", sort=False):
        rows.append(
            {
                "cluster": cluster,
                "x": float(np.nanmedian(subset["umap_1"].to_numpy(dtype=float))),
                "y": float(np.nanmedian(subset["umap_2"].to_numpy(dtype=float))),
                "n_cells": int(subset.shape[0]),
            }
        )
    return pd.DataFrame(rows)


def plot_cluster_grid(data: pd.DataFrame, output_prefix: Path, title: str, max_label_clusters: int = 80) -> None:
    study_order = (
        data[["study_order", "study_id", "study_label"]]
        .drop_duplicates()
        .sort_values(["study_order", "study_id"])
    )
    study_keys = list(study_order[["study_id", "study_label"]].itertuples(index=False, name=None))
    n = len(study_keys)
    ncols = min(3, n)
    nrows = int(np.ceil(n / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(5.0 * ncols, 4.4 * nrows), squeeze=False)
    for ax in axes.ravel()[n:]:
        ax.axis("off")

    for ax, (study_id, study_label) in zip(axes.ravel(), study_keys):
        subset = data.loc[(data["study_id"] == study_id) & (data["study_label"] == study_label)].copy()
        subset = subset[np.isfinite(subset["umap_1"]) & np.isfinite(subset["umap_2"])]
        clusters = sorted(subset["mapped_cluster_label"].astype(str).unique(), key=natural_sort_key)
        colors = cluster_palette(clusters)
        for cluster in clusters:
            cluster_data = subset.loc[subset["mapped_cluster_label"] == cluster]
            ax.scatter(
                cluster_data["umap_1"],
                cluster_data["umap_2"],
                s=0.35,
                c=[colors[cluster]],
                linewidths=0,
                rasterized=True,
            )
        positions = label_positions(subset)
        if positions.shape[0] <= max_label_clusters:
            for _, row in positions.iterrows():
                ax.text(
                    row["x"],
                    row["y"],
                    str(row["cluster"]).split(" - ", 1)[0],
                    ha="center",
                    va="center",
                    fontsize=6.5,
                    fontweight="bold",
                    color="black",
                    bbox={"boxstyle": "round,pad=0.12", "facecolor": "white", "edgecolor": "black", "linewidth": 0.25, "alpha": 0.78},
                )
        ax.set_title(f"{study_label}\n{len(clusters)} clusters, n={subset.shape[0]:,}", fontsize=10)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_aspect("equal", adjustable="box")
        for spine in ax.spines.values():
            spine.set_visible(False)
        legend_handles = [Line2D([0], [0], marker="o", color="none", markerfacecolor=colors[c], markersize=4, label=c) for c in clusters[:32]]
        if legend_handles:
            ax.legend(
                handles=legend_handles,
                title="Cluster",
                loc="center left",
                bbox_to_anchor=(1.01, 0.5),
                fontsize=5,
                title_fontsize=6,
                frameon=False,
                ncol=1,
                borderaxespad=0,
            )
    fig.suptitle(title, fontsize=14)
    fig.tight_layout(rect=(0, 0, 1, 0.965), w_pad=2.5, h_pad=1.2)
    output_prefix.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_prefix.with_suffix(".png"), dpi=450, bbox_inches="tight")
    fig.savefig(output_prefix.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(output_prefix.with_suffix(".svg"), bbox_inches="tight")
    plt.close(fig)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, default=DEFAULT_PROJECT_ROOT)
    parser.add_argument("--run-label", default=DEFAULT_RUN_LABEL)
    parser.add_argument("--outdir", type=Path, default=None)
    parser.add_argument("--exclude-study-id", action="append", default=["bershteyn_2025"])
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    outdir = args.outdir or run_dir(args.project_root, args.run_label) / "cluster_qc"
    outdir.mkdir(parents=True, exist_ok=True)
    studies = load_study_table(args.project_root, args.run_label)
    data, source_manifest = load_cluster_tables(args.project_root, args.run_label, studies)
    source_manifest.to_csv(outdir / "cluster_umap_qc_source_manifest.tsv", sep="\t", index=False)

    exclude_ids = parse_study_ids(args.exclude_study_id)
    figure_data, filter_summary = apply_figure_filters(data, exclude_ids)
    filter_summary.to_csv(outdir / "cluster_umap_qc_figure_filter_summary.tsv", sep="\t", index=False)

    pd.concat(
        [
            study_summary(data, "all_prepared_cells"),
            study_summary(figure_data, "figure_default_filtered"),
        ],
        ignore_index=True,
    ).to_csv(outdir / "cluster_umap_qc_study_summary.tsv", sep="\t", index=False)
    pd.concat(
        [
            cluster_counts(data, "all_prepared_cells"),
            cluster_counts(figure_data, "figure_default_filtered"),
        ],
        ignore_index=True,
    ).to_csv(outdir / "cluster_umap_qc_cluster_counts.tsv", sep="\t", index=False)
    pd.concat(
        [
            sample_cluster_counts(data, "all_prepared_cells"),
            sample_cluster_counts(figure_data, "figure_default_filtered"),
        ],
        ignore_index=True,
    ).to_csv(outdir / "cluster_umap_qc_sample_cluster_counts.tsv", sep="\t", index=False)

    plot_cluster_grid(
        figure_data,
        outdir / "cross_study_marker_expression_v12_cluster_umap_qc_figure_default",
        "Cross-study cluster UMAP QC - v12 figure-default studies",
    )
    plot_cluster_grid(
        data,
        outdir / "cross_study_marker_expression_v12_cluster_umap_qc_all_prepared_cells",
        "Cross-study cluster UMAP QC - all prepared marker tables",
    )
    print(f"Wrote cluster QC outputs to: {outdir}", flush=True)


if __name__ == "__main__":
    main()
