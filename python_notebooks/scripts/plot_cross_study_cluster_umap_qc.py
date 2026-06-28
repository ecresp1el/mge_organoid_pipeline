#!/usr/bin/env python3
"""Plot cross-study UMAP cluster QC grids from prepared marker-expression tables."""

from __future__ import annotations

import argparse
import logging
from pathlib import Path
import re
from typing import Iterable

import matplotlib

matplotlib.use("Agg")
matplotlib.rcParams.update(
    {
        "font.family": "Arial",
        "font.sans-serif": ["Arial", "Nimbus Sans", "Liberation Sans", "DejaVu Sans"],
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    }
)
logging.getLogger("matplotlib.font_manager").setLevel(logging.ERROR)

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd


DEFAULT_PROJECT_ROOT = Path("/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder")
DEFAULT_RUN_LABEL = "cross_study_marker_expression_v12"
BASE_COLUMNS = ["cell_id", "study_id", "study_label", "sample", "cluster", "umap_1", "umap_2"]
UMAP_POINT_SIZE = 1.6
FIGURE_EXPORT_DPI = 600

DIV30_PAPER_CLUSTER_MAP = {
    "0": ("1", "Radial glia"),
    "3": ("1", "Radial glia"),
    "7": ("1", "Radial glia"),
    "6": ("2", "Inhibitory progenitors"),
    "1": ("3", "SST+ cIN"),
    "4": ("4", "PV neuron precursor"),
    "2": ("5", "MGE subpallial neurons"),
}

DIV90_PUBLISHED_RECODE = {
    "0": {"id": "3", "label": "MGE Striatal/GP fated", "keep": True, "exclude_reason": ""},
    "1": {"id": "7", "label": "PV Precursors/Migrating cells/Cortical fated", "keep": True, "exclude_reason": ""},
    "2": {"id": "2", "label": "CRABP1+/PV Precursors", "keep": True, "exclude_reason": ""},
    "3": {"id": "1", "label": "SST+, NPY+ Cortical fated", "keep": True, "exclude_reason": ""},
    "4": {"id": "8", "label": "Pre-Astrocytes/Astrocytes", "keep": True, "exclude_reason": ""},
    "5": {"id": "4", "label": "LHX8+ vMGE GABAergic Striatal/GP fated 1", "keep": True, "exclude_reason": ""},
    "6": {"id": "", "label": "", "keep": False, "exclude_reason": "Stressed Cells"},
    "7": {"id": "", "label": "", "keep": False, "exclude_reason": "Stressed Cells"},
    "8": {"id": "5", "label": "LHX8+ vMGE GABAergic Striatal/GP fated 2", "keep": True, "exclude_reason": ""},
    "9": {"id": "10", "label": "Pre-OPCs/OPCs", "keep": True, "exclude_reason": ""},
    "10": {"id": "8", "label": "Pre-Astrocytes/Astrocytes", "keep": True, "exclude_reason": ""},
    "11": {"id": "6", "label": "PV Precursors", "keep": True, "exclude_reason": ""},
    "12": {"id": "9", "label": "Dividing cells", "keep": True, "exclude_reason": ""},
}

DIV90_PUBLISHED_ORDER = [
    "1. SST+, NPY+ Cortical fated",
    "2. CRABP1+/PV Precursors",
    "3. MGE Striatal/GP fated",
    "4. LHX8+ vMGE GABAergic Striatal/GP fated 1",
    "5. LHX8+ vMGE GABAergic Striatal/GP fated 2",
    "6. PV Precursors",
    "7. PV Precursors/Migrating cells/Cortical fated",
    "8. Pre-Astrocytes/Astrocytes",
    "9. Dividing cells",
    "10. Pre-OPCs/OPCs",
]


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


def apply_mapped_cluster_labels(data: pd.DataFrame, study_id: str, project_root: Path, div90_mapping: pd.DataFrame) -> pd.DataFrame:
    out = data.copy()
    out["raw_cluster"] = out["cluster"].astype(str)
    out["mapped_cluster_id"] = out["raw_cluster"]
    out["mapped_cluster_name"] = out["raw_cluster"]
    out["mapped_cluster_label"] = out["raw_cluster"]
    out["cluster_mapping_source"] = "raw_cluster"
    out["published_cluster_id"] = out["mapped_cluster_id"]
    out["published_cluster_name"] = out["mapped_cluster_name"]
    out["published_cluster_label"] = out["mapped_cluster_label"]
    out["published_keep"] = True
    out["published_exclude_reason"] = ""
    out["published_recode_source"] = "same_as_mapped_cluster"

    if study_id == "varela_div30":
        mapped = out["raw_cluster"].map(DIV30_PAPER_CLUSTER_MAP)
        if mapped.isna().any():
            missing = sorted(out.loc[mapped.isna(), "raw_cluster"].unique(), key=natural_sort_key)
            raise ValueError(f"Unmapped DIV30 clusters in marker table: {missing}")
        out["mapped_cluster_id"] = mapped.map(lambda value: value[0])
        out["mapped_cluster_name"] = mapped.map(lambda value: value[1])
        out["mapped_cluster_label"] = out["mapped_cluster_id"] + " - " + out["mapped_cluster_name"]
        out["cluster_mapping_source"] = "div30_paper_cluster_annotation_mapping"
        out["published_cluster_id"] = out["mapped_cluster_id"]
        out["published_cluster_name"] = out["mapped_cluster_name"]
        out["published_cluster_label"] = out["mapped_cluster_label"]
        out["published_recode_source"] = "div30_paper_cluster_annotation_mapping"
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
        out["cluster_mapping_source"] = str(div90_mapping_path(project_root))
        out = out.drop(columns=[col for col in out.columns if col.endswith("_div90")])
        recoded = out["raw_cluster"].map(DIV90_PUBLISHED_RECODE)
        if recoded.isna().any():
            missing = sorted(out.loc[recoded.isna(), "raw_cluster"].unique(), key=natural_sort_key)
            raise ValueError(f"Unmapped DIV90 published recode clusters in marker table: {missing}")
        out["published_cluster_id"] = recoded.map(lambda value: value["id"])
        out["published_cluster_name"] = recoded.map(lambda value: value["label"])
        out["published_keep"] = recoded.map(lambda value: bool(value["keep"]))
        out["published_exclude_reason"] = recoded.map(lambda value: value["exclude_reason"])
        out["published_cluster_label"] = np.where(
            out["published_keep"],
            out["published_cluster_id"].astype(str) + ". " + out["published_cluster_name"].astype(str),
            "EXCLUDED - " + out["published_exclude_reason"].astype(str),
        )
        out["published_recode_source"] = "div90_published_fig_d_10_class_recode"
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
        data = apply_mapped_cluster_labels(data, study_id, project_root, div90_mapping)
        manifest["n_cells_loaded"] = int(data.shape[0])
        manifest["mapping_source"] = str(data["cluster_mapping_source"].iloc[0])
        manifest["published_recode_source"] = str(data["published_recode_source"].iloc[0])
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
        if not kept.empty and "published_keep" in kept.columns:
            before_published_filter = kept.shape[0]
            kept = kept.loc[kept["published_keep"].astype(bool)].copy()
            if kept.shape[0] < before_published_filter:
                removed_reasons = sorted(group.loc[~group["published_keep"].astype(bool), "published_exclude_reason"].dropna().astype(str).unique())
                rule = f"{rule};published_recode_keep_filter" if rule != "none" else "published_recode_keep_filter"
                reason_text = ",".join(reason for reason in removed_reasons if reason)
                filter_label = f"{filter_label}; Published recode excludes {reason_text}".strip("; ")
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


def visualization_data(data: pd.DataFrame) -> pd.DataFrame:
    """Remove cells flagged out of published-style plotting while preserving tables."""
    if "published_keep" not in data.columns:
        return data.copy()
    return data.loc[data["published_keep"].astype(bool)].copy()


def cluster_counts(data: pd.DataFrame, scope: str) -> pd.DataFrame:
    counts = (
        data.groupby(
            [
                "study_id",
                "study_label",
                "raw_cluster",
                "mapped_cluster_id",
                "mapped_cluster_name",
                "mapped_cluster_label",
                "cluster_mapping_source",
                "published_cluster_id",
                "published_cluster_name",
                "published_cluster_label",
                "published_keep",
                "published_exclude_reason",
                "published_recode_source",
            ],
            sort=False,
        )
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
        data.groupby(
            [
                "study_id",
                "study_label",
                "sample",
                "raw_cluster",
                "mapped_cluster_id",
                "mapped_cluster_name",
                "mapped_cluster_label",
                "published_cluster_id",
                "published_cluster_name",
                "published_cluster_label",
                "published_keep",
                "published_exclude_reason",
            ],
            sort=False,
        )
        .size()
        .reset_index(name="n_cells")
    )
    out["scope"] = scope
    return out


def study_summary(data: pd.DataFrame, scope: str) -> pd.DataFrame:
    rows = []
    for (study_id, study_label), group in data.groupby(["study_id", "study_label"], sort=False):
        clusters = sorted(group["published_cluster_label"].astype(str).unique(), key=natural_sort_key)
        rows.append(
            {
                "scope": scope,
                "study_id": study_id,
                "study_label": study_label,
                "n_cells": int(group.shape[0]),
                "n_clusters": len(clusters),
                "cluster_labels": ";".join(clusters),
                "raw_cluster_labels": ";".join(sorted(group["raw_cluster"].astype(str).unique(), key=natural_sort_key)),
                "excluded_raw_cluster_labels": ";".join(sorted(group.loc[~group["published_keep"].astype(bool), "raw_cluster"].astype(str).unique(), key=natural_sort_key)),
                "excluded_reasons": ";".join(sorted(reason for reason in group["published_exclude_reason"].dropna().astype(str).unique() if reason)),
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
    return label_positions_for_coordinates(group, "umap_1", "umap_2")


def label_positions_for_coordinates(group: pd.DataFrame, x_col: str, y_col: str) -> pd.DataFrame:
    rows = []
    for cluster, subset in group.groupby("published_cluster_label", sort=False, observed=True):
        rows.append(
            {
                "cluster": cluster,
                "cluster_id": str(subset["published_cluster_id"].iloc[0]),
                "x": float(np.nanmedian(subset[x_col].to_numpy(dtype=float))),
                "y": float(np.nanmedian(subset[y_col].to_numpy(dtype=float))),
                "n_cells": int(subset.shape[0]),
            }
        )
    return pd.DataFrame(rows)


def add_plot_coordinates(subset: pd.DataFrame, study_id: str) -> pd.DataFrame:
    out = subset.copy()
    out["UMAP1_plot"] = out["umap_1"]
    out["UMAP2_plot"] = out["umap_2"]
    if study_id == "varela_div90":
        out["UMAP2_plot"] = -1.0 * out["umap_2"]
    return out


def plot_cluster_grid(data: pd.DataFrame, output_prefix: Path, title: str, max_label_clusters: int = 80) -> None:
    study_order = (
        data[["study_order", "study_id", "study_label"]]
        .drop_duplicates()
        .sort_values(["study_order", "study_id"])
    )
    study_keys = list(study_order[["study_id", "study_label"]].itertuples(index=False, name=None))
    n = len(study_keys)
    ncols = n
    nrows = 1
    fig, axes = plt.subplots(nrows, ncols, figsize=(4.6 * ncols, 4.4), squeeze=False)
    for ax in axes.ravel()[n:]:
        ax.axis("off")

    for ax, (study_id, study_label) in zip(axes.ravel(), study_keys):
        subset = data.loc[(data["study_id"] == study_id) & (data["study_label"] == study_label)].copy()
        subset = subset[np.isfinite(subset["umap_1"]) & np.isfinite(subset["umap_2"])]
        subset = add_plot_coordinates(subset, study_id)
        clusters = sorted(subset["published_cluster_label"].astype(str).unique(), key=natural_sort_key)
        colors = cluster_palette(clusters)
        for cluster in clusters:
            cluster_data = subset.loc[subset["published_cluster_label"] == cluster]
            ax.scatter(
                cluster_data["UMAP1_plot"],
                cluster_data["UMAP2_plot"],
                s=UMAP_POINT_SIZE,
                c=[colors[cluster]],
                linewidths=0,
                rasterized=True,
            )
        positions = label_positions_for_coordinates(subset, "UMAP1_plot", "UMAP2_plot")
        if positions.shape[0] <= max_label_clusters:
            for _, row in positions.iterrows():
                ax.text(
                    row["x"],
                    row["y"],
                    str(row["cluster_id"]),
                    ha="center",
                    va="center",
                    fontsize=6.5,
                    fontweight="bold",
                    color="black",
                    bbox={"boxstyle": "round,pad=0.12", "facecolor": "white", "edgecolor": "black", "linewidth": 0.25, "alpha": 0.78},
                )
        ax.set_title(f"{study_label}\n{len(clusters)} clusters, n={subset.shape[0]:,}", fontsize=10)
        ax.set_aspect("equal", adjustable="box")
        ax.set_axis_off()
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
    fig.savefig(output_prefix.with_suffix(".png"), dpi=FIGURE_EXPORT_DPI, bbox_inches="tight")
    fig.savefig(output_prefix.with_suffix(".pdf"), dpi=FIGURE_EXPORT_DPI, bbox_inches="tight")
    fig.savefig(output_prefix.with_suffix(".svg"), dpi=FIGURE_EXPORT_DPI, bbox_inches="tight")
    plt.close(fig)


def div90_published_data(data: pd.DataFrame) -> pd.DataFrame:
    div90 = data.loc[(data["study_id"] == "varela_div90") & data["published_keep"].astype(bool)].copy()
    if div90.empty:
        raise ValueError("No DIV90 cells remain after published recode keep filter.")
    div90["published_cluster_label"] = pd.Categorical(
        div90["published_cluster_label"],
        categories=DIV90_PUBLISHED_ORDER,
        ordered=True,
    )
    return div90


def div90_published_plot_data(data: pd.DataFrame) -> pd.DataFrame:
    div90 = div90_published_data(data)
    div90 = div90[np.isfinite(div90["umap_1"]) & np.isfinite(div90["umap_2"])].copy()
    if div90.empty:
        raise ValueError("No finite DIV90 UMAP coordinates remain after published recode keep filter.")
    # Plot-only orientation transform; original umap_1/umap_2 remain unchanged.
    div90["UMAP1_published"] = div90["umap_1"]
    div90["UMAP2_published"] = -1.0 * div90["umap_2"]
    return div90


def write_div90_published_recode_tables(data: pd.DataFrame, outdir: Path) -> None:
    div90_all = data.loc[data["study_id"] == "varela_div90"].copy()
    recode = (
        div90_all[
            [
                "raw_cluster",
                "mapped_cluster_label",
                "published_cluster_id",
                "published_cluster_name",
                "published_cluster_label",
                "published_keep",
                "published_exclude_reason",
                "published_recode_source",
            ]
        ]
        .drop_duplicates()
        .sort_values(["raw_cluster"], key=lambda col: col.map(lambda value: natural_sort_key(value) if col.name == "raw_cluster" else value))
    )
    recode.to_csv(outdir / "div90_published_fig_d_10_class_recode.tsv", sep="\t", index=False)

    div90 = div90_published_data(data)
    counts = (
        div90.groupby(["published_cluster_id", "published_cluster_name", "published_cluster_label"], observed=True)
        .size()
        .reset_index(name="n_cells")
    )
    counts["fraction_cells"] = counts["n_cells"] / counts["n_cells"].sum()
    counts["published_cluster_label"] = pd.Categorical(
        counts["published_cluster_label"],
        categories=DIV90_PUBLISHED_ORDER,
        ordered=True,
    )
    counts = counts.sort_values("published_cluster_label")
    counts.to_csv(outdir / "div90_published_fig_d_cluster_counts.tsv", sep="\t", index=False)

    sample_counts = (
        div90.groupby(["sample", "published_cluster_id", "published_cluster_name", "published_cluster_label"], observed=True)
        .size()
        .reset_index(name="n_cells")
    )
    sample_totals = sample_counts.groupby("sample")["n_cells"].transform("sum")
    sample_counts["fraction_cells"] = sample_counts["n_cells"] / sample_totals
    sample_counts["published_cluster_label"] = pd.Categorical(
        sample_counts["published_cluster_label"],
        categories=DIV90_PUBLISHED_ORDER,
        ordered=True,
    )
    sample_counts = sample_counts.sort_values(["sample", "published_cluster_label"])
    sample_counts.to_csv(outdir / "div90_published_fig_d_sample_composition.tsv", sep="\t", index=False)


def plot_div90_published_umap(data: pd.DataFrame, output_prefix: Path) -> None:
    div90 = div90_published_plot_data(data)
    labels = [label for label in DIV90_PUBLISHED_ORDER if label in set(div90["published_cluster_label"].astype(str))]
    colors = cluster_palette(labels)
    fig, ax = plt.subplots(figsize=(7.4, 5.8))
    for label in labels:
        subset = div90.loc[div90["published_cluster_label"].astype(str) == label]
        ax.scatter(
            subset["UMAP1_published"],
            subset["UMAP2_published"],
            s=UMAP_POINT_SIZE,
            c=[colors[label]],
            linewidths=0,
            rasterized=True,
            label=label,
        )
    positions = label_positions_for_coordinates(div90, "UMAP1_published", "UMAP2_published")
    for _, row in positions.iterrows():
        ax.text(
            row["x"],
            row["y"],
            str(row["cluster_id"]),
            ha="center",
            va="center",
            fontsize=7,
            fontweight="bold",
            color="black",
            bbox={"boxstyle": "round,pad=0.12", "facecolor": "white", "edgecolor": "black", "linewidth": 0.25, "alpha": 0.82},
        )
    ax.set_title(f"DIV90 published Fig. D classes\n{len(labels)} clusters, n={div90.shape[0]:,}", fontsize=11)
    ax.set_aspect("equal", adjustable="box")
    ax.set_axis_off()
    handles = [Line2D([0], [0], marker="o", color="none", markerfacecolor=colors[label], markersize=5, label=label) for label in labels]
    ax.legend(handles=handles, title="Published class", loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=False, fontsize=7, title_fontsize=8)
    fig.tight_layout()
    fig.savefig(output_prefix.with_suffix(".png"), dpi=FIGURE_EXPORT_DPI, bbox_inches="tight")
    fig.savefig(output_prefix.with_suffix(".pdf"), dpi=FIGURE_EXPORT_DPI, bbox_inches="tight")
    fig.savefig(output_prefix.with_suffix(".svg"), dpi=FIGURE_EXPORT_DPI, bbox_inches="tight")
    plt.close(fig)


def plot_div90_sample_composition(data: pd.DataFrame, output_prefix: Path) -> None:
    div90 = div90_published_data(data)
    labels = [label for label in DIV90_PUBLISHED_ORDER if label in set(div90["published_cluster_label"].astype(str))]
    colors = cluster_palette(labels)
    counts = (
        div90.groupby(["sample", "published_cluster_label"], observed=True)
        .size()
        .reset_index(name="n_cells")
    )
    pivot = counts.pivot_table(index="sample", columns="published_cluster_label", values="n_cells", fill_value=0, observed=True)
    pivot = pivot.reindex(columns=labels, fill_value=0)
    fractions = pivot.div(pivot.sum(axis=1), axis=0)
    fig, ax = plt.subplots(figsize=(8.2, 4.8))
    bottom = np.zeros(fractions.shape[0])
    x = np.arange(fractions.shape[0])
    for label in labels:
        values = fractions[label].to_numpy(dtype=float)
        ax.bar(x, values, bottom=bottom, color=colors[label], width=0.78, label=label)
        bottom += values
    ax.set_xticks(x)
    ax.set_xticklabels(fractions.index.astype(str), rotation=35, ha="right", fontsize=8)
    ax.set_ylim(0, 1)
    ax.set_ylabel("Fraction of cells")
    ax.set_title("DIV90 sample composition by published Fig. D class", fontsize=11)
    ax.legend(title="Published class", loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=False, fontsize=7, title_fontsize=8)
    fig.tight_layout()
    fig.savefig(output_prefix.with_suffix(".png"), dpi=FIGURE_EXPORT_DPI, bbox_inches="tight")
    fig.savefig(output_prefix.with_suffix(".pdf"), dpi=FIGURE_EXPORT_DPI, bbox_inches="tight")
    fig.savefig(output_prefix.with_suffix(".svg"), dpi=FIGURE_EXPORT_DPI, bbox_inches="tight")
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
    write_div90_published_recode_tables(data, outdir)

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
        visualization_data(data),
        outdir / "cross_study_marker_expression_v12_cluster_umap_qc_all_prepared_cells",
        "Cross-study cluster UMAP QC - all prepared marker tables",
    )
    plot_div90_published_umap(
        data,
        outdir / "div90_published_fig_d_10_class_umap",
    )
    plot_div90_sample_composition(
        data,
        outdir / "div90_published_fig_d_sample_composition",
    )
    print(f"Wrote cluster QC outputs to: {outdir}", flush=True)


if __name__ == "__main__":
    main()
