#!/usr/bin/env python3
"""Render final-candidate restricted Siletti river plots."""

from __future__ import annotations

import argparse
import gzip
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from plot_siletti_div90_jia_figure import PALETTE, draw_river, find_one


PALLIAL_ORDER = ["Pallial/cortical", "Subpallial"]
MAJOR_SUBTYPE_ORDER = [
    "Pallial/cortical Pvalb",
    "Pallial/cortical Chandelier",
    "Pallial/cortical Sst",
    "Pallial/cortical Lamp5 Lhx6",
    "Pallial/cortical Vip",
    "Pallial/cortical Lamp5",
    "Pallial/cortical Sncg",
    "Pallial/cortical Pax6",
    "Pallial/cortical Medium spiny neuron",
    "Subpallial Pvalb",
    "Subpallial Chandelier",
    "Subpallial Sst",
    "Subpallial Lamp5 Lhx6",
    "Subpallial Vip",
    "Subpallial Lamp5",
    "Subpallial Medium spiny neuron",
    "Subpallial Eccentric medium spiny neuron",
    "Subpallial Cholinergic neurons",
]
FINE_SUBTYPE_ORDER = [
    "Cortical PV+ basket neurons",
    "Cortical PV+ Chandelier neurons",
    "Cortical SST+ LRP neurons",
    "Cortical SST+ Mt neurons",
    "Cortical SST+ nMt neurons",
    "Pallial/cortical Medium spiny neuron",
    "Subpallial PV+ neurons",
    "Subpallial SST+ LRP neurons",
    "Subpallial SST+ neurons",
    "Subpallial Cholinergic neurons",
    "Subpallial Medium spiny neuron",
    "Subpallial Eccentric medium spiny neuron",
]
OTHER_SELECTED_REFERENCE = "Other selected reference"

PALETTE.update(
    {
        "Pallial/cortical": "#4c78a8",
        "Subpallial": "#f58518",
        "Pvalb": "#e45756",
        "Pallial/cortical Pvalb": "#e45756",
        "Subpallial Pvalb": "#f28e8b",
        "Chandelier": "#ff9da6",
        "Pallial/cortical Chandelier": "#ff9da6",
        "Subpallial Chandelier": "#ffc1ca",
        "Sst": "#54a24b",
        "Pallial/cortical Sst": "#54a24b",
        "Subpallial Sst": "#8cd17d",
        "Lamp5 Lhx6": "#b279a2",
        "Pallial/cortical Lamp5 Lhx6": "#b279a2",
        "Subpallial Lamp5 Lhx6": "#d4a6c8",
        "Vip": "#72b7b2",
        "Pallial/cortical Vip": "#72b7b2",
        "Subpallial Vip": "#a1d5d1",
        "Lamp5": "#9ecae9",
        "Pallial/cortical Lamp5": "#9ecae9",
        "Subpallial Lamp5": "#c6dbef",
        "Sncg": "#f2cf5b",
        "Pallial/cortical Sncg": "#f2cf5b",
        "Pax6": "#b5cf6b",
        "Pallial/cortical Pax6": "#b5cf6b",
        "Medium spiny neuron": "#8c6d31",
        "Eccentric medium spiny neuron": "#ad494a",
        "Pallial/cortical Medium spiny neuron": "#c7a76c",
        "Subpallial Medium spiny neuron": "#8c6d31",
        "Subpallial Eccentric medium spiny neuron": "#ad494a",
    }
)


def read_tsv_gz(path: Path) -> pd.DataFrame:
    with gzip.open(path, "rt") as handle:
        return pd.read_csv(handle, sep="\t")


def read_transfer_obs(transfer_dir: Path) -> pd.DataFrame:
    path = find_one("tables/*_query_obs_with_predictions.tsv.gz", transfer_dir)
    obs = read_tsv_gz(path)
    needed = ["seurat_cell_id", "div90_broad_class", "predicted.id", "prediction.score.max"]
    missing = [col for col in needed if col not in obs.columns]
    if missing:
        raise ValueError(f"{path} missing columns: {missing}")
    return obs


def read_bridge_reference(bridge_dir: Path) -> pd.DataFrame:
    path = bridge_dir.parent / "tables/siletti_reference_cell_metadata.tsv.gz"
    if not path.exists():
        raise FileNotFoundError(path)
    return read_tsv_gz(path)


def make_edges(obs: pd.DataFrame, source_col: str = "div90_broad_class") -> pd.DataFrame:
    edges = (
        obs.groupby([source_col, "predicted.id"], dropna=False)
        .size()
        .reset_index(name="n_cells")
        .rename(columns={source_col: "div90_class", "predicted.id": "adult_subtype"})
    )
    return edges


def complete_order(order: list[str] | None, labels: pd.Series) -> list[str] | None:
    if order is None:
        return None
    seen = set(order)
    extras = sorted(str(label) for label in labels.astype(str).unique() if str(label) not in seen)
    return list(order) + extras


def label_group(label: str) -> str:
    label = str(label)
    if label.startswith(("Pallial/cortical", "Cortical")):
        return "Pallial/cortical"
    if label.startswith("Subpallial"):
        return "Subpallial"
    if label in {"Medium spiny neuron", "Eccentric medium spiny neuron"}:
        return "Subpallial"
    return "Other"


def write_one(
    obs: pd.DataFrame,
    out_prefix: Path,
    source_col: str,
    right_order: list[str] | None,
    right_title: str,
    plot_title: str,
) -> pd.DataFrame:
    edges = make_edges(obs, source_col)
    edges.to_csv(out_prefix.with_suffix(".tsv"), sep="\t", index=False)
    right_order = complete_order(right_order, edges["adult_subtype"])
    draw_river(
        edges,
        out_prefix,
        right_order=right_order,
        right_title=right_title,
        plot_title=plot_title,
        right_group_getter=label_group,
        include_zero_right=right_order is not None,
    )
    return edges


def draw_combined_rivers(
    pallial_edges: pd.DataFrame,
    subtype_edges: pd.DataFrame,
    fine_edges: pd.DataFrame,
    out_prefix: Path,
) -> None:
    panels = [
        (
            pallial_edges,
            PALLIAL_ORDER,
            "Adult Siletti ROI bin",
            "1. Pallial vs subpallial",
        ),
        (
            subtype_edges,
            MAJOR_SUBTYPE_ORDER,
            "Adult Siletti major subtype",
            "2. Major cell types",
        ),
        (
            fine_edges,
            FINE_SUBTYPE_ORDER,
            "Adult Siletti fine subtype",
            "3. Finer cell types",
        ),
    ]
    fig, axes = plt.subplots(1, 3, figsize=(37, 13), constrained_layout=False)
    for ax, (edges, order, right_title, plot_title) in zip(axes, panels):
        order = complete_order(order, edges["adult_subtype"])
        draw_river(
            edges,
            None,
            right_order=order,
            right_title=right_title,
            plot_title=plot_title,
            right_group_getter=label_group,
            include_zero_right=True,
            ax=ax,
        )
    fig.subplots_adjust(left=0.04, right=0.985, top=0.9, bottom=0.04, wspace=0.55)
    fig.savefig(out_prefix.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(out_prefix.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def write_stacked_bar(counts: pd.DataFrame, index_col: str, column_col: str, value_col: str, out_prefix: Path, title: str) -> None:
    pivot = counts.pivot_table(index=index_col, columns=column_col, values=value_col, fill_value=0, aggfunc="sum")
    pivot = pivot.loc[pivot.sum(axis=1).sort_values(ascending=True).index]
    cols = sorted(pivot.columns.astype(str))
    pivot = pivot[cols]
    colors = dict(zip(cols, plt.get_cmap("tab20").colors[: len(cols)]))
    fig, ax = plt.subplots(figsize=(11, max(5, 0.32 * len(pivot))))
    left = np.zeros(pivot.shape[0])
    y = np.arange(pivot.shape[0])
    for col in cols:
        vals = pivot[col].to_numpy()
        ax.barh(y, vals, left=left, label=col, color=colors[col], edgecolor="white", linewidth=0.35)
        left += vals
    ax.set_yticks(y)
    ax.set_yticklabels(pivot.index.astype(str), fontsize=8)
    ax.set_xlabel("Cells")
    ax.set_title(title)
    ax.legend(frameon=False, bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=8)
    fig.savefig(out_prefix.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(out_prefix.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def write_simple_bar(counts: pd.DataFrame, label_col: str, value_col: str, out_prefix: Path, title: str) -> None:
    work = counts.sort_values(value_col, ascending=True)
    fig, ax = plt.subplots(figsize=(10, max(4, 0.32 * len(work))))
    ax.barh(np.arange(len(work)), work[value_col].to_numpy(), color="#4c78a8", edgecolor="white", linewidth=0.4)
    ax.set_yticks(np.arange(len(work)))
    ax.set_yticklabels(work[label_col].astype(str), fontsize=8)
    ax.set_xlabel("Cells")
    ax.set_title(title)
    fig.savefig(out_prefix.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(out_prefix.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def write_score_histograms(
    pallial_obs: pd.DataFrame,
    subtype_obs: pd.DataFrame,
    out_prefix: Path,
    fine_obs: pd.DataFrame | None = None,
) -> None:
    panels = [
        (pallial_obs, "Pallial/subpallial transfer scores", "#4c78a8"),
        (subtype_obs, "Major subtype transfer scores", "#f58518"),
    ]
    if fine_obs is not None:
        panels.append((fine_obs, "Fine subtype transfer scores", "#54a24b"))
    fig, axes = plt.subplots(1, len(panels), figsize=(6 * len(panels), 4), sharey=True)
    if len(panels) == 1:
        axes = [axes]
    for ax, (obs, title, color) in zip(axes, panels):
        scores = pd.to_numeric(obs["prediction.score.max"], errors="coerce").dropna()
        ax.hist(scores, bins=np.linspace(0, 1, 41), color=color, edgecolor="white", linewidth=0.35)
        ax.set_xlabel("prediction.score.max")
        ax.set_title(title)
    axes[0].set_ylabel("DIV90 cells")
    fig.savefig(out_prefix.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(out_prefix.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def write_audit_assets(
    ref: pd.DataFrame,
    pallial_obs: pd.DataFrame,
    subtype_obs: pd.DataFrame,
    plots_dir: Path,
    tables_dir: Path,
    source_col: str,
    fine_obs: pd.DataFrame | None = None,
) -> None:
    major_col = "major_interneuron_subtype_roi" if "major_interneuron_subtype_roi" in ref.columns else "major_interneuron_subtype"
    ref_super_bin = (
        ref.groupby(["source_supercluster", "roi_pallial_subpallial_bin"], dropna=False)
        .size()
        .reset_index(name="n_reference_cells")
    )
    ref_subtype_bin = (
        ref.groupby([major_col, "roi_pallial_subpallial_bin"], dropna=False)
        .size()
        .reset_index(name="n_reference_cells")
        .rename(columns={major_col: "major_interneuron_subtype"})
    )
    ref_subtype_total = (
        ref.groupby(major_col, dropna=False)
        .size()
        .reset_index(name="n_reference_cells")
        .rename(columns={major_col: "major_interneuron_subtype"})
    )
    ref_fine_bin = None
    ref_fine_total = None
    if "final_fine_subtype" in ref.columns:
        ref_fine_bin = (
            ref.groupby(["final_fine_subtype", "roi_pallial_subpallial_bin"], dropna=False)
            .size()
            .reset_index(name="n_reference_cells")
        )
        ref_fine_total = (
            ref.groupby("final_fine_subtype", dropna=False)
            .size()
            .reset_index(name="n_reference_cells")
        )
    query_pallial = (
        pallial_obs.groupby([source_col, "predicted.id"], dropna=False)
        .size()
        .reset_index(name="n_query_cells")
        .rename(columns={source_col: "div90_class", "predicted.id": "adult_pallial_subpallial_bin"})
    )
    query_subtype = (
        subtype_obs.groupby([source_col, "predicted.id"], dropna=False)
        .size()
        .reset_index(name="n_query_cells")
        .rename(columns={source_col: "div90_class", "predicted.id": "adult_major_interneuron_subtype"})
    )
    query_fine = None
    if fine_obs is not None:
        query_fine = (
            fine_obs.groupby([source_col, "predicted.id"], dropna=False)
            .size()
            .reset_index(name="n_query_cells")
            .rename(columns={source_col: "div90_class", "predicted.id": "adult_final_fine_subtype"})
        )
    totals = pd.DataFrame(
        [
            {"metric": "reference_cells_total", "n_cells": int(ref.shape[0])},
            {
                "metric": "reference_cells_used_for_subtype_transfer",
                "n_cells": int(ref.loc[~ref[major_col].eq(OTHER_SELECTED_REFERENCE)].shape[0]),
            },
            {
                "metric": "reference_cells_excluded_from_subtype_transfer_other_selected_reference",
                "n_cells": int(ref[major_col].eq(OTHER_SELECTED_REFERENCE).sum()),
            },
            {
                "metric": "reference_cells_used_for_final_fine_subtype_transfer",
                "n_cells": int(ref.loc[~ref["final_fine_subtype"].eq(OTHER_SELECTED_REFERENCE)].shape[0])
                if "final_fine_subtype" in ref.columns
                else 0,
            },
            {
                "metric": "reference_cells_excluded_from_final_fine_subtype_transfer_other_selected_reference",
                "n_cells": int(ref["final_fine_subtype"].eq(OTHER_SELECTED_REFERENCE).sum()) if "final_fine_subtype" in ref.columns else 0,
            },
            {"metric": "div90_cells_pallial_transfer", "n_cells": int(pallial_obs.shape[0])},
            {"metric": "div90_cells_major_subtype_transfer", "n_cells": int(subtype_obs.shape[0])},
            {
                "metric": "div90_cells_final_fine_subtype_transfer",
                "n_cells": int(fine_obs.shape[0]) if fine_obs is not None else 0,
            },
            {"metric": "river_pallial_edges_total", "n_cells": int(query_pallial["n_query_cells"].sum())},
            {"metric": "river_major_subtype_edges_total", "n_cells": int(query_subtype["n_query_cells"].sum())},
            {
                "metric": "river_final_fine_subtype_edges_total",
                "n_cells": int(query_fine["n_query_cells"].sum()) if query_fine is not None else 0,
            },
        ]
    )

    ref_super_bin.to_csv(tables_dir / "audit_reference_source_supercluster_by_pallial_subpallial_bin.tsv", sep="\t", index=False)
    ref_subtype_bin.to_csv(tables_dir / "audit_reference_major_subtype_by_pallial_subpallial_bin.tsv", sep="\t", index=False)
    ref_subtype_total.to_csv(tables_dir / "audit_reference_major_subtype_totals.tsv", sep="\t", index=False)
    query_pallial.to_csv(tables_dir / "audit_div90_class_to_pallial_subpallial_counts.tsv", sep="\t", index=False)
    query_subtype.to_csv(tables_dir / "audit_div90_class_to_major_subtype_counts.tsv", sep="\t", index=False)
    if ref_fine_bin is not None:
        ref_fine_bin.to_csv(tables_dir / "audit_reference_final_fine_subtype_by_pallial_subpallial_bin.tsv", sep="\t", index=False)
        ref_fine_total.to_csv(tables_dir / "audit_reference_final_fine_subtype_totals.tsv", sep="\t", index=False)
    if query_fine is not None:
        query_fine.to_csv(tables_dir / "audit_div90_class_to_final_fine_subtype_counts.tsv", sep="\t", index=False)
    totals.to_csv(tables_dir / "audit_accounting_totals.tsv", sep="\t", index=False)

    write_stacked_bar(
        ref_super_bin,
        "source_supercluster",
        "roi_pallial_subpallial_bin",
        "n_reference_cells",
        plots_dir / "audit_reference_supercluster_by_pallial_subpallial_bin",
        "Reference accounting: selected Siletti superclusters by pallial/subpallial bin",
    )
    write_stacked_bar(
        ref_subtype_bin,
        "major_interneuron_subtype",
        "roi_pallial_subpallial_bin",
        "n_reference_cells",
        plots_dir / "audit_reference_major_subtype_by_pallial_subpallial_bin",
        "Reference accounting: major subtype labels by pallial/subpallial bin",
    )
    write_simple_bar(
        ref_subtype_total,
        "major_interneuron_subtype",
        "n_reference_cells",
        plots_dir / "audit_reference_major_subtype_totals_including_excluded",
        "Reference accounting: major subtype labels including excluded pool",
    )
    if ref_fine_bin is not None:
        write_stacked_bar(
            ref_fine_bin,
            "final_fine_subtype",
            "roi_pallial_subpallial_bin",
            "n_reference_cells",
            plots_dir / "audit_reference_final_fine_subtype_by_pallial_subpallial_bin",
            "Reference accounting: fine final subtype labels by pallial/subpallial bin",
        )
        write_simple_bar(
            ref_fine_total,
            "final_fine_subtype",
            "n_reference_cells",
            plots_dir / "audit_reference_final_fine_subtype_totals_including_excluded",
            "Reference accounting: fine final subtype labels including excluded pool",
        )
    write_stacked_bar(
        query_pallial,
        "div90_class",
        "adult_pallial_subpallial_bin",
        "n_query_cells",
        plots_dir / "audit_div90_class_to_pallial_subpallial_counts",
        "DIV90 accounting: pallial/subpallial assignments by class",
    )
    write_stacked_bar(
        query_subtype,
        "div90_class",
        "adult_major_interneuron_subtype",
        "n_query_cells",
        plots_dir / "audit_div90_class_to_major_subtype_counts",
        "DIV90 accounting: major subtype assignments by class",
    )
    if query_fine is not None:
        write_stacked_bar(
            query_fine,
            "div90_class",
            "adult_final_fine_subtype",
            "n_query_cells",
            plots_dir / "audit_div90_class_to_final_fine_subtype_counts",
            "DIV90 accounting: fine final subtype assignments by class",
        )
    write_score_histograms(pallial_obs, subtype_obs, plots_dir / "audit_transfer_score_distributions", fine_obs=fine_obs)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pallial-transfer-dir", required=True, type=Path)
    parser.add_argument("--subtype-transfer-dir", required=True, type=Path)
    parser.add_argument("--fine-transfer-dir", type=Path, default=None)
    parser.add_argument("--bridge-dir", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--source-class-col", default="div90_broad_class")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    plots_dir = args.outdir / "plots"
    tables_dir = args.outdir / "tables"
    reports_dir = args.outdir / "reports"
    for directory in [plots_dir, tables_dir, reports_dir]:
        directory.mkdir(parents=True, exist_ok=True)

    pallial_obs = read_transfer_obs(args.pallial_transfer_dir)
    subtype_obs = read_transfer_obs(args.subtype_transfer_dir)
    fine_obs = read_transfer_obs(args.fine_transfer_dir) if args.fine_transfer_dir is not None else None
    ref = read_bridge_reference(args.bridge_dir)

    pallial_edges = write_one(
        pallial_obs,
        plots_dir / "river_div90_class_to_adult_pallial_subpallial_bin",
        args.source_class_col,
        PALLIAL_ORDER,
        "Adult Siletti ROI bin",
        "DIV90 classes mapped to adult Siletti pallial/subpallial ROI bin",
    )
    subtype_edges = write_one(
        subtype_obs,
        plots_dir / "river_div90_class_to_adult_major_interneuron_subtypes",
        args.source_class_col,
        MAJOR_SUBTYPE_ORDER,
        "Adult Siletti major subtype",
        "DIV90 classes mapped to adult Siletti major subtypes",
    )
    fine_edges = None
    if fine_obs is not None:
        fine_edges = write_one(
            fine_obs,
            plots_dir / "river_div90_class_to_adult_final_fine_subtypes",
            args.source_class_col,
            FINE_SUBTYPE_ORDER,
            "Adult Siletti fine subtype",
            "DIV90 classes mapped to adult Siletti fine final subtypes",
        )
        draw_combined_rivers(
            pallial_edges,
            subtype_edges,
            fine_edges,
            plots_dir / "river_div90_class_to_adult_combined_1x3_pallial_major_fine",
        )
    pallial_edges.to_csv(tables_dir / "river_div90_class_to_adult_pallial_subpallial_bin_edges.tsv", sep="\t", index=False)
    subtype_edges.to_csv(tables_dir / "river_div90_class_to_adult_major_interneuron_subtypes_edges.tsv", sep="\t", index=False)
    if fine_edges is not None:
        fine_edges.to_csv(tables_dir / "river_div90_class_to_adult_final_fine_subtypes_edges.tsv", sep="\t", index=False)

    pallial_obs.to_csv(tables_dir / "div90_query_with_pallial_subpallial_assignments.tsv.gz", sep="\t", index=False, compression="gzip")
    subtype_obs.to_csv(tables_dir / "div90_query_with_major_interneuron_subtype_assignments.tsv.gz", sep="\t", index=False, compression="gzip")
    if fine_obs is not None:
        fine_obs.to_csv(tables_dir / "div90_query_with_final_fine_subtype_assignments.tsv.gz", sep="\t", index=False, compression="gzip")
    write_audit_assets(ref, pallial_obs, subtype_obs, plots_dir, tables_dir, args.source_class_col, fine_obs=fine_obs)

    config = {
        "bridge_dir": str(args.bridge_dir),
        "pallial_transfer_dir": str(args.pallial_transfer_dir),
        "subtype_transfer_dir": str(args.subtype_transfer_dir),
        "fine_transfer_dir": str(args.fine_transfer_dir) if args.fine_transfer_dir is not None else None,
        "outdir": str(args.outdir),
        "source_class_col": args.source_class_col,
        "plots": [
            "river_div90_class_to_adult_pallial_subpallial_bin.png/pdf",
            "river_div90_class_to_adult_major_interneuron_subtypes.png/pdf",
            "river_div90_class_to_adult_final_fine_subtypes.png/pdf",
            "river_div90_class_to_adult_combined_1x3_pallial_major_fine.png/pdf",
        ],
        "audit_plots": [
            "audit_reference_supercluster_by_pallial_subpallial_bin.png/pdf",
            "audit_reference_major_subtype_by_pallial_subpallial_bin.png/pdf",
            "audit_reference_major_subtype_totals_including_excluded.png/pdf",
            "audit_reference_final_fine_subtype_by_pallial_subpallial_bin.png/pdf",
            "audit_reference_final_fine_subtype_totals_including_excluded.png/pdf",
            "audit_div90_class_to_pallial_subpallial_counts.png/pdf",
            "audit_div90_class_to_major_subtype_counts.png/pdf",
            "audit_div90_class_to_final_fine_subtype_counts.png/pdf",
            "audit_transfer_score_distributions.png/pdf",
        ],
    }
    (tables_dir / "plot_config.json").write_text(json.dumps(config, indent=2, sort_keys=True) + "\n")
    (reports_dir / "README_final_candidate_restricted_rivers.md").write_text(
        "# Final-Candidate Restricted Siletti River Plots\n\n"
        "Reference scope: MGE, CGE, LAMP5-LHX6/chandelier, medium spiny neuron, "
        "eccentric medium spiny neuron, and Splatter restricted to subpallial "
        "cholinergic/CHAT-like interneuron rows.\n\n"
        "Pallial/subpallial rule: CerebralCortex + Hippocampus = Pallial/cortical; "
        "all other primary ROI groups = Subpallial.\n\n"
        "Important caveat: this pallial/subpallial bin is based on adult Siletti "
        "primary ROI group, not developmental birthplace. Under this rule, many "
        "adult MGE/CGE/LAMP5-LHX6/chandelier interneurons are Pallial/cortical "
        "because they are sampled from CerebralCortex or Hippocampus. MSN, "
        "eccentric MSN, and Splatter-CHAT/cholinergic rows are mostly "
        "Subpallial.\n\n"
        "Fine subtype options in this run: cortical PV basket, cortical PV "
        "chandelier, cortical SST LRP, cortical SST Mt, cortical SST nMt, "
        "subpallial PV, subpallial SST LRP, subpallial SST, subpallial "
        "cholinergic, pallial/cortical MSN, subpallial MSN, and subpallial "
        "eccentric MSN. `Other selected reference` is retained for reference "
        "accounting and excluded from the fine-subtype transfer. Cortical PV+ "
        "Chandelier neurons are present in the adult reference; a zero label in "
        "the river means no DIV90 cells were assigned to that class in this "
        "transfer.\n\n"
        "Main river plots:\n\n"
        "- `plots/river_div90_class_to_adult_pallial_subpallial_bin.png`\n"
        "- `plots/river_div90_class_to_adult_major_interneuron_subtypes.png`\n"
        "- `plots/river_div90_class_to_adult_final_fine_subtypes.png`\n"
        "- `plots/river_div90_class_to_adult_combined_1x3_pallial_major_fine.png`\n"
    )

    print(json.dumps(config, indent=2, sort_keys=True))
    for path in sorted(args.outdir.rglob("*")):
        if path.is_file():
            print(path)


if __name__ == "__main__":
    main()
