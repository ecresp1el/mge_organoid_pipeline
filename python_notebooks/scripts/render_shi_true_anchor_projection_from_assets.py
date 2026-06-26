#!/usr/bin/env python3
"""Render true Seurat anchor projection plots from saved anchor assets."""

from __future__ import annotations

import argparse
import re
import shutil
import subprocess
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.collections import LineCollection


PROJECT_ROOT = Path("/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder")
DEFAULT_RUN_DIR = (
    PROJECT_ROOT
    / "results"
    / "cross_study_shi_seurat_label_transfer"
    / "cross_study_shi_seurat_anchor_projection_v1"
)
DEFAULT_PREDICTION_TABLE = (
    PROJECT_ROOT
    / "results"
    / "cross_study_shi_seurat_label_transfer"
    / "cross_study_shi_seurat_label_transfer_v1"
    / "tables"
    / "cross_study_shi_seurat_label_transfer_obs.tsv.gz"
)
DEFAULT_STUDIES = ("varela_div30", "varela_div90")
BACKGROUND_GREY = "#d0d0d0"
QUERY_ANCHOR = "#222222"
LINE_GREY = "#555555"
EXPORT_DPI = 600
CLUSTER_CMAP = "tab20"

MAJOR_COLORS = {
    "MGE": "#16697a",
    "LGE": "#4895ef",
    "CGE": "#52b788",
    "progenitor": "#f4a261",
    "Excitatory IPC": "#e76f51",
    "Excitatory neuron": "#d62828",
    "Thalamic neurons": "#7b2cbf",
    "Microglia": "#6c757d",
    "OPC": "#a7c957",
    "Endothelial": "#2a9d8f",
}

GW_COLORS = {
    "GW09": "#231611",
    "GW12": "#3F1C6A",
    "GW13": "#A02E6B",
    "GW16": "#EB5840",
    "GW18": "#FCC031",
}

LABEL_TEXT_OFFSETS = {
    "MGE": (-0.04, 0.00),
    "CGE": (-0.05, -0.05),
    "LGE": (0.02, -0.08),
    "progenitor": (0.02, 0.08),
    "Excitatory IPC": (0.13, 0.03),
    "Excitatory neuron": (0.14, -0.05),
    "Thalamic neurons": (0.12, 0.08),
    "OPC": (0.04, 0.11),
    "Microglia": (0.00, 0.20),
    "Endothelial": (0.02, 0.28),
    "GW09": (-0.06, -0.02),
    "GW12": (0.06, 0.04),
    "GW13": (0.10, -0.02),
    "GW16": (0.08, 0.10),
    "GW18": (0.00, 0.14),
}


def configure_matplotlib() -> None:
    mpl.rcParams.update(
        {
            "font.family": "Arial",
            "font.sans-serif": ["Arial", "Nimbus Sans", "Liberation Sans", "DejaVu Sans"],
            "svg.fonttype": "none",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "axes.linewidth": 0.8,
            "figure.dpi": 150,
            "savefig.dpi": EXPORT_DPI,
        }
    )


def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", low_memory=False)


def normalize_bool(series: pd.Series) -> pd.Series:
    return series.astype(str).str.lower().isin({"true", "t", "1", "yes"})


def canonical_gw_label(value: object) -> str:
    match = re.search(r"GW\s*0*([0-9]+)", str(value), flags=re.IGNORECASE)
    if not match:
        return "unlabeled"
    return f"GW{int(match.group(1)):02d}"


def scaled_xy(df: pd.DataFrame, x_col: str = "coord_1", y_col: str = "coord_2") -> np.ndarray:
    coords = df[[x_col, y_col]].astype(float).to_numpy()
    out = coords.copy()
    finite = np.isfinite(out).all(axis=1)
    if not finite.any():
        return np.zeros_like(out)
    center = np.nanmedian(out[finite], axis=0)
    out = out - center
    span = np.nanmax(out[finite], axis=0) - np.nanmin(out[finite], axis=0)
    scale = np.nanmax(span)
    if not np.isfinite(scale) or scale <= 0:
        scale = 1.0
    out = out / scale
    return out


def normalize_reference_ids(links: pd.DataFrame, ref_ids: set[str]) -> pd.DataFrame:
    links = links.copy()
    original = links["reference_cell_id"].astype(str)
    stripped = original.str.replace("_reference$", "", regex=True)
    links["reference_cell_id_original"] = original
    links["reference_cell_id_plot"] = np.where(original.isin(ref_ids), original, stripped)
    return links


def load_prediction_table(path: Path) -> pd.DataFrame:
    cols = [
        "cell_id",
        "study_id",
        "shi_seurat_full_predicted_shi_label",
        "shi_seurat_full_predicted_shi_week_label",
        "shi_seurat_full_prediction_score",
        "shi_seurat_full_week_prediction_score",
    ]
    return pd.read_csv(path, sep="\t", usecols=cols, dtype={"cell_id": str, "study_id": str}, low_memory=False)


def attach_query_predictions(query: pd.DataFrame, predictions: pd.DataFrame, study_id: str) -> pd.DataFrame:
    pred = predictions.loc[predictions["study_id"].astype(str) == study_id].copy()
    pred = pred.drop_duplicates(subset=["cell_id"])
    keep_cols = [
        "cell_id",
        "shi_seurat_full_predicted_shi_label",
        "shi_seurat_full_predicted_shi_week_label",
        "shi_seurat_full_prediction_score",
        "shi_seurat_full_week_prediction_score",
    ]
    out = query.merge(pred[keep_cols], on="cell_id", how="left")
    out["shi_seurat_full_predicted_shi_week_label"] = out["shi_seurat_full_predicted_shi_week_label"].map(canonical_gw_label)
    return out


def load_study(run_dir: Path, study_id: str, predictions: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    study_dir = run_dir / study_id
    tables = study_dir / "tables"
    links = read_tsv(tables / f"{study_id}_shi_full_anchor_pairs.tsv.gz")
    query = read_tsv(tables / f"{study_id}_query_coordinates_for_anchor_plot.tsv.gz")
    reference = read_tsv(tables / f"{study_id}_shi_reference_coordinates_for_anchor_plot.tsv.gz")
    query["cell_id"] = query["cell_id"].astype(str)
    reference["cell_id"] = reference["cell_id"].astype(str)
    if "shi_week_label" in reference.columns:
        reference["shi_week_label"] = reference["shi_week_label"].map(canonical_gw_label)

    if "plot_include" in links.columns:
        links = links.loc[normalize_bool(links["plot_include"])].copy()
    if "plot_include" in query.columns:
        query = query.loc[normalize_bool(query["plot_include"])].copy()

    if study_id == "varela_div90":
        query = query.copy()
        query["coord_2"] = -1.0 * query["coord_2"].astype(float)

    query = attach_query_predictions(query, predictions, study_id)
    links = normalize_reference_ids(links, set(reference["cell_id"].astype(str)))
    return links, query, reference


def join_links(links: pd.DataFrame, query: pd.DataFrame, reference: pd.DataFrame) -> pd.DataFrame:
    q = query[["cell_id", "coord_1", "coord_2"]].rename(
        columns={"cell_id": "query_cell_id", "coord_1": "query_coord_1", "coord_2": "query_coord_2"}
    )
    r_cols = ["cell_id", "coord_1", "coord_2"]
    for col in ("shi_label", "shi_week_label"):
        if col in reference.columns:
            r_cols.append(col)
    r = reference[r_cols].rename(
        columns={
            "cell_id": "reference_cell_id_plot",
            "coord_1": "reference_coord_1",
            "coord_2": "reference_coord_2",
        }
    )
    merged = links.merge(q, on="query_cell_id", how="inner").merge(r, on="reference_cell_id_plot", how="inner")
    return merged


def diagnostics_for(study_id: str, links: pd.DataFrame, query: pd.DataFrame, reference: pd.DataFrame, merged: pd.DataFrame) -> dict:
    original = links["reference_cell_id_original"].astype(str)
    stripped = links["reference_cell_id_plot"].astype(str)
    ref_ids = set(reference["cell_id"].astype(str))
    query_ids = set(query["cell_id"].astype(str))
    return {
        "study_id": study_id,
        "n_query_cells_plotted_no_downsampling": len(query),
        "n_reference_cells_plotted_no_downsampling": len(reference),
        "n_anchor_links_after_link_flag_filter": len(links),
        "n_query_cells_with_wta_major_label": int(query["shi_seurat_full_predicted_shi_label"].notna().sum()),
        "n_query_cells_with_wta_gw_label": int(query["shi_seurat_full_predicted_shi_week_label"].notna().sum()),
        "n_reference_ids_exact_match_before_suffix_fix": int(original.isin(ref_ids).sum()),
        "n_reference_ids_match_after_suffix_fix": int(stripped.isin(ref_ids).sum()),
        "n_query_ids_match_coordinates": int(links["query_cell_id"].astype(str).isin(query_ids).sum()),
        "n_links_after_coordinate_merge": len(merged),
        "n_anchor_links_removed_by_coordinate_visual_filter": len(links) - len(merged),
        "n_unique_query_anchor_cells_after_merge": int(merged["query_cell_id"].nunique()),
        "n_unique_reference_anchor_cells_after_merge": int(merged["reference_cell_id_plot"].nunique()),
    }


def colors_for(values: pd.Series, mode: str) -> tuple[list[str], dict[str, str]]:
    palette = MAJOR_COLORS if mode == "major_class" else GW_COLORS
    labels = values.fillna("unlabeled").astype(str)
    colors = [palette.get(v, "#777777") for v in labels]
    used = {k: palette[k] for k in palette if k in set(labels)}
    if any(v not in palette for v in labels):
        used["unlabeled"] = "#777777"
    return colors, used


def label_col_for_mode(mode: str, source: str) -> str:
    if source == "query":
        return "shi_seurat_full_predicted_shi_label" if mode == "major_class" else "shi_seurat_full_predicted_shi_week_label"
    return "shi_label" if mode == "major_class" else "shi_week_label"


def cluster_colors(values: pd.Series) -> tuple[list[str], dict[str, tuple[float, float, float, float]]]:
    labels = values.fillna("NA").astype(str)
    unique_labels = sorted(labels.unique(), key=lambda x: (not x.lstrip("-").isdigit(), int(x) if x.lstrip("-").isdigit() else x))
    cmap = mpl.colormaps[CLUSTER_CMAP]
    palette = {label: cmap(i % cmap.N) for i, label in enumerate(unique_labels)}
    return [palette[v] for v in labels], palette


def add_category_centroid_labels(ax: plt.Axes, df: pd.DataFrame, category_col: str, fontsize: float = 6.8) -> None:
    if category_col not in df.columns:
        return
    tmp = df[["x_plot", "y_plot", category_col]].copy()
    tmp[category_col] = tmp[category_col].fillna("unlabeled").astype(str)
    for label, group in tmp.groupby(category_col, sort=False):
        if not label or label.lower() in {"nan", "unlabeled"} or group.empty:
            continue
        dx, dy = LABEL_TEXT_OFFSETS.get(label, (0.0, 0.0))
        ax.text(
            float(group["x_plot"].median()) + dx,
            float(group["y_plot"].median()) + dy,
            label,
            ha="center",
            va="center",
            fontsize=fontsize,
            color="black",
            zorder=7,
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.78, "pad": 0.7},
        )


def add_cluster_centroid_labels(ax: plt.Axes, df: pd.DataFrame) -> None:
    if "query_cluster" not in df.columns:
        return
    tmp = df[["x_plot", "y_plot", "query_cluster"]].copy()
    tmp["query_cluster"] = tmp["query_cluster"].fillna("NA").astype(str)
    for cluster, group in tmp.groupby("query_cluster", sort=False):
        if group.empty:
            continue
        ax.text(
            float(group["x_plot"].median()),
            float(group["y_plot"].median()),
            cluster,
            ha="center",
            va="center",
            fontsize=7,
            color="black",
            zorder=6,
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.72, "pad": 0.6},
        )


def draw_reference_subpanel(
    ax: plt.Axes,
    study_id: str,
    query: pd.DataFrame,
    reference: pd.DataFrame,
    mode: str,
    title: str,
) -> dict[str, str]:
    query_xy = scaled_xy(query)
    ref_xy = scaled_xy(reference)
    query_plot = query.copy()
    ref_plot = reference.copy()
    query_plot["x_plot"] = query_xy[:, 0] - 0.7
    query_plot["y_plot"] = query_xy[:, 1]
    ref_plot["x_plot"] = ref_xy[:, 0] + 0.7
    ref_plot["y_plot"] = ref_xy[:, 1]

    query_label_col = label_col_for_mode(mode, "query")
    ref_label_col = label_col_for_mode(mode, "reference")
    query_colors, _ = colors_for(query_plot.get(query_label_col, pd.Series(index=query_plot.index, dtype=object)), mode)
    ref_colors, legend_colors = colors_for(ref_plot.get(ref_label_col, pd.Series(index=ref_plot.index, dtype=object)), mode)

    ax.scatter(query_plot["x_plot"], query_plot["y_plot"], s=2.2, c=query_colors, linewidths=0, rasterized=True, zorder=2)
    ax.scatter(ref_plot["x_plot"], ref_plot["y_plot"], s=2.2, c=ref_colors, linewidths=0, rasterized=True, zorder=2)
    add_cluster_centroid_labels(ax, query_plot)
    add_category_centroid_labels(ax, ref_plot, ref_label_col, fontsize=6.6 if mode == "major_class" else 7.0)

    ax.text(-0.7, 0.68, f"{title} WTA labels + clusters", ha="center", va="bottom", fontsize=8.5)
    ax.text(
        0.7,
        0.68,
        "Shi cell types" if mode == "major_class" else "Shi GW/stage",
        ha="center",
        va="bottom",
        fontsize=8.5,
    )
    ax.set_xlim(-1.32, 1.32)
    ax.set_ylim(-0.70, 0.78)
    ax.set_aspect("equal")
    ax.axis("off")
    return legend_colors


def draw_study_panel(
    ax: plt.Axes,
    study_id: str,
    query: pd.DataFrame,
    reference: pd.DataFrame,
    merged: pd.DataFrame,
    mode: str,
    title: str,
) -> dict[str, str]:
    query_xy = scaled_xy(query)
    ref_xy = scaled_xy(reference)
    query_plot = query.copy()
    ref_plot = reference.copy()
    query_plot["x_plot"] = query_xy[:, 0] - 0.7
    query_plot["y_plot"] = query_xy[:, 1]
    ref_plot["x_plot"] = ref_xy[:, 0] + 0.7
    ref_plot["y_plot"] = ref_xy[:, 1]

    query_lookup = query_plot.set_index("cell_id")[["x_plot", "y_plot"]]
    ref_lookup = ref_plot.set_index("cell_id")[["x_plot", "y_plot"]]
    line_df = merged.join(query_lookup, on="query_cell_id", rsuffix="_query_plot")
    line_df = line_df.rename(columns={"x_plot": "query_x_plot", "y_plot": "query_y_plot"})
    line_df = line_df.join(ref_lookup, on="reference_cell_id_plot")
    line_df = line_df.rename(columns={"x_plot": "reference_x_plot", "y_plot": "reference_y_plot"})
    line_df = line_df.dropna(subset=["query_x_plot", "query_y_plot", "reference_x_plot", "reference_y_plot"])

    query_label_col = label_col_for_mode(mode, "query")
    ref_label_col = label_col_for_mode(mode, "reference")
    query_colors, _ = colors_for(query_plot.get(query_label_col, pd.Series(index=query_plot.index, dtype=object)), mode)
    ref_colors, legend_colors = colors_for(ref_plot.get(ref_label_col, pd.Series(index=ref_plot.index, dtype=object)), mode)
    line_values = line_df.get(ref_label_col, pd.Series(index=line_df.index, dtype=object))
    line_colors, _ = colors_for(line_values, mode)
    scores = pd.to_numeric(line_df.get("score", 0.5), errors="coerce").fillna(0.5).to_numpy()
    widths = 0.15 + 1.6 * np.clip(scores, 0, 1)
    alphas = 0.08 + 0.28 * np.clip(scores, 0, 1)
    segments = np.stack(
        [
            line_df[["query_x_plot", "query_y_plot"]].to_numpy(),
            line_df[["reference_x_plot", "reference_y_plot"]].to_numpy(),
        ],
        axis=1,
    )
    rgba = [mpl.colors.to_rgba(color, alpha=float(alpha)) for color, alpha in zip(line_colors, alphas)]

    if len(segments):
        ax.add_collection(LineCollection(segments, colors=rgba, linewidths=widths, zorder=1))
    ax.scatter(query_plot["x_plot"], query_plot["y_plot"], s=2.0, c=query_colors, alpha=0.45, linewidths=0, rasterized=True, zorder=2)
    ax.scatter(ref_plot["x_plot"], ref_plot["y_plot"], s=2.0, c=ref_colors, alpha=0.45, linewidths=0, rasterized=True, zorder=2)

    anchor_query = query_plot.loc[query_plot["cell_id"].isin(set(line_df["query_cell_id"]))]
    anchor_ref = ref_plot.loc[ref_plot["cell_id"].isin(set(line_df["reference_cell_id_plot"]))]
    anchor_query_colors, _ = colors_for(anchor_query.get(query_label_col, pd.Series(index=anchor_query.index, dtype=object)), mode)
    anchor_ref_colors, _ = colors_for(anchor_ref.get(ref_label_col, pd.Series(index=anchor_ref.index, dtype=object)), mode)
    ax.scatter(anchor_query["x_plot"], anchor_query["y_plot"], s=10.0, c=anchor_query_colors, edgecolors="#111111", linewidths=0.15, zorder=4)
    ax.scatter(anchor_ref["x_plot"], anchor_ref["y_plot"], s=10.0, c=anchor_ref_colors, edgecolors="#111111", linewidths=0.12, zorder=5)

    title_box = {"facecolor": "white", "edgecolor": "none", "alpha": 0.82, "pad": 1.5}
    ax.text(-0.7, 0.70, f"{title} query UMAP", ha="center", va="bottom", fontsize=9, zorder=10, bbox=title_box)
    ax.text(0.7, 0.70, "Shi et al. reference UMAP", ha="center", va="bottom", fontsize=9, zorder=10, bbox=title_box)
    ax.text(
        0,
        -0.62,
        f"{len(line_df):,} true anchor links after visual filters",
        ha="center",
        va="top",
        fontsize=8,
        color="#333333",
    )
    ax.set_xlim(-1.32, 1.32)
    ax.set_ylim(-0.70, 0.78)
    ax.set_aspect("equal")
    ax.axis("off")
    return legend_colors


def save_figure(fig: plt.Figure, out_prefix: Path) -> None:
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "pdf", "svg"):
        fig.savefig(out_prefix.with_suffix(f".{ext}"), dpi=EXPORT_DPI, bbox_inches="tight")


def render(run_dir: Path, out_dir: Path, studies: list[str], prediction_table: Path) -> pd.DataFrame:
    predictions = load_prediction_table(prediction_table)
    all_data = {}
    diagnostics = []
    for study_id in studies:
        links, query, reference = load_study(run_dir, study_id, predictions)
        merged = join_links(links, query, reference)
        diagnostics.append(diagnostics_for(study_id, links, query, reference, merged))
        all_data[study_id] = (links, query, reference, merged)

    diag = pd.DataFrame(diagnostics)
    diag.to_csv(out_dir / "tables" / "shi_true_anchor_projection_render_diagnostics.tsv", sep="\t", index=False)

    title_map = {"varela_div30": "DIV30", "varela_div90": "DIV90"}
    for mode in ("major_class", "gw_stage"):
        height_ratios = []
        for _ in studies:
            height_ratios.extend([3.35, 2.05])
        fig, axes = plt.subplots(
            len(studies) * 2,
            1,
            figsize=(9.2, 5.85 * len(studies)),
            constrained_layout=True,
            gridspec_kw={"height_ratios": height_ratios},
        )
        axes = np.atleast_1d(axes)
        legend_colors_all: dict[str, str] = {}
        for idx, study_id in enumerate(studies):
            ax_projection = axes[idx * 2]
            ax_reference = axes[idx * 2 + 1]
            _, query, reference, merged = all_data[study_id]
            legend_colors_all.update(
                draw_study_panel(ax_projection, study_id, query, reference, merged, mode, title_map.get(study_id, study_id))
            )
            legend_colors_all.update(
                draw_reference_subpanel(ax_reference, study_id, query, reference, mode, title_map.get(study_id, study_id))
            )
        handles = [
            mpl.lines.Line2D([0], [0], marker="o", color="none", markerfacecolor=color, markersize=5, label=label)
            for label, color in legend_colors_all.items()
        ]
        handles.insert(
            0,
            mpl.lines.Line2D([0], [0], marker="o", color="#111111", markerfacecolor="white", markersize=5, label="anchor-linked cell"),
        )
        fig.legend(
            handles=handles,
            loc="upper center",
            bbox_to_anchor=(0.5, -0.01),
            ncol=min(6, max(1, len(handles))),
            frameon=False,
            fontsize=7,
        )
        save_figure(fig, out_dir / "figures" / f"div30_div90_shi_true_anchor_projection_{mode}_side_by_side_grid")
        plt.close(fig)

        for study_id in studies:
            fig, axes = plt.subplots(
                2,
                1,
                figsize=(9.2, 5.85),
                constrained_layout=True,
                gridspec_kw={"height_ratios": [3.35, 2.05]},
            )
            _, query, reference, merged = all_data[study_id]
            legend_colors = draw_study_panel(axes[0], study_id, query, reference, merged, mode, title_map.get(study_id, study_id))
            legend_colors.update(draw_reference_subpanel(axes[1], study_id, query, reference, mode, title_map.get(study_id, study_id)))
            handles = [
                mpl.lines.Line2D([0], [0], marker="o", color="none", markerfacecolor=color, markersize=5, label=label)
                for label, color in legend_colors.items()
            ]
            handles.insert(
                0,
                mpl.lines.Line2D([0], [0], marker="o", color="#111111", markerfacecolor="white", markersize=5, label="anchor-linked cell"),
            )
            fig.legend(
                handles=handles,
                loc="upper center",
                bbox_to_anchor=(0.5, -0.01),
                ncol=min(6, max(1, len(handles))),
                frameon=False,
                fontsize=7,
            )
            save_figure(fig, out_dir / "figures" / f"{study_id}_shi_true_anchor_projection_{mode}_side_by_side")
            plt.close(fig)

    return diag


def write_readme(out_dir: Path, run_dir: Path, prediction_table: Path, diag: pd.DataFrame) -> None:
    diag_text = diag.to_csv(sep="\t", index=False)
    lines = [
        "# Shi True Anchor Projection Plots",
        "",
        "These plots are rendered from saved Seurat anchor assets only. They do not rerun FindTransferAnchors.",
        "",
        "Each panel shows the query UMAP next to the corresponding Shi et al. reference UMAP.",
        "Lines connect all saved true Seurat anchor links from query cells to Shi reference cells after the visualization filters.",
        "Line width and opacity scale with the Seurat anchor score. DIV90 uses the same visualization-only stressed-cell removal and vertical UMAP flip used in finalized Shi plots.",
        "Query UMAPs are colored by the saved winner-take-all Seurat TransferData prediction labels: major class for the major-class view and GW/stage for the GW view.",
        "Shi reference UMAPs and anchor-link colors use the same palette as the query winner-take-all labels for the selected view.",
        "Each study also has a lower reference subpanel: the query UMAP is colored by winner-take-all label with query cluster numbers overlaid, and the Shi UMAP is colored and directly labeled by the plotted Shi label mode.",
        "UMAP backgrounds and reference subpanels use all cells present in the saved coordinate tables; no point downsampling is applied.",
        "The term filter in the diagnostics refers to visualization-only coordinate/cluster exclusions, not winner-take-all label assignment.",
        "",
        "The reference-cell ID suffix `_reference` was removed only for plotting joins, because the saved anchor table carries suffixed Shi IDs while the Shi coordinate table stores unsuffixed cell IDs.",
        "",
        f"Input asset directory: `{run_dir}`",
        f"Winner-take-all prediction table: `{prediction_table}`",
        "",
        "Diagnostics:",
        "",
        "```text",
        diag_text.strip(),
        "```",
        "",
    ]
    (out_dir / "README.md").write_text("\n".join(lines), encoding="utf-8")


def copy_provenance(out_dir: Path) -> None:
    prov = out_dir / "provenance"
    prov.mkdir(parents=True, exist_ok=True)
    script_path = Path(__file__)
    shutil.copy2(script_path, prov / script_path.name)
    try:
        status = subprocess.run(["git", "status", "--short"], check=False, capture_output=True, text=True)
        (prov / "git_status_short.txt").write_text(status.stdout, encoding="utf-8")
    except Exception as exc:
        (prov / "git_status_short.txt").write_text(f"git status failed: {exc}\n", encoding="utf-8")


def write_checksums(out_dir: Path) -> None:
    paths = sorted(p for p in out_dir.rglob("*") if p.is_file() and p.name != "sha256sums.txt")
    with (out_dir / "provenance" / "sha256sums.txt").open("w", encoding="utf-8") as handle:
        for path in paths:
            try:
                digest = subprocess.check_output(["sha256sum", str(path)], text=True).strip()
            except subprocess.CalledProcessError:
                continue
            handle.write(digest + "\n")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", type=Path, default=DEFAULT_RUN_DIR)
    parser.add_argument("--prediction-table", type=Path, default=DEFAULT_PREDICTION_TABLE)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_RUN_DIR / "true_anchor_projection_plots")
    parser.add_argument("--studies", nargs="+", default=list(DEFAULT_STUDIES))
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    configure_matplotlib()
    for subdir in ("figures", "tables", "provenance"):
        (args.out_dir / subdir).mkdir(parents=True, exist_ok=True)
    diag = render(args.run_dir, args.out_dir, args.studies, args.prediction_table)
    write_readme(args.out_dir, args.run_dir, args.prediction_table, diag)
    copy_provenance(args.out_dir)
    write_checksums(args.out_dir)
    print(diag.to_csv(sep="\t", index=False))
    print(f"Rendered true anchor projection plots to {args.out_dir}")


if __name__ == "__main__":
    main()
