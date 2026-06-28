#!/usr/bin/env python3
"""Static routing diagnostic from Varela DIV30/DIV90 cells to Shi predictions."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

from mge_organoid_python.cross_study_shi_prediction_plots import (
    FIGURE_EXPORT_DPI,
    LABEL_COLORS,
    SHI_LABEL_ORDER,
    add_plot_coordinates,
    apply_internal_umap_plot_filters,
    load_combined_table,
)


PROJECT_ROOT_DEFAULT = Path("/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder")
RUN_LABEL_DEFAULT = "cross_study_shi_seurat_label_transfer_v1"
OUTDIR_DEFAULT = (
    PROJECT_ROOT_DEFAULT
    / "results"
    / "cross_study_shi_seurat_label_transfer"
    / RUN_LABEL_DEFAULT
    / "plots"
    / "validation"
    / "shi_query_projection_routing_static"
)

STUDIES = ["varela_div30", "varela_div90"]
STUDY_TITLES = {
    "varela_div30": "This Study, DIV 30",
    "varela_div90": "This Study, DIV 90",
}

GW_ORDER = ["GW09", "GW12", "GW13", "GW16", "GW18"]
GW_COLORS = {
    "GW09": "#231611",
    "GW12": "#3F1C6A",
    "GW13": "#A02E6B",
    "GW16": "#EB5840",
    "GW18": "#FCC031",
}

DIV30_RECODE = {
    "3": ("1", "DIV30 class 1"),
    "0": ("1", "DIV30 class 1"),
    "7": ("1", "DIV30 class 1"),
    "6": ("2", "DIV30 class 2"),
    "1": ("3", "DIV30 class 3"),
    "4": ("4", "DIV30 class 4"),
    "2": ("5", "DIV30 class 5"),
}

DIV90_RECODE = {
    "0": ("3", "MGE Striatal/GP fated"),
    "1": ("7", "PV Precursors/Migrating cells/Cortical fated"),
    "2": ("2", "CRABP1+/PV Precursors"),
    "3": ("1", "SST+, NPY+ Cortical fated"),
    "4": ("8", "Pre-Astrocytes/Astrocytes"),
    "5": ("4", "LHX8+ vMGE GABAergic Striatal/GP fated 1"),
    "8": ("5", "LHX8+ vMGE GABAergic Striatal/GP fated 2"),
    "9": ("10", "Pre-OPCs/OPCs"),
    "10": ("8", "Pre-Astrocytes/Astrocytes"),
    "11": ("6", "PV Precursors"),
    "12": ("9", "Dividing cells"),
}

plt.rcParams.update(
    {
        "font.family": "Arial",
        "font.sans-serif": ["Arial", "Nimbus Sans", "Liberation Sans", "DejaVu Sans"],
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    }
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, default=PROJECT_ROOT_DEFAULT)
    parser.add_argument("--run-label", default=RUN_LABEL_DEFAULT)
    parser.add_argument("--outdir", type=Path, default=OUTDIR_DEFAULT)
    parser.add_argument("--max-cell-points-per-study", type=int, default=6500)
    parser.add_argument("--min-route-fraction", type=float, default=0.025)
    parser.add_argument("--top-routes-per-source", type=int, default=3)
    parser.add_argument("--seed", type=int, default=17)
    return parser.parse_args()


def clean_cluster(value: object) -> str:
    text = str(value)
    if text.endswith(".0"):
        text = text[:-2]
    return text


def add_source_groups(data: pd.DataFrame) -> pd.DataFrame:
    out = data.copy()
    raw_cluster = out["cluster"].map(clean_cluster)
    source_id = []
    source_name = []
    source_label = []
    for study_id, cluster in zip(out["study_id"].astype(str), raw_cluster):
        if study_id == "varela_div30":
            recoded, name = DIV30_RECODE.get(cluster, (f"raw {cluster}", f"DIV30 raw cluster {cluster}"))
            label = recoded
        elif study_id == "varela_div90":
            recoded, name = DIV90_RECODE.get(cluster, (f"raw {cluster}", f"DIV90 raw cluster {cluster}"))
            label = recoded
        else:
            recoded, name, label = cluster, f"cluster {cluster}", cluster
        source_id.append(recoded)
        source_name.append(name)
        source_label.append(label)
    out["raw_cluster"] = raw_cluster
    out["source_group_id"] = source_id
    out["source_group_name"] = source_name
    out["source_plot_label"] = source_label
    return out


def load_varela_plot_data(project_root: Path, run_label: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    data = load_combined_table(project_root, run_label)
    data, filter_summary = apply_internal_umap_plot_filters(data)
    data = add_plot_coordinates(data)
    data = data.loc[data["study_id"].astype(str).isin(STUDIES)].copy()
    finite = np.isfinite(pd.to_numeric(data["UMAP1_plot"], errors="coerce")) & np.isfinite(
        pd.to_numeric(data["UMAP2_plot"], errors="coerce")
    )
    data = data.loc[finite].copy()
    data = add_source_groups(data)
    return data, filter_summary.loc[filter_summary["study_id"].astype(str).isin(STUDIES)].copy()


def stratified_sample(data: pd.DataFrame, target_col: str, max_points: int, seed: int) -> pd.DataFrame:
    if data.shape[0] <= max_points:
        return data.copy()
    rng = np.random.default_rng(seed)
    groups = list(data.groupby(["source_group_id", target_col], observed=True, sort=False))
    per_group = max(8, int(np.ceil(max_points / max(len(groups), 1))))
    parts = []
    for _, group in groups:
        if group.shape[0] <= per_group:
            parts.append(group)
        else:
            idx = rng.choice(group.index.to_numpy(), size=per_group, replace=False)
            parts.append(group.loc[idx])
    sampled = pd.concat(parts, ignore_index=False)
    if sampled.shape[0] > max_points:
        idx = rng.choice(sampled.index.to_numpy(), size=max_points, replace=False)
        sampled = sampled.loc[idx]
    return sampled.copy()


def route_table(
    data: pd.DataFrame,
    *,
    target_col: str,
    score_col: str,
    target_order: list[str],
) -> pd.DataFrame:
    work = data.loc[data[target_col].astype(str).isin(target_order)].copy()
    work[target_col] = work[target_col].astype(str)
    work[score_col] = pd.to_numeric(work[score_col], errors="coerce")
    group_cols = [
        "study_id",
        "source_group_id",
        "source_group_name",
        "source_plot_label",
        target_col,
    ]
    routes = (
        work.groupby(group_cols, observed=True)
        .agg(
            n_cells=("cell_id", "size"),
            mean_prediction_score=(score_col, "mean"),
            source_umap1=("UMAP1_plot", "mean"),
            source_umap2=("UMAP2_plot", "mean"),
        )
        .reset_index()
    )
    totals = (
        work.groupby(["study_id", "source_group_id"], observed=True)
        .size()
        .rename("n_source_cells")
        .reset_index()
    )
    routes = routes.merge(totals, on=["study_id", "source_group_id"], how="left", validate="many_to_one")
    routes["fraction_of_source"] = routes["n_cells"] / routes["n_source_cells"]
    routes["target_order"] = routes[target_col].map({label: idx for idx, label in enumerate(target_order)})
    routes = routes.sort_values(["study_id", "source_group_id", "fraction_of_source"], ascending=[True, True, False])
    routes["route_rank_in_source"] = routes.groupby(["study_id", "source_group_id"], observed=True).cumcount() + 1
    return routes


def selected_routes(routes: pd.DataFrame, min_fraction: float, top_n: int) -> pd.DataFrame:
    keep = (routes["fraction_of_source"] >= min_fraction) | (routes["route_rank_in_source"] <= int(top_n))
    return routes.loc[keep].copy()


def target_positions(study_data: pd.DataFrame, target_order: list[str]) -> dict[str, tuple[float, float, float]]:
    x = pd.to_numeric(study_data["UMAP1_plot"], errors="coerce")
    y = pd.to_numeric(study_data["UMAP2_plot"], errors="coerce")
    xmin, xmax = float(x.min()), float(x.max())
    ymin, ymax = float(y.min()), float(y.max())
    x_pad = max((xmax - xmin) * 0.08, 1.0)
    y_pad = max((ymax - ymin) * 0.20, 2.0)
    if len(target_order) <= 6:
        xs = np.linspace(xmin + x_pad, xmax - x_pad, len(target_order))
        return {label: (float(xs[i]), ymax + y_pad, 1.35) for i, label in enumerate(target_order)}

    n_first_row = int(np.ceil(len(target_order) / 2))
    out: dict[str, tuple[float, float, float]] = {}
    for i, label in enumerate(target_order):
        row = 0 if i < n_first_row else 1
        within = i if row == 0 else i - n_first_row
        n_in_row = n_first_row if row == 0 else len(target_order) - n_first_row
        xs = np.linspace(xmin + x_pad, xmax - x_pad, n_in_row)
        out[label] = (float(xs[within]), ymax + y_pad * (1.0 + 0.85 * row), 1.35)
    return out


def set_axes_equalish(ax, data: pd.DataFrame, target_pos: dict[str, tuple[float, float, float]]) -> None:
    x = pd.to_numeric(data["UMAP1_plot"], errors="coerce")
    y = pd.to_numeric(data["UMAP2_plot"], errors="coerce")
    tx = [p[0] for p in target_pos.values()]
    ty = [p[1] for p in target_pos.values()]
    xmin, xmax = min(float(x.min()), min(tx)), max(float(x.max()), max(tx))
    ymin, ymax = min(float(y.min()), min(ty)), max(float(y.max()), max(ty))
    ax.set_xlim(xmin - 0.7, xmax + 0.7)
    ax.set_ylim(ymin - 0.7, ymax + 0.9)
    ax.set_zlim(-0.12, 1.55)


def target_label_text(label: str) -> str:
    replacements = {
        "Excitatory IPC": "Excitatory\nIPC",
        "Excitatory neuron": "Excitatory\nneuron",
        "Thalamic neurons": "Thalamic\nneurons",
    }
    return replacements.get(label, label)


def draw_panel(
    ax,
    *,
    study_id: str,
    data: pd.DataFrame,
    routes: pd.DataFrame,
    target_col: str,
    target_order: list[str],
    colors: dict[str, str],
    max_points: int,
    seed: int,
) -> None:
    study_data = data.loc[data["study_id"].astype(str) == study_id].copy()
    study_routes = routes.loc[routes["study_id"].astype(str) == study_id].copy()
    sampled = stratified_sample(study_data, target_col, max_points, seed)
    target_pos = target_positions(study_data, target_order)

    for label in target_order:
        sub = sampled.loc[sampled[target_col].astype(str) == label]
        if sub.empty:
            continue
        ax.scatter(
            pd.to_numeric(sub["UMAP1_plot"], errors="coerce"),
            pd.to_numeric(sub["UMAP2_plot"], errors="coerce"),
            np.zeros(sub.shape[0]),
            s=3.0,
            c=colors.get(label, "#999999"),
            alpha=0.24,
            depthshade=False,
            rasterized=True,
        )

    centroid_cols = ["source_group_id", "source_plot_label", "source_group_name"]
    centroids = (
        study_data.groupby(centroid_cols, observed=True)
        .agg(source_umap1=("UMAP1_plot", "mean"), source_umap2=("UMAP2_plot", "mean"), n_cells=("cell_id", "size"))
        .reset_index()
    )
    ax.scatter(
        centroids["source_umap1"],
        centroids["source_umap2"],
        np.zeros(centroids.shape[0]),
        s=np.clip(18 + np.sqrt(centroids["n_cells"]) * 1.4, 45, 140),
        c="#111111",
        alpha=0.95,
        depthshade=False,
    )
    for row in centroids.itertuples(index=False):
        ax.text(row.source_umap1, row.source_umap2, 0.08, str(row.source_plot_label), fontsize=8.5, ha="center", va="center")

    for label, (tx, ty, tz) in target_pos.items():
        ax.scatter([tx], [ty], [tz], s=85, c=colors.get(label, "#999999"), edgecolors="#111111", linewidths=0.45, depthshade=False)
        ax.text(tx, ty, tz + 0.08, target_label_text(label), fontsize=7.0, ha="center", va="bottom")

    for row in study_routes.itertuples(index=False):
        target = getattr(row, target_col)
        tx, ty, tz = target_pos[str(target)]
        sx = float(row.source_umap1)
        sy = float(row.source_umap2)
        width = 0.25 + 7.0 * float(row.fraction_of_source)
        alpha = float(np.clip(0.18 + 0.58 * row.mean_prediction_score, 0.20, 0.78))
        ax.plot(
            [sx, tx],
            [sy, ty],
            [0.04, tz - 0.04],
            color=colors.get(str(target), "#999999"),
            alpha=alpha,
            linewidth=width,
            solid_capstyle="round",
        )

    set_axes_equalish(ax, study_data, target_pos)
    ax.view_init(elev=23, azim=-62)
    ax.set_title(f"{STUDY_TITLES.get(study_id, study_id)}\n{study_data.shape[0]:,} cells after visualization filters", fontsize=10, pad=12)
    ax.set_xlabel("UMAP1", labelpad=2)
    ax.set_ylabel("UMAP2", labelpad=2)
    ax.set_zlabel("routing plane", labelpad=2)
    ax.set_zticks([0, 1.35])
    ax.set_zticklabels(["query UMAP", "Shi target"], fontsize=7)
    ax.tick_params(axis="both", labelsize=7, pad=0)
    ax.grid(True, alpha=0.23)
    ax.xaxis.pane.set_alpha(0.02)
    ax.yaxis.pane.set_alpha(0.02)
    ax.zaxis.pane.set_alpha(0.02)


def save_figure(fig: plt.Figure, outdir: Path, stem: str) -> None:
    for ext in ["png", "pdf", "svg"]:
        path = outdir / f"{stem}.{ext}"
        fig.savefig(path, dpi=FIGURE_EXPORT_DPI, bbox_inches="tight")


def render_one(
    *,
    data: pd.DataFrame,
    outdir: Path,
    target_col: str,
    score_col: str,
    target_order: list[str],
    colors: dict[str, str],
    stem: str,
    title: str,
    max_points: int,
    min_route_fraction: float,
    top_routes_per_source: int,
    seed: int,
) -> pd.DataFrame:
    routes_all = route_table(data, target_col=target_col, score_col=score_col, target_order=target_order)
    routes_plot = selected_routes(routes_all, min_fraction=min_route_fraction, top_n=top_routes_per_source)
    routes_all["plotted_route"] = False
    key_cols = ["study_id", "source_group_id", target_col]
    plot_keys = routes_plot[key_cols].drop_duplicates()
    routes_all = routes_all.merge(plot_keys.assign(plotted_route=True), on=key_cols, how="left", suffixes=("", "_selected"))
    routes_all["plotted_route"] = routes_all["plotted_route_selected"].fillna(False).astype(bool)
    routes_all = routes_all.drop(columns=["plotted_route_selected"])

    fig = plt.figure(figsize=(14.4, 7.8))
    fig.suptitle(title, fontsize=13, y=0.98)
    for idx, study_id in enumerate(STUDIES, start=1):
        ax = fig.add_subplot(1, 2, idx, projection="3d")
        draw_panel(
            ax,
            study_id=study_id,
            data=data,
            routes=routes_plot,
            target_col=target_col,
            target_order=target_order,
            colors=colors,
            max_points=max_points,
            seed=seed + idx,
        )
    handles = [
        plt.Line2D([0], [0], color=colors.get(label, "#999999"), lw=4, label=label)
        for label in target_order
        if (data[target_col].astype(str) == label).any()
    ]
    fig.legend(handles=handles, loc="lower center", ncol=min(5, len(handles)), frameon=False, fontsize=8, bbox_to_anchor=(0.5, 0.035))
    fig.text(
        0.5,
        0.105,
        "Source points are recoded DIV30/DIV90 classes; target nodes are Shi TransferData winners. Line width = fraction of source class; opacity = mean prediction score.",
        ha="center",
        va="center",
        fontsize=8,
    )
    fig.subplots_adjust(left=0.02, right=0.985, top=0.88, bottom=0.20, wspace=0.12)
    save_figure(fig, outdir, stem)
    plt.close(fig)
    return routes_all


def main() -> None:
    args = parse_args()
    outdir = args.outdir
    outdir.mkdir(parents=True, exist_ok=True)
    data, filter_summary = load_varela_plot_data(args.project_root, args.run_label)
    filter_summary.to_csv(outdir / "shi_query_projection_routing_filter_summary.tsv", sep="\t", index=False)

    stage_routes = render_one(
        data=data,
        outdir=outdir,
        target_col="shi_seurat_full_predicted_shi_week_label",
        score_col="shi_seurat_full_week_prediction_score",
        target_order=GW_ORDER,
        colors=GW_COLORS,
        stem="shi_query_projection_routing_predicted_stage_static",
        title="Static routing diagnostic: DIV30/DIV90 UMAP regions to Shi predicted stage",
        max_points=args.max_cell_points_per_study,
        min_route_fraction=args.min_route_fraction,
        top_routes_per_source=args.top_routes_per_source,
        seed=args.seed,
    )
    stage_routes.to_csv(outdir / "shi_query_projection_routing_predicted_stage_routes.tsv", sep="\t", index=False)

    major_routes = render_one(
        data=data,
        outdir=outdir,
        target_col="shi_seurat_full_predicted_shi_label",
        score_col="shi_seurat_full_prediction_score",
        target_order=SHI_LABEL_ORDER,
        colors=LABEL_COLORS,
        stem="shi_query_projection_routing_major_class_static",
        title="Static routing diagnostic: DIV30/DIV90 UMAP regions to Shi predicted major class",
        max_points=args.max_cell_points_per_study,
        min_route_fraction=args.min_route_fraction,
        top_routes_per_source=args.top_routes_per_source,
        seed=args.seed + 100,
    )
    major_routes.to_csv(outdir / "shi_query_projection_routing_major_class_routes.tsv", sep="\t", index=False)

    manifest = pd.DataFrame(
        [
            {"output": str(path), "exists": path.exists(), "size_bytes": path.stat().st_size if path.exists() else 0}
            for path in sorted(outdir.glob("*"))
        ]
    )
    manifest.to_csv(outdir / "shi_query_projection_routing_static_manifest.tsv", sep="\t", index=False)
    print(manifest.to_string(index=False))


if __name__ == "__main__":
    main()
