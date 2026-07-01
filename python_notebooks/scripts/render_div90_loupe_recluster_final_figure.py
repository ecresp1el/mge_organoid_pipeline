#!/usr/bin/env python3
"""Render DIV90 Loupe-export recluster UMAPs as a plot-only final figure.

The Loupe CSVs are treated as the figure source of truth. The DIV90 AnnData
cache is used only to add existing metadata/context and canonical UMAP
coordinates for subset-location panels.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import shutil
from pathlib import Path

import anndata as ad
import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


PROJECT_ROOT_DEFAULT = "/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder"
CSV_DIR_DEFAULT = (
    "/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/"
    "data/raw/loupe_exports/div90_recluster_annotations"
)
H5AD_DEFAULT = (
    "/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/"
    "results/python_anndata/varela_div90.h5ad"
)
LOUPE_TRACKS_DEFAULT = (
    "/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/"
    "data/raw/loupe_exports/div90_recluster_annotations/recovered_from_cloupe/"
    "loupe_recovered_selected_celltracks.tsv.gz"
)
FIGURE_ID_DEFAULT = "fig_div90_loupe_recluster_annotations_v1"

RECLUSTER_SPECS = [
    {
        "set_id": "cortical_only",
        "title": "Cortical-only recluster",
        "header_title": "Cortical",
        "track_name": "Cortical only recluster",
        "background_pattern": r"^UMAP-Projection_corticalonly\.csv$",
        "recluster_pattern": r"^Cortical only recluster-UMAP-Projection\.csv$",
        "color": "#08519C",
    },
    {
        "set_id": "progenitor_astro_glia",
        "title": "Progenitor and astro/glia recluster",
        "header_title": "Progenitors, Astro/Glia",
        "track_name": "Subcluster Progen & Astro",
        "background_pattern": r"^UMAP-Projection_progen_and_glia\.csv$",
        "recluster_pattern": r"^Subcluster Progen & Astro-UMAP-Projection\.csv$",
        "color": "#00796B",
    },
    {
        "set_id": "subpallial_only",
        "title": "Subpallial-only recluster",
        "header_title": "Subpallial",
        "track_name": "Subpallial only Recluster",
        "background_pattern": r"^UMAP-Projection_suppallialonly\.csv$",
        "recluster_pattern": r"^Subpallial only Recluster-UMAP-Projection\.csv$",
        "color": "#AD1457",
    },
]

PARENT_PALETTE = [
    "#1f77b4",
    "#2ca02c",
    "#d62728",
    "#9467bd",
    "#8c564b",
    "#e377c2",
    "#7f7f7f",
    "#bcbd22",
    "#17becf",
    "#ff7f0e",
    "#4daf4a",
    "#984ea3",
    "#a65628",
]

PANEL_LABEL_COLORS = {
    "cortical_only": {
        "Cluster 1": "#0B0B0B",
        "2 PV chand": "#08306B",
        "Cluster 3": "#005F99",
        "6+9 SST/NPY": "#00A3A3",
        "4 PLXNA4/LMO4": "#7A869A",
        "5 RND3/EFS": "#9AA7B2",
        "7 PDK1/PFKFB3": "#0077CC",
        "8 SST/CHODL": "#4D908E",
        "10 PCDH9 Basket": "#505A64",
        "11 VIM/SPARCL1": "#00A6D6",
        "Cluster 6": "#F2F2F2",
    },
    "progenitor_astro_glia": {
        "Cluster 1": "#00441B",
        "Cluster 2": "#2E7D32",
        "Cluster 3": "#A6D96A",
        "4-Sub pal fated": "#8C510A",
        "Cluster 5": "#BF812D",
        "6-Div Progen": "#01665E",
        "7-SST/NPY unknown": "#D9EF8B",
    },
    "subpallial_only": {
        "1+3 GP projection": "#542788",
        "7+9 striatal?": "#E6AB02",
        "2 LCN9": "#D95F02",
        "4 MSRA": "#F1A340",
        "5 NPTN": "#B2182B",
        "6 PCP4/ISL1": "#C51B7D",
        "Cluster 8": "#F768A1",
        "10 meis2/QKI/NR2F1/LCN9": "#FB9A99",
        "11 SOX14/SOX3": "#8E0152",
        "Cluster 12": "#FDD0A2",
    },
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", default=os.environ.get("PROJECT_ROOT", PROJECT_ROOT_DEFAULT))
    parser.add_argument("--csv-dir", default=CSV_DIR_DEFAULT)
    parser.add_argument("--h5ad", default=H5AD_DEFAULT)
    parser.add_argument("--loupe-tracks", default=LOUPE_TRACKS_DEFAULT)
    parser.add_argument("--figure-id", default=FIGURE_ID_DEFAULT)
    parser.add_argument("--outdir", default=None)
    parser.add_argument("--dpi", type=int, default=600)
    return parser.parse_args()


def safe_label(value: str) -> str:
    value = re.sub(r"^[0-9]+\s*-\s*", "", str(value))
    value = value.replace("GABergic", "GABAergic")
    return value


def figure_label(value: str) -> str:
    label = safe_label(value)
    replacements = [
        ("SST+, NPY +, Cortical Fated", "SST+/NPY+ cortical-fated"),
        ("PV precursors/Migrating cells/Cortical-fated", "PV/migrating cortical-fated"),
        ("Pre-Astrocytes/Astrocytes", "Pre-astrocytes/astrocytes"),
        ("LHX8+ vMGE GABAergic Striatal/GP fated", "LHX8+ striatal/GP-fated"),
        ("MGE Striatal/GP Fated", "MGE striatal/GP-fated"),
    ]
    for old, new in replacements:
        label = label.replace(old, new)
    return label


def natural_cluster_key(value: object) -> tuple[int, str]:
    text = str(value)
    try:
        return (0, f"{int(float(text)):04d}")
    except ValueError:
        return (1, text)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def find_one(csv_dir: Path, pattern: str) -> Path:
    hits = [p for p in csv_dir.glob("*.csv") if re.match(pattern, p.name)]
    if len(hits) != 1:
        raise FileNotFoundError(f"Expected one CSV matching {pattern!r}; found {len(hits)} in {csv_dir}")
    return hits[0]


def read_loupe_csv(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    required = ["Barcode", "X Coordinate", "Y Coordinate"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"{path} is missing required Loupe columns: {', '.join(missing)}")
    out = df[required].copy()
    out.columns = ["cell_id", "loupe_x", "loupe_y"]
    out["cell_id"] = out["cell_id"].astype(str)
    out["loupe_x"] = pd.to_numeric(out["loupe_x"], errors="coerce")
    out["loupe_y"] = pd.to_numeric(out["loupe_y"], errors="coerce")
    return out.dropna(subset=["cell_id", "loupe_x", "loupe_y"])


def load_obs(h5ad_path: Path) -> tuple[pd.DataFrame, str]:
    a = ad.read_h5ad(h5ad_path, backed="r")
    obs = a.obs.copy()
    obs.index.name = None
    obs["cell_id"] = obs.index.astype(str)
    if "cell_id" in obs.columns:
        obs["cell_id"] = obs["cell_id"].astype(str)
    reduction_key = "X_umap_seurat" if "X_umap_seurat" in a.obsm.keys() else "X_umap"
    coords = pd.DataFrame(a.obsm[reduction_key], index=a.obs_names, columns=["div90_umap_1", "div90_umap_2"])
    coords.index.name = None
    coords["cell_id"] = coords.index.astype(str)
    coords = coords.reset_index(drop=True)
    obs = obs.merge(coords, on="cell_id", how="left")
    return obs, reduction_key


def setup_dirs(outdir: Path) -> dict[str, Path]:
    dirs = {
        "root": outdir,
        "png": outdir / "figures" / "png",
        "pdf": outdir / "figures" / "pdf",
        "svg": outdir / "figures" / "svg",
        "tables": outdir / "tables",
        "logs": outdir / "logs",
        "provenance": outdir / "provenance",
        "source_csv": outdir / "tables" / "source_loupe_csv",
    }
    for path in dirs.values():
        path.mkdir(parents=True, exist_ok=True)
    return dirs


def parent_palette(obs: pd.DataFrame) -> dict[str, str]:
    clusters = sorted(obs["cluster_id"].dropna().astype(str).unique(), key=natural_cluster_key)
    return {cluster: PARENT_PALETTE[i % len(PARENT_PALETTE)] for i, cluster in enumerate(clusters)}


def apply_plot_style() -> None:
    matplotlib.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "DejaVu Sans", "Liberation Sans"],
            "svg.fonttype": "none",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "axes.linewidth": 0.7,
            "axes.labelsize": 15.0,
            "xtick.labelsize": 11.0,
            "ytick.labelsize": 11.0,
        }
    )


def format_axes(ax: plt.Axes, xlabel: str = "UMAP 1", ylabel: str = "UMAP 2") -> None:
    ax.set_aspect("equal", adjustable="datalim")
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.tick_params(length=0, width=0.0, labelbottom=False, labelleft=False)


def format_inset_axes(ax: plt.Axes) -> None:
    ax.set_aspect("equal", adjustable="datalim")
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_facecolor("white")


def save_figure(fig: plt.Figure, dirs: dict[str, Path], stem: str, dpi: int) -> None:
    fig.savefig(dirs["png"] / f"{stem}.png", dpi=dpi, bbox_inches="tight", facecolor="white")
    fig.savefig(dirs["pdf"] / f"{stem}.pdf", dpi=dpi, bbox_inches="tight", facecolor="white")
    fig.savefig(dirs["svg"] / f"{stem}.svg", dpi=dpi, bbox_inches="tight", facecolor="white")


def annotate_cluster_centers(ax: plt.Axes, df: pd.DataFrame, xcol: str, ycol: str) -> None:
    if "cluster_id" not in df.columns:
        return
    centers = (
        df.dropna(subset=["cluster_id", xcol, ycol])
        .assign(cluster_id=lambda x: x["cluster_id"].astype(str))
        .loc[lambda x: x["cluster_id"] != "unmatched_div90_metadata"]
        .groupby("cluster_id", observed=False)[[xcol, ycol]]
        .median()
        .reset_index()
    )
    for _, row in centers.iterrows():
        ax.text(
            row[xcol],
            row[ycol],
            row["cluster_id"],
            ha="center",
            va="center",
            fontsize=6.7,
            fontweight="bold",
            color="black",
            bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="0.55", lw=0.4, alpha=0.9),
        )


def short_loupe_label(label: str, index: object) -> str:
    label = str(label)
    if re.match(r"^\d+", label):
        return re.match(r"^\d+", label).group(0)
    if "+" in label and re.match(r"^\d+\+\d+", label):
        return re.match(r"^\d+\+\d+", label).group(0)
    try:
        return str(int(float(index)) + 1)
    except Exception:
        return label[:5]


def annotate_loupe_centers(ax: plt.Axes, df: pd.DataFrame, xcol: str, ycol: str) -> None:
    centers = (
        df.dropna(subset=["loupe_label", xcol, ycol])
        .groupby(["loupe_label", "loupe_label_index"], observed=False)[[xcol, ycol]]
        .median()
        .reset_index()
    )
    for _, row in centers.iterrows():
        ax.text(
            row[xcol],
            row[ycol],
            short_loupe_label(row["loupe_label"], row["loupe_label_index"]),
            ha="center",
            va="center",
            fontsize=13.0,
            fontweight="bold",
            color="black",
            bbox=dict(boxstyle="round,pad=0.20", fc="white", ec="0.50", lw=0.6, alpha=0.94),
        )


def cluster_key_rows(local: pd.DataFrame) -> list[dict[str, str]]:
    comp = (
        local.assign(
            loupe_label=lambda x: x["loupe_label"].astype(str),
            loupe_color=lambda x: x["loupe_color"].fillna("#777777"),
        )
        .groupby(["loupe_label", "loupe_color"], dropna=False, observed=False)
        .size()
        .reset_index(name="n_cells")
        .sort_values("n_cells", ascending=False)
    )
    comp["fraction"] = comp["n_cells"] / max(1, local.shape[0])
    rows = []
    for _, row in comp.iterrows():
        pct = 100 * float(row["fraction"])
        rows.append(
            {
                "label": str(row["loupe_label"]),
                "color": str(row["loupe_color"]),
                "text": f"{row['loupe_label']} ({int(row['n_cells']):,}; {pct:.1f}%)",
            }
        )
    return rows


def draw_cluster_key(ax: plt.Axes, local: pd.DataFrame) -> None:
    ax.axis("off")
    rows = cluster_key_rows(local)
    y0 = 0.96
    row_step = 0.083
    for i, row in enumerate(rows):
        x = 0.0
        y = y0 - i * row_step
        ax.add_patch(
            plt.Circle((x + 0.018, y - 0.001), 0.022, transform=ax.transAxes, color=row["color"], clip_on=False)
        )
        ax.text(
            x + 0.058,
            y,
            row["text"],
            ha="left",
            va="top",
            fontsize=12.6,
            linespacing=1.04,
        )


def load_loupe_tracks(path: Path) -> pd.DataFrame:
    tracks = pd.read_csv(path, sep="\t")
    required = {"cell_id", "track_name", "label_index", "label", "color_hex"}
    missing = required.difference(tracks.columns)
    if missing:
        raise ValueError(f"{path} is missing required columns: {', '.join(sorted(missing))}")
    tracks = tracks.loc[tracks["label"].notna()].copy()
    tracks["cell_id"] = tracks["cell_id"].astype(str)
    tracks["track_name"] = tracks["track_name"].astype(str)
    tracks["loupe_label"] = tracks["label"].astype(str)
    tracks["loupe_label_index"] = tracks["label_index"]
    tracks["loupe_color_original"] = tracks["color_hex"].astype(str)
    return tracks[["cell_id", "track_name", "loupe_label_index", "loupe_label", "loupe_color_original"]]


def rgb_distance(hex_a: str, hex_b: str) -> float:
    a = tuple(int(hex_a.lstrip("#")[i : i + 2], 16) for i in (0, 2, 4))
    b = tuple(int(hex_b.lstrip("#")[i : i + 2], 16) for i in (0, 2, 4))
    return sum((x - y) ** 2 for x, y in zip(a, b)) ** 0.5


def assert_palette_contrast(min_distance: float = 45.0) -> None:
    entries = [
        (set_id, label, color)
        for set_id, color_map in PANEL_LABEL_COLORS.items()
        for label, color in color_map.items()
    ]
    too_close = []
    for i, left in enumerate(entries):
        for right in entries[i + 1 :]:
            dist = rgb_distance(left[2], right[2])
            if dist < min_distance:
                too_close.append((left, right, dist))
    assert not too_close, "Palette colors too similar: " + "; ".join(
        f"{a[0]}:{a[1]} {a[2]} vs {b[0]}:{b[1]} {b[2]} dist={dist:.1f}"
        for a, b, dist in too_close
    )


def assign_panel_palette(local: pd.DataFrame, set_id: str) -> pd.DataFrame:
    color_map = PANEL_LABEL_COLORS[set_id]
    labels_in_data = set(local["loupe_label"].dropna().astype(str).unique())
    missing = labels_in_data.difference(color_map)
    assert not missing, f"Missing color assignments for {set_id}: {sorted(missing)}"
    local["loupe_color"] = local["loupe_label"].astype(str).map(color_map).fillna("#BDBDBD")
    return local


def build_tables(
    specs: list[dict[str, object]],
    obs: pd.DataFrame,
    dirs: dict[str, Path],
    csv_dir: Path,
    h5ad_path: Path,
    loupe_tracks: pd.DataFrame,
    loupe_tracks_path: Path,
    reduction_key: str,
) -> tuple[dict[str, pd.DataFrame], pd.DataFrame, pd.DataFrame]:
    data: dict[str, pd.DataFrame] = {}
    inventory_rows = []
    composition_rows = []
    membership_frames = []

    for spec in specs:
        bg_path = find_one(csv_dir, str(spec["background_pattern"]))
        recluster_path = find_one(csv_dir, str(spec["recluster_pattern"]))
        shutil.copy2(bg_path, dirs["source_csv"] / bg_path.name)
        shutil.copy2(recluster_path, dirs["source_csv"] / recluster_path.name)

        bg = read_loupe_csv(bg_path).merge(obs, on="cell_id", how="left")
        local = read_loupe_csv(recluster_path).merge(obs, on="cell_id", how="left")
        track_name = str(spec["track_name"])
        track_values = loupe_tracks.loc[loupe_tracks["track_name"] == track_name].copy()
        local = local.merge(
            track_values.drop(columns=["track_name"]),
            on="cell_id",
            how="left",
        )
        missing_labels = int(local["loupe_label"].isna().sum())
        if missing_labels:
            local["loupe_label"] = local["loupe_label"].fillna("Unassigned in recovered Loupe track")
            local["loupe_color_original"] = local["loupe_color_original"].fillna("#BDBDBD")
        for frame in (bg, local):
            if "cluster_id" in frame.columns:
                frame["cluster_id"] = frame["cluster_id"].astype(object).where(
                    frame["cluster_id"].notna(),
                    "unmatched_div90_metadata",
                )
            if "cluster_number_name" in frame.columns:
                frame["cluster_number_name"] = frame["cluster_number_name"].astype(object).where(
                    frame["cluster_number_name"].notna(),
                    "Unmatched DIV90 metadata",
                )
            if "orig.ident" in frame.columns:
                frame["orig.ident"] = frame["orig.ident"].astype(object).where(
                    frame["orig.ident"].notna(),
                    "unmatched",
                )
        set_id = str(spec["set_id"])
        local = assign_panel_palette(local, set_id)
        bg["set_id"] = set_id
        local["set_id"] = set_id
        local["set_title"] = str(spec["title"])
        data[set_id] = {"background": bg, "local": local, "spec": spec}

        bg_barcodes = set(bg["cell_id"])
        local_barcodes = set(local["cell_id"])
        matched_meta = int(local["cluster_id"].notna().sum()) if "cluster_id" in local else 0
        inventory_rows.extend(
            [
                {
                    "set_id": set_id,
                    "title": spec["title"],
                    "csv_role": "loupe_projection_background",
                    "source_csv": str(bg_path),
                    "copied_csv": str(dirs["source_csv"] / bg_path.name),
                    "n_rows": bg.shape[0],
                    "n_unique_barcodes": bg["cell_id"].nunique(),
                    "sha256": sha256(bg_path),
                },
                {
                    "set_id": set_id,
                    "title": spec["title"],
                    "csv_role": "loupe_recluster_layout",
                    "source_csv": str(recluster_path),
                    "copied_csv": str(dirs["source_csv"] / recluster_path.name),
                    "recovered_loupe_track": track_name,
                    "n_rows": local.shape[0],
                    "n_unique_barcodes": local["cell_id"].nunique(),
                    "n_barcodes_in_background": len(local_barcodes & bg_barcodes),
                    "n_barcodes_with_div90_metadata": matched_meta,
                    "n_barcodes_with_recovered_loupe_annotation": int((local["loupe_label"] != "Unassigned in recovered Loupe track").sum()),
                    "sha256": sha256(recluster_path),
                },
            ]
        )

        if "loupe_label" in local.columns:
            comp = (
                local.assign(
                    loupe_label=lambda x: x["loupe_label"].astype(str),
                    loupe_color=lambda x: x["loupe_color"].astype(str),
                )
                .groupby(["loupe_label", "loupe_color"], dropna=False, observed=False)
                .size()
                .reset_index(name="n_cells")
            )
            comp["set_id"] = set_id
            comp["title"] = spec["title"]
            comp["fraction_of_recluster"] = comp["n_cells"] / max(1, local.shape[0])
            composition_rows.append(comp)

        membership_frames.append(
            local[
                [
                    "cell_id",
                    "set_id",
                    "set_title",
                    "cluster_id",
                    "cluster_number_name",
                    "loupe_label_index",
                    "loupe_label",
                    "loupe_color",
                    "loupe_color_original",
                    "orig.ident",
                    "div90_umap_1",
                    "div90_umap_2",
                ]
            ].copy()
        )

    inventory = pd.DataFrame(inventory_rows)
    composition = pd.concat(composition_rows, ignore_index=True) if composition_rows else pd.DataFrame()
    membership = pd.concat(membership_frames, ignore_index=True)

    inventory.to_csv(dirs["tables"] / "div90_loupe_recluster_source_inventory.tsv", sep="\t", index=False)
    composition.to_csv(dirs["tables"] / "div90_loupe_recluster_recovered_loupe_annotation_composition.tsv", sep="\t", index=False)
    membership.to_csv(dirs["tables"] / "div90_loupe_recluster_membership.tsv.gz", sep="\t", index=False)
    shutil.copy2(loupe_tracks_path, dirs["tables"] / Path(loupe_tracks_path).name)

    provenance = {
        "figure_id": dirs["root"].name,
        "csv_dir": str(csv_dir),
        "h5ad": str(h5ad_path),
        "recovered_loupe_tracks": str(loupe_tracks_path),
        "h5ad_reduction_used_for_context": reduction_key,
        "render_mode": "plot-only; Loupe CSV coordinates plus recovered .cloupe cell-track labels are source of truth",
        "notes": [
            "No clustering or expression analysis is performed.",
            "The Loupe full-projection context panel colors selected cells by the same recovered Loupe reannotation labels as the main recluster map.",
        ],
    }
    (dirs["provenance"] / "render_provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
    return data, inventory, composition


def draw_loupe_projection_context(ax: plt.Axes, bg: pd.DataFrame, local: pd.DataFrame) -> None:
    loupe_projection = bg.merge(
        local[["cell_id", "loupe_label", "loupe_color"]].drop_duplicates("cell_id"),
        on="cell_id",
        how="left",
    )
    background = loupe_projection["loupe_label"].isna()
    ax.scatter(
        loupe_projection.loc[background, "loupe_x"],
        loupe_projection.loc[background, "loupe_y"],
        s=4.5,
        c="#D7D7D7",
        lw=0,
        alpha=0.43,
    )
    for _, sub in loupe_projection.loc[~background].groupby("loupe_label", observed=False):
        ax.scatter(
            sub["loupe_x"],
            sub["loupe_y"],
            s=12.0,
            c=str(sub["loupe_color"].iloc[0]),
            lw=0,
            alpha=0.93,
        )
    format_inset_axes(ax)


def plot_main_panel(data: dict[str, dict[str, object]], obs: pd.DataFrame, palette: dict[str, str], dirs: dict[str, Path], dpi: int) -> None:
    fig = plt.figure(figsize=(19.2, 15.8), constrained_layout=True)
    grid = fig.add_gridspec(nrows=3, ncols=3, height_ratios=[1.45, 4.65, 3.45], hspace=0.045, wspace=0.10)
    for i, set_id in enumerate([str(s["set_id"]) for s in RECLUSTER_SPECS]):
        entry = data[set_id]
        spec = entry["spec"]
        bg = entry["background"]
        local = entry["local"]
        color = str(spec["color"])

        loupe_ax = fig.add_subplot(grid[0, i])
        draw_loupe_projection_context(loupe_ax, bg, local)

        ax = fig.add_subplot(grid[1, i])
        ax.set_facecolor("#FCFCFC")
        for label, sub in local.groupby("loupe_label", observed=False):
            ax.scatter(
                sub["loupe_x"],
                sub["loupe_y"],
                s=28.8,
                c=str(sub["loupe_color"].iloc[0]),
                lw=0,
                alpha=0.86,
            )
        annotate_loupe_centers(ax, local, "loupe_x", "loupe_y")
        ax.set_title(
            f"{spec['header_title']}",
            fontsize=23.5,
            fontweight="bold",
            pad=12,
            color=color,
        )
        format_axes(ax, "Loupe recluster UMAP 1", "Loupe recluster UMAP 2")
        ax.text(
            0.015,
            0.015,
            f"n={local.shape[0]:,}",
            transform=ax.transAxes,
            fontsize=18.0,
            fontweight="bold",
            ha="left",
            va="bottom",
            bbox=dict(fc="white", ec="0.72", lw=0.65, alpha=0.94),
        )

        key_ax = fig.add_subplot(grid[2, i])
        draw_cluster_key(key_ax, local)

    fig.text(
        0.5,
        -0.01,
        "Colors are grouped by recluster family; shades distinguish Loupe-reassigned subclusters within each family.",
        ha="center",
        va="top",
        fontsize=13.0,
        color="#404040",
    )
    save_figure(fig, dirs, "div90_loupe_recluster_final_panel", dpi)
    plt.close(fig)


def spaced_values(top: float, bottom: float, n: int) -> list[float]:
    if n <= 1:
        return [(top + bottom) / 2]
    step = (top - bottom) / (n - 1)
    return [top - i * step for i in range(n)]


def callout_side(row: pd.Series, x_mid: float, y_mid: float, x_range: float, y_range: float) -> str:
    dx = (float(row["loupe_x"]) - x_mid) / x_range
    dy = (float(row["loupe_y"]) - y_mid) / y_range
    if abs(dy) > abs(dx) * 1.15:
        return "top" if dy >= 0 else "bottom"
    return "right" if dx >= 0 else "left"


def elbow_points(x0: float, y0: float, x1: float, y1: float, side: str, x_range: float, y_range: float) -> tuple[list[float], list[float]]:
    if side in {"left", "right"}:
        elbow_x = x0 + (x1 - x0) * 0.45
        return [x0, elbow_x, elbow_x, x1], [y0, y0, y1, y1]
    elbow_y = y0 + (y1 - y0) * 0.45
    return [x0, x0, x1, x1], [y0, elbow_y, elbow_y, y1]


def draw_callout_umap(ax: plt.Axes, local: pd.DataFrame, spec: dict[str, object]) -> None:
    ax.set_facecolor("#FCFCFC")
    for _, sub in local.groupby("loupe_label", observed=False):
        ax.scatter(
            sub["loupe_x"],
            sub["loupe_y"],
            s=28.8,
            c=str(sub["loupe_color"].iloc[0]),
            lw=0,
            alpha=0.86,
            zorder=2,
        )

    centers = (
        local.dropna(subset=["loupe_label", "loupe_x", "loupe_y"])
        .groupby("loupe_label", observed=False)
        .agg(loupe_x=("loupe_x", "median"), loupe_y=("loupe_y", "median"), loupe_color=("loupe_color", "first"))
        .reset_index()
    )
    rows = pd.DataFrame(cluster_key_rows(local))[["label", "text"]]
    centers = centers.merge(rows, left_on="loupe_label", right_on="label", how="left")

    xmin, xmax = float(local["loupe_x"].min()), float(local["loupe_x"].max())
    ymin, ymax = float(local["loupe_y"].min()), float(local["loupe_y"].max())
    x_range = max(xmax - xmin, 1.0)
    y_range = max(ymax - ymin, 1.0)
    x_pad = x_range * 0.86
    y_pad = y_range * 0.32
    ax.set_xlim(xmin - x_pad, xmax + x_pad)
    ax.set_ylim(ymin - y_pad, ymax + y_pad)

    x_mid = (xmin + xmax) / 2
    y_mid = (ymin + ymax) / 2
    centers["side"] = centers.apply(callout_side, axis=1, args=(x_mid, y_mid, x_range, y_range))
    side_order = {"top": "loupe_x", "bottom": "loupe_x", "left": "loupe_y", "right": "loupe_y"}
    for side in ["top", "right", "bottom", "left"]:
        frame = centers.loc[centers["side"] == side].copy()
        if frame.empty:
            continue
        if side in {"top", "bottom"}:
            frame = frame.sort_values("loupe_x").reset_index(drop=True)
            x_positions = np.linspace(xmin - x_pad * 0.18, xmax + x_pad * 0.18, frame.shape[0])
            text_y = ymax + y_pad * 0.58 if side == "top" else ymin - y_pad * 0.58
            va = "bottom" if side == "top" else "top"
            ha = "center"
            positions = [(float(x), text_y, ha, va) for x in x_positions]
        else:
            frame = frame.sort_values("loupe_y", ascending=False).reset_index(drop=True)
            y_positions = spaced_values(ymax + y_pad * 0.20, ymin - y_pad * 0.20, frame.shape[0])
            text_x = xmin - x_pad * 0.50 if side == "left" else xmax + x_pad * 0.16
            ha = "right" if side == "left" else "left"
            positions = [(text_x, float(y), ha, "center") for y in y_positions]

        for (_, row), (text_x, text_y, ha, va) in zip(frame.iterrows(), positions):
            color = str(row["loupe_color"])
            line_x, line_y = elbow_points(
                float(row["loupe_x"]),
                float(row["loupe_y"]),
                float(text_x),
                float(text_y),
                side,
                x_range,
                y_range,
            )
            ax.plot(line_x, line_y, color=color, lw=2.4, alpha=0.98, zorder=4)
            ax.text(
                text_x,
                text_y,
                str(row["text"]),
                ha=ha,
                va=va,
                fontsize=17.0,
                color="#151515",
                bbox=dict(boxstyle="round,pad=0.26", fc="white", ec=color, lw=2.0, alpha=0.96),
                zorder=5,
            )

    ax.set_title(
        f"{spec['header_title']}",
        fontsize=25.0,
        fontweight="bold",
        pad=10,
        color=str(spec["color"]),
    )
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.text(
        0.50,
        0.02,
        f"n={local.shape[0]:,}",
        transform=ax.transAxes,
        fontsize=16.0,
        fontweight="bold",
        ha="center",
        va="bottom",
        bbox=dict(fc="white", ec="0.72", lw=0.65, alpha=0.94),
        zorder=6,
    )


def plot_callout_panel(data: dict[str, dict[str, object]], dirs: dict[str, Path], dpi: int) -> None:
    fig = plt.figure(figsize=(20.0, 18.5), constrained_layout=True)
    grid = fig.add_gridspec(nrows=3, ncols=2, width_ratios=[1.1, 6.0], hspace=0.12, wspace=0.05)
    for i, set_id in enumerate([str(s["set_id"]) for s in RECLUSTER_SPECS]):
        entry = data[set_id]
        spec = entry["spec"]
        bg = entry["background"]
        local = entry["local"]

        loupe_ax = fig.add_subplot(grid[i, 0])
        draw_loupe_projection_context(loupe_ax, bg, local)

        ax = fig.add_subplot(grid[i, 1])
        draw_callout_umap(ax, local, spec)

    fig.text(
        0.5,
        -0.012,
        "Colors are grouped by recluster family; shades distinguish Loupe-reassigned subclusters within each family.",
        ha="center",
        va="top",
        fontsize=13.0,
        color="#404040",
    )
    save_figure(fig, dirs, "div90_loupe_recluster_final_panel_callouts", dpi)
    plt.close(fig)


def plot_loupe_background_panel(data: dict[str, dict[str, object]], dirs: dict[str, Path], dpi: int) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(11.2, 3.5), constrained_layout=True)
    for ax, spec in zip(axes, RECLUSTER_SPECS):
        entry = data[str(spec["set_id"])]
        bg = entry["background"]
        local = entry["local"]
        selected = bg["cell_id"].isin(set(local["cell_id"]))
        ax.scatter(bg.loc[~selected, "loupe_x"], bg.loc[~selected, "loupe_y"], s=6.0, c="#D1D1D1", lw=0, alpha=0.4)
        ax.scatter(bg.loc[selected, "loupe_x"], bg.loc[selected, "loupe_y"], s=10.8, c=str(spec["color"]), lw=0, alpha=0.88)
        ax.set_title(f"{spec['title']}\non Loupe exported projection", fontsize=9, fontweight="bold")
        format_axes(ax, "Loupe projection 1", "Loupe projection 2")
    save_figure(fig, dirs, "div90_loupe_recluster_membership_on_loupe_projection", dpi)
    plt.close(fig)


def plot_local_grid(data: dict[str, dict[str, object]], palette: dict[str, str], dirs: dict[str, Path], dpi: int) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(11.2, 3.5), constrained_layout=True)
    for ax, spec in zip(axes, RECLUSTER_SPECS):
        local = data[str(spec["set_id"])]["local"]
        for label, sub in local.groupby("loupe_label", observed=False):
            ax.scatter(sub["loupe_x"], sub["loupe_y"], s=14.4, c=str(sub["loupe_color"].iloc[0]), lw=0, alpha=0.86)
        ax.set_title(f"{spec['title']}\nLoupe recluster layout", fontsize=9, fontweight="bold")
        format_axes(ax, "Loupe UMAP 1", "Loupe UMAP 2")
    save_figure(fig, dirs, "div90_loupe_recluster_local_layouts_by_loupe_annotation", dpi)
    plt.close(fig)


def write_readme(dirs: dict[str, Path], inventory: pd.DataFrame, composition: pd.DataFrame) -> None:
    top_comp = composition.sort_values(["set_id", "n_cells"], ascending=[True, False]).copy()
    top_comp["label"] = top_comp["loupe_label"].astype(str)
    inventory_md = ["| set_id | csv_role | n_rows | source_csv |", "| --- | --- | ---: | --- |"]
    for _, row in inventory[["set_id", "csv_role", "n_rows", "source_csv"]].iterrows():
        inventory_md.append(f"| {row['set_id']} | {row['csv_role']} | {int(row['n_rows']):,} | `{row['source_csv']}` |")
    lines = [
        "# DIV90 Loupe Recluster Annotation Final Figure",
        "",
        "Plot-only render from Loupe Browser CSV coordinate exports plus recovered Loupe `.cloupe` cell-track annotations. No reclustering, differential expression, or expression scoring is performed here.",
        "",
        "The displayed cluster labels are the Loupe reannotation tracks recovered from `D90 combined SO multi reclustered.cloupe`. Each column places the Loupe-projection context above the main recluster UMAP, so no context panel covers the subcluster map.",
        "",
        "## Main Outputs",
        "",
        "- `figures/png/div90_loupe_recluster_final_panel.png`",
        "- `figures/pdf/div90_loupe_recluster_final_panel.pdf`",
        "- `figures/svg/div90_loupe_recluster_final_panel.svg`",
        "- `figures/png/div90_loupe_recluster_final_panel_callouts.png`",
        "- `figures/pdf/div90_loupe_recluster_final_panel_callouts.pdf`",
        "- `figures/svg/div90_loupe_recluster_final_panel_callouts.svg`",
        "- `figures/png/div90_loupe_recluster_membership_on_loupe_projection.png`",
        "- `figures/png/div90_loupe_recluster_local_layouts_by_loupe_annotation.png`",
        "- `tables/div90_loupe_recluster_recovered_loupe_annotation_composition.tsv`",
        "",
        "## Source CSV Inventory",
        "",
        "\n".join(inventory_md),
        "",
        "## Recovered Loupe Annotation Composition",
        "",
    ]
    for set_id, sub in top_comp.groupby("set_id", sort=False):
        lines.append(f"### {set_id}")
        for _, row in sub.iterrows():
            lines.append(f"- {row['label']}: {int(row['n_cells']):,} cells ({row['fraction_of_recluster']:.1%})")
        lines.append("")
    (dirs["root"] / "README.md").write_text("\n".join(lines) + "\n")


def main() -> None:
    opt = parse_args()
    project_root = Path(opt.project_root)
    csv_dir = Path(opt.csv_dir)
    h5ad_path = Path(opt.h5ad)
    loupe_tracks_path = Path(opt.loupe_tracks)
    outdir = Path(opt.outdir) if opt.outdir else project_root / "final_figures" / opt.figure_id

    apply_plot_style()
    assert_palette_contrast()
    dirs = setup_dirs(outdir)
    obs, reduction_key = load_obs(h5ad_path)
    loupe_tracks = load_loupe_tracks(loupe_tracks_path)
    palette = parent_palette(obs)
    data, inventory, composition = build_tables(
        RECLUSTER_SPECS,
        obs,
        dirs,
        csv_dir,
        h5ad_path,
        loupe_tracks,
        loupe_tracks_path,
        reduction_key,
    )
    plot_main_panel(data, obs, palette, dirs, opt.dpi)
    plot_callout_panel(data, dirs, opt.dpi)
    plot_loupe_background_panel(data, dirs, opt.dpi)
    plot_local_grid(data, palette, dirs, opt.dpi)
    write_readme(dirs, inventory, composition)
    print(f"Done: {outdir}")


if __name__ == "__main__":
    main()
