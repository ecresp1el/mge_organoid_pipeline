#!/usr/bin/env python3
"""Plot Velmeshev 2023 interneuron UMAPs from staged UCSC metadata."""

import gzip
import os
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


PROJECT_ROOT_DEFAULT = "/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder"
RUN_LABEL_DEFAULT = "velmeshev_2023_metadata_inventory_v1"
RESULTS_DIRNAME = "velmeshev_2023_metadata_inventory"
VIEW_ID = "in_rna"

CELL_TYPE_COLORS = {
    "Progenitors": "#f4a261",
    "INT": "#6c757d",
    "SST": "#2a9d8f",
    "PV": "#7b2cbf",
    "PV_MP": "#9d4edd",
    "CALB2": "#4895ef",
    "SST_RELN": "#52b788",
    "VIP": "#ef476f",
    "CCK": "#ffb703",
    "SV2C": "#8ecae6",
    "RELN": "#fb8500",
    "NOS": "#264653",
}

AGE_COLORS = {
    "2nd trimester": "#2166ac",
    "3rd trimester": "#67a9cf",
    "0-1 years": "#d1e5f0",
    "1-2 years": "#fddbc7",
    "2-4 years": "#f4a582",
    "4-10 years": "#d6604d",
    "10-20 years": "#b2182b",
    "Adult": "#67001f",
}


def log(message):
    print("[Velmeshev2023UMAP] {}".format(message), flush=True)


def run_dir():
    root = Path(os.environ.get("PROJECT_ROOT", PROJECT_ROOT_DEFAULT)).expanduser().resolve()
    label = os.environ.get("VELMESHEV_2023_RUN_LABEL", RUN_LABEL_DEFAULT)
    return root / "results" / RESULTS_DIRNAME / label


def read_inputs(root):
    meta_path = root / "metadata" / "{}.meta.tsv.gz".format(VIEW_ID)
    umap_path = root / "umap" / "{}.UMAP.coords.tsv.gz".format(VIEW_ID)
    crosswalk_path = root / "tables" / "velmeshev_2023_cluster_label_crosswalk.tsv"
    if not meta_path.exists():
        raise FileNotFoundError("Missing metadata: {}".format(meta_path))
    if not umap_path.exists():
        raise FileNotFoundError("Missing UMAP coordinates: {}".format(umap_path))
    if not crosswalk_path.exists():
        raise FileNotFoundError("Missing cluster crosswalk: {}".format(crosswalk_path))

    log("Reading metadata: {}".format(meta_path))
    meta = pd.read_csv(meta_path, sep="\t")
    id_col = "cellId" if "cellId" in meta.columns else "Cell_ID"
    meta = meta.rename(columns={id_col: "cell_id"})
    meta["cell_id"] = meta["cell_id"].astype(str)

    log("Reading UMAP coordinates: {}".format(umap_path))
    umap = pd.read_csv(
        umap_path,
        sep="\t",
        header=None,
        names=["cell_id", "umap_1", "umap_2"],
    )
    umap["cell_id"] = umap["cell_id"].astype(str)

    log("Joining metadata to UMAP coordinates")
    data = umap.merge(meta, on="cell_id", how="left", sort=False)
    if data["Seurat_Clusters"].isna().any():
        missing = int(data["Seurat_Clusters"].isna().sum())
        raise ValueError("UMAP join left {} cells without Seurat_Clusters".format(missing))

    crosswalk = pd.read_csv(crosswalk_path, sep="\t")
    crosswalk = crosswalk.loc[crosswalk["view_id"] == VIEW_ID].copy()
    return data, crosswalk


def cluster_sort_key(value):
    parts = []
    for chunk in str(value).replace("_", ".").split("."):
        if chunk.isdigit():
            parts.append((0, int(chunk)))
        else:
            parts.append((1, chunk))
    return parts


def prep_axis(ax):
    ax.set_aspect("equal", adjustable="box")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel("UMAP 1")
    ax.set_ylabel("UMAP 2")
    for spine in ax.spines.values():
        spine.set_visible(False)


def save(fig, png_path):
    png_path.parent.mkdir(parents=True, exist_ok=True)
    pdf_path = png_path.with_suffix(".pdf")
    fig.savefig(png_path, dpi=300, bbox_inches="tight")
    fig.savefig(pdf_path, bbox_inches="tight")
    plt.close(fig)
    return [png_path, pdf_path]


def plot_clusters(data, crosswalk, out_dir):
    cluster_counts = data["Seurat_Clusters"].astype(str).value_counts().to_dict()
    clusters = sorted(cluster_counts, key=cluster_sort_key)
    label_by_cluster = crosswalk.set_index("cluster_id")["dominant_label"].astype(str).to_dict()
    color_by_cluster = {}
    for cluster in clusters:
        label = label_by_cluster.get(cluster, "unknown")
        color_by_cluster[cluster] = CELL_TYPE_COLORS.get(label, "#999999")

    fig, ax = plt.subplots(figsize=(10.5, 8.5))
    ax.scatter(data["umap_1"], data["umap_2"], s=0.8, c="#e7e7e7", linewidths=0, rasterized=True)
    handles = []
    for cluster in clusters:
        subset = data.loc[data["Seurat_Clusters"].astype(str) == cluster]
        label = label_by_cluster.get(cluster, "unknown")
        color = color_by_cluster[cluster]
        artist = ax.scatter(
            subset["umap_1"],
            subset["umap_2"],
            s=1.2,
            c=color,
            linewidths=0,
            rasterized=True,
            label="cluster {}: {} (n={:,})".format(cluster, label, subset.shape[0]),
        )
        handles.append(artist)

        x = float(np.median(subset["umap_1"]))
        y = float(np.median(subset["umap_2"]))
        ax.text(
            x,
            y,
            "{}\n{}".format(cluster, label),
            ha="center",
            va="center",
            fontsize=6.4,
            color="black",
            bbox={"boxstyle": "round,pad=0.18", "fc": "white", "ec": color, "lw": 0.7, "alpha": 0.82},
        )

    prep_axis(ax)
    ax.set_title("Velmeshev 2023 interneuron UMAP: Seurat clusters and Cell_Type labels")
    ax.legend(
        handles=handles,
        loc="center left",
        bbox_to_anchor=(1.02, 0.5),
        frameon=False,
        markerscale=5,
        fontsize=8,
        title="Cluster: Cell_Type",
    )
    return save(fig, out_dir / "velmeshev_2023_in_rna_umap_clusters_cell_type_labels.png")


def plot_age(data, out_dir):
    age_order = [
        "2nd trimester",
        "3rd trimester",
        "0-1 years",
        "1-2 years",
        "2-4 years",
        "4-10 years",
        "10-20 years",
        "Adult",
    ]
    counts = data["Age_Range"].astype(str).value_counts().to_dict()
    fig, ax = plt.subplots(figsize=(8.4, 7.2))
    handles = []
    for age in age_order:
        subset = data.loc[data["Age_Range"].astype(str) == age]
        if subset.empty:
            continue
        artist = ax.scatter(
            subset["umap_1"],
            subset["umap_2"],
            s=1.0,
            c=AGE_COLORS.get(age, "#999999"),
            linewidths=0,
            rasterized=True,
            label="{} (n={:,})".format(age, counts.get(age, 0)),
        )
        handles.append(artist)
    prep_axis(ax)
    ax.set_title("Velmeshev 2023 interneuron UMAP: Age_Range")
    ax.legend(
        handles=handles,
        loc="center left",
        bbox_to_anchor=(1.02, 0.5),
        frameon=False,
        markerscale=5,
        fontsize=8,
        title="Age_Range",
    )
    return save(fig, out_dir / "velmeshev_2023_in_rna_umap_age_range.png")


def write_manifest(root, paths, data):
    table_dir = root / "tables"
    rows = []
    for path in paths:
        rows.append(
            {
                "plot_path": str(path),
                "view_id": VIEW_ID,
                "n_cells": data.shape[0],
                "cluster_field": "Seurat_Clusters",
                "label_field": "Cell_Type",
                "age_field": "Age_Range",
            }
        )
    manifest = table_dir / "velmeshev_2023_in_rna_umap_plot_manifest.tsv"
    pd.DataFrame(rows).to_csv(manifest, sep="\t", index=False)
    return manifest


def main():
    root = run_dir()
    data, crosswalk = read_inputs(root)
    out_dir = root / "plots" / "in_rna_umap"
    plot_paths = []
    log("Plotting cluster/Cell_Type UMAP")
    plot_paths.extend(plot_clusters(data, crosswalk, out_dir))
    log("Plotting Age_Range UMAP")
    plot_paths.extend(plot_age(data, out_dir))
    manifest = write_manifest(root, plot_paths, data)
    log("Wrote manifest: {}".format(manifest))
    log("Wrote plots to: {}".format(out_dir))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
