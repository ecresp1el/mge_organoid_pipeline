#!/usr/bin/env python
"""Generate Jia-style Siletti adult-reference plots for DIV90 neurons.

This script adapts the Jia fetal/adult comparison figure logic to the current
DIV90 organoid setting:

* fetal hMGE cells -> DIV90 neuron-lineage query cells
* adult MGE-derived reference -> Siletti MGE + LAMP5-LHX6/chandelier cells
* fetal classes -> DIV90 cluster/broad-class labels
* developmental stages -> DIV90 samples

The input transfer is the fast KNN Siletti run. By default, the UMAP is fit on
the adult Siletti reference and DIV90 cells are projected into that adult
manifold. This keeps the workflow Python-only while making Panel B an explicit
adult-reference overlay rather than a de novo combined UMAP.
"""

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
import seaborn as sns
import umap
from matplotlib.patches import PathPatch, Rectangle
from matplotlib.path import Path as MplPath
from scipy import sparse
from sklearn.decomposition import TruncatedSVD
from sklearn.preprocessing import normalize

from siletti_div90_fast_knn_label_transfer import (
    log_normalize,
    read_bridge_counts,
    safe_token,
    sparse_variance,
)


PROJECT_ROOT_DEFAULT = Path("/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder")
DEFAULT_BRIDGE_DIR = (
    PROJECT_ROOT_DEFAULT
    / "results/siletti_2023_whb_reference_label_transfer"
    / "siletti_div90_seurat_bridge_v1/mge_llc/seurat_bridge"
)
DEFAULT_TRANSFER_DIR = (
    PROJECT_ROOT_DEFAULT
    / "results/siletti_2023_whb_reference_label_transfer"
    / "siletti_div90_fast_knn_transfer_full_mge_llc_v1/mge_llc/svd50_k50_ref0_query0"
)
DEFAULT_OUTDIR = (
    PROJECT_ROOT_DEFAULT
    / "results/siletti_2023_whb_reference_label_transfer"
    / "siletti_div90_jia_style_figure"
)

ADULT_ORDER = [
    "Cortical SST+ LRP neurons",
    "Cortical SST+ nMt neurons",
    "Cortical SST+ Mt neurons",
    "Cortical PV+ basket neurons",
    "Cortical PV+ Chandelier neurons",
    "Subpallial SST+ LRP neurons",
    "Subpallial SST+ neurons",
    "Subpallial PV+ neurons",
    "Subpallial Cholinergic neurons",
]

BROAD_CLASS_MAP = {
    "Cortical SST+ LRP neurons": "Cortical SST+",
    "Cortical SST+ nMt neurons": "Cortical SST+",
    "Cortical SST+ Mt neurons": "Cortical SST+",
    "Cortical PV+ basket neurons": "Cortical PV+",
    "Cortical PV+ Chandelier neurons": "Cortical PV+",
    "Subpallial SST+ LRP neurons": "Subpallial SST+",
    "Subpallial SST+ neurons": "Subpallial SST+",
    "Subpallial PV+ neurons": "Subpallial PV+",
    "Subpallial Cholinergic neurons": "Subpallial Cholinergic",
}

BROAD_ORDER = [
    "Cortical SST+",
    "Cortical PV+",
    "Subpallial SST+",
    "Subpallial PV+",
    "Subpallial Cholinergic",
    "Other adult label",
    "Unassigned",
]

PALETTE = {
    "Cortical SST+ LRP neurons": "#4c78a8",
    "Cortical SST+ nMt neurons": "#72b7b2",
    "Cortical SST+ Mt neurons": "#54a24b",
    "Cortical PV+ basket neurons": "#e45756",
    "Cortical PV+ Chandelier neurons": "#ff9da6",
    "Subpallial SST+ LRP neurons": "#f58518",
    "Subpallial SST+ neurons": "#b279a2",
    "Subpallial PV+ neurons": "#9d755d",
    "Subpallial Cholinergic neurons": "#bab0ac",
    "Unassigned": "#8f8f8f",
    "Other adult label": "#8f8f8f",
}


def read_tsv_gz(path: Path) -> pd.DataFrame:
    with gzip.open(path, "rt") as handle:
        return pd.read_csv(handle, sep="\t")


def write_tsv_gz(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt") as handle:
        df.to_csv(handle, sep="\t", index=False)


def find_one(path_glob: str, root: Path) -> Path:
    matches = sorted(root.glob(path_glob))
    if len(matches) != 1:
        raise FileNotFoundError(f"Expected exactly one match for {root / path_glob}, found {len(matches)}")
    return matches[0]


def classify_side(label: str) -> str:
    label = str(label)
    if label.startswith("Cortical"):
        return "Cortical"
    if label.startswith("Subpallial"):
        return "Subpallial"
    return "Unassigned"


def classify_broad(label: str) -> str:
    return BROAD_CLASS_MAP.get(str(label), "Other adult label")


def prepare_embedding(
    bridge_dir: Path,
    transfer_dir: Path,
    label_column: str,
    exclude_labels: list[str],
    nfeatures: int,
    n_components: int,
    seed: int,
    umap_mode: str,
) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    x_ref, genes_ref, meta_ref = read_bridge_counts("reference", bridge_dir)
    x_query, genes_query, meta_query = read_bridge_counts("query", bridge_dir)
    if genes_ref != genes_query:
        raise ValueError("Reference and query bridge gene orders differ.")

    if label_column not in meta_ref.columns:
        raise ValueError(f"Reference metadata missing label column: {label_column}")
    meta_ref[label_column] = meta_ref[label_column].fillna("unlabeled_or_na").astype(str)
    if exclude_labels:
        keep_ref = ~meta_ref[label_column].isin(exclude_labels)
        x_ref = x_ref[keep_ref.to_numpy(), :]
        meta_ref = meta_ref.loc[keep_ref].reset_index(drop=True)

    features_path = transfer_dir / "fast_knn/selected_transfer_features.tsv"
    if features_path.exists():
        features = pd.read_csv(features_path, sep="\t")["gene"].astype(str).tolist()
        gene_to_idx = {gene: i for i, gene in enumerate(genes_ref)}
        feature_idx = np.array([gene_to_idx[g] for g in features if g in gene_to_idx], dtype=int)
    else:
        x_tmp = sparse.vstack([log_normalize(x_ref), log_normalize(x_query)], format="csr")
        var = sparse_variance(x_tmp)
        feature_idx = np.argsort(var)[::-1][: min(nfeatures, x_ref.shape[1])]
        feature_idx.sort()
        features = np.asarray(genes_ref, dtype=object)[feature_idx].astype(str).tolist()

    x_ref_norm = log_normalize(x_ref)[:, feature_idx]
    x_query_norm = log_normalize(x_query)[:, feature_idx]
    n_components_use = min(n_components, len(feature_idx) - 1, x_ref_norm.shape[0] - 1, x_query_norm.shape[0] - 1)

    reducer = umap.UMAP(
        n_neighbors=30,
        min_dist=0.3,
        metric="cosine",
        random_state=seed,
        low_memory=True,
    )

    if umap_mode == "reference_project":
        svd = TruncatedSVD(n_components=n_components_use, random_state=seed)
        ref_pcs = normalize(svd.fit_transform(x_ref_norm))
        query_pcs = normalize(svd.transform(x_query_norm))
        ref_coords = reducer.fit_transform(ref_pcs)
        query_coords = reducer.transform(query_pcs)
        fit_description = "SVD and UMAP fit on adult Siletti reference; DIV90 projected with the fitted models"
    elif umap_mode == "combined":
        svd = TruncatedSVD(n_components=n_components_use, random_state=seed)
        ref_pcs = svd.fit_transform(x_ref_norm)
        query_pcs = svd.transform(x_query_norm)
        pcs = np.vstack([normalize(ref_pcs), normalize(query_pcs)])
        coords = reducer.fit_transform(pcs)
        ref_coords = coords[: ref_pcs.shape[0], :]
        query_coords = coords[ref_pcs.shape[0] :, :]
        fit_description = "SVD fit on adult Siletti reference; UMAP fit on combined adult reference plus DIV90 query"
    else:
        raise ValueError(f"Unknown umap_mode: {umap_mode}")

    meta_ref = meta_ref.copy()
    meta_ref["dataset"] = "Siletti adult reference"
    meta_ref["umap_1"] = ref_coords[:, 0]
    meta_ref["umap_2"] = ref_coords[:, 1]

    meta_query = meta_query.copy()
    meta_query["dataset"] = "DIV90 query"
    meta_query["umap_1"] = query_coords[:, 0]
    meta_query["umap_2"] = query_coords[:, 1]

    diag = {
        "n_reference_cells_umap": int(meta_ref.shape[0]),
        "n_query_cells_umap": int(meta_query.shape[0]),
        "n_features_umap": int(len(feature_idx)),
        "n_components_umap": int(n_components_use),
        "svd_explained_variance_ratio_sum": float(svd.explained_variance_ratio_.sum()),
        "umap_n_neighbors": 30,
        "umap_min_dist": 0.3,
        "umap_metric": "cosine",
        "umap_mode": umap_mode,
        "embedding_fit_description": fit_description,
    }
    return meta_ref, meta_query, diag


def attach_predictions(meta_query: pd.DataFrame, transfer_dir: Path, min_score: float) -> pd.DataFrame:
    obs_path = find_one("tables/*_query_obs_with_predictions.tsv.gz", transfer_dir)
    obs = read_tsv_gz(obs_path)
    needed = ["seurat_cell_id", "predicted.id", "prediction.score.max"]
    missing = [col for col in needed if col not in obs.columns]
    if missing:
        raise ValueError(f"Prediction obs table missing columns: {missing}")

    keep_cols = [col for col in obs.columns if col not in meta_query.columns or col == "seurat_cell_id"]
    if "predicted.id" not in keep_cols:
        keep_cols.append("predicted.id")
    if "prediction.score.max" not in keep_cols:
        keep_cols.append("prediction.score.max")
    merged = meta_query.merge(obs[keep_cols], on="seurat_cell_id", how="left", validate="one_to_one")
    merged["predicted.id"] = merged["predicted.id"].fillna("Unassigned").astype(str)
    merged["prediction.score.max"] = pd.to_numeric(merged["prediction.score.max"], errors="coerce").fillna(0.0)
    merged["predicted_side"] = merged["predicted.id"].map(classify_side)
    merged["predicted_broad_class"] = merged["predicted.id"].map(classify_broad)
    assigned = (merged["predicted_side"] != "Unassigned") & (merged["prediction.score.max"] >= min_score)
    merged["assignment_status"] = np.where(assigned, "assigned", "unassigned_or_below_threshold")
    merged.loc[~assigned, "predicted_side"] = "Unassigned"
    merged.loc[~assigned, "predicted_broad_class"] = "Unassigned"
    return merged


def adult_label_order(labels: list[str]) -> list[str]:
    known = [label for label in ADULT_ORDER if label in labels]
    extra = sorted(label for label in labels if label not in known)
    return known + extra


def save_umap_plots(ref: pd.DataFrame, query: pd.DataFrame, plots_dir: Path, label_column: str) -> None:
    plots_dir.mkdir(parents=True, exist_ok=True)
    labels = adult_label_order(sorted(ref[label_column].dropna().astype(str).unique()))

    fig, ax = plt.subplots(figsize=(8, 7))
    for label in labels:
        sub = ref.loc[ref[label_column].astype(str) == label]
        ax.scatter(
            sub["umap_1"],
            sub["umap_2"],
            s=2.2,
            c=PALETTE.get(label, "#777777"),
            label=label,
            alpha=0.55,
            linewidths=0,
            rasterized=True,
        )
    ax.scatter(query["umap_1"], query["umap_2"], s=3, c="#111111", alpha=0.35, linewidths=0, rasterized=True, label="DIV90 query")
    ax.set_xlabel("UMAP 1")
    ax.set_ylabel("UMAP 2")
    ax.set_title("B. Siletti adult MGE/LLC reference with DIV90 overlay")
    ax.legend(markerscale=4, fontsize=7, frameon=False, loc="center left", bbox_to_anchor=(1.02, 0.5))
    fig.savefig(plots_dir / "figure_B_shared_umap_adult_subtypes_div90_overlay.png", dpi=300, bbox_inches="tight")
    fig.savefig(plots_dir / "figure_B_shared_umap_adult_subtypes_div90_overlay.pdf", bbox_inches="tight")
    plt.close(fig)

    assigned = query.loc[query["assignment_status"] == "assigned"].copy()
    omitted = query.loc[query["assignment_status"] != "assigned"].copy()
    side_palette = {"Cortical": "#2f6db3", "Subpallial": "#d4782a", "Unassigned": "#999999"}
    fig, ax = plt.subplots(figsize=(8, 7))
    ax.scatter(ref["umap_1"], ref["umap_2"], s=1.2, c="#d7d7d7", alpha=0.45, linewidths=0, rasterized=True)
    for side in ["Cortical", "Subpallial"]:
        sub = assigned.loc[assigned["predicted_side"] == side]
        if sub.empty:
            continue
        ax.scatter(
            sub["umap_1"],
            sub["umap_2"],
            s=7,
            c=side_palette[side],
            label=f"DIV90 assigned {side}",
            alpha=0.78,
            linewidths=0,
            rasterized=True,
        )
    if not omitted.empty:
        ax.scatter(
            omitted["umap_1"],
            omitted["umap_2"],
            s=4,
            c=side_palette["Unassigned"],
            label="DIV90 omitted/unassigned",
            alpha=0.18,
            linewidths=0,
            rasterized=True,
        )
    ax.set_xlabel("UMAP 1")
    ax.set_ylabel("UMAP 2")
    ax.set_title("B. DIV90 cortical/subpallial assignments in adult reference space")
    ax.legend(markerscale=3, fontsize=8, frameon=False, loc="best")
    fig.savefig(plots_dir / "figure_B_shared_umap_div90_cortical_subpallial_assignments.png", dpi=300, bbox_inches="tight")
    fig.savefig(plots_dir / "figure_B_shared_umap_div90_cortical_subpallial_assignments.pdf", bbox_inches="tight")
    plt.close(fig)


def make_river_edges(query: pd.DataFrame, source_col: str) -> pd.DataFrame:
    assigned = query.loc[query["assignment_status"] == "assigned"].copy()
    if source_col not in assigned.columns:
        raise ValueError(f"Query table missing source class column: {source_col}")
    edges = (
        assigned.groupby([source_col, "predicted.id"], dropna=False)
        .size()
        .reset_index(name="n_cells")
        .rename(columns={source_col: "div90_class", "predicted.id": "adult_subtype"})
    )
    return edges


def draw_river(edges: pd.DataFrame, outpath: Path) -> None:
    if edges.empty:
        raise ValueError("No assigned cells available for river plot.")
    left_totals = edges.groupby("div90_class")["n_cells"].sum().sort_values(ascending=False)
    right_order = adult_label_order(edges["adult_subtype"].astype(str).unique().tolist())
    right_totals = edges.groupby("adult_subtype")["n_cells"].sum().reindex(right_order).dropna()
    total = float(edges["n_cells"].sum())

    def positions(totals: pd.Series) -> dict[str, tuple[float, float]]:
        gap = 0.012
        usable = 1.0 - gap * max(0, len(totals) - 1)
        y = 1.0
        pos = {}
        for label, value in totals.items():
            h = usable * float(value) / total
            pos[str(label)] = (y - h, y)
            y -= h + gap
        return pos

    left_pos = positions(left_totals)
    right_pos = positions(right_totals)
    left_cursor = {k: v[1] for k, v in left_pos.items()}
    right_cursor = {k: v[1] for k, v in right_pos.items()}

    fig, ax = plt.subplots(figsize=(12, max(7, 0.38 * max(len(left_totals), len(right_totals)))))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    left_x0, left_x1 = 0.05, 0.12
    right_x0, right_x1 = 0.88, 0.95
    for label, (y0, y1) in left_pos.items():
        ax.add_patch(Rectangle((left_x0, y0), left_x1 - left_x0, y1 - y0, color="#b8c7d9", ec="white", lw=0.6))
        ax.text(left_x0 - 0.015, (y0 + y1) / 2, f"{label}\n{int(left_totals[label]):,}", ha="right", va="center", fontsize=8)
    for label, (y0, y1) in right_pos.items():
        color = PALETTE.get(label, "#777777")
        ax.add_patch(Rectangle((right_x0, y0), right_x1 - right_x0, y1 - y0, color=color, ec="white", lw=0.6))
        ax.text(right_x1 + 0.015, (y0 + y1) / 2, f"{label}\n{int(right_totals[label]):,}", ha="left", va="center", fontsize=8)

    edges_sorted = edges.sort_values(["div90_class", "adult_subtype"])
    for _, row in edges_sorted.iterrows():
        left = str(row["div90_class"])
        right = str(row["adult_subtype"])
        value = float(row["n_cells"])
        if value <= 0 or left not in left_pos or right not in right_pos:
            continue
        h = value / total
        ly1 = left_cursor[left]
        ly0 = ly1 - h
        left_cursor[left] = ly0
        ry1 = right_cursor[right]
        ry0 = ry1 - h
        right_cursor[right] = ry0
        verts = [
            (left_x1, ly0),
            (0.38, ly0),
            (0.62, ry0),
            (right_x0, ry0),
            (right_x0, ry1),
            (0.62, ry1),
            (0.38, ly1),
            (left_x1, ly1),
            (left_x1, ly0),
        ]
        codes = [
            MplPath.MOVETO,
            MplPath.CURVE4,
            MplPath.CURVE4,
            MplPath.CURVE4,
            MplPath.LINETO,
            MplPath.CURVE4,
            MplPath.CURVE4,
            MplPath.CURVE4,
            MplPath.CLOSEPOLY,
        ]
        patch = PathPatch(
            MplPath(verts, codes),
            facecolor=PALETTE.get(right, "#777777"),
            alpha=0.28,
            edgecolor="none",
        )
        ax.add_patch(patch)

    ax.text((left_x0 + left_x1) / 2, 1.03, "DIV90 class", ha="center", va="bottom", fontsize=11, fontweight="bold")
    ax.text((right_x0 + right_x1) / 2, 1.03, "Adult Siletti subtype", ha="center", va="bottom", fontsize=11, fontweight="bold")
    ax.set_title("C. DIV90 classes mapped to adult inhibitory-neuron subtypes", fontsize=13, pad=20)
    fig.savefig(outpath.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(outpath.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def sample_proportions(query: pd.DataFrame, sample_col: str) -> pd.DataFrame:
    assigned = query.loc[query["assignment_status"] == "assigned"].copy()
    if sample_col not in assigned.columns:
        raise ValueError(f"Query table missing sample column: {sample_col}")
    counts = (
        assigned.groupby([sample_col, "predicted_broad_class"], dropna=False)
        .size()
        .reset_index(name="n_cells")
    )
    totals = counts.groupby(sample_col)["n_cells"].transform("sum")
    counts["fraction"] = counts["n_cells"] / totals
    counts = counts.rename(columns={sample_col: "sample"})
    return counts


def draw_sample_barplot(props: pd.DataFrame, outpath: Path) -> None:
    if props.empty:
        raise ValueError("No assigned cells available for sample bar plot.")
    pivot = props.pivot_table(index="sample", columns="predicted_broad_class", values="fraction", fill_value=0)
    cols = [col for col in BROAD_ORDER if col in pivot.columns] + sorted(col for col in pivot.columns if col not in BROAD_ORDER)
    pivot = pivot[cols]
    colors = {
        "Cortical SST+": "#4c78a8",
        "Cortical PV+": "#e45756",
        "Subpallial SST+": "#f58518",
        "Subpallial PV+": "#9d755d",
        "Subpallial Cholinergic": "#bab0ac",
        "Other adult label": "#8f8f8f",
        "Unassigned": "#8f8f8f",
    }

    fig, ax = plt.subplots(figsize=(max(8, 0.7 * pivot.shape[0]), 5.5))
    bottom = np.zeros(pivot.shape[0])
    x = np.arange(pivot.shape[0])
    for col in pivot.columns:
        vals = pivot[col].to_numpy()
        ax.bar(x, vals, bottom=bottom, label=col, color=colors.get(col, "#777777"), edgecolor="white", linewidth=0.4)
        bottom += vals
    ax.set_xticks(x)
    ax.set_xticklabels(pivot.index.astype(str), rotation=45, ha="right")
    ax.set_ylim(0, 1)
    ax.set_ylabel("Fraction of assigned DIV90 neurons")
    ax.set_title("D. Adult-like broad inhibitory classes across DIV90 samples")
    ax.legend(frameon=False, bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=8)
    fig.savefig(outpath.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(outpath.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bridge-dir", type=Path, default=DEFAULT_BRIDGE_DIR)
    parser.add_argument("--transfer-dir", type=Path, default=DEFAULT_TRANSFER_DIR)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument("--label-column", default="candidate_jia_group")
    parser.add_argument("--exclude-label", default="Excluded / not assigned to Jia-style 9 groups")
    parser.add_argument("--exclude-labels", default=None, help="Optional '||'-separated labels to remove from the reference UMAP.")
    parser.add_argument("--source-class-col", default="div90_broad_class")
    parser.add_argument("--sample-col", default="orig.ident")
    parser.add_argument("--min-score-for-assigned", type=float, default=0.0)
    parser.add_argument("--nfeatures", type=int, default=3000)
    parser.add_argument("--n-components", type=int, default=50)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument(
        "--umap-mode",
        choices=["reference_project", "combined"],
        default="reference_project",
        help="reference_project fits UMAP on adult reference and projects DIV90; combined fits UMAP on both datasets.",
    )
    return parser.parse_args()


def excluded_labels(args: argparse.Namespace) -> list[str]:
    raw = args.exclude_labels if args.exclude_labels is not None else args.exclude_label
    if not raw or str(raw).upper() == "NONE":
        return []
    return [label.strip() for label in str(raw).split("||") if label.strip()]


def main() -> None:
    args = parse_args()
    labels_to_exclude = excluded_labels(args)
    plots_dir = args.outdir / "plots"
    tables_dir = args.outdir / "tables"
    reports_dir = args.outdir / "reports"
    for directory in [plots_dir, tables_dir, reports_dir]:
        directory.mkdir(parents=True, exist_ok=True)

    ref, query_umap, diag = prepare_embedding(
        args.bridge_dir,
        args.transfer_dir,
        args.label_column,
        labels_to_exclude,
        args.nfeatures,
        args.n_components,
        args.seed,
        args.umap_mode,
    )
    query = attach_predictions(query_umap, args.transfer_dir, args.min_score_for_assigned)

    save_umap_plots(ref, query, plots_dir, args.label_column)

    edges = make_river_edges(query, args.source_class_col)
    edges.to_csv(tables_dir / "figure_C_div90_class_to_adult_subtype_edges.tsv", sep="\t", index=False)
    draw_river(edges, plots_dir / "figure_C_div90_class_to_adult_subtype_river")

    props = sample_proportions(query, args.sample_col)
    props.to_csv(tables_dir / "figure_D_sample_adult_broad_class_proportions.tsv", sep="\t", index=False)
    draw_sample_barplot(props, plots_dir / "figure_D_sample_adult_broad_class_proportions")

    assignment_tracking = (
        query.groupby(["assignment_status", "predicted_side", "predicted_broad_class", "predicted.id"], dropna=False)
        .size()
        .reset_index(name="n_cells")
        .sort_values(["assignment_status", "predicted_side", "predicted_broad_class", "predicted.id"])
    )
    assignment_tracking.to_csv(tables_dir / "assignment_tracking.tsv", sep="\t", index=False)
    write_tsv_gz(query, tables_dir / "div90_query_cells_with_siletti_umap_and_assignments.tsv.gz")
    ref[["seurat_cell_id", args.label_column, "transferred_mtg_label", "jia_side", "anatomy_bin", "umap_1", "umap_2"]].to_csv(
        tables_dir / "siletti_reference_cells_with_umap.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )

    config = {
        "bridge_dir": str(args.bridge_dir),
        "transfer_dir": str(args.transfer_dir),
        "outdir": str(args.outdir),
        "label_column": args.label_column,
        "exclude_label": args.exclude_label,
        "exclude_labels": labels_to_exclude,
        "source_class_col": args.source_class_col,
        "sample_col": args.sample_col,
        "min_score_for_assigned": args.min_score_for_assigned,
        "nfeatures_requested": args.nfeatures,
        "n_components_requested": args.n_components,
        "seed": args.seed,
        "umap_mode": args.umap_mode,
        **diag,
    }
    (tables_dir / "plot_config.json").write_text(json.dumps(config, indent=2, sort_keys=True) + "\n")

    readme = [
        "# Siletti DIV90 Jia-Style Figure",
        "",
        "This run adapts the Jia fetal/adult comparison plots to DIV90 neurons.",
        "",
        "- Panel B: shared UMAP of Siletti adult MGE/LLC reference and DIV90 query cells.",
        "- Panel C: river plot from DIV90 class to predicted adult Siletti subtype.",
        "- Panel D: sample-wise proportions of predicted broad adult inhibitory classes.",
        "",
        "`transferred_mtg_label` is an adult-reference annotation from the Siletti/Linnarsson workbook,",
        "not a label transferred from Siletti to DIV90 in this run.",
        "",
        f"UMAP mode: {args.umap_mode}.",
        "",
        f"Assigned-cell threshold: prediction.score.max >= {args.min_score_for_assigned}",
        f"Reference labels excluded from UMAP: {', '.join(labels_to_exclude) if labels_to_exclude else 'none'}",
        "",
        "Generated tables record omitted/unassigned cells so the plots can omit them without losing accounting.",
    ]
    (reports_dir / "README_siletti_div90_jia_style_figure.md").write_text("\n".join(readme) + "\n")

    print(json.dumps(config, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
