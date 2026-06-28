#!/usr/bin/env python3
"""Audit GAD1/GAD2 expression in the Splatter cholinergic candidate rows."""

from __future__ import annotations

import argparse
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import sparse


PROJECT_ROOT_DEFAULT = Path("/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder")
H5AD_DEFAULT = (
    PROJECT_ROOT_DEFAULT
    / "results/siletti_2023_whb_reference_label_transfer/source_cellxgene_superclusters/h5ad/siletti_whb_splatter.h5ad"
)
OUTDIR_DEFAULT = (
    PROJECT_ROOT_DEFAULT
    / "results/siletti_2023_whb_reference_label_transfer/siletti_splatter_candidate_chol_gad_feature_audit_v1"
)

MARKER_GENES = ("GAD1", "GAD2", "CHAT", "SLC5A7", "LHX8", "ISL1")
CANDIDATES = (
    {
        "version": "v8",
        "candidate_set": "v8 cluster 400 NT-CHOL NT-VGLUT3",
        "short_label": "v8 cluster 400",
        "cluster_id": "400",
        "subcluster_ids": {"1634", "1635", "1636", "1637", "1638", "1640", "1641", "1642"},
        "transfer_label": "Subpallial Cholinergic neurons",
    },
    {
        "version": "v9b",
        "candidate_set": "v9b cluster 392/subcluster 1639 NT-CHOL NT-GABA NT-VGLUT3",
        "short_label": "v9b 392/1639",
        "cluster_id": "392",
        "subcluster_ids": {"1639"},
        "transfer_label": "Subpallial Cholinergic-GABA neurons",
    },
    {
        "version": "v10",
        "candidate_set": "v10 cluster 398/subcluster 1532 pure NT-CHOL",
        "short_label": "v10 398/1532",
        "cluster_id": "398",
        "subcluster_ids": {"1532"},
        "transfer_label": "Subpallial Pure cholinergic neurons",
    },
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--h5ad", type=Path, default=H5AD_DEFAULT)
    parser.add_argument("--outdir", type=Path, default=OUTDIR_DEFAULT)
    parser.add_argument("--background-cells", type=int, default=25000)
    parser.add_argument("--seed", type=int, default=0)
    return parser.parse_args()


def dense_vector(x) -> np.ndarray:
    if sparse.issparse(x):
        return np.asarray(x.toarray()).ravel()
    return np.asarray(x).ravel()


def gene_indices(var: pd.DataFrame, genes: tuple[str, ...]) -> dict[str, int]:
    gene_names = var["Gene"].astype(str) if "Gene" in var.columns else pd.Series(var.index.astype(str), index=var.index)
    out: dict[str, int] = {}
    for gene in genes:
        hits = np.flatnonzero(gene_names.to_numpy() == gene)
        if hits.size == 0:
            raise KeyError(f"Gene not found in Siletti Splatter var table: {gene}")
        out[gene] = int(hits[0])
    return out


def candidate_metadata(obs: pd.DataFrame) -> pd.DataFrame:
    meta = obs.copy()
    meta["cluster_id"] = meta["cluster_id"].astype(str)
    meta["subcluster_id"] = meta["subcluster_id"].astype(str)
    meta["candidate_version"] = "not_selected"
    meta["candidate_set"] = "not_selected"
    meta["candidate_short_label"] = "not_selected"
    meta["transfer_label"] = "not_selected"
    for candidate in CANDIDATES:
        mask = meta["cluster_id"].eq(candidate["cluster_id"]) & meta["subcluster_id"].isin(candidate["subcluster_ids"])
        meta.loc[mask, "candidate_version"] = candidate["version"]
        meta.loc[mask, "candidate_set"] = candidate["candidate_set"]
        meta.loc[mask, "candidate_short_label"] = candidate["short_label"]
        meta.loc[mask, "transfer_label"] = candidate["transfer_label"]
    return meta


def write_summary_tables(cell_df: pd.DataFrame, outdir: Path) -> None:
    tables = outdir / "tables"
    tables.mkdir(parents=True, exist_ok=True)
    cell_df.to_csv(tables / "siletti_splatter_candidate_chol_clusters_gad_feature_cells.tsv.gz", sep="\t", index=False)

    group_cols = ["candidate_version", "candidate_set", "cluster_id", "transfer_label"]
    rows = []
    for group_keys, group in cell_df.groupby(group_cols, sort=False):
        base = dict(zip(group_cols, group_keys, strict=False))
        base["n_cells"] = int(group.shape[0])
        for gene in MARKER_GENES:
            values = group[gene].to_numpy()
            base[f"{gene}_pct_detected"] = float((values > 0).mean() * 100)
            base[f"{gene}_mean"] = float(np.mean(values))
            base[f"{gene}_median"] = float(np.median(values))
            base[f"{gene}_max"] = float(np.max(values))
        rows.append(base)
    pd.DataFrame(rows).to_csv(tables / "siletti_splatter_candidate_chol_clusters_gad_feature_summary.tsv", sep="\t", index=False)

    sub_rows = []
    sub_cols = ["candidate_version", "candidate_set", "cluster_id", "subcluster_id", "transfer_label"]
    for group_keys, group in cell_df.groupby(sub_cols, sort=True):
        base = dict(zip(sub_cols, group_keys, strict=False))
        base["n_cells"] = int(group.shape[0])
        for gene in MARKER_GENES:
            values = group[gene].to_numpy()
            base[f"{gene}_pct_detected"] = float((values > 0).mean() * 100)
            base[f"{gene}_mean"] = float(np.mean(values))
            base[f"{gene}_median"] = float(np.median(values))
            base[f"{gene}_max"] = float(np.max(values))
        sub_rows.append(base)
    pd.DataFrame(sub_rows).to_csv(
        tables / "siletti_splatter_candidate_chol_clusters_gad_feature_by_subcluster.tsv",
        sep="\t",
        index=False,
    )


def plot_umap_feature(background_umap: np.ndarray, cell_df: pd.DataFrame, gene: str, plots_dir: Path) -> None:
    versions = ["v8", "v9b", "v10"]
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.2), sharex=True, sharey=True)
    vmax = np.percentile(cell_df[gene].to_numpy(), 99) if np.any(cell_df[gene].to_numpy() > 0) else 1.0
    for ax, version in zip(axes, versions, strict=True):
        sub = cell_df[cell_df["candidate_version"].eq(version)]
        ax.scatter(background_umap[:, 0], background_umap[:, 1], s=1, c="#d4d4d4", alpha=0.18, linewidths=0)
        scatter = ax.scatter(
            sub["UMAP_1"],
            sub["UMAP_2"],
            c=sub[gene],
            s=12,
            cmap="viridis",
            vmin=0,
            vmax=max(vmax, 1e-9),
            linewidths=0,
        )
        ax.set_title(f"{version}: {sub['candidate_short_label'].iloc[0]}\nn={sub.shape[0]}", fontsize=10)
        ax.set_xlabel("UMAP 1")
        ax.set_ylabel("UMAP 2")
        ax.set_aspect("equal", adjustable="box")
    cbar = fig.colorbar(scatter, ax=axes.ravel().tolist(), shrink=0.74)
    cbar.set_label(f"{gene} expression")
    fig.suptitle(f"Siletti Splatter candidate cholinergic clusters: {gene}", fontsize=12)
    fig.tight_layout(rect=(0, 0, 0.95, 0.93))
    for suffix in ("png", "pdf"):
        fig.savefig(plots_dir / f"siletti_splatter_candidate_chol_clusters_feature_{gene}.{suffix}", dpi=300)
    plt.close(fig)


def plot_feature_grid(background_umap: np.ndarray, cell_df: pd.DataFrame, plots_dir: Path) -> None:
    versions = ["v8", "v9b", "v10"]
    fig, axes = plt.subplots(
        len(MARKER_GENES),
        len(versions),
        figsize=(11.5, 17.0),
        sharex=True,
        sharey=True,
        constrained_layout=True,
    )
    for row, gene in enumerate(MARKER_GENES):
        values = cell_df[gene].to_numpy()
        vmax = np.percentile(values, 99) if np.any(values > 0) else 1.0
        for col, version in enumerate(versions):
            ax = axes[row, col]
            sub = cell_df[cell_df["candidate_version"].eq(version)]
            ax.scatter(background_umap[:, 0], background_umap[:, 1], s=0.7, c="#d8d8d8", alpha=0.16, linewidths=0)
            scatter = ax.scatter(
                sub["UMAP_1"],
                sub["UMAP_2"],
                c=sub[gene],
                s=9,
                cmap="viridis",
                vmin=0,
                vmax=max(vmax, 1e-9),
                linewidths=0,
            )
            if row == 0:
                ax.set_title(f"{sub['candidate_short_label'].iloc[0]}\nn={sub.shape[0]}", fontsize=10)
            if col == 0:
                ax.set_ylabel(gene, fontsize=11)
            ax.set_aspect("equal", adjustable="box")
            ax.set_xticks([])
            ax.set_yticks([])
        cbar = fig.colorbar(scatter, ax=axes[row, :].ravel().tolist(), shrink=0.72, pad=0.01)
        cbar.ax.tick_params(labelsize=7)
    fig.suptitle("Siletti Splatter candidate cholinergic clusters: FeaturePlot view", fontsize=13)
    for suffix in ("png", "pdf"):
        fig.savefig(plots_dir / f"siletti_splatter_candidate_chol_clusters_featureplot_view_all_markers.{suffix}", dpi=300)
    plt.close(fig)


def plot_marker_violins(cell_df: pd.DataFrame, plots_dir: Path) -> None:
    labels = ["v8 cluster 400", "v9b 392/1639", "v10 398/1532"]
    fig, axes = plt.subplots(2, 3, figsize=(12, 7), sharex=True)
    for ax, gene in zip(axes.ravel(), MARKER_GENES, strict=True):
        values = [cell_df.loc[cell_df["candidate_short_label"].eq(label), gene].to_numpy() for label in labels]
        parts = ax.violinplot(values, showmedians=True, showextrema=False)
        for body in parts["bodies"]:
            body.set_facecolor("#5d8cc1")
            body.set_edgecolor("#2c4a66")
            body.set_alpha(0.85)
        parts["cmedians"].set_color("#111111")
        ax.set_title(gene)
        ax.set_xticks(np.arange(1, len(labels) + 1), labels, rotation=25, ha="right")
        ax.set_ylabel("Expression")
    fig.tight_layout()
    for suffix in ("png", "pdf"):
        fig.savefig(plots_dir / f"siletti_splatter_candidate_chol_clusters_marker_violins.{suffix}", dpi=300)
    plt.close(fig)


def plot_gad_pct(cell_df: pd.DataFrame, plots_dir: Path) -> None:
    labels = ["v8 cluster 400", "v9b 392/1639", "v10 398/1532"]
    x = np.arange(len(labels))
    width = 0.34
    pct = {
        gene: [
            float((cell_df.loc[cell_df["candidate_short_label"].eq(label), gene].to_numpy() > 0).mean() * 100)
            for label in labels
        ]
        for gene in ("GAD1", "GAD2")
    }
    fig, ax = plt.subplots(figsize=(7.5, 4.2))
    ax.bar(x - width / 2, pct["GAD1"], width, label="GAD1", color="#4f7cac")
    ax.bar(x + width / 2, pct["GAD2"], width, label="GAD2", color="#b56a3a")
    ax.set_xticks(x, labels, rotation=20, ha="right")
    ax.set_ylabel("Detected cells (%)")
    ax.set_ylim(0, 100)
    ax.legend(frameon=False)
    ax.set_title("GAD1/GAD2 detection in Siletti Splatter candidate rows")
    fig.tight_layout()
    for suffix in ("png", "pdf"):
        fig.savefig(plots_dir / f"siletti_splatter_candidate_chol_clusters_gad_pct_detected.{suffix}", dpi=300)
    plt.close(fig)


def write_readme(outdir: Path, h5ad: Path) -> None:
    text = f"""# Siletti Splatter candidate cholinergic GAD feature audit

Input H5AD:
{h5ad}

This audit is restricted to the Splatter rows used in the recent cholinergic
candidate versions, not MSN/eccentric MSN:

- v8: cluster 400, subclusters 1634, 1635, 1636, 1637, 1638, 1640, 1641, 1642
- v9b: cluster 392, subcluster 1639
- v10: cluster 398, subcluster 1532

Feature genes plotted/summarized:
GAD1, GAD2, CHAT, SLC5A7, LHX8, ISL1
"""
    reports = outdir / "reports"
    reports.mkdir(parents=True, exist_ok=True)
    (reports / "README_siletti_splatter_candidate_chol_gad_feature_audit.md").write_text(text)


def main() -> None:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    plots_dir = args.outdir / "plots"
    plots_dir.mkdir(parents=True, exist_ok=True)

    rng = np.random.default_rng(args.seed)
    a = ad.read_h5ad(args.h5ad, backed="r")
    meta = candidate_metadata(a.obs)
    selected = meta["candidate_version"].ne("not_selected").to_numpy()
    selected_idx = np.flatnonzero(selected)
    if selected_idx.size == 0:
        raise RuntimeError("No selected candidate cells found in Splatter H5AD.")

    gidx = gene_indices(a.var, MARKER_GENES)
    sub = a[selected_idx, list(gidx.values())].to_memory()
    expr = pd.DataFrame(
        {gene: dense_vector(sub.X[:, i]) for i, gene in enumerate(gidx.keys())},
        index=meta.index[selected_idx],
    )

    umap = np.asarray(a.obsm["X_UMAP"][selected_idx, :])
    cell_df = meta.iloc[selected_idx].reset_index().rename(columns={"index": "cell_barcode"})
    cell_df["UMAP_1"] = umap[:, 0]
    cell_df["UMAP_2"] = umap[:, 1]
    cell_df = pd.concat([cell_df.reset_index(drop=True), expr.reset_index(drop=True)], axis=1)

    n_background = min(args.background_cells, a.n_obs)
    background_idx = rng.choice(np.arange(a.n_obs), size=n_background, replace=False)
    background_umap = np.asarray(a.obsm["X_UMAP"][background_idx, :])

    write_summary_tables(cell_df, args.outdir)
    for gene in MARKER_GENES:
        plot_umap_feature(background_umap, cell_df, gene, plots_dir)
    plot_feature_grid(background_umap, cell_df, plots_dir)
    plot_marker_violins(cell_df, plots_dir)
    plot_gad_pct(cell_df, plots_dir)
    write_readme(args.outdir, args.h5ad)

    summary = pd.read_csv(
        args.outdir / "tables/siletti_splatter_candidate_chol_clusters_gad_feature_summary.tsv",
        sep="\t",
    )
    print(summary.to_string(index=False))
    print(f"Wrote audit to: {args.outdir}")


if __name__ == "__main__":
    main()
