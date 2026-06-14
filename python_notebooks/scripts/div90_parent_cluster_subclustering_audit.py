#!/usr/bin/env python3
"""Audit whether DIV90 parent clusters support biologically meaningful subclusters.

This is a local-only analysis from the cached DIV90 H5AD. It does not submit jobs
and does not load the source Seurat RDS.
"""

from __future__ import annotations

import argparse
import gzip
import os
from pathlib import Path
from typing import Iterable

import anndata as ad
import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import chi2_contingency
from sklearn.metrics import adjusted_rand_score


PROJECT_ROOT_DEFAULT = "/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder"
NUMERIC_COVARIATES = ["nCount_RNA", "nFeature_RNA", "percent.mt", "S.Score", "G2M.Score", "CC.Difference"]
CATEGORICAL_COVARIATES = ["orig.ident", "cell_line", "Phase"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", default=os.environ.get("PROJECT_ROOT", PROJECT_ROOT_DEFAULT))
    parser.add_argument("--h5ad", default=None)
    parser.add_argument(
        "--sample-map",
        default="metadata/div30_div90_sample_id_to_biolabel_map.tsv",
        help="Sample map used to derive cell_line from biological_label.",
    )
    parser.add_argument(
        "--outdir",
        default=None,
        help="Default: PROJECT_ROOT/results/div90_parent_cluster_subclustering_audit/div90_parent_cluster_subclustering_audit_v1",
    )
    parser.add_argument("--exclude-clusters", default="12", help="Comma-separated parent cluster IDs to skip.")
    parser.add_argument(
        "--parent-clusters",
        default="all",
        help="Comma-separated parent cluster IDs to audit, or 'all'. Applied after --exclude-clusters.",
    )
    parser.add_argument("--resolutions", default="0.1,0.2,0.3,0.4,0.6,0.8,1.0")
    parser.add_argument("--n-top-genes", type=int, default=2000)
    parser.add_argument("--n-pcs", type=int, default=30)
    parser.add_argument("--n-neighbors", type=int, default=30)
    parser.add_argument("--min-parent-cells", type=int, default=100)
    parser.add_argument("--min-subcluster-cells", type=int, default=40)
    parser.add_argument("--marker-logfc", type=float, default=0.25)
    parser.add_argument("--marker-padj", type=float, default=0.05)
    parser.add_argument("--min-markers-per-subcluster", type=int, default=2)
    parser.add_argument("--max-cramers-v", type=float, default=0.65)
    parser.add_argument("--max-eta-squared", type=float, default=0.25)
    parser.add_argument("--max-dominant-sample-fraction", type=float, default=0.80)
    return parser.parse_args()


def split_csv(value: str) -> list[str]:
    return [x.strip() for x in value.split(",") if x.strip()]


def natural_cluster_order(values: Iterable[str]) -> list[str]:
    def key(value: str) -> tuple[int, str]:
        try:
            return (0, f"{int(float(str(value))):05d}")
        except ValueError:
            return (1, str(value))

    return sorted({str(v) for v in values}, key=key)


def to_numeric_obs(adata: ad.AnnData) -> None:
    for col in NUMERIC_COVARIATES:
        if col in adata.obs:
            adata.obs[col] = pd.to_numeric(adata.obs[col].astype(str), errors="coerce")


def add_cell_line(adata: ad.AnnData, sample_map_path: Path) -> None:
    if not sample_map_path.exists() or "orig.ident" not in adata.obs:
        adata.obs["cell_line"] = "unknown"
        return
    sample_map = pd.read_csv(sample_map_path, sep="\t")
    div90 = sample_map.loc[sample_map["DIV"].astype(str).str.upper() == "DIV90"].copy()
    lookup = div90.set_index("run_sample_id")["biological_label"].astype(str).to_dict()
    labels = adata.obs["orig.ident"].astype(str).map(lookup).fillna("unknown")
    adata.obs["cell_line"] = labels.str.replace(r"_rep\d+$", "", regex=True)


def cramers_v(labels: pd.Series, covariate: pd.Series) -> float:
    df = pd.DataFrame({"label": labels.astype(str), "covariate": covariate.astype(str)})
    table = pd.crosstab(df["label"], df["covariate"])
    if table.shape[0] < 2 or table.shape[1] < 2:
        return 0.0
    chi2 = chi2_contingency(table, correction=False)[0]
    n = table.to_numpy().sum()
    denom = n * max(1, min(table.shape[0] - 1, table.shape[1] - 1))
    return float(np.sqrt(chi2 / denom)) if denom else 0.0


def dominant_fraction(labels: pd.Series, covariate: pd.Series) -> float:
    df = pd.DataFrame({"label": labels.astype(str), "covariate": covariate.astype(str)})
    max_frac = 0.0
    for _, group in df.groupby("label"):
        counts = group["covariate"].value_counts(normalize=True, dropna=False)
        if not counts.empty:
            max_frac = max(max_frac, float(counts.iloc[0]))
    return max_frac


def eta_squared(labels: pd.Series, values: pd.Series) -> float:
    numeric = pd.to_numeric(values, errors="coerce")
    df = pd.DataFrame({"label": labels.astype(str), "value": numeric}).dropna()
    if df["label"].nunique() < 2 or df.shape[0] < 3:
        return 0.0
    grand = df["value"].mean()
    ss_total = float(((df["value"] - grand) ** 2).sum())
    if ss_total <= 0:
        return 0.0
    ss_between = 0.0
    for _, group in df.groupby("label"):
        ss_between += float(group.shape[0] * (group["value"].mean() - grand) ** 2)
    return ss_between / ss_total


def marker_counts(marker_df: pd.DataFrame, logfc_min: float, padj_max: float) -> pd.DataFrame:
    if marker_df.empty:
        return pd.DataFrame(columns=["subcluster", "n_enriched_markers", "top_marker", "top_marker_logfc", "top_marker_padj"])
    df = marker_df.copy()
    df["logfoldchanges"] = pd.to_numeric(df.get("logfoldchanges", np.nan), errors="coerce")
    df["pvals_adj"] = pd.to_numeric(df.get("pvals_adj", np.nan), errors="coerce")
    keep = df.loc[(df["logfoldchanges"] >= logfc_min) & (df["pvals_adj"] <= padj_max)].copy()
    rows = []
    for group, sub in df.groupby("group", sort=False):
        enriched = keep.loc[keep["group"].astype(str) == str(group)].sort_values(["pvals_adj", "logfoldchanges"], ascending=[True, False])
        top = enriched.iloc[0] if not enriched.empty else sub.sort_values("pvals_adj").iloc[0]
        rows.append(
            {
                "subcluster": str(group),
                "n_enriched_markers": int(enriched.shape[0]),
                "top_marker": str(top["names"]),
                "top_marker_logfc": float(top["logfoldchanges"]) if pd.notna(top["logfoldchanges"]) else np.nan,
                "top_marker_padj": float(top["pvals_adj"]) if pd.notna(top["pvals_adj"]) else np.nan,
            }
        )
    return pd.DataFrame(rows)


def save_umap(adata: ad.AnnData, color: list[str], path: Path, title: str) -> None:
    n = len(color)
    fig, axes = plt.subplots(1, n, figsize=(4.2 * n, 4), squeeze=False)
    for ax, col in zip(axes[0], color):
        sc.pl.umap(adata, color=col, ax=ax, show=False, frameon=False, size=10, title=col)
    fig.suptitle(title)
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=220)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def write_tsv_gz(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt") as handle:
        df.to_csv(handle, sep="\t", index=False)


def recommendation_for_parent(summary: pd.DataFrame, args: argparse.Namespace) -> dict[str, object]:
    parent = summary.iloc[0]["parent_cluster_id"]
    viable = summary.loc[
        (summary["n_subclusters"] >= 2)
        & (summary["min_subcluster_cells"] >= args.min_subcluster_cells)
        & (summary["marker_supported_fraction"] >= 0.75)
        & (~summary["technical_confounding_flag"])
    ].copy()
    if viable.empty:
        best = summary.sort_values(
            ["marker_supported_fraction", "n_subclusters", "max_technical_effect"],
            ascending=[False, True, True],
        ).iloc[0]
        return {
            "parent_cluster_id": parent,
            "recommendation": "keep_parent_only",
            "recommended_resolution": "",
            "recommended_n_subclusters": 1,
            "reason": (
                "No swept resolution had enough marker-supported subclusters while staying below "
                "technical/sample/cell-line/cell-cycle confounding thresholds."
            ),
            "best_tested_resolution": best["resolution"],
            "best_tested_n_subclusters": best["n_subclusters"],
            "best_marker_supported_fraction": best["marker_supported_fraction"],
            "best_max_technical_effect": best["max_technical_effect"],
        }
    chosen = viable.sort_values(["n_subclusters", "resolution"]).iloc[0]
    return {
        "parent_cluster_id": parent,
        "recommendation": "candidate_subcluster",
        "recommended_resolution": chosen["resolution"],
        "recommended_n_subclusters": chosen["n_subclusters"],
        "reason": "Lowest resolution passing marker-support and technical-confounding filters.",
        "best_tested_resolution": chosen["resolution"],
        "best_tested_n_subclusters": chosen["n_subclusters"],
        "best_marker_supported_fraction": chosen["marker_supported_fraction"],
        "best_max_technical_effect": chosen["max_technical_effect"],
    }


def run_parent(parent: str, adata: ad.AnnData, args: argparse.Namespace, outdir: Path) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    parent_mask = adata.obs["cluster_id"].astype(str) == str(parent)
    sub = adata[parent_mask].copy()
    sub.layers["log_expr"] = sub.X.copy()
    parent_name = str(sub.obs["cluster_number_name"].iloc[0]) if "cluster_number_name" in sub.obs else str(parent)
    if sub.n_obs < args.min_parent_cells:
        row = {
            "parent_cluster_id": parent,
            "parent_cluster_name": parent_name,
            "n_parent_cells": sub.n_obs,
            "resolution": "",
            "n_subclusters": 1,
            "min_subcluster_cells": sub.n_obs,
            "marker_supported_fraction": 0.0,
            "max_technical_effect": 0.0,
            "technical_confounding_flag": False,
            "skip_reason": f"parent has fewer than {args.min_parent_cells} cells",
        }
        return pd.DataFrame([row]), pd.DataFrame(), pd.DataFrame()

    sc.pp.highly_variable_genes(sub, n_top_genes=min(args.n_top_genes, sub.n_vars), flavor="seurat")
    sub = sub[:, sub.var["highly_variable"].to_numpy()].copy()
    sc.pp.scale(sub, max_value=10)
    sc.tl.pca(sub, n_comps=min(args.n_pcs, sub.n_obs - 2, sub.n_vars - 1), svd_solver="arpack")
    sc.pp.neighbors(sub, n_neighbors=min(args.n_neighbors, max(5, sub.n_obs - 1)), n_pcs=min(args.n_pcs, sub.obsm["X_pca"].shape[1]))
    sc.tl.umap(sub, random_state=0)

    resolutions = [float(x) for x in split_csv(args.resolutions)]
    summary_rows = []
    marker_rows = []
    top_marker_rows = []
    previous_labels = None
    previous_resolution = None

    for res in resolutions:
        key = f"subcluster_r{str(res).replace('.', '_')}"
        sc.tl.leiden(sub, resolution=res, key_added=key, random_state=0)
        labels = sub.obs[key].astype(str)
        counts = labels.value_counts()
        n_subclusters = int(counts.shape[0])
        min_cells = int(counts.min())
        adjacent_ari = np.nan if previous_labels is None else adjusted_rand_score(previous_labels, labels)

        if n_subclusters >= 2:
            sc.tl.rank_genes_groups(sub, groupby=key, method="wilcoxon", pts=True, layer="log_expr")
            markers = sc.get.rank_genes_groups_df(sub, group=None)
            markers["parent_cluster_id"] = parent
            markers["parent_cluster_name"] = parent_name
            markers["resolution"] = res
            markers["subcluster_key"] = key
            marker_count = marker_counts(markers, args.marker_logfc, args.marker_padj)
        else:
            markers = pd.DataFrame()
            marker_count = pd.DataFrame(
                [{"subcluster": "0", "n_enriched_markers": 0, "top_marker": "", "top_marker_logfc": np.nan, "top_marker_padj": np.nan}]
            )

        marker_count["parent_cluster_id"] = parent
        marker_count["parent_cluster_name"] = parent_name
        marker_count["resolution"] = res
        marker_count["subcluster_key"] = key
        marker_count["subcluster_cells"] = marker_count["subcluster"].map(counts.astype(int).to_dict()).fillna(0).astype(int)
        marker_rows.append(marker_count)
        if not markers.empty:
            top_marker_rows.append(markers.groupby("group", group_keys=False).head(25))

        supported = marker_count["n_enriched_markers"] >= args.min_markers_per_subcluster
        marker_supported_fraction = float(supported.mean()) if marker_count.shape[0] else 0.0

        categorical_effects = {}
        dominant_effects = {}
        for cov in CATEGORICAL_COVARIATES:
            if cov in sub.obs:
                categorical_effects[f"{cov}_cramers_v"] = cramers_v(labels, sub.obs[cov])
                dominant_effects[f"{cov}_max_dominant_fraction"] = dominant_fraction(labels, sub.obs[cov])
        numeric_effects = {}
        for cov in NUMERIC_COVARIATES:
            if cov in sub.obs:
                numeric_effects[f"{cov}_eta_squared"] = eta_squared(labels, sub.obs[cov])

        max_cramer = max(categorical_effects.values()) if categorical_effects else 0.0
        max_eta = max(numeric_effects.values()) if numeric_effects else 0.0
        max_dom_sample = dominant_effects.get("orig.ident_max_dominant_fraction", 0.0)
        max_dom_cell_line = dominant_effects.get("cell_line_max_dominant_fraction", 0.0)
        technical_flag = (
            max_cramer >= args.max_cramers_v
            or max_eta >= args.max_eta_squared
            or max_dom_sample >= args.max_dominant_sample_fraction
            or max_dom_cell_line >= args.max_dominant_sample_fraction
        )
        max_technical_effect = max(max_cramer, max_eta, max_dom_sample, max_dom_cell_line)

        summary_rows.append(
            {
                "parent_cluster_id": parent,
                "parent_cluster_name": parent_name,
                "n_parent_cells": sub.n_obs,
                "resolution": res,
                "subcluster_key": key,
                "n_subclusters": n_subclusters,
                "min_subcluster_cells": min_cells,
                "marker_supported_subclusters": int(supported.sum()),
                "marker_supported_fraction": marker_supported_fraction,
                "adjacent_previous_resolution": previous_resolution if previous_resolution is not None else "",
                "adjacent_ari": adjacent_ari,
                "max_cramers_v": max_cramer,
                "max_eta_squared": max_eta,
                "max_dominant_sample_fraction": max_dom_sample,
                "max_dominant_cell_line_fraction": max_dom_cell_line,
                "max_technical_effect": max_technical_effect,
                "technical_confounding_flag": bool(technical_flag),
                "skip_reason": "",
                **categorical_effects,
                **dominant_effects,
                **numeric_effects,
            }
        )
        previous_labels = labels.copy()
        previous_resolution = res

    summary = pd.DataFrame(summary_rows)
    rec = recommendation_for_parent(summary, args)
    chosen_key = None
    if rec["recommendation"] == "candidate_subcluster":
        chosen_res = float(rec["recommended_resolution"])
        chosen_key = f"subcluster_r{str(chosen_res).replace('.', '_')}"
    else:
        best_res = summary.sort_values(["marker_supported_fraction", "n_subclusters"], ascending=[False, True]).iloc[0]["resolution"]
        chosen_key = f"subcluster_r{str(float(best_res)).replace('.', '_')}"

    plot_cols = [chosen_key]
    for cov in ["orig.ident", "cell_line", "Phase", "nCount_RNA", "nFeature_RNA", "percent.mt", "S.Score", "G2M.Score"]:
        if cov in sub.obs:
            plot_cols.append(cov)
    safe_parent = str(parent).replace("/", "_")
    save_umap(
        sub,
        color=plot_cols[:6],
        path=outdir / "plots" / f"parent_cluster_{safe_parent}_local_umap_audit.png",
        title=f"Parent cluster {parent}: {parent_name}",
    )

    marker_table = pd.concat(marker_rows, ignore_index=True) if marker_rows else pd.DataFrame()
    top_markers = pd.concat(top_marker_rows, ignore_index=True) if top_marker_rows else pd.DataFrame()
    rec_df = pd.DataFrame([rec])
    return summary, marker_table, top_markers.assign(parent_cluster_id=parent, parent_cluster_name=parent_name) if not top_markers.empty else top_markers


def main() -> int:
    args = parse_args()
    project_root = Path(args.project_root).expanduser().resolve()
    h5ad_path = Path(args.h5ad).expanduser().resolve() if args.h5ad else project_root / "results/python_anndata/varela_div90.h5ad"
    outdir = (
        Path(args.outdir).expanduser().resolve()
        if args.outdir
        else project_root / "results/div90_parent_cluster_subclustering_audit/div90_parent_cluster_subclustering_audit_v1"
    )
    table_dir = outdir / "tables"
    table_dir.mkdir(parents=True, exist_ok=True)
    sample_map = Path(args.sample_map)
    if not sample_map.is_absolute():
        sample_map = Path.cwd() / sample_map

    print(f"Reading {h5ad_path}")
    adata = sc.read_h5ad(h5ad_path)
    to_numeric_obs(adata)
    add_cell_line(adata, sample_map)
    adata.obs["cluster_id"] = adata.obs["cluster_id"].astype(str)
    exclude = set(split_csv(args.exclude_clusters))
    parents = [cluster for cluster in natural_cluster_order(adata.obs["cluster_id"]) if cluster not in exclude]
    if args.parent_clusters.lower() != "all":
        requested = set(split_csv(args.parent_clusters))
        parents = [cluster for cluster in parents if cluster in requested]
    if not parents:
        raise SystemExit("No parent clusters selected after applying --parent-clusters and --exclude-clusters")

    all_summary = []
    all_marker_counts = []
    all_top_markers = []
    all_recs = []
    for parent in parents:
        print(f"Auditing parent cluster {parent}")
        summary, marker_counts_df, top_markers = run_parent(parent, adata, args, outdir)
        rec = recommendation_for_parent(summary, args)
        all_summary.append(summary)
        all_recs.append(pd.DataFrame([rec]))
        if not marker_counts_df.empty:
            all_marker_counts.append(marker_counts_df)
        if not top_markers.empty:
            all_top_markers.append(top_markers)

    summary_df = pd.concat(all_summary, ignore_index=True)
    rec_df = pd.concat(all_recs, ignore_index=True)
    marker_counts_df = pd.concat(all_marker_counts, ignore_index=True) if all_marker_counts else pd.DataFrame()
    top_markers_df = pd.concat(all_top_markers, ignore_index=True) if all_top_markers else pd.DataFrame()

    summary_df.to_csv(table_dir / "div90_parent_cluster_subclustering_resolution_summary.tsv", sep="\t", index=False)
    rec_df.to_csv(table_dir / "div90_parent_cluster_subclustering_recommendations.tsv", sep="\t", index=False)
    marker_counts_df.to_csv(table_dir / "div90_parent_cluster_subclustering_marker_support.tsv", sep="\t", index=False)
    write_tsv_gz(top_markers_df, table_dir / "div90_parent_cluster_subclustering_top_markers.tsv.gz")

    report = [
        "# DIV90 Parent Cluster Subclustering Audit",
        "",
        "Local run from cached DIV90 H5AD. Dividing cells were excluded from the parent-cluster audit.",
        "",
        "## Recommendation Counts",
        "",
        rec_df["recommendation"].value_counts().to_string(),
        "",
        "## Candidate Subcluster Resolutions",
        "",
    ]
    candidates = rec_df.loc[rec_df["recommendation"] == "candidate_subcluster"].copy()
    if candidates.empty:
        report.append("No parent cluster passed the marker-support and confounding filters.")
    else:
        for row in candidates.itertuples():
            report.append(
                f"- Parent cluster `{row.parent_cluster_id}`: resolution `{row.recommended_resolution}`, "
                f"{row.recommended_n_subclusters} subclusters; marker-supported fraction "
                f"{row.best_marker_supported_fraction:.2f}; max technical effect {row.best_max_technical_effect:.2f}."
            )
    report.extend(["", "## Keep Parent Only", ""])
    for row in rec_df.loc[rec_df["recommendation"] == "keep_parent_only"].itertuples():
        report.append(f"- Parent cluster `{row.parent_cluster_id}`: {row.reason}")
    (outdir / "div90_parent_cluster_subclustering_audit_report.md").write_text("\n".join(report) + "\n")
    print(f"Wrote: {outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
