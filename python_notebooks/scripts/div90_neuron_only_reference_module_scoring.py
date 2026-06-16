#!/usr/bin/env python3
"""Score Siletti/Data S9 MGE divergence modules in DIV90 neuron-only clusters."""

from __future__ import annotations

import argparse
import os
from pathlib import Path

import anndata as ad
import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from scipy import sparse


PROJECT_ROOT_DEFAULT = "/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder"
DEFAULT_EXCEL = "reference/science.adw1803_data_s9.xlsx"
NEURON_PARENT_CLUSTERS = ["0", "1", "2", "3", "5", "8", "11"]
NEURON_RESOLUTION = 0.6
NEURON_KEY = "neuron_leiden_r0_6"

ANCHOR_CLASSES = {
    "class1_EPHA5_MEF2C": ["EPHA5", "MEF2C"],
    "class2_LHX6_NFIA": ["LHX6", "NFIA"],
    "class3_CRABP1_ANGPT2": ["CRABP1", "ANGPT2"],
    "class4_NR2F1_NR2F2": ["NR2F1", "NR2F2"],
    "class5_LHX8_ISL1": ["LHX8", "ISL1"],
}

ANCHOR_LABELS = {
    "class1_EPHA5_MEF2C": "EPHA5/MEF2C",
    "class2_LHX6_NFIA": "LHX6/NFIA",
    "class3_CRABP1_ANGPT2": "CRABP1/ANGPT2",
    "class4_NR2F1_NR2F2": "NR2F1/NR2F2",
    "class5_LHX8_ISL1": "LHX8/ISL1",
}

NUMERIC_OBS = ["nCount_RNA", "nFeature_RNA", "percent.mt", "S.Score", "G2M.Score", "CC.Difference"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", default=os.environ.get("PROJECT_ROOT", PROJECT_ROOT_DEFAULT))
    parser.add_argument("--h5ad", default=None)
    parser.add_argument("--excel", default=DEFAULT_EXCEL)
    parser.add_argument("--sample-map", default="metadata/div30_div90_sample_id_to_biolabel_map.tsv")
    parser.add_argument(
        "--outdir",
        default=None,
        help="Default: PROJECT_ROOT/results/div90_neuron_only_reference_modules/div90_neuron_only_reference_modules_v1",
    )
    parser.add_argument("--n-top-genes", type=int, default=2500)
    parser.add_argument("--n-pcs", type=int, default=40)
    parser.add_argument("--n-neighbors", type=int, default=30)
    return parser.parse_args()


def to_numeric_obs(adata: ad.AnnData) -> None:
    for col in NUMERIC_OBS:
        if col in adata.obs:
            adata.obs[col] = pd.to_numeric(adata.obs[col].astype(str), errors="coerce")


def add_cell_line(adata: ad.AnnData, sample_map_path: Path) -> None:
    if "cell_line" in adata.obs:
        return
    if not sample_map_path.exists() or "orig.ident" not in adata.obs:
        adata.obs["cell_line"] = "unknown"
        return
    sample_map = pd.read_csv(sample_map_path, sep="\t")
    div90 = sample_map.loc[sample_map["DIV"].astype(str).str.upper() == "DIV90"].copy()
    lookup = div90.set_index("run_sample_id")["biological_label"].astype(str).to_dict()
    labels = adata.obs["orig.ident"].astype(str).map(lookup).fillna("unknown")
    adata.obs["cell_line"] = labels.str.replace(r"_rep\d+$", "", regex=True)


def resolution_key(resolution: float) -> str:
    return f"neuron_leiden_r{str(resolution).replace('.', '_')}"


def read_excel_modules(path: Path) -> dict[str, list[str]]:
    xl = pd.ExcelFile(path)
    modules: dict[str, list[str]] = {}
    for sheet in xl.sheet_names:
        df = pd.read_excel(path, sheet_name=sheet)
        if "gene" not in df.columns:
            continue
        genes = (
            df["gene"]
            .dropna()
            .astype(str)
            .str.strip()
            .loc[lambda s: s.ne("")]
            .drop_duplicates()
            .tolist()
        )
        if not genes:
            continue
        safe = sheet.replace("（", "_").replace("）", "").replace(" ", "_").lower()
        safe = safe.replace("__", "_")
        modules[f"s9_{safe}"] = genes
    return modules


def matrix_values(adata: ad.AnnData, genes: list[str]) -> tuple[np.ndarray, list[str]]:
    present = [gene for gene in genes if gene in adata.var_names]
    if not present:
        return np.zeros((adata.n_obs, 0), dtype=float), []
    mat = adata[:, present].X
    if sparse.issparse(mat):
        mat = mat.toarray()
    return np.asarray(mat, dtype=float), present


def add_module_score(adata: ad.AnnData, genes: list[str], name: str) -> list[str]:
    values, present = matrix_values(adata, genes)
    adata.obs[name] = np.nan if values.shape[1] == 0 else np.nanmean(values, axis=1)
    return present


def zscore_columns(df: pd.DataFrame, cols: list[str]) -> pd.DataFrame:
    out = df.copy()
    for col in cols:
        vals = pd.to_numeric(out[col], errors="coerce")
        sd = vals.std(ddof=0)
        out[col] = 0.0 if not np.isfinite(sd) or sd == 0 else (vals - vals.mean()) / sd
    return out


def run_neuron_reclustering(adata: ad.AnnData, args: argparse.Namespace) -> ad.AnnData:
    neuron = adata[adata.obs["cluster_id"].astype(str).isin(NEURON_PARENT_CLUSTERS)].copy()
    neuron.obs["parent_cluster"] = neuron.obs["cluster_id"].astype(str) + " - " + neuron.obs["cluster_number_name"].astype(str)
    work = neuron.copy()
    sc.pp.highly_variable_genes(work, n_top_genes=min(args.n_top_genes, work.n_vars), flavor="seurat")
    work = work[:, work.var["highly_variable"].to_numpy()].copy()
    sc.pp.scale(work, max_value=10)
    sc.tl.pca(work, n_comps=min(args.n_pcs, work.n_obs - 2, work.n_vars - 1), svd_solver="arpack")
    sc.pp.neighbors(work, n_neighbors=min(args.n_neighbors, work.n_obs - 1), n_pcs=min(args.n_pcs, work.obsm["X_pca"].shape[1]))
    sc.tl.umap(work, random_state=0)
    sc.tl.leiden(work, resolution=NEURON_RESOLUTION, key_added=NEURON_KEY, random_state=0)
    neuron.obsm["X_umap"] = work.obsm["X_umap"].copy()
    neuron.obs[NEURON_KEY] = work.obs[NEURON_KEY].astype(str).to_numpy()
    return neuron


def gene_set_table(adata: ad.AnnData, modules: dict[str, list[str]]) -> pd.DataFrame:
    rows = []
    var_names = set(map(str, adata.var_names))
    for module, genes in modules.items():
        present = [gene for gene in genes if gene in var_names]
        missing = [gene for gene in genes if gene not in var_names]
        rows.append(
            {
                "module": module,
                "n_genes_total": len(genes),
                "n_genes_present": len(present),
                "fraction_present": len(present) / len(genes) if genes else np.nan,
                "present_genes": ",".join(present),
                "missing_genes": ",".join(missing),
            }
        )
    return pd.DataFrame(rows)


def summarize_scores(adata: ad.AnnData, group_col: str, score_cols: list[str], zscore_suffix: str = "_z") -> pd.DataFrame:
    rows = []
    obs = adata.obs[[group_col] + score_cols].copy()
    for group, df in obs.groupby(group_col, observed=False):
        row = {group_col: str(group), "n_cells": int(df.shape[0])}
        for col in score_cols:
            row[f"mean_{col}"] = float(pd.to_numeric(df[col], errors="coerce").mean())
            row[f"median_{col}"] = float(pd.to_numeric(df[col], errors="coerce").median())
        rows.append(row)
    summary = pd.DataFrame(rows)
    mean_cols = [col for col in summary.columns if col.startswith("mean_")]
    z = zscore_columns(summary, mean_cols)
    for col in mean_cols:
        summary[f"{col}{zscore_suffix}"] = z[col]
    return summary


def best_module_calls(summary: pd.DataFrame, group_col: str, raw_score_cols: list[str], prefix: str) -> pd.DataFrame:
    rows = []
    z_cols = [f"mean_{col}_z" for col in raw_score_cols]
    for _, row in summary.iterrows():
        vals = pd.Series({col: row[col] for col in z_cols}, dtype=float).sort_values(ascending=False)
        best = vals.index[0]
        second = vals.index[1] if len(vals) > 1 else None
        best_module = best.removeprefix("mean_").removesuffix("_z")
        second_module = second.removeprefix("mean_").removesuffix("_z") if second else ""
        rows.append(
            {
                group_col: row[group_col],
                f"best_{prefix}_module": best_module,
                f"best_{prefix}_module_label": ANCHOR_LABELS.get(best_module, best_module),
                f"best_{prefix}_z": float(vals.iloc[0]),
                f"second_{prefix}_module": second_module,
                f"second_{prefix}_module_label": ANCHOR_LABELS.get(second_module, second_module),
                f"second_{prefix}_z": float(vals.iloc[1]) if len(vals) > 1 else np.nan,
                f"{prefix}_margin_z": float(vals.iloc[0] - vals.iloc[1]) if len(vals) > 1 else np.nan,
            }
        )
    return pd.DataFrame(rows)


def plot_heatmap(summary: pd.DataFrame, group_col: str, score_cols: list[str], path: Path, title: str) -> None:
    z_cols = [f"mean_{col}_z" for col in score_cols]
    mat = summary.set_index(group_col)[z_cols].copy()
    mat.columns = [ANCHOR_LABELS.get(col.removeprefix("mean_").removesuffix("_z"), col.removeprefix("mean_").removesuffix("_z")) for col in mat.columns]
    mat = mat.sort_index(key=lambda idx: idx.astype(str).map(lambda x: int(x) if x.isdigit() else x))
    fig, ax = plt.subplots(figsize=(max(7, 0.6 * mat.shape[1] + 3), max(4, 0.42 * mat.shape[0] + 2)))
    sns.heatmap(mat, cmap="vlag", center=0, linewidths=0.4, linecolor="white", ax=ax)
    ax.set_title(title)
    ax.set_xlabel("")
    ax.set_ylabel(group_col)
    fig.tight_layout()
    fig.savefig(path, dpi=240)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def save_umap(adata: ad.AnnData, colors: list[str], path: Path, title: str, ncols: int = 3) -> None:
    fig = sc.pl.umap(adata, color=colors, ncols=ncols, show=False, frameon=False, size=8, return_fig=True)
    fig.suptitle(title, y=1.01)
    fig.savefig(path, dpi=240, bbox_inches="tight")
    fig.savefig(path.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def write_report(
    outdir: Path,
    gene_sets: pd.DataFrame,
    cluster_summary: pd.DataFrame,
    anchor_calls: pd.DataFrame,
    excel_calls: pd.DataFrame,
) -> None:
    lines = [
        "# DIV90 Neuron-Only S9 Reference Module Scoring",
        "",
        "This run scores the Science Data S9 five MGE divergence modules and TF-only module sheets",
        "against the DIV90 neuron-only reclustering at resolution 0.6.",
        "",
        "It also scores the five anchor classes provided by EC:",
        "",
        "```text",
        "EPHA5/MEF2C",
        "LHX6/NFIA",
        "CRABP1/ANGPT2",
        "NR2F1/NR2F2",
        "LHX8/ISL1",
        "```",
        "",
        "## Gene Set Coverage",
        "",
        simple_markdown_table(gene_sets[["module", "n_genes_total", "n_genes_present", "fraction_present"]]),
        "",
        "## Best Anchor Class By Neuron-Only Cluster",
        "",
        simple_markdown_table(anchor_calls),
        "",
        "## Best Excel Full Module By Neuron-Only Cluster",
        "",
        simple_markdown_table(excel_calls),
        "",
        "## Main Interpretation",
        "",
        "Use `best_anchor_module_label` as the cleanest five-class working call.",
        "Use the S9 full-module and TF-module heatmaps to inspect whether each call is driven by the",
        "whole module or mainly by the small anchor pair.",
        "",
        "Strong confidence means a high positive best z-score and a large margin over the second-best class.",
        "Low margin means the neuron-only cluster is mixed or intermediate between S9 module classes.",
    ]
    (outdir / "div90_neuron_only_s9_reference_module_report.md").write_text("\n".join(lines) + "\n")


def simple_markdown_table(df: pd.DataFrame) -> str:
    if df.empty:
        return "_No rows._"
    display = df.copy()
    for col in display.columns:
        if pd.api.types.is_float_dtype(display[col]):
            display[col] = display[col].map(lambda x: "" if pd.isna(x) else f"{x:.3g}")
        else:
            display[col] = display[col].map(lambda x: "" if pd.isna(x) else str(x))
    cols = [str(col) for col in display.columns]
    lines = ["| " + " | ".join(cols) + " |", "| " + " | ".join(["---"] * len(cols)) + " |"]
    for row in display.astype(str).values.tolist():
        lines.append("| " + " | ".join(row) + " |")
    return "\n".join(lines)


def main() -> int:
    args = parse_args()
    project_root = Path(args.project_root).expanduser().resolve()
    h5ad = Path(args.h5ad).expanduser().resolve() if args.h5ad else project_root / "results/python_anndata/varela_div90.h5ad"
    excel = Path(args.excel).expanduser()
    if not excel.is_absolute():
        excel = project_root / excel
    sample_map = Path(args.sample_map).expanduser()
    if not sample_map.is_absolute():
        sample_map = Path.cwd() / sample_map
    outdir = (
        Path(args.outdir).expanduser().resolve()
        if args.outdir
        else project_root / "results/div90_neuron_only_reference_modules/div90_neuron_only_reference_modules_v1"
    )
    outdir.mkdir(parents=True, exist_ok=True)

    print(f"Reading {h5ad}")
    adata = sc.read_h5ad(h5ad)
    adata.obs["cluster_id"] = adata.obs["cluster_id"].astype(str)
    to_numeric_obs(adata)
    add_cell_line(adata, sample_map)

    print("Recomputing neuron-only resolution 0.6")
    neuron = run_neuron_reclustering(adata, args)

    print(f"Reading modules from {excel}")
    excel_modules = read_excel_modules(excel)
    full_modules = {k: v for k, v in excel_modules.items() if "_tfs" not in k}
    tf_modules = {k: v for k, v in excel_modules.items() if "_tfs" in k}
    all_modules = {**ANCHOR_CLASSES, **full_modules, **tf_modules}

    gene_sets = gene_set_table(neuron, all_modules)
    gene_sets.to_csv(outdir / "div90_neuron_only_s9_gene_set_coverage.tsv", sep="\t", index=False)

    score_cols = []
    for name, genes in all_modules.items():
        present = add_module_score(neuron, genes, name)
        score_cols.append(name)
        print(f"{name}: {len(present)}/{len(genes)} genes present")

    cluster_summary = summarize_scores(neuron, NEURON_KEY, score_cols)
    parent_summary = summarize_scores(neuron, "parent_cluster", score_cols)
    cell_line_summary = summarize_scores(neuron, "cell_line", score_cols)
    sample_summary = summarize_scores(neuron, "orig.ident", score_cols)

    anchor_cols = list(ANCHOR_CLASSES)
    full_cols = list(full_modules)
    tf_cols = list(tf_modules)
    anchor_calls = best_module_calls(cluster_summary, NEURON_KEY, anchor_cols, "anchor")
    excel_calls = best_module_calls(cluster_summary, NEURON_KEY, full_cols, "excel_full")
    tf_calls = best_module_calls(cluster_summary, NEURON_KEY, tf_cols, "excel_tf")
    combined_calls = anchor_calls.merge(excel_calls, on=NEURON_KEY, how="left").merge(tf_calls, on=NEURON_KEY, how="left")

    cluster_summary.to_csv(outdir / "div90_neuron_only_s9_module_scores_by_neuron_cluster.tsv", sep="\t", index=False)
    parent_summary.to_csv(outdir / "div90_neuron_only_s9_module_scores_by_parent_cluster.tsv", sep="\t", index=False)
    cell_line_summary.to_csv(outdir / "div90_neuron_only_s9_module_scores_by_cell_line.tsv", sep="\t", index=False)
    sample_summary.to_csv(outdir / "div90_neuron_only_s9_module_scores_by_sample.tsv", sep="\t", index=False)
    combined_calls.to_csv(outdir / "div90_neuron_only_s9_best_module_calls.tsv", sep="\t", index=False)

    plot_heatmap(cluster_summary, NEURON_KEY, anchor_cols, outdir / "div90_neuron_only_anchor_class_heatmap.png", "Anchor class scores by neuron-only cluster")
    plot_heatmap(cluster_summary, NEURON_KEY, full_cols, outdir / "div90_neuron_only_s9_full_module_heatmap.png", "S9 full module scores by neuron-only cluster")
    plot_heatmap(cluster_summary, NEURON_KEY, tf_cols, outdir / "div90_neuron_only_s9_tf_module_heatmap.png", "S9 TF-module scores by neuron-only cluster")
    plot_heatmap(parent_summary, "parent_cluster", anchor_cols, outdir / "div90_parent_cluster_anchor_class_heatmap.png", "Anchor class scores by original parent cluster")
    plot_heatmap(cell_line_summary, "cell_line", anchor_cols, outdir / "div90_cell_line_anchor_class_heatmap.png", "Anchor class scores by cell line")

    save_umap(neuron, [NEURON_KEY, "parent_cluster", "cell_line"], outdir / "div90_neuron_only_reference_context_umap.png", "DIV90 neuron-only context")
    save_umap(neuron, anchor_cols, outdir / "div90_neuron_only_anchor_class_umaps.png", "Anchor class module scores", ncols=3)
    save_umap(neuron, full_cols, outdir / "div90_neuron_only_s9_full_module_umaps.png", "S9 full module scores", ncols=3)
    save_umap(neuron, tf_cols, outdir / "div90_neuron_only_s9_tf_module_umaps.png", "S9 TF-module scores", ncols=3)

    write_report(outdir, gene_sets, cluster_summary, anchor_calls, excel_calls)
    print(f"Wrote S9 reference module outputs to {outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
