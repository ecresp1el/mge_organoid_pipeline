#!/usr/bin/env python
"""Export focused MGEO RGC/IPC inputs for a Monocle3 Slurm stage."""

import argparse
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmwrite


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Export focused MGEO RGC/IPC matrix and metadata for Monocle3."
    )
    parser.add_argument("--project-root", required=True, help="Runtime PROJECT_ROOT.")
    parser.add_argument("--repo-root", required=True, help="Repository root.")
    parser.add_argument(
        "--outdir",
        default=None,
        help="Output directory. Defaults to PROJECT_ROOT/results/mgeo_rgc_ipc_monocle3/inputs.",
    )
    parser.add_argument(
        "--activity-quantile",
        type=float,
        default=0.30,
        help="Keep cells with max(RGC, IPC) score >= this quantile.",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Rewrite existing exported inputs.",
    )
    return parser.parse_args()


def clean_gene_name(gene):
    if pd.isna(gene):
        return None
    gene = str(gene).strip()
    if gene == "" or gene.lower() in {"nan", "none", "na"}:
        return None
    if "00:00:00" in gene:
        return None
    return gene


def resolve_genes_to_anndata(adata, genes):
    var_names = pd.Index(adata.var_names.astype(str))
    upper_to_actual = {}
    for gene in var_names:
        upper_to_actual.setdefault(gene.upper(), gene)

    found = []
    missing = []
    for gene in genes:
        if gene in var_names:
            found.append(gene)
        elif gene.upper() in upper_to_actual:
            found.append(upper_to_actual[gene.upper()])
        else:
            missing.append(gene)
    return list(dict.fromkeys(found)), missing


def score_gene_set_mean_expression(adata, genes, score_name):
    if len(genes) == 0:
        raise ValueError(f"No genes available for score {score_name}")
    x = adata[:, genes].X
    if sparse.issparse(x):
        scores = np.asarray(x.mean(axis=1)).ravel()
    else:
        scores = np.asarray(x).mean(axis=1)
    adata.obs[score_name] = scores.astype(np.float32)


def add_div_metadata(adata, div_label):
    out = adata.copy()
    out.obs["DIV"] = div_label
    out.obs["source_h5ad"] = div_label
    out.obs["original_cell_id"] = out.obs_names.astype(str)
    return out


def main() -> int:
    args = parse_args()
    repo_root = Path(args.repo_root).expanduser().resolve()
    project_root = Path(args.project_root).expanduser().resolve()
    src = repo_root / "python_notebooks" / "src"
    sys.path.insert(0, str(src))

    from mge_organoid_python import default_studies, load_cached_anndata

    outdir = (
        Path(args.outdir).expanduser().resolve()
        if args.outdir
        else project_root / "results" / "mgeo_rgc_ipc_monocle3" / "inputs"
    )
    outdir.mkdir(parents=True, exist_ok=True)

    matrix_path = outdir / "mgeo_rgc_ipc_expression_genes_by_cells.mtx"
    cell_metadata_path = outdir / "mgeo_rgc_ipc_cell_metadata.csv"
    gene_metadata_path = outdir / "mgeo_rgc_ipc_gene_metadata.csv"
    full_metadata_path = outdir / "mgeo_full_scored_metadata.csv"
    manifest_path = outdir / "mgeo_rgc_ipc_input_manifest.tsv"

    needed = [matrix_path, cell_metadata_path, gene_metadata_path, full_metadata_path]
    if all(path.exists() for path in needed) and not args.force:
        print("All Monocle3 input files already exist. Use --force to rewrite.")
        for path in needed:
            print(path)
        return 0

    studies = {study.study_id: study for study in default_studies()}
    div30, _ = load_cached_anndata(studies["varela_div30"], project_root=project_root, backed=None)
    div90, _ = load_cached_anndata(studies["varela_div90"], project_root=project_root, backed=None)

    div30 = add_div_metadata(div30, "DIV30")
    div90 = add_div_metadata(div90, "DIV90")
    mgeo = ad.concat(
        [div30, div90],
        axis=0,
        join="inner",
        label="DIV_source",
        keys=["DIV30", "DIV90"],
        index_unique="-",
    )
    if sparse.issparse(mgeo.X):
        mgeo.X = mgeo.X.astype(np.float32)
    else:
        mgeo.X = np.asarray(mgeo.X, dtype=np.float32)
    mgeo.obs["DIV"] = pd.Categorical(
        mgeo.obs["DIV"].astype(str),
        categories=["DIV30", "DIV90"],
        ordered=True,
    )

    table_s5_path = (
        project_root
        / "reference"
        / "shi_2021_tables_s2_to_s9"
        / "science.abj6641_table_s5.xlsx"
    )
    if not table_s5_path.exists():
        raise FileNotFoundError(f"Missing Shi Table S5: {table_s5_path}")

    s5 = pd.read_excel(table_s5_path, sheet_name="RGC IPC DEGs", header=1)
    s5.columns = ["p_val", "avg_logFC", "pct.1", "pct.2", "p_val_adj", "cluster", "gene"]
    s5["gene"] = s5["gene"].map(clean_gene_name)
    s5 = s5.dropna(subset=["gene", "cluster"]).copy()
    s5["avg_logFC"] = pd.to_numeric(s5["avg_logFC"], errors="coerce")
    s5["p_val_adj"] = pd.to_numeric(s5["p_val_adj"], errors="coerce")
    s5_sig = s5[(s5["p_val_adj"] < 0.05) & (s5["avg_logFC"] > 0)].copy()

    rgc_degs_published = (
        s5_sig[s5_sig["cluster"] == "RGC"]
        .sort_values("avg_logFC", ascending=False)
        .drop_duplicates("gene")["gene"]
        .tolist()
    )
    ipc_degs_published = (
        s5_sig[s5_sig["cluster"] == "IPC"]
        .sort_values("avg_logFC", ascending=False)
        .drop_duplicates("gene")["gene"]
        .tolist()
    )

    rgc_degs, missing_rgc = resolve_genes_to_anndata(mgeo, rgc_degs_published)
    ipc_degs, missing_ipc = resolve_genes_to_anndata(mgeo, ipc_degs_published)
    score_gene_set_mean_expression(mgeo, rgc_degs, "shi_s5_RGC_score")
    score_gene_set_mean_expression(mgeo, ipc_degs, "shi_s5_IPC_score")
    mgeo.obs["shi_s5_IPC_minus_RGC_score"] = (
        mgeo.obs["shi_s5_IPC_score"] - mgeo.obs["shi_s5_RGC_score"]
    )

    rgc_ipc_activity = mgeo.obs[["shi_s5_RGC_score", "shi_s5_IPC_score"]].max(axis=1)
    activity_cutoff = rgc_ipc_activity.quantile(args.activity_quantile)
    subset_mask = (rgc_ipc_activity >= activity_cutoff) & rgc_ipc_activity.notna()
    mgeo_rgc_ipc = mgeo[subset_mask.to_numpy(), :].copy()

    root_df = mgeo_rgc_ipc.obs[
        ["DIV", "shi_s5_RGC_score", "shi_s5_IPC_score", "shi_s5_IPC_minus_RGC_score"]
    ].copy()
    root_candidates = root_df[
        (root_df["DIV"].astype(str) == "DIV30")
        & (root_df["shi_s5_RGC_score"] >= root_df["shi_s5_RGC_score"].quantile(0.90))
        & (root_df["shi_s5_IPC_score"] <= root_df["shi_s5_IPC_score"].quantile(0.50))
    ].copy()
    if len(root_candidates) < 25:
        root_candidates = root_df[
            (root_df["DIV"].astype(str) == "DIV30")
            & (root_df["shi_s5_RGC_score"] >= root_df["shi_s5_RGC_score"].quantile(0.80))
        ].copy()
    if len(root_candidates) == 0:
        root_candidates = root_df[
            (root_df["shi_s5_RGC_score"] >= root_df["shi_s5_RGC_score"].quantile(0.95))
            & (root_df["shi_s5_IPC_score"] <= root_df["shi_s5_IPC_score"].quantile(0.50))
        ].copy()
    if len(root_candidates) == 0:
        raise RuntimeError("No RGC-root candidates found.")

    root_candidates["rgc_root_score"] = (
        root_candidates["shi_s5_RGC_score"] - root_candidates["shi_s5_IPC_score"]
    )
    root_seed_cell = root_candidates["rgc_root_score"].idxmax()
    mgeo_rgc_ipc.obs["monocle_root_candidate"] = False
    mgeo_rgc_ipc.obs.loc[root_candidates.index, "monocle_root_candidate"] = True
    mgeo_rgc_ipc.obs["monocle_root_seed_cell"] = mgeo_rgc_ipc.obs_names == root_seed_cell

    expression_genes_by_cells = (
        mgeo_rgc_ipc.X.T.tocsc()
        if sparse.issparse(mgeo_rgc_ipc.X)
        else sparse.csc_matrix(np.asarray(mgeo_rgc_ipc.X).T)
    )
    mmwrite(matrix_path, expression_genes_by_cells)

    cell_metadata = mgeo_rgc_ipc.obs.copy()
    cell_metadata.index.name = "cell_id"
    for col in cell_metadata.columns:
        if isinstance(cell_metadata[col].dtype, pd.CategoricalDtype):
            cell_metadata[col] = cell_metadata[col].astype(str)
    cell_metadata.to_csv(cell_metadata_path)

    gene_metadata = pd.DataFrame(index=mgeo_rgc_ipc.var_names.astype(str))
    gene_metadata.index.name = "gene_id"
    gene_metadata["gene_short_name"] = gene_metadata.index.astype(str)
    gene_metadata.to_csv(gene_metadata_path)

    full_metadata = mgeo.obs[
        [
            "DIV",
            "original_cell_id",
            "shi_s5_RGC_score",
            "shi_s5_IPC_score",
            "shi_s5_IPC_minus_RGC_score",
        ]
    ].copy()
    full_metadata["in_monocle3_rgc_ipc_subset"] = full_metadata.index.isin(mgeo_rgc_ipc.obs_names)
    full_metadata.index.name = "cell_id"
    full_metadata.to_csv(full_metadata_path)

    manifest = pd.DataFrame(
        [
            ("matrix_genes_by_cells", matrix_path),
            ("cell_metadata", cell_metadata_path),
            ("gene_metadata", gene_metadata_path),
            ("full_scored_metadata", full_metadata_path),
            ("shi_table_s5", table_s5_path),
            ("root_seed_cell", root_seed_cell),
            ("activity_quantile", args.activity_quantile),
            ("activity_cutoff", activity_cutoff),
            ("n_full_cells", mgeo.n_obs),
            ("n_subset_cells", mgeo_rgc_ipc.n_obs),
            ("n_genes", mgeo_rgc_ipc.n_vars),
            ("n_rgc_genes_found", len(rgc_degs)),
            ("n_ipc_genes_found", len(ipc_degs)),
            ("n_rgc_genes_missing", len(missing_rgc)),
            ("n_ipc_genes_missing", len(missing_ipc)),
        ],
        columns=["key", "value"],
    )
    manifest.to_csv(manifest_path, sep="\t", index=False)

    print("Wrote Monocle3 inputs:")
    print("  matrix:", matrix_path)
    print("  cell metadata:", cell_metadata_path)
    print("  gene metadata:", gene_metadata_path)
    print("  full metadata:", full_metadata_path)
    print("  manifest:", manifest_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
