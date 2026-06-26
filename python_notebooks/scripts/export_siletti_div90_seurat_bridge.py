#!/usr/bin/env python3
"""Export Siletti reference and DIV90 query sparse matrices for Seurat transfer."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import io, sparse


PROJECT_ROOT_DEFAULT = "/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder"
ALL_SUPERCLUSTERS = (
    "Upper-layer intratelencephalic",
    "Deep-layer intratelencephalic",
    "Deep-layer near-projecting",
    "Deep-layer corticothalamic and 6b",
    "MGE interneuron",
    "CGE interneuron",
    "LAMP5-LHX6 and Chandelier",
    "Miscellaneous",
    "Hippocampal CA1-3",
    "Hippocampal CA4",
    "Hippocampal dentate gyrus",
    "Amygdala excitatory",
    "Medium spiny neuron",
    "Eccentric medium spiny neuron",
    "Splatter",
    "Mammillary body",
    "Thalamic excitatory",
    "Midbrain-derived inhibitory",
    "Upper rhombic lip",
    "Cerebellar inhibitory",
    "Lower rhombic lip",
    "Oligodendrocyte",
    "Committed oligodendrocyte precursor",
    "Oligodendrocyte precursor",
    "Astrocyte",
    "Ependymal",
    "Microglia",
    "Vascular",
    "Bergmann glia",
    "Fibroblast",
    "Choroid plexus",
)
SCOPES = {
    "all_superclusters": ALL_SUPERCLUSTERS,
    "mge_llc": ("MGE interneuron", "LAMP5-LHX6 and Chandelier"),
    "mge_llc_cholinergic": ("MGE interneuron", "LAMP5-LHX6 and Chandelier", "Splatter"),
    "mge_cge_llc": ("MGE interneuron", "CGE interneuron", "LAMP5-LHX6 and Chandelier"),
    "mge_cge_llc_cholinergic": ("MGE interneuron", "CGE interneuron", "LAMP5-LHX6 and Chandelier", "Splatter"),
    "mge_cge_llc_splatter": ("MGE interneuron", "CGE interneuron", "LAMP5-LHX6 and Chandelier", "Splatter"),
}
H5AD_BY_SUPERCLUSTER = {
    "MGE interneuron": "siletti_whb_mge_interneuron.h5ad",
    "CGE interneuron": "siletti_whb_cge_interneuron.h5ad",
    "LAMP5-LHX6 and Chandelier": "siletti_whb_lamp5_lhx6_and_chandelier.h5ad",
    "Splatter": "siletti_whb_splatter.h5ad",
}
CHOLINERGIC_SPLATTER_CLUSTER_ID = "400"
CHOLINERGIC_SPLATTER_SUBCLUSTER_IDS = {"1634", "1635", "1636", "1637", "1638", "1640", "1641", "1642"}
CHOLINERGIC_JIA_GROUP = "Subpallial Cholinergic neurons"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", default=os.environ.get("PROJECT_ROOT", PROJECT_ROOT_DEFAULT))
    parser.add_argument("--scope", choices=sorted(SCOPES), required=True)
    parser.add_argument("--run-label", default="siletti_div90_seurat_bridge_v1")
    parser.add_argument("--outdir", default=None)
    parser.add_argument("--max-ref-cells-per-subcluster", type=int, default=100)
    parser.add_argument("--max-ref-cells-total", type=int, default=60000)
    parser.add_argument("--seed", type=int, default=0)
    return parser.parse_args()


def safe_token(value: str) -> str:
    return "".join(ch if ch.isalnum() else "_" for ch in str(value)).strip("_").lower()


def h5ad_filename(supercluster: str) -> str:
    return H5AD_BY_SUPERCLUSTER.get(supercluster, f"siletti_whb_{safe_token(supercluster)}.h5ad")


def unique_gene_index(var: pd.DataFrame, preferred_col: str) -> tuple[pd.Index, pd.DataFrame]:
    genes = var[preferred_col].astype(str)
    valid = genes.notna() & genes.ne("") & genes.ne("nan")
    counts = genes[valid].value_counts()
    unique = valid & genes.map(counts).fillna(0).eq(1)
    dropped = pd.DataFrame({"gene": genes, "kept_unique": unique.to_numpy()})
    return pd.Index(genes[unique]), dropped


def read_jia_mapping(project_root: Path) -> pd.DataFrame:
    path = (
        project_root
        / "results/siletti_2023_whb_reference_metadata/siletti_jia9_fetal_marker_pair_validation_v1/tables/row_level_fetal_marker_pair_scores_all_rows.csv"
    )
    mapping = pd.read_csv(path)
    mapping["cluster_id"] = pd.to_numeric(mapping["Cluster"], errors="coerce").astype("Int64").astype(str)
    mapping["subcluster_id"] = pd.to_numeric(mapping["Subcluster"], errors="coerce").astype("Int64").astype(str)
    keep_cols = [
        "cluster_id",
        "subcluster_id",
        "Supercluster",
        "Transferred MTG Label (Transferred from cluster level)",
        "candidate_jia_group",
        "jia_side",
        "anatomy_bin",
        "best_fetal_pair",
        "best_fetal_pair_score",
    ]
    return mapping[keep_cols].drop_duplicates()


def load_reference_metadata_and_indices(
    project_root: Path,
    source_dir: Path,
    scope: str,
    max_per_subcluster: int,
    max_total: int,
    seed: int,
) -> tuple[list[ad.AnnData], pd.DataFrame, pd.DataFrame]:
    rng = np.random.default_rng(seed)
    jia = read_jia_mapping(project_root)
    adatas = []
    metadata_rows = []
    global_offset = 0
    for supercluster in SCOPES[scope]:
        path = source_dir / h5ad_filename(supercluster)
        if not path.exists():
            raise FileNotFoundError(path)
        a = ad.read_h5ad(path, backed="r")
        adatas.append(a)
        obs = a.obs.copy()
        obs["source_h5ad"] = str(path)
        obs["source_supercluster"] = supercluster
        obs["source_row_index"] = np.arange(a.n_obs)
        obs["global_row_index"] = global_offset + np.arange(a.n_obs)
        global_offset += a.n_obs
        metadata_rows.append(obs)
    meta = pd.concat(metadata_rows, axis=0, ignore_index=True)
    meta["cluster_id"] = meta["cluster_id"].astype(str)
    meta["subcluster_id"] = meta["subcluster_id"].astype(str)
    meta = meta.merge(jia, on=["cluster_id", "subcluster_id"], how="left")
    meta["siletti_supercluster_label"] = meta["supercluster_term"].astype(str)
    meta["siletti_cluster_label"] = meta["source_supercluster"].astype(str) + "_cluster_" + meta["cluster_id"].astype(str)
    meta["siletti_subcluster_label"] = meta["source_supercluster"].astype(str) + "_subcluster_" + meta["subcluster_id"].astype(str)
    meta["transferred_mtg_label"] = meta["Transferred MTG Label (Transferred from cluster level)"].fillna("unlabeled_or_na").astype(str)
    meta["candidate_jia_group"] = meta["candidate_jia_group"].fillna("unassigned_jia_group").astype(str)
    meta["best_fetal_pair"] = meta["best_fetal_pair"].fillna("No fetal marker-pair hit").astype(str)

    scope_audit = []
    for supercluster, sub in meta.groupby("source_supercluster", sort=False):
        scope_audit.append(
            {
                "stage": "loaded_before_scope_filter",
                "source_supercluster": supercluster,
                "candidate_jia_group": "ALL",
                "n_cells": int(sub.shape[0]),
            }
        )

    if scope.endswith("_cholinergic"):
        is_splatter = meta["source_supercluster"].eq("Splatter")
        is_cholinergic_splatter = (
            is_splatter
            & meta["cluster_id"].eq(CHOLINERGIC_SPLATTER_CLUSTER_ID)
            & meta["subcluster_id"].isin(CHOLINERGIC_SPLATTER_SUBCLUSTER_IDS)
            & meta["candidate_jia_group"].eq(CHOLINERGIC_JIA_GROUP)
        )
        keep = (~is_splatter) | is_cholinergic_splatter
        scope_audit.append(
            {
                "stage": "splatter_cholinergic_filter",
                "source_supercluster": "Splatter",
                "candidate_jia_group": CHOLINERGIC_JIA_GROUP,
                "n_cells": int(is_cholinergic_splatter.sum()),
            }
        )
        scope_audit.append(
            {
                "stage": "splatter_excluded_by_scope_filter",
                "source_supercluster": "Splatter",
                "candidate_jia_group": "non_cholinergic_or_not_cluster_400",
                "n_cells": int((is_splatter & ~is_cholinergic_splatter).sum()),
            }
        )
        meta = meta.loc[keep].copy()

    selected_indices = []
    for _, sub in meta.groupby("siletti_subcluster_label", sort=False):
        idx = sub.index.to_numpy()
        if max_per_subcluster > 0 and idx.size > max_per_subcluster:
            idx = rng.choice(idx, size=max_per_subcluster, replace=False)
        selected_indices.extend(idx.tolist())
    selected_indices = np.array(selected_indices, dtype=int)
    if max_total > 0 and selected_indices.size > max_total:
        selected_indices = rng.choice(selected_indices, size=max_total, replace=False)
    selected_indices.sort()
    selected_meta = meta.loc[selected_indices].copy()
    selection_stage = "exported_all_after_scope_filter" if max_per_subcluster <= 0 and max_total <= 0 else "exported_after_subsampling"
    for (supercluster, group), sub in selected_meta.groupby(["source_supercluster", "candidate_jia_group"], sort=False):
        scope_audit.append(
            {
                "stage": selection_stage,
                "source_supercluster": supercluster,
                "candidate_jia_group": group,
                "n_cells": int(sub.shape[0]),
            }
        )
    return adatas, selected_meta, pd.DataFrame(scope_audit)


def collect_reference_matrix(adatas: list[ad.AnnData], meta: pd.DataFrame, genes: list[str], ref_gene_to_idx: dict[str, int]) -> sparse.csr_matrix:
    blocks = []
    for a in adatas:
        source_path = str(Path(a.filename).resolve())
        rows = meta.loc[meta["source_h5ad"].map(lambda x: str(Path(x).resolve())).eq(source_path), "source_row_index"].to_numpy(dtype=int)
        if rows.size == 0:
            continue
        col_idx = np.array([ref_gene_to_idx[g] for g in genes], dtype=int)
        mat = a.X[rows, :][:, col_idx]
        if not sparse.issparse(mat):
            mat = sparse.csr_matrix(mat)
        blocks.append(mat.tocsr())
    if not blocks:
        raise RuntimeError("No reference matrix blocks were collected; check source_h5ad path matching.")
    return sparse.vstack(blocks, format="csr")


def write_bridge(outdir: Path, prefix: str, mat_cells_by_genes: sparse.spmatrix, metadata: pd.DataFrame, genes: list[str]) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    mat = mat_cells_by_genes.T.tocoo()
    io.mmwrite(outdir / f"{prefix}_counts.mtx", mat)
    pd.Series(genes, name="gene").to_csv(outdir / f"{prefix}_genes.tsv", sep="\t", index=False)
    metadata.to_csv(outdir / f"{prefix}_metadata.tsv.gz", sep="\t", index=False)
    pd.Series(metadata["seurat_cell_id"].astype(str), name="cell_id").to_csv(outdir / f"{prefix}_barcodes.tsv", sep="\t", index=False)


def main() -> None:
    args = parse_args()
    project_root = Path(args.project_root)
    outdir = (
        Path(args.outdir)
        if args.outdir
        else project_root / "results/siletti_2023_whb_reference_label_transfer" / args.run_label / args.scope
    )
    table_dir = outdir / "tables"
    bridge_dir = outdir / "seurat_bridge"
    report_dir = outdir / "reports"
    for directory in [table_dir, bridge_dir, report_dir]:
        directory.mkdir(parents=True, exist_ok=True)

    source_dir = project_root / "results/siletti_2023_whb_reference_label_transfer/source_cellxgene_superclusters/h5ad"
    query_h5ad = project_root / "results/python_anndata/varela_div90.h5ad"
    query_cells_path = (
        project_root
        / "results/siletti_2023_whb_reference_label_transfer/siletti_div90_neuron_prep_v1/tables/div90_siletti_query_neuron_cells.tsv"
    )

    ref_adatas, ref_meta, scope_audit = load_reference_metadata_and_indices(
        project_root,
        source_dir,
        args.scope,
        args.max_ref_cells_per_subcluster,
        args.max_ref_cells_total,
        args.seed,
    )
    ref_var = ref_adatas[0].var.copy()
    ref_genes, ref_gene_audit = unique_gene_index(ref_var, "Gene")

    query = ad.read_h5ad(query_h5ad)
    query_gene_series = pd.Series(query.var_names.astype(str), index=query.var_names.astype(str))
    query_gene_counts = query_gene_series.value_counts()
    query_unique_genes = pd.Index(query_gene_series[query_gene_series.map(query_gene_counts).eq(1)])
    shared_genes = sorted(set(ref_genes).intersection(set(query_unique_genes)))
    if len(shared_genes) < 500:
        raise RuntimeError(f"Only {len(shared_genes)} shared unique genes; refusing export.")

    ref_gene_to_idx = {str(gene): int(idx) for idx, gene in enumerate(ref_var["Gene"].astype(str)) if str(gene) in shared_genes}
    query_gene_to_idx = {str(gene): int(idx) for idx, gene in enumerate(query.var_names.astype(str)) if str(gene) in shared_genes}

    query_cells = pd.read_csv(query_cells_path, sep="\t")
    query_obs_names = query_cells["obs_name"].astype(str)
    query_row_index = pd.Index(query.obs_names.astype(str)).get_indexer(query_obs_names)
    if (query_row_index < 0).any():
        raise RuntimeError("Some query obs_names were not found in DIV90 H5AD.")

    ref_meta = ref_meta.reset_index(drop=True)
    ref_meta["seurat_cell_id"] = [f"SILETTI_{args.scope}_{i:06d}" for i in range(ref_meta.shape[0])]
    ref_mat = collect_reference_matrix(ref_adatas, ref_meta, shared_genes, ref_gene_to_idx)
    for a in ref_adatas:
        a.file.close()

    q_col_idx = np.array([query_gene_to_idx[g] for g in shared_genes], dtype=int)
    source = query.layers["counts"] if "counts" in query.layers else query.X
    query_mat = source[query_row_index, :][:, q_col_idx]
    if not sparse.issparse(query_mat):
        query_mat = sparse.csr_matrix(query_mat)
    query_meta = query.obs.iloc[query_row_index].copy().reset_index(names="obs_name")
    query_meta = query_meta.merge(query_cells, on="obs_name", how="left", suffixes=("", "_manifest"))
    query_meta["seurat_cell_id"] = ["DIV90_QUERY_%06d" % i for i in range(query_meta.shape[0])]

    write_bridge(bridge_dir, "reference", ref_mat, ref_meta, shared_genes)
    write_bridge(bridge_dir, "query", query_mat, query_meta, shared_genes)

    ref_meta.to_csv(table_dir / "siletti_reference_cell_metadata.tsv.gz", sep="\t", index=False)
    query_meta.to_csv(table_dir / "div90_query_cell_metadata.tsv.gz", sep="\t", index=False)
    pd.DataFrame({"gene": shared_genes}).to_csv(table_dir / "siletti_div90_shared_genes.tsv", sep="\t", index=False)
    ref_gene_audit.to_csv(table_dir / "siletti_reference_gene_uniqueness_audit.tsv.gz", sep="\t", index=False)
    scope_audit.to_csv(table_dir / "siletti_reference_scope_and_subsampling_audit.tsv", sep="\t", index=False)

    label_counts = []
    for col in ["siletti_supercluster_label", "cell_type", "siletti_cluster_label", "siletti_subcluster_label", "transferred_mtg_label", "candidate_jia_group", "best_fetal_pair"]:
        vc = ref_meta[col].astype(str).value_counts()
        label_counts.extend({"label_column": col, "label": label, "n_reference_cells": int(n)} for label, n in vc.items())
    pd.DataFrame(label_counts).to_csv(table_dir / "siletti_reference_label_counts.tsv", sep="\t", index=False)

    config = {
        "scope": args.scope,
        "superclusters": SCOPES[args.scope],
        "cholinergic_splatter_cluster_id": CHOLINERGIC_SPLATTER_CLUSTER_ID if args.scope.endswith("_cholinergic") else None,
        "cholinergic_splatter_subcluster_ids": sorted(CHOLINERGIC_SPLATTER_SUBCLUSTER_IDS) if args.scope.endswith("_cholinergic") else None,
        "cholinergic_jia_group": CHOLINERGIC_JIA_GROUP if args.scope.endswith("_cholinergic") else None,
        "max_ref_cells_per_subcluster": args.max_ref_cells_per_subcluster,
        "max_ref_cells_total": args.max_ref_cells_total,
        "seed": args.seed,
        "n_reference_cells_exported": int(ref_meta.shape[0]),
        "n_query_cells_exported": int(query_meta.shape[0]),
        "n_shared_genes": int(len(shared_genes)),
        "query_h5ad": str(query_h5ad),
        "query_cells_path": str(query_cells_path),
    }
    (table_dir / "siletti_div90_seurat_bridge_config.json").write_text(json.dumps(config, indent=2, sort_keys=True) + "\n")
    report = [
        "# Siletti DIV90 Seurat Bridge Export",
        "",
        f"Scope: `{args.scope}`",
        f"Reference cells exported: {ref_meta.shape[0]:,}",
        f"Query cells exported: {query_meta.shape[0]:,}",
        f"Shared unique genes: {len(shared_genes):,}",
        "",
        "Bridge directory:",
        "",
        f"`{bridge_dir}`",
    ]
    (report_dir / "siletti_div90_seurat_bridge_report.md").write_text("\n".join(report) + "\n")
    print(json.dumps(config, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
