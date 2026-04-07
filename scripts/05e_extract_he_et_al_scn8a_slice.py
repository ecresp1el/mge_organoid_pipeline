#!/usr/bin/env python3
"""
Extract a lightweight gene slice + UMAP from He et al hnoca_allmeta.h5ad.

Outputs under --outdir:
  - scn8a_counts.mtx            (genes x cells MatrixMarket; filename kept for compatibility)
  - genes.tsv                   (one row per extracted gene)
  - barcodes.tsv                (cell IDs)
  - umap.tsv.gz                 (cell_id, UMAP_1, UMAP_2)
  - obs_metadata.tsv.gz         (subset of metadata columns + cell_id)
  - slice_summary.tsv           (key-value summary)
"""

import argparse
import gzip
import os
import sys
from typing import List


def log(msg: str) -> None:
    sys.stderr.write(msg + "\n")
    sys.stderr.flush()


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Extract gene slice + UMAP from He et al h5ad")
    p.add_argument("--h5ad", required=True, help="Path to hnoca_allmeta.h5ad")
    p.add_argument("--outdir", required=True, help="Output directory for extracted slice")
    p.add_argument(
        "--gene",
        default="SCN8A",
        help="Single target gene symbol (default: SCN8A). Ignored when --genes is provided.",
    )
    p.add_argument(
        "--genes",
        default="",
        help="Comma-separated target genes (for example: SCN8A,LHX6,NKX2.1)",
    )
    p.add_argument(
        "--umap-key",
        default="X_umap",
        help="Preferred obsm key for UMAP coordinates (default: X_umap)",
    )
    return p.parse_args()


def normalize_series(values) -> List[str]:
    out: List[str] = []
    for v in values:
        if v is None:
            out.append("")
        else:
            out.append(str(v))
    return out


def canonical_gene_name(gene: str) -> str:
    return gene.strip().upper().replace(".", "-")


def gene_aliases(gene: str) -> List[str]:
    g = gene.strip().upper()
    if not g:
        return []
    aliases = [g, g.replace(".", "-"), g.replace("-", ".")]
    # preserve order while unique
    seen = set()
    out: List[str] = []
    for a in aliases:
        if a not in seen:
            seen.add(a)
            out.append(a)
    return out


def find_gene_column(var_df, query_gene: str) -> dict:
    candidate_columns = [
        None,
        "feature_name",
        "gene_symbol",
        "gene_symbols",
        "symbol",
        "gene_name",
    ]

    aliases = gene_aliases(query_gene)
    if not aliases:
        raise ValueError("Empty gene query")

    for col in candidate_columns:
        if col is None:
            vals_raw = normalize_series(var_df.index.tolist())
            source = "var_index"
        else:
            if col not in var_df.columns:
                continue
            vals_raw = normalize_series(var_df[col].tolist())
            source = col

        vals = [v.strip().upper() for v in vals_raw]
        for alias in aliases:
            for i, v in enumerate(vals):
                if v == alias:
                    return {
                        "index": i,
                        "matched_token": vals_raw[i],
                        "matched_alias": alias,
                        "source": source,
                    }

    raise ValueError(
        f"Gene '{query_gene}' not found in h5ad var index/columns. "
        f"Checked index + {candidate_columns[1:]}"
    )


def select_umap_key(adata, preferred_key: str) -> str:
    keys = list(adata.obsm.keys())
    if preferred_key in keys:
        return preferred_key
    fallback = ["X_umap", "umap", "X_umap_harmony"]
    for k in fallback:
        if k in keys:
            return k
    raise ValueError(f"No UMAP-like embedding found in obsm keys: {keys}")


def parse_gene_targets(args: argparse.Namespace) -> List[str]:
    if args.genes.strip():
        raw = [x.strip() for x in args.genes.split(",")]
        genes = [x for x in raw if x]
    else:
        genes = [args.gene.strip()]

    if not genes:
        raise ValueError("No target genes provided")

    # canonicalize, preserve order, unique
    out: List[str] = []
    seen = set()
    for g in genes:
        cg = canonical_gene_name(g)
        if cg and cg not in seen:
            seen.add(cg)
            out.append(cg)
    if not out:
        raise ValueError("No valid target genes after normalization")
    return out


def main() -> int:
    args = parse_args()

    try:
        import anndata as ad
        import numpy as np
        import pandas as pd
        import scipy.io
        import scipy.sparse
    except Exception as exc:  # pragma: no cover
        log("ERROR: missing required Python packages (anndata, numpy, pandas, scipy).")
        log(f"Import error: {exc}")
        log("Install in your analysis env, for example:")
        log("  pip install anndata pandas scipy numpy")
        return 2

    try:
        target_genes = parse_gene_targets(args)
    except Exception as exc:
        log(f"ERROR: {exc}")
        return 2

    h5ad_path = os.path.abspath(args.h5ad)
    outdir = os.path.abspath(args.outdir)

    if not os.path.exists(h5ad_path):
        log(f"ERROR: h5ad not found: {h5ad_path}")
        return 2

    os.makedirs(outdir, exist_ok=True)

    log(f"Opening h5ad (backed): {h5ad_path}")
    adata = ad.read_h5ad(h5ad_path, backed="r")

    n_cells = int(adata.n_obs)
    n_features = int(adata.n_vars)
    log(f"Cells: {n_cells:,}; features: {n_features:,}")

    gene_hits = []
    for g in target_genes:
        try:
            hit = find_gene_column(adata.var, g)
        except Exception as exc:
            log(f"ERROR: {exc}")
            return 2
        hit["requested_gene"] = g
        gene_hits.append(hit)
        log(
            f"Resolved {g} -> token '{hit['matched_token']}' "
            f"(alias '{hit['matched_alias']}', source={hit['source']}, idx={hit['index']})"
        )

    gene_indices = [h["index"] for h in gene_hits]
    out_genes = [h["requested_gene"] for h in gene_hits]

    log("Extracting target expression columns from X")
    x = adata[:, gene_indices].X  # n_cells x n_genes

    if scipy.sparse.issparse(x):
        x_csc = x.tocsc()
    else:
        x_arr = np.asarray(x)
        x_csc = scipy.sparse.csc_matrix(x_arr)

    # Convert to genes x cells for Seurat-style readMM downstream.
    counts_gene_by_cell = x_csc.T.tocoo()

    counts_path = os.path.join(outdir, "scn8a_counts.mtx")
    genes_path = os.path.join(outdir, "genes.tsv")
    barcodes_path = os.path.join(outdir, "barcodes.tsv")
    umap_path = os.path.join(outdir, "umap.tsv.gz")
    meta_path = os.path.join(outdir, "obs_metadata.tsv.gz")
    summary_path = os.path.join(outdir, "slice_summary.tsv")

    log(f"Writing MatrixMarket: {counts_path}")
    scipy.io.mmwrite(counts_path, counts_gene_by_cell)

    cells = [str(x) for x in adata.obs_names.tolist()]

    with open(genes_path, "w", encoding="utf-8") as fh:
        for g in out_genes:
            fh.write(g + "\n")

    with open(barcodes_path, "w", encoding="utf-8") as fh:
        for cell in cells:
            fh.write(cell + "\n")

    umap_key = select_umap_key(adata, args.umap_key)
    log(f"Using UMAP embedding key: {umap_key}")
    umap = adata.obsm[umap_key]
    umap_arr = np.asarray(umap)
    if umap_arr.ndim != 2 or umap_arr.shape[0] != len(cells) or umap_arr.shape[1] < 2:
        log(
            "ERROR: invalid UMAP matrix shape. "
            f"Expected (n_cells, >=2), got {umap_arr.shape}"
        )
        return 2

    with gzip.open(umap_path, "wt", encoding="utf-8") as fh:
        fh.write("cell_id\tUMAP_1\tUMAP_2\n")
        for i, cell in enumerate(cells):
            fh.write(f"{cell}\t{umap_arr[i, 0]}\t{umap_arr[i, 1]}\n")

    # Keep metadata small and reproducible; avoid copying all columns.
    obs = adata.obs
    keep_obs_cols = [
        "dataset_id",
        "batch",
        "sample",
        "cell_type",
        "tissue",
        "assay",
        "development_stage",
        "donor_id",
    ]
    keep_obs_cols = [c for c in keep_obs_cols if c in obs.columns]
    if keep_obs_cols:
        obs_small = obs[keep_obs_cols].copy()
    else:
        obs_small = pd.DataFrame(index=obs.index)
    obs_small.insert(0, "cell_id", cells)
    obs_small.to_csv(meta_path, sep="\t", index=False, compression="gzip")

    with open(summary_path, "w", encoding="utf-8") as fh:
        fh.write("key\tvalue\n")
        fh.write(f"h5ad_path\t{h5ad_path}\n")
        fh.write(f"outdir\t{outdir}\n")
        fh.write(f"genes\t{','.join(out_genes)}\n")
        fh.write(f"n_cells\t{n_cells}\n")
        fh.write(f"n_features\t{n_features}\n")
        fh.write(f"umap_key\t{umap_key}\n")
        fh.write(f"counts_mtx\t{counts_path}\n")
        fh.write(f"genes_tsv\t{genes_path}\n")
        fh.write(f"barcodes_tsv\t{barcodes_path}\n")
        fh.write(f"umap_tsv_gz\t{umap_path}\n")
        fh.write(f"obs_metadata_tsv_gz\t{meta_path}\n")

    log("Done.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
