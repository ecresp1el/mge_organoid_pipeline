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
from typing import Dict, List


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
        default="X_umap_scpoli",
        help="Preferred obsm key for UMAP coordinates (default: X_umap_scpoli)",
    )
    p.add_argument(
        "--row-chunk",
        type=int,
        default=4000,
        help="Number of cell rows to process per X chunk (default: 4000)",
    )
    p.add_argument(
        "--io-chunk",
        type=int,
        default=200000,
        help="Rows per chunk when writing text outputs (default: 200000)",
    )
    return p.parse_args()


def decode_value(v) -> str:
    if isinstance(v, bytes):
        return v.decode("utf-8", errors="replace")
    return str(v)


def normalize_series(values) -> List[str]:
    out: List[str] = []
    for v in values:
        if v is None:
            out.append("")
        else:
            out.append(decode_value(v))
    return out


def canonical_gene_name(gene: str) -> str:
    return gene.strip().upper().replace(".", "-")


def gene_aliases(gene: str) -> List[str]:
    g = gene.strip().upper()
    if not g:
        return []
    aliases = [g, g.replace(".", "-"), g.replace("-", ".")]
    seen = set()
    out: List[str] = []
    for a in aliases:
        if a not in seen:
            seen.add(a)
            out.append(a)
    return out


def parse_gene_targets(args: argparse.Namespace) -> List[str]:
    if args.genes.strip():
        raw = [x.strip() for x in args.genes.split(",")]
        genes = [x for x in raw if x]
    else:
        genes = [args.gene.strip()]

    if not genes:
        raise ValueError("No target genes provided")

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


def load_var_strings(var_group, key: str) -> List[str]:
    obj = var_group[key]
    if hasattr(obj, "keys") and "categories" in obj and "codes" in obj:
        cats = normalize_series(obj["categories"][()])
        codes = obj["codes"][()]
        out = []
        n_cats = len(cats)
        for code in codes:
            c = int(code)
            if c < 0 or c >= n_cats:
                out.append("")
            else:
                out.append(cats[c])
        return out
    return normalize_series(obj[()])


def find_gene_index(var_group, query_gene: str) -> Dict[str, object]:
    candidate_keys = ["_index", "feature_name", "gene_symbol", "gene_symbols", "symbol", "gene_name"]
    aliases = gene_aliases(query_gene)
    if not aliases:
        raise ValueError("Empty gene query")

    for key in candidate_keys:
        if key not in var_group:
            continue
        vals_raw = load_var_strings(var_group, key)
        vals = [v.strip().upper() for v in vals_raw]
        for alias in aliases:
            for i, v in enumerate(vals):
                if v == alias:
                    return {
                        "index": i,
                        "matched_token": vals_raw[i],
                        "matched_alias": alias,
                        "source": key,
                    }

    raise ValueError(
        f"Gene '{query_gene}' not found in h5ad var index/columns. "
        f"Checked keys: {candidate_keys}"
    )


def choose_umap_key(obsm_group, preferred_key: str) -> str:
    keys = list(obsm_group.keys())
    if preferred_key in keys:
        return preferred_key

    priority = [
        "X_umap",
        "umap",
        "X_umap_scpoli",
        "X_umap_pca_unintegrated",
        "X_umap_pca_rss",
    ]
    for k in priority:
        if k in keys:
            return k

    umap_like = [k for k in keys if k.lower().startswith("x_umap") or "umap" in k.lower()]
    if umap_like:
        return sorted(umap_like)[0]

    raise ValueError(f"No UMAP-like embedding found in obsm keys: {keys}")


def sanitize_tsv(v: str) -> str:
    return v.replace("\t", " ").replace("\n", " ").replace("\r", " ")


def read_obs_column_chunk(obs_group, col: str, start: int, end: int, cache: Dict[str, List[str]]) -> List[str]:
    obj = obs_group[col]
    if hasattr(obj, "keys") and "categories" in obj and "codes" in obj:
        if col not in cache:
            cache[col] = normalize_series(obj["categories"][()])
        cats = cache[col]
        n_cats = len(cats)
        codes = obj["codes"][start:end]
        out = []
        for code in codes:
            c = int(code)
            if c < 0 or c >= n_cats:
                out.append("")
            else:
                out.append(cats[c])
        return out
    return normalize_series(obj[start:end])


def main() -> int:
    args = parse_args()

    try:
        import h5py
        import numpy as np
        import scipy.io
        import scipy.sparse
    except Exception as exc:  # pragma: no cover
        log("ERROR: missing required Python packages (h5py, numpy, scipy).")
        log(f"Import error: {exc}")
        log("Install in your analysis env, for example:")
        log("  pip install h5py numpy scipy")
        return 2

    try:
        target_genes = parse_gene_targets(args)
    except Exception as exc:
        log(f"ERROR: {exc}")
        return 2

    h5ad_path = os.path.abspath(args.h5ad)
    outdir = os.path.abspath(args.outdir)
    row_chunk = max(1, int(args.row_chunk))
    io_chunk = max(1, int(args.io_chunk))

    if not os.path.exists(h5ad_path):
        log(f"ERROR: h5ad not found: {h5ad_path}")
        return 2

    os.makedirs(outdir, exist_ok=True)

    counts_path = os.path.join(outdir, "scn8a_counts.mtx")
    genes_path = os.path.join(outdir, "genes.tsv")
    barcodes_path = os.path.join(outdir, "barcodes.tsv")
    umap_path = os.path.join(outdir, "umap.tsv.gz")
    meta_path = os.path.join(outdir, "obs_metadata.tsv.gz")
    summary_path = os.path.join(outdir, "slice_summary.tsv")

    log(f"Opening h5ad (h5py): {h5ad_path}")
    with h5py.File(h5ad_path, "r") as f:
        if "X" not in f or "indptr" not in f["X"] or "indices" not in f["X"] or "data" not in f["X"]:
            log("ERROR: unsupported X layout; expected CSR-like group with data/indices/indptr")
            return 2

        x_group = f["X"]
        indptr = x_group["indptr"][()]
        indices_ds = x_group["indices"]
        data_ds = x_group["data"]

        n_cells = int(indptr.shape[0] - 1)
        n_features = int(f["var"]["_index"].shape[0]) if "var" in f and "_index" in f["var"] else -1
        log(f"Cells: {n_cells:,}; features: {n_features:,}")

        if "var" not in f:
            log("ERROR: missing var table in h5ad")
            return 2
        var_group = f["var"]

        gene_hits = []
        for g in target_genes:
            try:
                hit = find_gene_index(var_group, g)
            except Exception as exc:
                log(f"ERROR: {exc}")
                return 2
            hit["requested_gene"] = g
            gene_hits.append(hit)
            log(
                f"Resolved {g} -> token '{hit['matched_token']}' "
                f"(alias '{hit['matched_alias']}', source={hit['source']}, idx={hit['index']})"
            )

        target_gene_indices = np.array([int(h["index"]) for h in gene_hits], dtype=np.int64)
        out_genes = [h["requested_gene"] for h in gene_hits]

        log("Extracting target expression columns from CSR X (chunked)")
        row_chunks = []
        col_chunks = []
        val_chunks = []
        total_hits = 0

        for r0 in range(0, n_cells, row_chunk):
            r1 = min(r0 + row_chunk, n_cells)
            start = int(indptr[r0])
            end = int(indptr[r1])
            if end <= start:
                continue

            idx_block = indices_ds[start:end]
            mask = np.isin(idx_block, target_gene_indices)
            if not np.any(mask):
                if r0 == 0 or (r0 // row_chunk) % 50 == 0:
                    log(f"Processed rows {r0:,}-{r1:,}; cumulative hits {total_hits:,}")
                continue

            pos = np.nonzero(mask)[0].astype(np.int64, copy=False)
            sel_idx = idx_block[pos]
            sel_vals = data_ds[start:end][pos]

            row_ptr = indptr[r0 : r1 + 1] - start
            cell_offsets = np.searchsorted(row_ptr[1:], pos, side="right")
            cell_ids = cell_offsets.astype(np.int64, copy=False) + r0

            gene_rows = np.full(pos.shape[0], -1, dtype=np.int16)
            for g_row, g_idx in enumerate(target_gene_indices):
                gene_rows[sel_idx == g_idx] = g_row

            valid = gene_rows >= 0
            if np.any(valid):
                rows = gene_rows[valid].astype(np.int32, copy=False)
                cols = cell_ids[valid].astype(np.int64, copy=False)
                vals = np.asarray(sel_vals[valid], dtype=np.float32)
                row_chunks.append(rows)
                col_chunks.append(cols)
                val_chunks.append(vals)
                total_hits += int(vals.shape[0])

            if r0 == 0 or (r0 // row_chunk) % 50 == 0:
                log(f"Processed rows {r0:,}-{r1:,}; cumulative hits {total_hits:,}")

        if total_hits > 0:
            rows_all = np.concatenate(row_chunks)
            cols_all = np.concatenate(col_chunks)
            vals_all = np.concatenate(val_chunks)
            counts_gene_by_cell = scipy.sparse.coo_matrix(
                (vals_all, (rows_all, cols_all)),
                shape=(len(out_genes), n_cells),
                dtype=np.float32,
            )
            counts_gene_by_cell.sum_duplicates()
        else:
            counts_gene_by_cell = scipy.sparse.coo_matrix((len(out_genes), n_cells), dtype=np.float32)

        log(f"Writing MatrixMarket: {counts_path}")
        scipy.io.mmwrite(counts_path, counts_gene_by_cell)

        with open(genes_path, "w", encoding="utf-8") as fh:
            for g in out_genes:
                fh.write(g + "\n")

        if "obs" not in f or "_index" not in f["obs"]:
            log("ERROR: missing obs/_index in h5ad")
            return 2
        obs_group = f["obs"]
        obs_index_ds = obs_group["_index"]

        log(f"Writing barcodes: {barcodes_path}")
        with open(barcodes_path, "w", encoding="utf-8") as fh:
            for i in range(0, n_cells, io_chunk):
                j = min(i + io_chunk, n_cells)
                cells = normalize_series(obs_index_ds[i:j])
                for cell in cells:
                    fh.write(cell + "\n")

        if "obsm" not in f:
            log("ERROR: missing obsm in h5ad")
            return 2
        umap_key = choose_umap_key(f["obsm"], args.umap_key)
        umap_ds = f["obsm"][umap_key]
        if len(umap_ds.shape) != 2 or umap_ds.shape[0] != n_cells or umap_ds.shape[1] < 2:
            log(
                "ERROR: invalid UMAP matrix shape. "
                f"Expected (n_cells, >=2), got {umap_ds.shape}"
            )
            return 2

        log(f"Using UMAP embedding key: {umap_key}")
        log(f"Writing UMAP table: {umap_path}")
        with gzip.open(umap_path, "wt", encoding="utf-8") as fh:
            fh.write("cell_id\tUMAP_1\tUMAP_2\n")
            for i in range(0, n_cells, io_chunk):
                j = min(i + io_chunk, n_cells)
                cells = normalize_series(obs_index_ds[i:j])
                coords = umap_ds[i:j, 0:2]
                for k, cell in enumerate(cells):
                    fh.write(f"{cell}\t{coords[k, 0]}\t{coords[k, 1]}\n")

        keep_obs_cols = ["dataset_id", "batch", "sample", "cell_type", "tissue", "assay", "development_stage", "donor_id"]
        keep_obs_cols = [c for c in keep_obs_cols if c in obs_group]
        cat_cache: Dict[str, List[str]] = {}

        log(f"Writing obs metadata: {meta_path}")
        with gzip.open(meta_path, "wt", encoding="utf-8") as fh:
            header = ["cell_id"] + keep_obs_cols
            fh.write("\t".join(header) + "\n")
            for i in range(0, n_cells, io_chunk):
                j = min(i + io_chunk, n_cells)
                cells = normalize_series(obs_index_ds[i:j])
                col_data = {
                    col: read_obs_column_chunk(obs_group, col, i, j, cat_cache)
                    for col in keep_obs_cols
                }
                chunk_len = j - i
                for k in range(chunk_len):
                    row = [sanitize_tsv(cells[k])]
                    for col in keep_obs_cols:
                        row.append(sanitize_tsv(col_data[col][k]))
                    fh.write("\t".join(row) + "\n")

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
