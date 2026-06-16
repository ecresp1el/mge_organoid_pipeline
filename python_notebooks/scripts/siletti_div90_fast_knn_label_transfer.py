#!/usr/bin/env python
"""Fast Siletti adult reference -> DIV90 label transfer from bridge matrices.

What this script transfers
--------------------------
The input is an already-exported "Seurat bridge" directory containing:

* ``reference_counts.mtx`` and ``reference_metadata.tsv.gz``: adult Siletti
  reference cells from a chosen scope, usually ``mge_llc``.
* ``query_counts.mtx`` and ``query_metadata.tsv.gz``: DIV90 neuron-lineage query
  cells from the Varela DIV90 object.

For the Jia-style MGE question, the intended label column is
``candidate_jia_group``. That column has one deliberately non-biological bucket,
``Excluded / not assigned to Jia-style 9 groups``. By default this script drops
that label before transfer, because otherwise the classifier would learn to
assign query cells to "not a Jia-style group" rather than one of the adult
MGE/LLC subtype candidates.

What the caps mean
------------------
``--max-reference-cells`` and ``--max-query-cells`` are compute/debug caps, not
biological filters. If a cap is >0, cells are uniformly sampled with the fixed
seed after label filtering. If a cap is 0, all available cells are used. The
current scale-up path is:

1. smoke: tiny caps to prove file/label logic.
2. pilot: 5,000 adult reference x 3,000 DIV90 query cells.
3. full mge_llc: no caps after filtering the excluded label, i.e. 18,459 adult
   reference cells x 16,206 DIV90 query cells.

Why not mge_cge_llc for candidate_jia_group?
--------------------------------------------
For ``candidate_jia_group``, the added CGE cells are almost entirely excluded
because Jia's comparison is MGE-derived inhibitory neuron biology. If the goal
changes to broader interneuron MTG labels such as Vip/Lamp5/Sncg/Pax6, then
``mge_cge_llc`` and another label column such as ``transferred_mtg_label`` would
be appropriate.

Algorithm
---------
This is a deliberately transparent transfer path for debugging/scaling:
counts -> log-normalized sparse matrix -> variable genes -> TruncatedSVD ->
cosine kNN label voting. It uses the same bridge files as the Seurat run, but
avoids the opaque Seurat ``FindTransferAnchors`` bottleneck.
"""

from __future__ import annotations

import argparse
import gzip
import json
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmread
from sklearn.decomposition import TruncatedSVD
from sklearn.neighbors import NearestNeighbors
from sklearn.preprocessing import normalize


def timestamp() -> str:
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def write_progress(path: Path, step: str, status: str, detail: str = "") -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    exists = path.exists()
    with path.open("a") as handle:
        if not exists:
            handle.write("timestamp\tstep\tstatus\tdetail\n")
        handle.write(f"{timestamp()}\t{step}\t{status}\t{detail}\n")
    print(f"[{timestamp()}] {step} [{status}] {detail}", flush=True)


def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t")


def read_bridge_counts(prefix: str, bridge_dir: Path) -> tuple[sparse.csr_matrix, list[str], pd.DataFrame]:
    counts_path = bridge_dir / f"{prefix}_counts.mtx"
    genes_path = bridge_dir / f"{prefix}_genes.tsv"
    barcodes_path = bridge_dir / f"{prefix}_barcodes.tsv"
    meta_path = bridge_dir / f"{prefix}_metadata.tsv.gz"

    mat = mmread(counts_path).tocsr()
    genes = read_tsv(genes_path)["gene"].astype(str).tolist()
    barcodes = read_tsv(barcodes_path)["cell_id"].astype(str).tolist()
    meta = read_tsv(meta_path)

    if mat.shape != (len(genes), len(barcodes)):
        raise ValueError(f"{prefix}: matrix shape {mat.shape} does not match genes/barcodes")
    if "seurat_cell_id" not in meta.columns:
        raise ValueError(f"{prefix}: metadata missing seurat_cell_id")

    # Work in cells x genes.
    x = mat.T.tocsr()
    meta = meta.set_index("seurat_cell_id").loc[barcodes].reset_index()
    return x, genes, meta


def log_normalize(x: sparse.csr_matrix, scale_factor: float = 1e4) -> sparse.csr_matrix:
    x = x.astype(np.float32).tocsr(copy=True)
    totals = np.asarray(x.sum(axis=1)).ravel()
    totals[totals == 0] = 1.0
    factors = scale_factor / totals
    x = sparse.diags(factors.astype(np.float32)).dot(x).tocsr()
    x.data = np.log1p(x.data)
    return x


def sparse_variance(x: sparse.csr_matrix) -> np.ndarray:
    mean = np.asarray(x.mean(axis=0)).ravel()
    mean_sq = np.asarray(x.power(2).mean(axis=0)).ravel()
    return mean_sq - np.square(mean)


def subset_rows(x: sparse.csr_matrix, meta: pd.DataFrame, max_cells: int, seed: int) -> tuple[sparse.csr_matrix, pd.DataFrame]:
    if max_cells <= 0 or x.shape[0] <= max_cells:
        return x, meta
    rng = np.random.default_rng(seed)
    idx = np.sort(rng.choice(x.shape[0], size=max_cells, replace=False))
    return x[idx, :], meta.iloc[idx].reset_index(drop=True)


def vote_labels(
    distances: np.ndarray,
    indices: np.ndarray,
    ref_labels: np.ndarray,
    labels: list[str],
) -> tuple[np.ndarray, np.ndarray, pd.DataFrame]:
    # Cosine distance is [0, 2]. Convert to positive similarity and normalize.
    weights = np.clip(1.0 - distances, a_min=0.0, a_max=None)
    row_sums = weights.sum(axis=1)
    zero_rows = row_sums == 0
    if np.any(zero_rows):
        weights[zero_rows, :] = 1.0
        row_sums = weights.sum(axis=1)
    weights = weights / row_sums[:, None]

    label_to_col = {label: i for i, label in enumerate(labels)}
    scores = np.zeros((indices.shape[0], len(labels)), dtype=np.float32)
    neighbor_labels = ref_labels[indices]
    for i in range(indices.shape[0]):
        for j in range(indices.shape[1]):
            scores[i, label_to_col[neighbor_labels[i, j]]] += weights[i, j]

    best_idx = scores.argmax(axis=1)
    predicted = np.asarray(labels, dtype=object)[best_idx]
    max_score = scores[np.arange(scores.shape[0]), best_idx]
    score_df = pd.DataFrame(scores, columns=[f"score_{safe_token(label)}" for label in labels])
    return predicted, max_score, score_df


def safe_token(value: str) -> str:
    token = "".join(ch if ch.isalnum() else "_" for ch in str(value).strip())
    while "__" in token:
        token = token.replace("__", "_")
    token = token.strip("_")
    return token or "value"


def write_tsv_gz(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt") as handle:
        df.to_csv(handle, sep="\t", index=False)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bridge-dir", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--label-column", default="candidate_jia_group", help="Reference metadata column to transfer onto DIV90 cells.")
    parser.add_argument("--exclude-label", default="Excluded / not assigned to Jia-style 9 groups", help="Reference label to remove before transfer; use NONE to keep all labels.")
    parser.add_argument("--max-reference-cells", type=int, default=0, help="Uniform random cap after label filtering; 0 means use all reference cells.")
    parser.add_argument("--max-query-cells", type=int, default=0, help="Uniform random cap for DIV90 query cells; 0 means use all query cells.")
    parser.add_argument("--nfeatures", type=int, default=3000, help="Number of high-variance shared genes used for SVD/KNN.")
    parser.add_argument("--n-components", type=int, default=50, help="Number of SVD dimensions used for cosine kNN.")
    parser.add_argument("--k", type=int, default=50, help="Number of adult reference neighbors voting on each DIV90 query cell.")
    parser.add_argument("--seed", type=int, default=0)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    tables_dir = args.outdir / "tables"
    seurat_dir = args.outdir / "fast_knn"
    tables_dir.mkdir(parents=True, exist_ok=True)
    seurat_dir.mkdir(parents=True, exist_ok=True)
    progress_path = seurat_dir / "run_progress.tsv"

    write_progress(progress_path, "run", "start", f"outdir={args.outdir}")
    x_ref, genes_ref, meta_ref = read_bridge_counts("reference", args.bridge_dir)
    x_query, genes_query, meta_query = read_bridge_counts("query", args.bridge_dir)
    if genes_ref != genes_query:
        raise ValueError("Reference and query gene order differ; rebuild bridge or add explicit gene intersection.")
    write_progress(progress_path, "load", "end", f"reference={x_ref.shape}; query={x_query.shape}; genes={len(genes_ref)}")

    if args.label_column not in meta_ref.columns:
        raise ValueError(f"Reference metadata missing label column: {args.label_column}")
    meta_ref[args.label_column] = meta_ref[args.label_column].fillna("unlabeled_or_na").astype(str)
    if args.exclude_label and args.exclude_label.upper() != "NONE":
        keep = meta_ref[args.label_column] != args.exclude_label
        before = x_ref.shape[0]
        x_ref = x_ref[keep.to_numpy(), :]
        meta_ref = meta_ref.loc[keep].reset_index(drop=True)
        write_progress(progress_path, "filter_reference_label", "end", f"before={before}; after={x_ref.shape[0]}; excluded={args.exclude_label}")

    x_ref, meta_ref = subset_rows(x_ref, meta_ref, args.max_reference_cells, args.seed)
    x_query, meta_query = subset_rows(x_query, meta_query, args.max_query_cells, args.seed + 1)
    write_progress(progress_path, "downsample", "end", f"reference={x_ref.shape[0]}; query={x_query.shape[0]}")

    x_ref = log_normalize(x_ref)
    x_query = log_normalize(x_query)
    write_progress(progress_path, "log_normalize", "end", f"reference_nnz={x_ref.nnz}; query_nnz={x_query.nnz}")

    var = sparse_variance(sparse.vstack([x_ref, x_query], format="csr"))
    nfeatures = min(args.nfeatures, x_ref.shape[1])
    feature_idx = np.argsort(var)[::-1][:nfeatures]
    feature_idx.sort()
    features = pd.DataFrame({"gene": np.asarray(genes_ref, dtype=object)[feature_idx], "variance": var[feature_idx]})
    features.to_csv(seurat_dir / "selected_transfer_features.tsv", sep="\t", index=False)
    x_ref = x_ref[:, feature_idx]
    x_query = x_query[:, feature_idx]
    write_progress(progress_path, "select_features", "end", f"nfeatures={nfeatures}")

    n_components = min(args.n_components, nfeatures - 1, x_ref.shape[0] - 1, x_query.shape[0] - 1)
    svd = TruncatedSVD(n_components=n_components, random_state=args.seed)
    ref_pcs = svd.fit_transform(x_ref)
    query_pcs = svd.transform(x_query)
    ref_pcs = normalize(ref_pcs)
    query_pcs = normalize(query_pcs)
    write_progress(progress_path, "svd", "end", f"n_components={n_components}; explained_variance={svd.explained_variance_ratio_.sum():.4f}")

    labels = sorted(meta_ref[args.label_column].dropna().astype(str).unique().tolist())
    k = min(args.k, x_ref.shape[0])
    nn = NearestNeighbors(n_neighbors=k, metric="cosine", algorithm="brute", n_jobs=-1)
    nn.fit(ref_pcs)
    distances, indices = nn.kneighbors(query_pcs)
    predicted, max_score, score_df = vote_labels(distances, indices, meta_ref[args.label_column].astype(str).to_numpy(), labels)
    write_progress(progress_path, "knn_transfer", "end", f"k={k}; labels={len(labels)}")

    cell_col = "seurat_cell_id" if "seurat_cell_id" in meta_query.columns else meta_query.columns[0]
    predictions = pd.DataFrame({
        "cell_id": meta_query[cell_col].astype(str).to_numpy(),
        "predicted.id": predicted,
        "prediction.score.max": max_score,
    })
    predictions = pd.concat([predictions, score_df], axis=1)

    prefix = f"siletti_div90__{safe_token(args.label_column)}__fast_knn__svd{n_components}__k{k}"
    write_tsv_gz(predictions, seurat_dir / f"{prefix}_predictions.tsv.gz")
    write_tsv_gz(predictions[["cell_id"] + [c for c in predictions.columns if c.startswith("score_")]], seurat_dir / f"{prefix}_prediction_scores.tsv.gz")

    obs = pd.concat([meta_query.reset_index(drop=True), predictions[["predicted.id", "prediction.score.max"]]], axis=1)
    write_tsv_gz(obs, tables_dir / f"{prefix}_query_obs_with_predictions.tsv.gz")

    cluster_col = "cluster_id" if "cluster_id" in obs.columns else "cluster_id_manifest" if "cluster_id_manifest" in obs.columns else None
    if cluster_col:
        cluster_counts = obs.groupby([cluster_col, "predicted.id"], dropna=False).size().reset_index(name="n_cells")
        cluster_counts.to_csv(tables_dir / f"{prefix}_cluster_label_counts.tsv", sep="\t", index=False)

    diag = {
        "method": "fast_knn_svd_cosine",
        "bridge_dir": str(args.bridge_dir),
        "label_column": args.label_column,
        "exclude_label": args.exclude_label,
        "n_reference_cells": int(x_ref.shape[0]),
        "n_query_cells": int(x_query.shape[0]),
        "n_features": int(nfeatures),
        "n_components": int(n_components),
        "k": int(k),
        "seed": int(args.seed),
        "labels": labels,
    }
    (seurat_dir / f"{prefix}_transfer_diagnostics.json").write_text(json.dumps(diag, indent=2) + "\n")
    write_progress(progress_path, "run", "end", f"prefix={prefix}")


if __name__ == "__main__":
    main()
