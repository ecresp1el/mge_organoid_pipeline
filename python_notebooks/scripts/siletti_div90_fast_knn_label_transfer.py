#!/usr/bin/env python
"""Fast Siletti adult reference -> DIV90 label transfer from bridge matrices.

This script is the transparent, sparse-matrix label-transfer implementation used
for the final Siletti adult-reference river figures. The final v8 figure used
this script with:

* ``--label-column unified_leaf_subtype``
* ``--exclude-label "Other selected reference"``
* ``--max-reference-cells 0`` and ``--max-query-cells 0`` (no subsampling)
* ``--nfeatures 3000``, ``--n-components 50``, and ``--k 50``

What this script is, and what it is not
---------------------------------------
Despite the historical directory name "Seurat bridge", this script does not run
Seurat and does not call ``FindTransferAnchors``. It also does not perform CCA,
RPCA, Harmony, or any anchor-based integration. The "bridge" is simply a pair
of matched sparse count matrices and metadata tables that were already exported
onto the same ordered set of shared unique genes.

The algorithm here is deliberately explicit:

1. Read reference and query count matrices from Matrix Market files.
2. Optionally remove unwanted reference labels before training. In the final v8
   run, ``Other selected reference`` was removed at this step.
3. Optionally downsample reference/query cells for smoke tests. A value of 0
   means no downsampling; the final v8 run used all exported cells.
4. Independently library-size normalize reference and query cells to 10,000
   counts per cell and apply ``log1p``.
5. Select high-variance genes jointly across the stacked normalized reference
   and query matrices. This makes the feature set responsive to both datasets,
   while the dimensional basis below is still fit on the adult reference.
6. Fit ``TruncatedSVD`` on the adult reference matrix and project DIV90 query
   cells with ``svd.transform`` into the same reference-derived latent basis.
7. L2-normalize those SVD coordinates and perform cosine nearest-neighbor
   search from each query cell to adult reference cells.
8. Convert cosine distances to positive similarities, normalize the 50 neighbor
   weights within each query cell, sum weights by reference label, and report
   the highest-scoring label.

Input bridge files
------------------
The input is an already-exported bridge directory containing:

* ``reference_counts.mtx`` and ``reference_metadata.tsv.gz``: adult Siletti
  reference cells from a chosen scope. For final v8 this was the restricted
  MGE/LAMP5-LHX6/CHAT scope.
* ``query_counts.mtx`` and ``query_metadata.tsv.gz``: DIV90 neuron-lineage query
  cells from the Varela DIV90 object.
* matching ``*_genes.tsv`` and ``*_barcodes.tsv`` files. The reference and query
  gene files must be identical and in the same order.

Label columns
-------------
For early Jia-style MGE tests, the label column was often ``candidate_jia_group``
and the default excluded label was ``Excluded / not assigned to Jia-style 9
groups``. For the final v8 Siletti adult-reference figure, the transferred label
was ``unified_leaf_subtype`` and the excluded label was ``Other selected
reference``. Broader plotting labels such as ``unified_major_subtype_roi`` and
``unified_pallial_subpallial_bin`` were derived downstream from the winning
``unified_leaf_subtype`` by the river-plotting script; they were not separate
independent transfers.

What the caps mean
------------------
``--max-reference-cells`` and ``--max-query-cells`` are compute/debug caps, not
biological filters. If a cap is >0, cells are uniformly sampled with the fixed
seed after reference-label filtering. If a cap is 0, all available cells are
used. The final v8 run used 0 for both caps.

Why the SVD is reference-fit
----------------------------
Feature selection is joint, but the SVD model is fit on the adult reference and
then applied to the DIV90 query. In code this is:

``ref_pcs = svd.fit_transform(x_ref)``
``query_pcs = svd.transform(x_query)``

This means adult identities are represented by a reference-derived
transcriptional basis, and query cells are classified by proximity to adult
reference cells in that basis.
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
    """Append an auditable progress row and echo it to stdout.

    The final v8 run kept this file as a compact execution trace. It records the
    matrix sizes after loading, the reference count after excluded-label removal,
    the number of selected features, SVD dimensionality, and kNN label count.
    """
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
    """Read one side of the bridge and return a cells x genes sparse matrix.

    Bridge counts are stored as Matrix Market features x cells, matching the
    export convention used elsewhere in the project. The transfer algorithm is
    easier to express as cells x genes, so this function transposes the matrix
    after validating that matrix shape agrees with the gene and barcode files.

    The barcode order defines matrix-column order. Metadata are reindexed to
    that same order with ``seurat_cell_id``. That column name is historical: the
    fast-kNN workflow uses the sparse matrices directly and does not create a
    Seurat object here.
    """
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
    """Library-size normalize each cell to ``scale_factor`` and apply log1p.

    Reference and query matrices are normalized separately. The transform is:

    ``log1p(raw_count / total_counts_in_cell * 10000)``

    This keeps the workflow close to common scRNA-seq log-normalized expression
    without borrowing information across reference/query cells during
    normalization.
    """
    x = x.astype(np.float32).tocsr(copy=True)
    totals = np.asarray(x.sum(axis=1)).ravel()
    totals[totals == 0] = 1.0
    factors = scale_factor / totals
    x = sparse.diags(factors.astype(np.float32)).dot(x).tocsr()
    x.data = np.log1p(x.data)
    return x


def sparse_variance(x: sparse.csr_matrix) -> np.ndarray:
    """Compute per-gene variance for a sparse cells x genes matrix."""
    mean = np.asarray(x.mean(axis=0)).ravel()
    mean_sq = np.asarray(x.power(2).mean(axis=0)).ravel()
    return mean_sq - np.square(mean)


def subset_rows(x: sparse.csr_matrix, meta: pd.DataFrame, max_cells: int, seed: int) -> tuple[sparse.csr_matrix, pd.DataFrame]:
    """Uniformly downsample cells for smoke/pilot runs.

    A ``max_cells`` value of 0 disables downsampling. For final figure runs this
    should remain 0 unless the figure is explicitly meant to show a sampled
    diagnostic.
    """
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
    """Convert neighbor distances into a normalized label vote.

    ``NearestNeighbors(metric="cosine")`` returns smaller values for closer
    cells. For each query cell, this function:

    1. Converts cosine distance to a nonnegative similarity with ``1 - d``.
       Negative similarities are clipped to 0.
    2. Normalizes the k neighbor similarities to sum to 1 within that query
       cell. If all similarities are 0, it falls back to an unweighted vote.
    3. Sums normalized weights by adult reference label.
    4. Reports the label with the largest summed weight and keeps the full score
       matrix for downstream audits.

    Therefore the prediction is not the single closest reference cell and not a
    plain majority vote. It is a cosine-ranked, similarity-weighted kNN vote.
    """
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
    parser.add_argument(
        "--exclude-labels",
        default=None,
        help="Optional '||'-separated labels to remove before transfer. Overrides --exclude-label. Use NONE to keep all labels.",
    )
    parser.add_argument("--max-reference-cells", type=int, default=0, help="Uniform random cap after label filtering; 0 means use all reference cells.")
    parser.add_argument("--max-query-cells", type=int, default=0, help="Uniform random cap for DIV90 query cells; 0 means use all query cells.")
    parser.add_argument("--nfeatures", type=int, default=3000, help="Number of high-variance shared genes used for SVD/KNN.")
    parser.add_argument("--n-components", type=int, default=50, help="Number of SVD dimensions used for cosine kNN.")
    parser.add_argument("--k", type=int, default=50, help="Number of adult reference neighbors voting on each DIV90 query cell.")
    parser.add_argument("--seed", type=int, default=0)
    return parser.parse_args()


def excluded_labels(args: argparse.Namespace) -> list[str]:
    """Parse labels to remove from the reference before fitting/neighbor search."""
    raw = args.exclude_labels if args.exclude_labels is not None else args.exclude_label
    if not raw or str(raw).upper() == "NONE":
        return []
    return [label.strip() for label in str(raw).split("||") if label.strip()]


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

    # Remove non-biological or accounting-only labels from the adult reference
    # before any dimensionality reduction or neighbor search. In the final v8
    # Siletti figure this removed "Other selected reference" from transfer
    # training, while keeping it available in upstream reference-accounting
    # tables.
    if args.label_column not in meta_ref.columns:
        raise ValueError(f"Reference metadata missing label column: {args.label_column}")
    meta_ref[args.label_column] = meta_ref[args.label_column].fillna("unlabeled_or_na").astype(str)
    labels_to_exclude = excluded_labels(args)
    if labels_to_exclude:
        keep = ~meta_ref[args.label_column].isin(labels_to_exclude)
        before = x_ref.shape[0]
        x_ref = x_ref[keep.to_numpy(), :]
        meta_ref = meta_ref.loc[keep].reset_index(drop=True)
        write_progress(
            progress_path,
            "filter_reference_label",
            "end",
            f"before={before}; after={x_ref.shape[0]}; excluded={' || '.join(labels_to_exclude)}",
        )

    # Optional compute/debug caps. These happen after excluded-label removal so
    # any pilot run samples from the actual training reference, not from labels
    # that would later be dropped.
    x_ref, meta_ref = subset_rows(x_ref, meta_ref, args.max_reference_cells, args.seed)
    x_query, meta_query = subset_rows(x_query, meta_query, args.max_query_cells, args.seed + 1)
    write_progress(progress_path, "downsample", "end", f"reference={x_ref.shape[0]}; query={x_query.shape[0]}")

    # Normalize reference and query independently, then use the transformed
    # matrices for both feature selection and SVD/kNN. No precomputed embedding,
    # integrated assay, anchor object, or reference UMAP is used here.
    x_ref = log_normalize(x_ref)
    x_query = log_normalize(x_query)
    write_progress(progress_path, "log_normalize", "end", f"reference_nnz={x_ref.nnz}; query_nnz={x_query.nnz}")

    # Joint feature selection: stack the normalized reference and query matrices
    # and select the highest-variance genes across both datasets. This is the
    # only step before SVD where reference and query are considered together.
    # The final v8 run selected 3,000 genes from 17,849 shared unique genes.
    var = sparse_variance(sparse.vstack([x_ref, x_query], format="csr"))
    nfeatures = min(args.nfeatures, x_ref.shape[1])
    feature_idx = np.argsort(var)[::-1][:nfeatures]
    feature_idx.sort()
    features = pd.DataFrame({"gene": np.asarray(genes_ref, dtype=object)[feature_idx], "variance": var[feature_idx]})
    features.to_csv(seurat_dir / "selected_transfer_features.tsv", sep="\t", index=False)
    x_ref = x_ref[:, feature_idx]
    x_query = x_query[:, feature_idx]
    write_progress(progress_path, "select_features", "end", f"nfeatures={nfeatures}")

    # Reference-fit latent basis: the SVD model is fit on adult reference cells,
    # then the DIV90 query is projected into that reference-derived basis. This
    # supports wording like "adult reference cells were embedded with
    # TruncatedSVD and DIV90 cells were projected into the same latent space."
    n_components = min(args.n_components, nfeatures - 1, x_ref.shape[0] - 1, x_query.shape[0] - 1)
    svd = TruncatedSVD(n_components=n_components, random_state=args.seed)
    ref_pcs = svd.fit_transform(x_ref)
    query_pcs = svd.transform(x_query)

    # L2-normalize coordinates before cosine nearest-neighbor search. Cosine
    # distance on normalized vectors is a scale-insensitive comparison of the
    # SVD expression profiles.
    ref_pcs = normalize(ref_pcs)
    query_pcs = normalize(query_pcs)
    write_progress(progress_path, "svd", "end", f"n_components={n_components}; explained_variance={svd.explained_variance_ratio_.sum():.4f}")

    # Find the k nearest adult reference cells for each DIV90 query cell and
    # assign the adult label with the largest normalized similarity-weighted
    # vote. The full per-label vote table is saved, which lets downstream audits
    # inspect score margins and second-best labels.
    labels = sorted(meta_ref[args.label_column].dropna().astype(str).unique().tolist())
    k = min(args.k, x_ref.shape[0])
    nn = NearestNeighbors(n_neighbors=k, metric="cosine", algorithm="brute", n_jobs=-1)
    nn.fit(ref_pcs)
    distances, indices = nn.kneighbors(query_pcs)
    predicted, max_score, score_df = vote_labels(distances, indices, meta_ref[args.label_column].astype(str).to_numpy(), labels)
    write_progress(progress_path, "knn_transfer", "end", f"k={k}; labels={len(labels)}")

    # Prediction outputs:
    # - fast_knn/*_predictions.tsv.gz: compact predictions + per-label scores.
    # - fast_knn/*_prediction_scores.tsv.gz: score matrix only.
    # - tables/*_query_obs_with_predictions.tsv.gz: original query metadata plus
    #   predicted label and maximum vote score.
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

    # Cluster-level table used by figure/audit code to summarize which adult
    # labels were assigned to each original DIV90 cluster.
    cluster_col = "cluster_id" if "cluster_id" in obs.columns else "cluster_id_manifest" if "cluster_id_manifest" in obs.columns else None
    if cluster_col:
        cluster_counts = obs.groupby([cluster_col, "predicted.id"], dropna=False).size().reset_index(name="n_cells")
        cluster_counts.to_csv(tables_dir / f"{prefix}_cluster_label_counts.tsv", sep="\t", index=False)

    # Store the exact transfer settings in a machine-readable sidecar. This is
    # the quickest way to confirm, after the fact, that the run used the intended
    # label column, reference/query counts, feature count, SVD dimensions, and k.
    diag = {
        "method": "fast_knn_svd_cosine",
        "bridge_dir": str(args.bridge_dir),
        "label_column": args.label_column,
        "exclude_label": args.exclude_label,
        "exclude_labels": labels_to_exclude,
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
