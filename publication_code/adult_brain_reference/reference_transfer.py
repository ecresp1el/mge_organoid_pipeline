"""Similarity-weighted label transfer from an adult single-cell reference."""
from __future__ import annotations
import numpy as np
import pandas as pd
from scipy import sparse
from sklearn.decomposition import TruncatedSVD
from sklearn.neighbors import NearestNeighbors
from sklearn.preprocessing import normalize

def log_normalize(x: sparse.csr_matrix, scale_factor: float = 1e4) -> sparse.csr_matrix:
    """Return float32 CSR counts normalized per cell to scale_factor, then log1p.
    
    Input: cells x genes raw counts. Totals use all supplied genes; zero-total
    cells remain zero. Reference and query cells are normalized separately."""
    x = x.astype(np.float32).tocsr(copy=True)
    totals = np.asarray(x.sum(axis=1)).ravel()
    totals[totals == 0] = 1.0
    factors = scale_factor / totals
    x = sparse.diags(factors.astype(np.float32)).dot(x).tocsr()
    x.data = np.log1p(x.data)
    return x


def sparse_variance(x: sparse.csr_matrix) -> np.ndarray:
    """Return per-gene population variance from sparse first and second moments.
    
    Input: normalized cells x genes matrix. Implicit zeros contribute to both moments."""
    mean = np.asarray(x.mean(axis=0)).ravel()
    mean_sq = np.asarray(x.power(2).mean(axis=0)).ravel()
    return mean_sq - np.square(mean)


def vote_labels(
    distances: np.ndarray,
    indices: np.ndarray,
    ref_labels: np.ndarray,
    labels: list[str],
) -> tuple[np.ndarray, np.ndarray, pd.DataFrame]:
    """Sum normalized cosine-similarity neighbor weights by reference label.
    
    Inputs: neighbor distances, reference indices, reference labels, ordered labels.
    Weights are max(1 - distance, 0), normalized within each query cell. When all
    weights are zero the vote is uniform across its neighbors. Returns winning
    labels, maximum vote scores, and all label scores. Scores are relative
    neighbor support, not calibrated probabilities."""
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
    """Convert a label to an alphanumeric/underscore score-column suffix."""
    token = "".join(ch if ch.isalnum() else "_" for ch in str(value).strip())
    while "__" in token:
        token = token.replace("__", "_")
    token = token.strip("_")
    return token or "value"

def transfer_labels(reference_counts, query_counts, reference_metadata, query_metadata,
                    genes, label_column="unified_leaf_subtype", nfeatures=3000,
                    n_components=50, k=50, seed=0, n_jobs=1,
                    exclude_labels=("Other selected reference",)):
    """Transfer reference labels using shared-gene expression and weighted kNN.

    Inputs: reference/query sparse raw counts (cells x shared genes), row-aligned
    metadata, and a single ordered gene list shared by both matrices. The label
    column defines reference classes; excluded labels are removed before fitting.
    No cell subsampling is performed. n_jobs controls the neighbor-search workers.

    Returns: dictionary with predictions (query metadata plus winning label and
    score), full label_scores, selected_features, and resolved parameters.
    Gene selection uses reference and query variance jointly; SVD is fitted only
    to reference cells, and query cells are projected into that reference basis.
    """
    x_ref = sparse.csr_matrix(reference_counts)
    x_query = sparse.csr_matrix(query_counts)
    meta_ref = reference_metadata.copy().reset_index(drop=True)
    meta_query = query_metadata.copy().reset_index(drop=True)
    if x_ref.shape[1] != x_query.shape[1] or x_ref.shape[1] != len(genes):
        raise ValueError("Count columns must match the shared gene list")
    if x_ref.shape[0] != len(meta_ref) or x_query.shape[0] != len(meta_query):
        raise ValueError("Metadata rows must match count rows")
    if len(set(genes)) != len(genes) or label_column not in meta_ref:
        raise ValueError("Supply unique genes and a reference label column")
    if meta_ref[label_column].isna().any():
        raise ValueError("Reference labels must not be missing")
    if min(nfeatures, n_components, k) < 1:
        raise ValueError("Feature, component, and neighbor counts must be positive")
    for matrix in (x_ref, x_query):
        if not np.isfinite(matrix.data).all() or (matrix.data < 0).any():
            raise ValueError("Counts must be finite and nonnegative")
    # Exclude accounting labels before normalization, feature selection, and fitting.
    keep = ~meta_ref[label_column].astype(str).isin(exclude_labels)
    x_ref = x_ref[keep.to_numpy()]
    meta_ref = meta_ref.loc[keep].reset_index(drop=True)
    if min(x_ref.shape[0], x_query.shape[0], x_ref.shape[1]) < 2:
        raise ValueError("Need at least two cells per dataset and two shared genes")
    x_ref = log_normalize(x_ref)
    x_query = log_normalize(x_query)
    # Rank shared genes jointly, retaining the original column order after selection.
    variance = sparse_variance(sparse.vstack([x_ref, x_query], format="csr"))
    feature_idx = np.argsort(variance)[::-1][:min(nfeatures, len(genes))]
    feature_idx.sort()
    features = pd.DataFrame({"gene": np.asarray(genes)[feature_idx], "variance": variance[feature_idx]})
    dimensions = min(n_components, len(feature_idx) - 1, x_ref.shape[0] - 1, x_query.shape[0] - 1)
    if dimensions < 1:
        raise ValueError("Not enough selected features for SVD")
    # Reference-fit SVD supplies a common basis; query cells are transformed only.
    svd = TruncatedSVD(n_components=dimensions, random_state=seed)
    ref_coords = normalize(svd.fit_transform(x_ref[:, feature_idx]))
    query_coords = normalize(svd.transform(x_query[:, feature_idx]))
    neighbors = min(k, x_ref.shape[0])
    model = NearestNeighbors(n_neighbors=neighbors, metric="cosine", algorithm="brute", n_jobs=n_jobs)
    model.fit(ref_coords)
    distances, indices = model.kneighbors(query_coords)
    labels = sorted(meta_ref[label_column].astype(str).unique())
    predicted, max_score, scores = vote_labels(distances, indices, meta_ref[label_column].astype(str).to_numpy(), labels)
    meta_query["predicted.id"] = predicted
    meta_query["prediction.score.max"] = max_score
    return {"predictions": meta_query, "label_scores": scores, "selected_features": features,
            "parameters": {"nfeatures": len(feature_idx), "n_components": dimensions,
                           "k": neighbors, "seed": seed, "label_column": label_column}}
