#!/usr/bin/env python3
"""Summarize marker-like genes for selected Shi-transfer labels from exported MTX."""

from __future__ import annotations

import argparse
import gzip
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.io import mmread
from scipy import sparse


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--matrix", required=True, type=Path, help="Feature x cell MTX, preferably log-normalized data.")
    parser.add_argument("--features", required=True, type=Path, help="features.tsv with one feature/gene column.")
    parser.add_argument("--obs", required=True, type=Path, help="Shi transfer obs TSV/TSV.GZ in the same cell order as the MTX.")
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--label-col", default="shi_seurat_full_predicted_shi_label")
    parser.add_argument("--score-col", default="shi_seurat_full_prediction_score")
    parser.add_argument(
        "--target-label",
        action="append",
        dest="target_labels",
        default=[],
        help="Label to summarize. Repeatable.",
    )
    parser.add_argument("--top-n", type=int, default=40)
    parser.add_argument("--min-pct-group", type=float, default=0.05)
    return parser.parse_args()


def read_features(path: Path) -> pd.Series:
    features = pd.read_csv(path, sep="\t")
    if features.empty:
        raise ValueError(f"No features found: {path}")
    return features.iloc[:, 0].astype(str)


def read_obs(path: Path, label_col: str, score_col: str) -> pd.DataFrame:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as handle:
        header = handle.readline().rstrip("\n").split("\t")
    usecols = [label_col]
    for col in ["cell_id", "cell_id_for_join", "orig.ident", "seurat_clusters", score_col]:
        if col in header and col not in usecols:
            usecols.append(col)
    return pd.read_csv(path, sep="\t", usecols=usecols)


def summarize_label(
    matrix: sparse.csc_matrix,
    total_sum: np.ndarray,
    total_nnz: np.ndarray,
    features: pd.Series,
    labels: np.ndarray,
    scores: pd.Series | None,
    target: str,
    top_n: int,
    min_pct_group: float,
) -> tuple[pd.DataFrame, dict[str, object]]:
    group_mask = labels == target
    n_group = int(group_mask.sum())
    n_cells = int(labels.shape[0])
    n_rest = n_cells - n_group
    if n_group == 0:
        return pd.DataFrame(), {"shi_label": target, "n_cells": 0, "mean_prediction_score": np.nan}

    group = matrix[:, np.flatnonzero(group_mask)]
    group_sum = np.asarray(group.sum(axis=1)).ravel()
    group_nnz = group.getnnz(axis=1)
    rest_sum = total_sum - group_sum
    rest_nnz = total_nnz - group_nnz

    mean_group = group_sum / n_group
    mean_rest = rest_sum / max(n_rest, 1)
    pct_group = group_nnz / n_group
    pct_rest = rest_nnz / max(n_rest, 1)
    mean_diff = mean_group - mean_rest
    pct_diff = pct_group - pct_rest

    keep = (pct_group >= min_pct_group) & (mean_diff > 0)
    stats = pd.DataFrame(
        {
            "shi_label": target,
            "gene": features.to_numpy(),
            "n_group_cells": n_group,
            "n_rest_cells": n_rest,
            "mean_log_expr_group": mean_group,
            "mean_log_expr_rest": mean_rest,
            "mean_log_expr_diff": mean_diff,
            "pct_expr_group": pct_group,
            "pct_expr_rest": pct_rest,
            "pct_expr_diff": pct_diff,
        }
    )
    stats = stats.loc[keep].sort_values(
        ["mean_log_expr_diff", "pct_expr_diff", "mean_log_expr_group"],
        ascending=[False, False, False],
    )
    stats.insert(2, "rank", np.arange(1, stats.shape[0] + 1))

    summary = {
        "shi_label": target,
        "n_cells": n_group,
        "mean_prediction_score": float(pd.to_numeric(scores[group_mask], errors="coerce").mean())
        if scores is not None
        else np.nan,
    }
    return stats.head(top_n), summary


def main() -> None:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    targets = args.target_labels or ["Thalamic neurons", "Excitatory neuron", "Excitatory IPC"]

    features = read_features(args.features)
    obs = read_obs(args.obs, args.label_col, args.score_col)
    labels = obs[args.label_col].astype(str).to_numpy()
    scores = obs[args.score_col] if args.score_col in obs.columns else None

    matrix = mmread(args.matrix)
    matrix = matrix.tocsc() if sparse.issparse(matrix) else sparse.csc_matrix(matrix)
    if matrix.shape[0] != features.shape[0]:
        raise ValueError(f"Feature count mismatch: matrix has {matrix.shape[0]}, features has {features.shape[0]}")
    if matrix.shape[1] != obs.shape[0]:
        raise ValueError(f"Cell count mismatch: matrix has {matrix.shape[1]}, obs has {obs.shape[0]}")

    total_sum = np.asarray(matrix.sum(axis=1)).ravel()
    total_nnz = matrix.getnnz(axis=1)

    marker_tables = []
    summaries = []
    for target in targets:
        markers, summary = summarize_label(
            matrix=matrix,
            total_sum=total_sum,
            total_nnz=total_nnz,
            features=features,
            labels=labels,
            scores=scores,
            target=target,
            top_n=args.top_n,
            min_pct_group=args.min_pct_group,
        )
        if not markers.empty:
            marker_tables.append(markers)
        summaries.append(summary)

    if marker_tables:
        pd.concat(marker_tables, ignore_index=True).to_csv(
            args.outdir / "shi_seurat_unexpected_label_marker_summary.tsv",
            sep="\t",
            index=False,
        )
    pd.DataFrame(summaries).to_csv(
        args.outdir / "shi_seurat_unexpected_label_marker_counts.tsv",
        sep="\t",
        index=False,
    )

    group_cols = [col for col in ["orig.ident", "seurat_clusters", args.label_col] if col in obs.columns]
    if len(group_cols) == 3:
        counts = obs.groupby(group_cols, observed=True).size().reset_index(name="n_cells")
        counts.to_csv(args.outdir / "shi_seurat_label_counts_by_sample_cluster.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
