"""Canonical Shi Seurat label-transfer column schema helpers."""

from __future__ import annotations

import re
from typing import Iterable

import numpy as np
import pandas as pd


PREDICTED_LABEL_COL = "shi_seurat_full_predicted_shi_label"
PREDICTION_SCORE_COL = "shi_seurat_full_prediction_score"
UNCERTAINTY_SCORE_COL = "shi_seurat_full_uncertainty_score"
LABEL_SCORE_PREFIX = "shi_seurat_full_prediction_score_"

PREDICTED_WEEK_LABEL_COL = "shi_seurat_full_predicted_shi_week_label"
WEEK_PREDICTION_SCORE_COL = "shi_seurat_full_week_prediction_score"
WEEK_UNCERTAINTY_SCORE_COL = "shi_seurat_full_week_uncertainty_score"
WEEK_SCORE_PREFIX = "shi_seurat_full_week_prediction_score_"


def sanitize_shi_label_token(value: object) -> str:
    """Return the stable suffix used by R and Python Shi score exports."""
    token = re.sub(r"[^A-Za-z0-9]+", "_", str(value).strip())
    token = re.sub(r"_+", "_", token).strip("_")
    return token or "value"


def canonical_score_column(label: object, prefix: str = LABEL_SCORE_PREFIX) -> str:
    return prefix + sanitize_shi_label_token(label)


def canonical_score_columns(columns: Iterable[str], prefix: str = LABEL_SCORE_PREFIX) -> list[str]:
    return [col for col in columns if col.startswith(prefix) and col != prefix.rstrip("_")]


def validate_canonical_prediction_scores(
    data: pd.DataFrame,
    score_cols: list[str],
    max_score_col: str = PREDICTION_SCORE_COL,
    uncertainty_col: str = UNCERTAINTY_SCORE_COL,
    context: str = "Shi Seurat prediction table",
    atol: float = 1e-8,
) -> None:
    """Assert max score and uncertainty agree with canonical per-label scores."""
    if not score_cols:
        raise ValueError(f"{context}: no canonical per-label score columns found.")
    missing = [col for col in [max_score_col, uncertainty_col, *score_cols] if col not in data.columns]
    if missing:
        raise ValueError(f"{context}: missing required score columns: {', '.join(missing)}")
    score_matrix = data[score_cols].apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float)
    max_score = pd.to_numeric(data[max_score_col], errors="coerce").to_numpy(dtype=float)
    uncertainty = pd.to_numeric(data[uncertainty_col], errors="coerce").to_numpy(dtype=float)
    if not np.isfinite(score_matrix).all() or not np.isfinite(max_score).all() or not np.isfinite(uncertainty).all():
        raise ValueError(f"{context}: score columns contain non-finite values.")
    bad_max = (max_score < -atol) | (max_score > 1.0 + atol)
    bad_scores = (score_matrix < -atol) | (score_matrix > 1.0 + atol)
    if bad_max.any() or bad_scores.any():
        raise ValueError(f"{context}: score columns contain values outside [0,1].")
    if not np.allclose(max_score, score_matrix.max(axis=1), atol=atol):
        raise ValueError(f"{context}: max score is not the row-wise maximum of per-label scores.")
    if not np.allclose(uncertainty, 1.0 - max_score, atol=atol):
        raise ValueError(f"{context}: uncertainty is not 1 - max score.")
