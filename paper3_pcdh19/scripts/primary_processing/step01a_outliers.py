"""Per-sample robust MAD boundaries and candidate-flag summaries for Step 01a."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

from .step01a_models import Step01aSettings


@dataclass(frozen=True)
class MADBoundary:
    """Describe one sample, metric, direction, and MAD-stringency boundary."""

    sample: str
    metric: str
    transformed_metric: str
    direction: str
    stringency_mad: float
    median_transformed: float
    raw_mad_transformed: float
    scaled_mad_transformed: float
    boundary_transformed: float
    boundary_original_scale: float


@dataclass(frozen=True)
class Step01aResult:
    """Collect all reproducible Step 01a tables produced from cell metadata."""

    boundaries: pd.DataFrame
    per_cell: pd.DataFrame
    flag_summary: pd.DataFrame
    overlap_summary: pd.DataFrame


class RobustMADBoundaryCalculator:
    """Calculate one-sided robust boundaries independently within each sample."""

    METRICS = (
        ("total_counts", "log1p_total_counts", "low"),
        ("n_genes_by_counts", "log1p_n_genes_by_counts", "low"),
        ("pct_counts_mt", "pct_counts_mt", "high"),
    )

    def __init__(self, settings: Step01aSettings):
        """Store the fixed stringencies and scaled-MAD convention."""

        self.settings = settings

    def calculate(self, obs: pd.DataFrame) -> Step01aResult:
        """Return boundaries and flags without modifying the source AnnData."""

        working = obs[[
            self.settings.sample_field,
            "total_counts",
            "n_genes_by_counts",
            "pct_counts_mt",
        ]].copy()
        working.insert(0, "cell_id", obs.index.astype(str))
        boundaries: list[dict[str, object]] = []
        for sample, sample_frame in working.groupby(self.settings.sample_field, sort=False, observed=True):
            index = sample_frame.index
            for metric, transformed_metric, direction in self.METRICS:
                original = sample_frame[metric].to_numpy(dtype=float)
                transformed = np.log1p(original) if transformed_metric.startswith("log1p_") else original
                median = float(np.median(transformed))
                raw_mad = float(np.median(np.abs(transformed - median)))
                scaled_mad = self.settings.mad_scale * raw_mad
                for stringency in self.settings.stringencies:
                    signed = -1.0 if direction == "low" else 1.0
                    boundary_transformed = median + signed * stringency * scaled_mad
                    boundary_original = (
                        max(0.0, float(np.expm1(boundary_transformed)))
                        if transformed_metric.startswith("log1p_")
                        else max(0.0, float(boundary_transformed))
                    )
                    suffix = self._suffix(stringency)
                    flag_name = self._flag_name(metric, suffix)
                    if direction == "low":
                        working.loc[index, flag_name] = original < boundary_original
                    else:
                        working.loc[index, flag_name] = original > boundary_original
                    boundary = MADBoundary(
                        sample=str(sample),
                        metric=metric,
                        transformed_metric=transformed_metric,
                        direction=direction,
                        stringency_mad=stringency,
                        median_transformed=median,
                        raw_mad_transformed=raw_mad,
                        scaled_mad_transformed=scaled_mad,
                        boundary_transformed=boundary_transformed,
                        boundary_original_scale=boundary_original,
                    )
                    boundaries.append(boundary.__dict__)
        for stringency in self.settings.stringencies:
            suffix = self._suffix(stringency)
            individual = self._individual_flag_columns(suffix)
            working[f"any_candidate_{suffix}"] = working[individual].any(axis=1)
            working[f"all_three_candidate_{suffix}"] = working[individual].all(axis=1)
            working[f"candidate_pattern_{suffix}"] = self._patterns(working, suffix)
        boolean_columns = [column for column in working if column.startswith(("low_", "high_", "any_candidate_", "all_three_candidate_"))]
        working[boolean_columns] = working[boolean_columns].astype(bool)
        boundary_frame = pd.DataFrame(boundaries)
        return Step01aResult(
            boundaries=boundary_frame,
            per_cell=working,
            flag_summary=self._flag_summary(working, boundary_frame),
            overlap_summary=self._overlap_summary(working),
        )

    def _flag_summary(self, per_cell: pd.DataFrame, boundaries: pd.DataFrame) -> pd.DataFrame:
        """Count individual and joint candidate flags per sample and stringency."""

        rows: list[dict[str, object]] = []
        for sample, frame in per_cell.groupby(self.settings.sample_field, sort=False, observed=True):
            n_cells = len(frame)
            for stringency in self.settings.stringencies:
                suffix = self._suffix(stringency)
                boundary_rows = boundaries.loc[
                    (boundaries["sample"] == str(sample))
                    & (boundaries["stringency_mad"] == stringency)
                ].set_index("metric")
                row: dict[str, object] = {
                    "technical_sample_id": sample,
                    "stringency_mad": stringency,
                    "n_cells": n_cells,
                    "low_total_counts_boundary": boundary_rows.loc["total_counts", "boundary_original_scale"],
                    "low_n_genes_boundary": boundary_rows.loc["n_genes_by_counts", "boundary_original_scale"],
                    "high_pct_counts_mt_boundary": boundary_rows.loc["pct_counts_mt", "boundary_original_scale"],
                }
                for label, column in (
                    ("low_total_counts", f"low_total_counts_{suffix}"),
                    ("low_n_genes", f"low_n_genes_by_counts_{suffix}"),
                    ("high_pct_counts_mt", f"high_pct_counts_mt_{suffix}"),
                    ("any_candidate", f"any_candidate_{suffix}"),
                    ("all_three_candidate", f"all_three_candidate_{suffix}"),
                ):
                    count = int(frame[column].sum())
                    row[f"{label}_n"] = count
                    row[f"{label}_pct"] = 100.0 * count / n_cells
                rows.append(row)
        return pd.DataFrame(rows)

    def _overlap_summary(self, per_cell: pd.DataFrame) -> pd.DataFrame:
        """Count exact three-criterion flag patterns per sample and stringency."""

        rows: list[dict[str, object]] = []
        for sample, frame in per_cell.groupby(self.settings.sample_field, sort=False, observed=True):
            n_cells = len(frame)
            for stringency in self.settings.stringencies:
                suffix = self._suffix(stringency)
                counts = frame[f"candidate_pattern_{suffix}"].value_counts()
                for pattern in ("none", "counts_only", "genes_only", "mt_only", "counts+genes", "counts+mt", "genes+mt", "all_three"):
                    count = int(counts.get(pattern, 0))
                    rows.append({
                        "technical_sample_id": sample,
                        "stringency_mad": stringency,
                        "candidate_pattern": pattern,
                        "n_cells": count,
                        "pct_cells": 100.0 * count / n_cells,
                    })
        return pd.DataFrame(rows)

    @staticmethod
    def _suffix(stringency: float) -> str:
        """Return a stable field-name suffix for one MAD multiplier."""

        return f"{stringency:g}mad"

    @staticmethod
    def _flag_name(metric: str, suffix: str) -> str:
        """Map a source metric to its one-sided candidate-flag field."""

        direction = "high" if metric == "pct_counts_mt" else "low"
        return f"{direction}_{metric}_{suffix}"

    @staticmethod
    def _individual_flag_columns(suffix: str) -> list[str]:
        """Return the three prespecified one-sided candidate flags."""

        return [
            f"low_total_counts_{suffix}",
            f"low_n_genes_by_counts_{suffix}",
            f"high_pct_counts_mt_{suffix}",
        ]

    @classmethod
    def _patterns(cls, frame: pd.DataFrame, suffix: str) -> pd.Series:
        """Encode exact overlap among the three candidate criteria."""

        count_flag, gene_flag, mt_flag = (frame[column].to_numpy(dtype=bool) for column in cls._individual_flag_columns(suffix))
        code = count_flag.astype(np.uint8) + 2 * gene_flag.astype(np.uint8) + 4 * mt_flag.astype(np.uint8)
        labels = np.array(["none", "counts_only", "genes_only", "counts+genes", "mt_only", "counts+mt", "genes+mt", "all_three"], dtype=object)
        return pd.Series(labels[code], index=frame.index, dtype="string")

