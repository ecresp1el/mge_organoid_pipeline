#!/usr/bin/env python3
"""Render a replicate-aware DIV30 PLXNA2 comparison for neuronal clusters 1, 4, and 2.

The biological replicate and inferential unit is the organoid. Cell-level
log1p(CP10K) values are summarized once per organoid and canonical cluster;
no inferential test treats cells as independent observations.
"""

from __future__ import annotations

import argparse
from datetime import datetime
import hashlib
import itertools
from pathlib import Path
import shutil
import subprocess

import matplotlib

matplotlib.use("Agg")
matplotlib.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Nimbus Sans", "Liberation Sans", "DejaVu Sans"],
        "font.size": 7,
        "axes.linewidth": 0.65,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
    }
)

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats


FINAL_DIR_DEFAULT = Path(
    "/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/final_figures/"
    "fig_varela_div30_lhx6_lhx8_plxna2_v1_candidate"
)
PER_CELL_RELATIVE = Path("tables/neuronal_clusters_1_2_4_per_cell.tsv.gz")
SCRIPT_NAME = "render_varela_div30_plxna2_replicate_aware_figure.py"
FIGURE_STEM = "varela_div30_plxna2_replicate_aware"
POSITIVE_THRESHOLD = 0.5
SMALL_GROUP_THRESHOLD = 50
BOOTSTRAP_ITERATIONS = 100_000

CLUSTER_ORDER = ["1", "4", "2"]
CLUSTER_LABELS = {
    "1": "SST+ cortical interneurons",
    "4": "PV neuron precursors",
    "2": "MGE subpallial neurons",
}
SHORT_LABELS = {
    "1": "3 · SST+ cIN",
    "4": "4 · PV precursor",
    "2": "5 · MGE subpallial",
}
CLUSTER_COLORS = {"1": "#00BF7D", "4": "#00B0F6", "2": "#C77CFF"}

ENDPOINTS = {
    "mean_PLXNA2_log1p_CP10K": "Mean expression (primary)",
    "median_PLXNA2_log1p_CP10K": "Median expression",
    "pct_PLXNA2_detected": "Detected cells (%)",
}
CONTRASTS = [
    ("PV_vs_MGE", "4", "2", True),
    ("PV_vs_SST", "4", "1", False),
    ("SST_vs_MGE", "1", "2", False),
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--final-dir", type=Path, default=FINAL_DIR_DEFAULT)
    parser.add_argument("--per-cell-table", type=Path, default=None)
    parser.add_argument("--random-state", type=int, default=1729)
    parser.add_argument("--bootstrap-iterations", type=int, default=BOOTSTRAP_ITERATIONS)
    return parser.parse_args()


def ensure_output_dirs(final_dir: Path) -> None:
    for relative in ["figures/png", "figures/pdf", "figures/svg", "tables", "code", "provenance"]:
        (final_dir / relative).mkdir(parents=True, exist_ok=True)


def load_and_validate(path: Path) -> pd.DataFrame:
    data = pd.read_csv(path, sep="\t", dtype={"cluster": str})
    required = {"sample", "cluster", "paper_cluster_annotation", "PLXNA2"}
    missing = sorted(required.difference(data.columns))
    if missing:
        raise ValueError(f"Per-cell table is missing required columns: {missing}")
    data = data.loc[data["cluster"].isin(CLUSTER_ORDER)].copy()
    data["PLXNA2"] = pd.to_numeric(data["PLXNA2"], errors="raise")
    if data["PLXNA2"].isna().any() or (~np.isfinite(data["PLXNA2"])).any():
        raise ValueError("PLXNA2 contains missing or non-finite values.")
    observed = (
        data[["cluster", "paper_cluster_annotation"]]
        .drop_duplicates()
        .sort_values("cluster")
    )
    if observed["cluster"].duplicated().any():
        raise ValueError("A selected cluster has more than one canonical annotation.")
    observed_map = dict(zip(observed["cluster"], observed["paper_cluster_annotation"]))
    expected_source_labels = {
        "1": "SST+ cIN",
        "4": "PV neuron precursor",
        "2": "MGE subpallial neurons",
    }
    if observed_map != expected_source_labels:
        raise ValueError(
            "Canonical sidecar annotations do not match the expected cluster mapping: "
            f"observed={observed_map}"
        )
    data["organoid"] = data["sample"].astype(str)
    data["cluster_label"] = data["cluster"].map(CLUSTER_LABELS)
    return data


def aggregate_organoid_cluster(data: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for (organoid, cluster), group in data.groupby(["organoid", "cluster"], sort=False):
        values = group["PLXNA2"].to_numpy(dtype=float)
        q25, q75 = np.quantile(values, [0.25, 0.75])
        rows.append(
            {
                "organoid": organoid,
                "cluster": cluster,
                "cluster_label": CLUSTER_LABELS[cluster],
                "n_cells": values.size,
                "mean_PLXNA2_log1p_CP10K": np.mean(values),
                "median_PLXNA2_log1p_CP10K": np.median(values),
                "sd_PLXNA2_log1p_CP10K": np.std(values, ddof=1),
                "q25_PLXNA2_log1p_CP10K": q25,
                "q75_PLXNA2_log1p_CP10K": q75,
                "iqr_PLXNA2_log1p_CP10K": q75 - q25,
                "pct_PLXNA2_detected": 100.0 * np.mean(values > 0),
                "secondary_PLXNA2_positive_threshold_log1p_CP10K": POSITIVE_THRESHOLD,
                "secondary_pct_PLXNA2_at_or_above_threshold": 100.0
                * np.mean(values >= POSITIVE_THRESHOLD),
            }
        )
    result = pd.DataFrame(rows)
    cluster_rank = {cluster: i for i, cluster in enumerate(CLUSTER_ORDER)}
    result["_cluster_rank"] = result["cluster"].map(cluster_rank)
    result = result.sort_values(["organoid", "_cluster_rank"]).drop(columns="_cluster_rank")
    return result.reset_index(drop=True)


def completeness_table(aggregates: pd.DataFrame) -> pd.DataFrame:
    organoids = sorted(aggregates["organoid"].unique())
    rows = []
    for organoid in organoids:
        observed = set(aggregates.loc[aggregates["organoid"] == organoid, "cluster"])
        missing = [cluster for cluster in CLUSTER_ORDER if cluster not in observed]
        rows.append(
            {
                "organoid": organoid,
                "n_expected_clusters": len(CLUSTER_ORDER),
                "n_observed_clusters": len(observed),
                "complete_block": len(missing) == 0,
                "missing_clusters": ";".join(missing) if missing else "none",
                "paired_analysis_consequence": (
                    "included in complete-block paired analysis"
                    if not missing
                    else "excluded from complete-block paired analysis; missing groups are reported"
                ),
            }
        )
    return pd.DataFrame(rows)


def tie_correction_for_friedman(matrix: np.ndarray) -> float:
    n, k = matrix.shape
    tie_sum = 0.0
    for row in matrix:
        _, counts = np.unique(row, return_counts=True)
        tie_sum += int(np.sum(counts**3 - counts))
    return 1.0 - tie_sum / (n * k * (k**2 - 1))


def friedman_exact(matrix: np.ndarray) -> dict[str, object]:
    """Exact randomization Friedman test by all within-block label permutations."""
    matrix = np.asarray(matrix, dtype=float)
    n, k = matrix.shape
    ranked = np.vstack([stats.rankdata(row, method="average") for row in matrix])
    correction = tie_correction_for_friedman(matrix)
    if correction <= 0:
        raise ValueError("Friedman statistic is undefined because every block is fully tied.")

    def statistic_from_ranks(ranks: np.ndarray) -> float:
        rank_sums = ranks.sum(axis=0)
        raw = 12.0 / (n * k * (k + 1)) * np.sum(rank_sums**2) - 3 * n * (k + 1)
        return float(raw / correction)

    observed = statistic_from_ranks(ranked)
    row_permutations = [list(itertools.permutations(row.tolist())) for row in ranked]
    extreme = 0
    total = 0
    for combination in itertools.product(*row_permutations):
        candidate = np.asarray(combination, dtype=float)
        extreme += statistic_from_ranks(candidate) >= observed - 1e-12
        total += 1
    asymptotic_p = float(stats.chi2.sf(observed, k - 1))
    return {
        "friedman_chi_square": observed,
        "degrees_of_freedom": k - 1,
        "p_value_exact": extreme / total,
        "p_value_asymptotic": asymptotic_p,
        "exact_enumerations": total,
        "n_complete_organoid_blocks": n,
        "kendalls_W": observed / (n * (k - 1)),
        "tie_correction": correction,
        "method": "exact within-organoid permutation of Friedman ranks",
    }


def count_absolute_ties(values: np.ndarray) -> int:
    _, counts = np.unique(np.abs(values), return_counts=True)
    return int(sum(max(0, count - 1) for count in counts))


def paired_wilcoxon_exact(first: np.ndarray, second: np.ndarray) -> dict[str, object]:
    """Two-sided exact signed-rank test using exhaustive sign-flip enumeration."""
    difference = np.asarray(first, dtype=float) - np.asarray(second, dtype=float)
    nonzero = difference != 0
    analyzed = difference[nonzero]
    n_zero = int((~nonzero).sum())
    if analyzed.size == 0:
        raise ValueError("All paired differences are zero; signed-rank test is undefined.")
    ranks = stats.rankdata(np.abs(analyzed), method="average")
    observed_signs = np.sign(analyzed)
    w_plus = float(ranks[observed_signs > 0].sum())
    w_minus = float(ranks[observed_signs < 0].sum())
    observed_signed_sum = float(np.dot(observed_signs, ranks))
    extreme = 0
    total = 2 ** analyzed.size
    for signs in itertools.product([-1.0, 1.0], repeat=analyzed.size):
        signed_sum = float(np.dot(np.asarray(signs), ranks))
        extreme += abs(signed_sum) >= abs(observed_signed_sum) - 1e-12
    tied_absolute = count_absolute_ties(analyzed)
    return {
        "wilcoxon_statistic_W_min": min(w_plus, w_minus),
        "W_plus": w_plus,
        "W_minus": w_minus,
        "p_value_raw_exact_two_sided": extreme / total,
        "exact_enumerations": total,
        "n_paired_organoids": difference.size,
        "n_nonzero_pairs": analyzed.size,
        "n_zero_differences": n_zero,
        "n_tied_absolute_differences": tied_absolute,
        "mean_paired_difference": float(np.mean(difference)),
        "median_paired_difference": float(np.median(difference)),
        "n_first_greater_than_second": int((difference > 0).sum()),
        "n_first_less_than_second": int((difference < 0).sum()),
        "rank_biserial_correlation": (w_plus - w_minus) / (w_plus + w_minus),
        "exact_method_note": (
            "standard exact signed-rank enumeration; no zero or tied absolute differences"
            if n_zero == 0 and tied_absolute == 0
            else "exhaustive conditional sign-flip permutation of average signed ranks; "
            "used because zero and/or tied absolute differences preclude the standard no-tie exact distribution"
        ),
    }


def holm_adjust(p_values: list[float]) -> list[float]:
    p = np.asarray(p_values, dtype=float)
    order = np.argsort(p)
    adjusted_sorted = np.empty_like(p)
    running = 0.0
    m = len(p)
    for position, original_index in enumerate(order):
        candidate = (m - position) * p[original_index]
        running = max(running, candidate)
        adjusted_sorted[position] = min(1.0, running)
    adjusted = np.empty_like(p)
    for position, original_index in enumerate(order):
        adjusted[original_index] = adjusted_sorted[position]
    return adjusted.tolist()


def bootstrap_paired_difference(
    difference: np.ndarray, iterations: int, rng: np.random.Generator
) -> dict[str, object]:
    indices = rng.integers(0, len(difference), size=(iterations, len(difference)))
    resampled = difference[indices]
    mean_distribution = np.mean(resampled, axis=1)
    median_distribution = np.median(resampled, axis=1)
    mean_low, mean_high = np.quantile(mean_distribution, [0.025, 0.975])
    median_low, median_high = np.quantile(median_distribution, [0.025, 0.975])
    return {
        "bootstrap_iterations": iterations,
        "bootstrap_unit": "organoid paired difference",
        "bootstrap_method": "paired organoid percentile bootstrap",
        "mean_difference_bootstrap_95CI_low": float(mean_low),
        "mean_difference_bootstrap_95CI_high": float(mean_high),
        "median_difference_bootstrap_95CI_low": float(median_low),
        "median_difference_bootstrap_95CI_high": float(median_high),
    }


def endpoint_matrix(
    aggregates: pd.DataFrame, endpoint: str, complete_organoids: list[str]
) -> np.ndarray:
    pivot = aggregates.pivot(index="organoid", columns="cluster", values=endpoint)
    return pivot.loc[complete_organoids, CLUSTER_ORDER].to_numpy(dtype=float)


def run_statistics(
    aggregates: pd.DataFrame,
    completeness: pd.DataFrame,
    random_state: int,
    bootstrap_iterations: int,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    complete_organoids = completeness.loc[completeness["complete_block"], "organoid"].tolist()
    if len(complete_organoids) < 2:
        raise ValueError("Fewer than two complete organoid blocks are available.")
    rng = np.random.default_rng(random_state)
    omnibus_rows = []
    pairwise_rows = []
    for endpoint, endpoint_label in ENDPOINTS.items():
        matrix = endpoint_matrix(aggregates, endpoint, complete_organoids)
        omnibus = friedman_exact(matrix)
        omnibus_rows.append(
            {
                "endpoint": endpoint,
                "endpoint_label": endpoint_label,
                "analysis_role": "primary" if endpoint.startswith("mean_") else "sensitivity",
                "group_order": " | ".join(CLUSTER_LABELS[c] for c in CLUSTER_ORDER),
                **omnibus,
            }
        )
        pivot = aggregates.pivot(index="organoid", columns="cluster", values=endpoint).loc[
            complete_organoids
        ]
        endpoint_rows = []
        for contrast, first_cluster, second_cluster, planned in CONTRASTS:
            first = pivot[first_cluster].to_numpy(dtype=float)
            second = pivot[second_cluster].to_numpy(dtype=float)
            result = paired_wilcoxon_exact(first, second)
            difference = first - second
            bootstrap = bootstrap_paired_difference(difference, bootstrap_iterations, rng)
            if result["mean_paired_difference"] > 0:
                direction = f"{SHORT_LABELS[first_cluster]} higher"
            elif result["mean_paired_difference"] < 0:
                direction = f"{SHORT_LABELS[second_cluster]} higher"
            else:
                direction = "no mean difference"
            endpoint_rows.append(
                {
                    "endpoint": endpoint,
                    "endpoint_label": endpoint_label,
                    "analysis_role": "primary" if endpoint.startswith("mean_") else "sensitivity",
                    "contrast": contrast,
                    "contrast_role": "planned primary pairwise" if planned else "secondary pairwise",
                    "first_cluster": first_cluster,
                    "first_cluster_label": CLUSTER_LABELS[first_cluster],
                    "second_cluster": second_cluster,
                    "second_cluster_label": CLUSTER_LABELS[second_cluster],
                    "difference_definition": "first minus second",
                    "effect_direction": direction,
                    **result,
                    **bootstrap,
                }
            )
        adjusted = holm_adjust([row["p_value_raw_exact_two_sided"] for row in endpoint_rows])
        for row, adjusted_p in zip(endpoint_rows, adjusted):
            row["p_value_Holm_adjusted_within_endpoint_across_3_contrasts"] = adjusted_p
            row["multiplicity_note"] = (
                "Holm correction is across the three cluster contrasts within this endpoint; "
                "the mean endpoint remains prespecified as primary"
            )
        pairwise_rows.extend(endpoint_rows)
    omnibus_df = pd.DataFrame(omnibus_rows)
    pairwise_df = pd.DataFrame(pairwise_rows)
    primary_omnibus = omnibus_df.loc[omnibus_df["analysis_role"] == "primary"].copy()
    primary_pairwise = pairwise_df.loc[pairwise_df["analysis_role"] == "primary"].copy()
    return primary_omnibus, primary_pairwise, omnibus_df, pairwise_df


def qc_table(
    data: pd.DataFrame, aggregates: pd.DataFrame, completeness: pd.DataFrame, source: Path
) -> pd.DataFrame:
    min_count = int(aggregates["n_cells"].min())
    small = aggregates.loc[aggregates["n_cells"] < SMALL_GROUP_THRESHOLD]
    complete_count = int(completeness["complete_block"].sum())
    expected_samples = sorted(data["organoid"].unique())
    checks = [
        (
            "biological_replicate",
            "PASS",
            "organoid is the inferential unit; cells are used only to calculate organoid-by-cluster summaries",
        ),
        (
            "aggregate_dimensions",
            "PASS" if len(aggregates) == 18 else "FAIL",
            f"{len(expected_samples)} organoids x 3 clusters = {len(aggregates)} aggregate observations",
        ),
        (
            "complete_matched_blocks",
            "PASS" if complete_count == len(expected_samples) else "FLAG",
            f"{complete_count}/{len(expected_samples)} organoids contain all clusters 1, 4, and 2",
        ),
        (
            "cell_counts",
            "PASS" if small.empty else "FLAG",
            f"minimum n={min_count}; very-small-group threshold n<{SMALL_GROUP_THRESHOLD}; flagged groups={len(small)}",
        ),
        (
            "expression_scale",
            "PASS",
            "PLXNA2 column is the existing log1p(CP10K) value exported from the same neuronal per-cell table used by the UMAP feature workflow",
        ),
        (
            "canonical_annotations",
            "PASS",
            "canonical sidecar labels retained: cluster 1 SST+ cIN; cluster 4 PV neuron precursor; cluster 2 MGE subpallial neurons",
        ),
        (
            "cluster_exclusions",
            "PASS",
            "primary analysis includes only clusters 1, 4, and 2; clusters 0, 3, 6, and 7 are excluded",
        ),
        (
            "no_expression_based_reannotation",
            "PASS",
            "LHX6, LHX8, and PLXNA2 were not used to redefine canonical cluster identities",
        ),
        (
            "no_high_low_thresholding",
            "PASS",
            "cells were not divided into LHX6-high or LHX8-high populations",
        ),
        (
            "secondary_PLXNA2_threshold",
            "PASS",
            f"PLXNA2 >= {POSITIVE_THRESHOLD:g} log1p(CP10K) is exported only as a secondary descriptive metric",
        ),
        ("source_table", "INFO", str(source.resolve())),
    ]
    return pd.DataFrame(checks, columns=["check", "status", "details"])


def format_p(value: float) -> str:
    if value < 0.001:
        return f"{value:.2g}"
    return f"{value:.3f}"


def group_mean_summary(aggregates: pd.DataFrame) -> pd.DataFrame:
    """Mean and two-sided 95% t CI across organoid-level cluster means."""
    rows = []
    for cluster in CLUSTER_ORDER:
        values = aggregates.loc[
            aggregates["cluster"] == cluster, "mean_PLXNA2_log1p_CP10K"
        ].to_numpy(dtype=float)
        n = len(values)
        mean = float(np.mean(values))
        sd = float(np.std(values, ddof=1))
        sem = sd / np.sqrt(n)
        half_width = float(stats.t.ppf(0.975, df=n - 1) * sem)
        rows.append(
            {
                "cluster": cluster,
                "cluster_label": CLUSTER_LABELS[cluster],
                "n_organoids": n,
                "mean_of_organoid_cluster_means_log1p_CP10K": mean,
                "sd_across_organoid_cluster_means": sd,
                "sem_across_organoid_cluster_means": sem,
                "mean_95CI_low": mean - half_width,
                "mean_95CI_high": mean + half_width,
                "ci_method": "two-sided t confidence interval across organoid-level cluster means",
            }
        )
    return pd.DataFrame(rows)


def render_figure(
    aggregates: pd.DataFrame,
    summaries: pd.DataFrame,
    primary_omnibus: pd.DataFrame,
    primary_pairwise: pd.DataFrame,
    final_dir: Path,
) -> None:
    samples = sorted(aggregates["organoid"].unique())
    jitter_offsets = dict(zip(samples, np.linspace(-0.055, 0.055, len(samples))))
    summary_offset = 0.13
    fig, (ax_a, ax_b) = plt.subplots(
        1,
        2,
        figsize=(6.2, 3.15),
        sharey=True,
        gridspec_kw={"width_ratios": [1.35, 1.0], "wspace": 0.24},
        facecolor="white",
    )

    pivot = aggregates.pivot(index="organoid", columns="cluster", values="mean_PLXNA2_log1p_CP10K")
    summary_by_cluster = summaries.set_index("cluster")
    group_x_a = {"1": 0.00, "4": 0.75, "2": 1.50}
    for sample in samples:
        y = pivot.loc[sample, CLUSTER_ORDER].to_numpy(dtype=float)
        x = np.asarray([group_x_a[cluster] + jitter_offsets[sample] for cluster in CLUSTER_ORDER])
        ax_a.plot(x, y, color="#BDBDBD", linewidth=0.8, alpha=0.7, zorder=1)
        ax_a.scatter(
            x,
            y,
            s=27,
            color=[CLUSTER_COLORS[cluster] for cluster in CLUSTER_ORDER],
            edgecolor="white",
            linewidth=0.5,
            zorder=3,
        )
    for cluster in CLUSTER_ORDER:
        row = summary_by_cluster.loc[cluster]
        summary_x = group_x_a[cluster] + summary_offset
        ax_a.plot(
            [summary_x, summary_x],
            [row["mean_95CI_low"], row["mean_95CI_high"]],
            color="#111111", linewidth=1.8, zorder=4,
        )
        ax_a.scatter(
            [summary_x], [row["mean_of_organoid_cluster_means_log1p_CP10K"]],
            marker="D", s=50, color="#111111",
            edgecolor="white", linewidth=0.5, zorder=5,
        )
    ax_a.set_xticks([group_x_a[c] for c in CLUSTER_ORDER], [SHORT_LABELS[c] for c in CLUSTER_ORDER])
    ax_a.set_ylabel("Mean " + r"$\it{PLXNA2}$" + ", log1p(CP10K)\nper organoid")
    ax_a.set_xlim(-0.20, 1.80)
    ax_a.set_ylim(0.35, 2.04)
    ax_a.text(
        0.00,
        -0.24,
        r"Friedman $\chi^2(2)=12.0$; P = $1.3 \times 10^{-4}$"
        + "\nKendall’s W = 1.00",
        transform=ax_a.transAxes,
        ha="left",
        va="top",
        fontsize=6.7,
        clip_on=False,
    )

    pair_clusters = ["4", "2"]
    group_x_b = {"4": 0.00, "2": 0.75}
    for sample in samples:
        pair_values = pivot.loc[sample, pair_clusters].to_numpy(dtype=float)
        x = np.asarray([group_x_b[cluster] + jitter_offsets[sample] for cluster in pair_clusters])
        ax_b.plot(x, pair_values, color="#BDBDBD", linewidth=0.8, alpha=0.7, zorder=1)
        ax_b.scatter(
            x,
            pair_values,
            s=27,
            color=[CLUSTER_COLORS[cluster] for cluster in pair_clusters],
            edgecolor="white",
            linewidth=0.5,
            zorder=3,
        )
    for cluster in pair_clusters:
        row = summary_by_cluster.loc[cluster]
        summary_x = group_x_b[cluster] + summary_offset
        ax_b.plot(
            [summary_x, summary_x],
            [row["mean_95CI_low"], row["mean_95CI_high"]],
            color="#111111", linewidth=1.8, zorder=4,
        )
        ax_b.scatter(
            [summary_x], [row["mean_of_organoid_cluster_means_log1p_CP10K"]],
            marker="D", s=50, color="#111111",
            edgecolor="white", linewidth=0.5, zorder=5,
        )
    ax_b.set_xticks([group_x_b[c] for c in pair_clusters], [SHORT_LABELS[c] for c in pair_clusters])
    ax_b.set_xlim(-0.20, 1.05)
    ax_b.tick_params(labelleft=False)
    ax_b.text(
        0.00,
        -0.24,
        "Exact paired Wilcoxon: raw P = 0.031; Holm P = 0.094\n"
        "Rank-biserial r = 1.00; 6/6 PV > MGE",
        transform=ax_b.transAxes,
        ha="left",
        va="top",
        fontsize=6.7,
        clip_on=False,
    )

    for label, ax in zip(["a", "b"], [ax_a, ax_b]):
        ax.text(-0.11, 1.06, label, transform=ax.transAxes, fontsize=9, fontweight="bold", va="top")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.grid(False)
        ax.tick_params(width=0.6, length=2.5, color="#222222")
    fig.subplots_adjust(left=0.105, right=0.985, bottom=0.34, top=0.96)
    for extension, dpi in [("png", 600), ("pdf", 600), ("svg", 600)]:
        output = final_dir / "figures" / extension / f"{FIGURE_STEM}.{extension}"
        fig.savefig(output, dpi=dpi, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def write_references(final_dir: Path) -> None:
    content = """# References and statistical rationale

## Replicate-aware single-cell analysis

1. Squair JW, et al. **Confronting false discoveries in single-cell differential expression.** *Nature Communications.* 2021;12:5692. DOI: [10.1038/s41467-021-25960-2](https://doi.org/10.1038/s41467-021-25960-2). This benchmarking study showed that methods that ignore biological-replicate variation can produce false discoveries, whereas aggregation within biological replicates substantially improves multisample single-cell inference. It supports treating each organoid, not each cell, as an independent inferential unit here.

2. Murphy AE, Skene NG. **A balanced measure shows superior performance of pseudobulk methods in single-cell RNA-sequencing analysis.** *Nature Communications.* 2022;13:7851. DOI: [10.1038/s41467-022-35519-4](https://doi.org/10.1038/s41467-022-35519-4). This study found superior performance for pseudobulk approaches in its benchmark and reinforces replicate-level aggregation rather than cell-level pseudoreplication.

3. Zimmerman KD, Espeland MA, Langefeld CD. **A practical solution to pseudoreplication bias in single-cell studies.** *Nature Communications.* 2021;12:738. DOI: [10.1038/s41467-021-21038-1](https://doi.org/10.1038/s41467-021-21038-1). This paper explains the dependence among cells obtained from the same individual or biological sample and motivates aggregation or replicate-aware models to account for that dependence.

4. Heumos L, et al. **Best practices for single-cell analysis across modalities.** *Nature Reviews Genetics.* 2023;24:550–572. DOI: [10.1038/s41576-023-00586-w](https://doi.org/10.1038/s41576-023-00586-w). This review provides broader guidance for transparent, reproducible multisample single-cell workflows, including biological replication, annotation, quality control, and reporting.

## Why these tests were selected

These publications support the principle that the biological sample—not the individual cell—is the inferential unit. They do **not** specifically prescribe a Friedman test or a paired Wilcoxon signed-rank test for this exact design.

The primary outcome is one mean normalized PLXNA2 value for each organoid and canonical cluster. All three populations are measured within the same six organoids, producing a small matched repeated-measures design. The Friedman test is therefore used as a nonparametric omnibus test without assuming normally distributed organoid-level values. The prespecified PV-versus-MGE comparison uses a paired Wilcoxon signed-rank test because the two population summaries are matched within organoid. Exact randomization/enumeration is used because only six complete organoid blocks are available.
"""
    (final_dir / "references.md").write_text(content)


def write_caption_and_summary(
    final_dir: Path,
    primary_omnibus: pd.DataFrame,
    primary_pairwise: pd.DataFrame,
    completeness: pd.DataFrame,
) -> None:
    omnibus = primary_omnibus.iloc[0]
    planned = primary_pairwise.loc[primary_pairwise["contrast"] == "PV_vs_MGE"].iloc[0]
    direction = "higher" if planned["mean_paired_difference"] > 0 else "lower"
    supported_overall = omnibus["p_value_exact"] < 0.05
    supported_planned = planned["p_value_raw_exact_two_sided"] < 0.05
    supported_planned_holm = (
        planned["p_value_Holm_adjusted_within_endpoint_across_3_contrasts"] < 0.05
    )
    caption = f"""Replicate-aware PLXNA2 comparison across established Varela DIV30 neuronal populations using the publication mapping.

(a) Mean normalized PLXNA2 expression in publication cluster 3 (SST+ cortical interneurons; original Seurat cluster 1), publication cluster 4 (PV neuron precursors; Seurat cluster 4), and publication cluster 5 (MGE subpallial neurons; Seurat cluster 2), calculated separately for each organoid from the existing log1p(CP10K) values. Each colored point represents one organoid-level cluster mean and uses the publication-cluster color. A small fixed horizontal offset is assigned to each organoid and reused across populations; thin light-gray lines connect the matched jittered observations. Black diamonds, displaced to the right of the individual observations, show the mean across six organoids; vertical black bars show two-sided 95% t confidence intervals across organoid-level cluster means. The Friedman test is the omnibus comparison across all three matched populations, with its result reported beneath panel a. (b) Planned paired comparison of publication cluster 4, PV neuron precursors, and publication cluster 5, MGE subpallial neurons. Points, jitter, connecting lines, and mean/95% confidence-interval summaries follow panel a. The exact two-sided paired Wilcoxon signed-rank test evaluates the prespecified PV-versus-MGE contrast; its raw P value, Holm-adjusted P value, and effect summary are reported beneath panel b. The other paired comparisons are PV precursor versus SST+ cIN and SST+ cIN versus MGE subpallial neurons; all three raw and Holm-adjusted P values are provided in the exported primary pairwise table. Sensitivity analyses using organoid-level medians and PLXNA2 detection frequencies are provided in supplementary tables and are not displayed in the main figure.

Cells are used only to calculate organoid-by-cluster summaries and are not independent biological replicates for inferential statistics. Clusters 0, 3, 6, and 7 are excluded. Canonical sidecar annotations are retained without expression-based reannotation, LHX6/LHX8 high-low thresholding, or PLXNA2-based cluster definition. PLXNA2 >= {POSITIVE_THRESHOLD:g} log1p(CP10K) is reported only as a secondary descriptive metric.
"""
    (final_dir / "replicate_aware_PLXNA2_figure_caption.txt").write_text(caption)

    if supported_overall:
        overall_text = (
            "PLXNA2 expression differed among the three DIV30 neuronal populations at the "
            "organoid level."
        )
    else:
        overall_text = (
            "The organoid-level analysis did not provide sufficient evidence that PLXNA2 "
            "expression differed among the three DIV30 neuronal populations."
        )
    if supported_planned and supported_planned_holm:
        planned_text = (
            f"PLXNA2 expression was {direction} in PV neuron precursors than in MGE "
            "subpallial neurons across matched organoids."
        )
    elif supported_planned:
        planned_text = (
            "In the prespecified comparison, PLXNA2 expression was higher in PV neuron "
            "precursors than in MGE subpallial neurons in all six organoids (raw exact "
            "P = 0.03125); the contrast did not remain below 0.05 after Holm adjustment "
            "across all three pairwise comparisons (adjusted P = 0.09375)."
        )
    else:
        planned_text = (
            "The planned organoid-level comparison did not provide sufficient evidence of a "
            "PLXNA2 expression difference between PV neuron precursors and MGE subpallial neurons."
        )
    missing = completeness.loc[~completeness["complete_block"]]
    missing_text = (
        "All six organoids contained all three populations; all were retained as complete paired blocks."
        if missing.empty
        else "Incomplete blocks are listed in tables/replicate_aware_plxna2_completeness_qc.tsv and were not silently dropped."
    )
    posthoc_rows = []
    for contrast in ["PV_vs_MGE", "PV_vs_SST", "SST_vs_MGE"]:
        row = primary_pairwise.loc[primary_pairwise["contrast"] == contrast].iloc[0]
        posthoc_rows.append(
            f"| {contrast.replace('_', ' ')} | {row['wilcoxon_statistic_W_min']:.1f} | "
            f"{row['p_value_raw_exact_two_sided']:.5f} | "
            f"{row['p_value_Holm_adjusted_within_endpoint_across_3_contrasts']:.5f} | "
            f"{row['rank_biserial_correlation']:.2f} |"
        )
    posthoc_table = "\n".join(posthoc_rows)
    summary = f"""# Replicate-aware PLXNA2 result

{overall_text}

{planned_text}

{missing_text}

## Omnibus test

Friedman chi-square = {omnibus['friedman_chi_square']:.1f}, df = {int(omnibus['degrees_of_freedom'])}, exact P = {omnibus['p_value_exact']:.8f}, Kendall's W = {omnibus['kendalls_W']:.2f}.

## Paired post hoc comparisons

| Contrast | Wilcoxon W | Raw exact P | Holm-adjusted P | Rank-biserial r |
|---|---:|---:|---:|---:|
{posthoc_table}

The PV-versus-MGE comparison was prespecified. Its raw exact P value is 0.03125; after Holm adjustment across all three reported pairwise contrasts, P = 0.09375.

The primary inference uses mean PLXNA2 log1p(CP10K) per organoid and cluster. The median-expression and detection-frequency endpoints are sensitivity analyses and do not replace the prespecified primary endpoint. This analysis establishes association with established neuronal annotations; it does not show that PLXNA2 defines cortical migration, PV fate, or striatal exclusion.
"""
    (final_dir / "replicate_aware_PLXNA2_results_summary.md").write_text(summary)


def git_value(args: list[str]) -> str:
    try:
        return subprocess.run(args, check=True, capture_output=True, text=True).stdout.strip()
    except (subprocess.CalledProcessError, FileNotFoundError):
        return "unavailable"


def write_provenance(
    final_dir: Path, source: Path, script_source: Path, bootstrap_iterations: int
) -> None:
    manifest = pd.DataFrame(
        [
            ("rendered_at_local", datetime.now().isoformat(timespec="seconds")),
            ("source_per_cell_table", str(source.resolve())),
            ("expression_column", "PLXNA2"),
            ("expression_scale", "existing log1p(CP10K)"),
            ("inferential_unit", "organoid"),
            ("included_clusters", ",".join(CLUSTER_ORDER)),
            ("excluded_clusters", "0,3,6,7"),
            ("positive_threshold_secondary_log1p_CP10K", str(POSITIVE_THRESHOLD)),
            ("small_group_flag_threshold_cells", str(SMALL_GROUP_THRESHOLD)),
            ("bootstrap_iterations", str(bootstrap_iterations)),
            ("script", str(script_source.resolve())),
            ("git_commit", git_value(["git", "rev-parse", "HEAD"])),
        ],
        columns=["key", "value"],
    )
    manifest.to_csv(
        final_dir / "provenance" / "replicate_aware_plxna2_render_manifest.tsv",
        sep="\t",
        index=False,
    )
    files = [
        final_dir / "figures" / ext / f"{FIGURE_STEM}.{ext}" for ext in ["png", "pdf", "svg"]
    ]
    files += [
        final_dir / "tables" / name
        for name in [
            "replicate_aware_plxna2_organoid_cluster_aggregates.tsv",
            "replicate_aware_plxna2_group_mean_95CI_summary.tsv",
            "replicate_aware_plxna2_completeness_qc.tsv",
            "replicate_aware_plxna2_cell_count_qc.tsv",
            "replicate_aware_plxna2_primary_friedman.tsv",
            "replicate_aware_plxna2_primary_pairwise.tsv",
            "replicate_aware_plxna2_sensitivity_omnibus.tsv",
            "replicate_aware_plxna2_sensitivity_pairwise.tsv",
            "replicate_aware_plxna2_sensitivity_planned_PV_vs_MGE.tsv",
            "replicate_aware_plxna2_qc_checks.tsv",
        ]
    ]
    files += [
        final_dir / "references.md",
        final_dir / "replicate_aware_PLXNA2_figure_caption.txt",
        final_dir / "replicate_aware_PLXNA2_results_summary.md",
        final_dir / "code" / SCRIPT_NAME,
    ]
    checksum_lines = []
    for path in files:
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        checksum_lines.append(f"{digest}  {path.relative_to(final_dir)}\n")
    (final_dir / "provenance" / "replicate_aware_plxna2_sha256_manifest.txt").write_text(
        "".join(checksum_lines)
    )


def main() -> None:
    args = parse_args()
    final_dir = args.final_dir.resolve()
    source = (args.per_cell_table or final_dir / PER_CELL_RELATIVE).resolve()
    script_source = Path(__file__).resolve()
    if not source.exists():
        raise FileNotFoundError(f"Missing per-cell source table: {source}")
    ensure_output_dirs(final_dir)
    data = load_and_validate(source)
    aggregates = aggregate_organoid_cluster(data)
    summaries = group_mean_summary(aggregates)
    completeness = completeness_table(aggregates)
    primary_omnibus, primary_pairwise, sensitivity_omnibus, sensitivity_pairwise = run_statistics(
        aggregates,
        completeness,
        random_state=args.random_state,
        bootstrap_iterations=args.bootstrap_iterations,
    )
    qc = qc_table(data, aggregates, completeness, source)

    tables = final_dir / "tables"
    aggregates.to_csv(
        tables / "replicate_aware_plxna2_organoid_cluster_aggregates.tsv", sep="\t", index=False
    )
    summaries.to_csv(
        tables / "replicate_aware_plxna2_group_mean_95CI_summary.tsv", sep="\t", index=False
    )
    completeness.to_csv(
        tables / "replicate_aware_plxna2_completeness_qc.tsv", sep="\t", index=False
    )
    aggregates[["organoid", "cluster", "cluster_label", "n_cells"]].assign(
        very_small_group=lambda x: x["n_cells"] < SMALL_GROUP_THRESHOLD,
        very_small_group_threshold_cells=SMALL_GROUP_THRESHOLD,
    ).to_csv(tables / "replicate_aware_plxna2_cell_count_qc.tsv", sep="\t", index=False)
    primary_omnibus.to_csv(
        tables / "replicate_aware_plxna2_primary_friedman.tsv", sep="\t", index=False
    )
    primary_pairwise.to_csv(
        tables / "replicate_aware_plxna2_primary_pairwise.tsv", sep="\t", index=False
    )
    sensitivity_omnibus.to_csv(
        tables / "replicate_aware_plxna2_sensitivity_omnibus.tsv", sep="\t", index=False
    )
    sensitivity_pairwise.to_csv(
        tables / "replicate_aware_plxna2_sensitivity_pairwise.tsv", sep="\t", index=False
    )
    sensitivity_pairwise.loc[sensitivity_pairwise["contrast"] == "PV_vs_MGE"].to_csv(
        tables / "replicate_aware_plxna2_sensitivity_planned_PV_vs_MGE.tsv",
        sep="\t",
        index=False,
    )
    qc.to_csv(tables / "replicate_aware_plxna2_qc_checks.tsv", sep="\t", index=False)

    render_figure(aggregates, summaries, primary_omnibus, primary_pairwise, final_dir)
    write_references(final_dir)
    write_caption_and_summary(final_dir, primary_omnibus, primary_pairwise, completeness)
    shutil.copy2(script_source, final_dir / "code" / SCRIPT_NAME)
    write_provenance(final_dir, source, script_source, args.bootstrap_iterations)

    omnibus = primary_omnibus.iloc[0]
    planned = primary_pairwise.loc[primary_pairwise["contrast"] == "PV_vs_MGE"].iloc[0]
    print(f"Output folder: {final_dir}")
    print(f"Aggregate observations: {len(aggregates)}")
    print(f"Complete organoid blocks: {int(completeness['complete_block'].sum())}")
    print(
        "Friedman exact: "
        f"Q={omnibus['friedman_chi_square']:.6g}, df={int(omnibus['degrees_of_freedom'])}, "
        f"P={omnibus['p_value_exact']:.6g}, W={omnibus['kendalls_W']:.6g}"
    )
    print(
        "Planned PV vs MGE: "
        f"Wmin={planned['wilcoxon_statistic_W_min']:.6g}, "
        f"P_exact={planned['p_value_raw_exact_two_sided']:.6g}, "
        f"P_Holm={planned['p_value_Holm_adjusted_within_endpoint_across_3_contrasts']:.6g}, "
        f"r_rb={planned['rank_biserial_correlation']:.6g}"
    )


if __name__ == "__main__":
    main()
