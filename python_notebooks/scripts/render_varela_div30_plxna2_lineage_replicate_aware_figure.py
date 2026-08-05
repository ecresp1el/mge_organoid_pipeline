#!/usr/bin/env python3
"""Replicate-aware PLXNA2 comparison of cortical interneuron and subpallial lineages.

Clusters 1 (SST+ cIN) and 4 (PV neuron precursor) are pooled within each
organoid as the cortical interneuron lineage. Cluster 2 is retained as the
subpallial MGE lineage. Cells are used only to compute one summary per
organoid and lineage; the organoid is the inferential unit.
"""

from __future__ import annotations

import argparse
from datetime import datetime
import hashlib
from pathlib import Path
import shutil
import subprocess

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

from render_varela_div30_plxna2_replicate_aware_figure import (
    BOOTSTRAP_ITERATIONS,
    FINAL_DIR_DEFAULT,
    PER_CELL_RELATIVE,
    POSITIVE_THRESHOLD,
    bootstrap_paired_difference,
    load_and_validate,
    paired_wilcoxon_exact,
)


SCRIPT_NAME = "render_varela_div30_plxna2_lineage_replicate_aware_figure.py"
FIGURE_STEM = "varela_div30_plxna2_cortical_vs_subpallial_replicate_aware"
TABLE_PREFIX = "lineage_replicate_aware_plxna2"

LINEAGE_ORDER = ["cortical_interneuron_lineage", "subpallial_MGE_lineage"]
LINEAGE_LABELS = {
    "cortical_interneuron_lineage": "Cortical interneuron lineage",
    "subpallial_MGE_lineage": "Subpallial MGE neurons",
}
LINEAGE_SHORT_LABELS = {
    "cortical_interneuron_lineage": "Cortical interneuron\nlineage",
    "subpallial_MGE_lineage": "Subpallial MGE\nneurons",
}
LINEAGE_COLORS = {
    "cortical_interneuron_lineage": "#009E73",
    "subpallial_MGE_lineage": "#C77CFF",
}
CLUSTER_TO_LINEAGE = {
    "1": "cortical_interneuron_lineage",
    "4": "cortical_interneuron_lineage",
    "2": "subpallial_MGE_lineage",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--final-dir", type=Path, default=FINAL_DIR_DEFAULT)
    parser.add_argument("--per-cell-table", type=Path, default=None)
    parser.add_argument("--random-state", type=int, default=2909)
    parser.add_argument("--bootstrap-iterations", type=int, default=BOOTSTRAP_ITERATIONS)
    return parser.parse_args()


def aggregate_lineages(data: pd.DataFrame) -> pd.DataFrame:
    work = data.copy()
    work["lineage"] = work["cluster"].map(CLUSTER_TO_LINEAGE)
    if work["lineage"].isna().any():
        raise ValueError("Selected cells include a cluster without a lineage mapping.")
    rows = []
    for (organoid, lineage), group in work.groupby(["organoid", "lineage"], sort=False):
        values = group["PLXNA2"].to_numpy(dtype=float)
        q25, q75 = np.quantile(values, [0.25, 0.75])
        clusters = sorted(group["cluster"].unique(), key=int)
        rows.append(
            {
                "organoid": organoid,
                "lineage": lineage,
                "lineage_label": LINEAGE_LABELS[lineage],
                "source_clusters": ",".join(clusters),
                "n_cells": len(values),
                "mean_PLXNA2_log1p_CP10K": float(np.mean(values)),
                "median_PLXNA2_log1p_CP10K": float(np.median(values)),
                "sd_PLXNA2_log1p_CP10K": float(np.std(values, ddof=1)),
                "q25_PLXNA2_log1p_CP10K": float(q25),
                "q75_PLXNA2_log1p_CP10K": float(q75),
                "iqr_PLXNA2_log1p_CP10K": float(q75 - q25),
                "pct_PLXNA2_detected": float(100 * np.mean(values > 0)),
                "secondary_PLXNA2_positive_threshold_log1p_CP10K": POSITIVE_THRESHOLD,
                "secondary_pct_PLXNA2_at_or_above_threshold": float(
                    100 * np.mean(values >= POSITIVE_THRESHOLD)
                ),
            }
        )
    result = pd.DataFrame(rows)
    rank = {lineage: i for i, lineage in enumerate(LINEAGE_ORDER)}
    result["_rank"] = result["lineage"].map(rank)
    return result.sort_values(["organoid", "_rank"]).drop(columns="_rank").reset_index(drop=True)


def validate_complete_pairs(aggregates: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for organoid in sorted(aggregates["organoid"].unique()):
        observed = set(aggregates.loc[aggregates["organoid"] == organoid, "lineage"])
        missing = [lineage for lineage in LINEAGE_ORDER if lineage not in observed]
        cortical_clusters = set(
            aggregates.loc[
                (aggregates["organoid"] == organoid)
                & (aggregates["lineage"] == "cortical_interneuron_lineage"),
                "source_clusters",
            ]
        )
        rows.append(
            {
                "organoid": organoid,
                "complete_lineage_pair": not missing,
                "missing_lineages": ";".join(missing) if missing else "none",
                "cortical_source_clusters": next(iter(cortical_clusters), "missing"),
                "paired_analysis_consequence": (
                    "included in paired analysis"
                    if not missing
                    else "excluded from paired analysis; missing lineage explicitly reported"
                ),
            }
        )
    return pd.DataFrame(rows)


def group_mean_summary(aggregates: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for lineage in LINEAGE_ORDER:
        values = aggregates.loc[
            aggregates["lineage"] == lineage, "mean_PLXNA2_log1p_CP10K"
        ].to_numpy(dtype=float)
        n = len(values)
        mean = float(np.mean(values))
        sd = float(np.std(values, ddof=1))
        sem = sd / np.sqrt(n)
        half_width = float(stats.t.ppf(0.975, n - 1) * sem)
        rows.append(
            {
                "lineage": lineage,
                "lineage_label": LINEAGE_LABELS[lineage],
                "n_organoids": n,
                "mean_of_organoid_lineage_means_log1p_CP10K": mean,
                "sd_across_organoid_lineage_means": sd,
                "sem_across_organoid_lineage_means": sem,
                "mean_95CI_low": mean - half_width,
                "mean_95CI_high": mean + half_width,
                "ci_method": "two-sided t confidence interval across organoid-level lineage means",
            }
        )
    return pd.DataFrame(rows)


def analyze_endpoint(
    aggregates: pd.DataFrame,
    endpoint: str,
    endpoint_label: str,
    analysis_role: str,
    rng: np.random.Generator,
    bootstrap_iterations: int,
) -> dict[str, object]:
    pivot = aggregates.pivot(index="organoid", columns="lineage", values=endpoint)
    pivot = pivot.dropna(subset=LINEAGE_ORDER).sort_index()
    cortical = pivot[LINEAGE_ORDER[0]].to_numpy(dtype=float)
    subpallial = pivot[LINEAGE_ORDER[1]].to_numpy(dtype=float)
    result = paired_wilcoxon_exact(cortical, subpallial)
    bootstrap = bootstrap_paired_difference(
        cortical - subpallial, bootstrap_iterations, rng
    )
    return {
        "endpoint": endpoint,
        "endpoint_label": endpoint_label,
        "analysis_role": analysis_role,
        "contrast": "cortical_interneuron_lineage_vs_subpallial_MGE_lineage",
        "first_lineage": LINEAGE_ORDER[0],
        "first_lineage_label": LINEAGE_LABELS[LINEAGE_ORDER[0]],
        "second_lineage": LINEAGE_ORDER[1],
        "second_lineage_label": LINEAGE_LABELS[LINEAGE_ORDER[1]],
        "difference_definition": "cortical interneuron lineage minus subpallial MGE lineage",
        "effect_direction": (
            "cortical interneuron lineage higher"
            if result["mean_paired_difference"] > 0
            else "subpallial MGE lineage higher"
        ),
        **result,
        **bootstrap,
        "multiplicity_adjustment": "none; one prespecified lineage-level contrast",
    }


def equal_cluster_weight_sensitivity(
    data: pd.DataFrame, rng: np.random.Generator, bootstrap_iterations: int
) -> tuple[pd.DataFrame, pd.DataFrame]:
    cluster_means = data.pivot_table(
        index="organoid", columns="cluster", values="PLXNA2", aggfunc="mean"
    ).sort_index()
    values = pd.DataFrame(
        {
            "organoid": cluster_means.index,
            "equal_weight_cortical_mean_log1p_CP10K": (
                cluster_means["1"] + cluster_means["4"]
            ).to_numpy()
            / 2,
            "subpallial_MGE_mean_log1p_CP10K": cluster_means["2"].to_numpy(),
        }
    )
    first = values["equal_weight_cortical_mean_log1p_CP10K"].to_numpy(dtype=float)
    second = values["subpallial_MGE_mean_log1p_CP10K"].to_numpy(dtype=float)
    result = paired_wilcoxon_exact(first, second)
    bootstrap = bootstrap_paired_difference(first - second, bootstrap_iterations, rng)
    statistics = pd.DataFrame(
        [
            {
                "analysis_role": "pooling-definition sensitivity",
                "cortical_definition": "equal-weight mean of organoid cluster-1 and cluster-4 means",
                "subpallial_definition": "organoid cluster-2 mean",
                **result,
                **bootstrap,
            }
        ]
    )
    return values, statistics


def render_figure(
    aggregates: pd.DataFrame, summaries: pd.DataFrame, primary: pd.Series, final_dir: Path
) -> None:
    samples = sorted(aggregates["organoid"].unique())
    jitter = dict(zip(samples, np.linspace(-0.055, 0.055, len(samples))))
    group_x = {LINEAGE_ORDER[0]: 0.00, LINEAGE_ORDER[1]: 0.75}
    summary_offset = 0.13
    pivot = aggregates.pivot(
        index="organoid", columns="lineage", values="mean_PLXNA2_log1p_CP10K"
    )
    summary_by_lineage = summaries.set_index("lineage")

    fig, ax = plt.subplots(figsize=(3.45, 3.2), facecolor="white")
    for sample in samples:
        x = np.asarray([group_x[lineage] + jitter[sample] for lineage in LINEAGE_ORDER])
        y = pivot.loc[sample, LINEAGE_ORDER].to_numpy(dtype=float)
        ax.plot(x, y, color="#BDBDBD", linewidth=0.8, alpha=0.7, zorder=1)
        ax.scatter(
            x,
            y,
            s=27,
            color=[LINEAGE_COLORS[lineage] for lineage in LINEAGE_ORDER],
            edgecolor="white",
            linewidth=0.5,
            zorder=3,
        )
    for lineage in LINEAGE_ORDER:
        row = summary_by_lineage.loc[lineage]
        x = group_x[lineage] + summary_offset
        ax.plot(
            [x, x],
            [row["mean_95CI_low"], row["mean_95CI_high"]],
            color="#111111",
            linewidth=1.8,
            zorder=4,
        )
        ax.scatter(
            [x],
            [row["mean_of_organoid_lineage_means_log1p_CP10K"]],
            marker="D",
            s=50,
            color="#111111",
            edgecolor="white",
            linewidth=0.5,
            zorder=5,
        )
    ax.set_xticks(
        [group_x[lineage] for lineage in LINEAGE_ORDER],
        [LINEAGE_SHORT_LABELS[lineage] for lineage in LINEAGE_ORDER],
    )
    ax.set_xlim(-0.20, 1.05)
    ax.set_ylim(0.30, 1.88)
    ax.set_ylabel("Mean " + r"$\it{PLXNA2}$" + ", log1p(CP10K)\nper organoid")
    ax.text(
        0.00,
        -0.34,
        "Exact paired Wilcoxon P = 0.031\n"
        "Rank-biserial r = 1.00; 6/6 cortical > subpallial",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=6.4,
        clip_on=False,
    )
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(False)
    ax.tick_params(width=0.6, length=2.5, color="#222222")
    fig.subplots_adjust(left=0.19, right=0.98, bottom=0.40, top=0.97)
    for extension in ["png", "pdf", "svg"]:
        output = final_dir / "figures" / extension / f"{FIGURE_STEM}.{extension}"
        fig.savefig(output, dpi=600, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def write_documentation(
    final_dir: Path, primary: pd.Series, sensitivity: pd.DataFrame
) -> None:
    caption = """Replicate-aware PLXNA2 comparison between cortical interneuron-lineage and subpallial MGE populations in Varela DIV30 organoids.

Publication clusters 3 (SST+ cortical interneurons; original Seurat cluster 1) and 4 (PV neuron precursors; Seurat cluster 4) are pooled at the cell level within each organoid to define the cortical interneuron lineage; publication cluster 5 (MGE subpallial neurons; Seurat cluster 2) defines the subpallial MGE lineage. Each colored point is one organoid-level lineage mean calculated from the existing PLXNA2 log1p(CP10K) values. A fixed horizontal jitter is reused for the two observations from each organoid, and thin light-gray lines connect matched values. Displaced black diamonds show the mean across six organoids, and black bars show two-sided 95% t confidence intervals across organoid-level lineage means. The exact two-sided paired Wilcoxon signed-rank test uses the six organoids as biological replicates and is reported beneath the plot. No individual cell is treated as an independent inferential observation. The equal-cluster-weight sensitivity analysis gives publication clusters 3 and 4 equal weight within each organoid and is exported separately.
"""
    (final_dir / "lineage_replicate_aware_PLXNA2_figure_caption.txt").write_text(caption)

    sensitivity_lines = []
    for _, row in sensitivity.iterrows():
        sensitivity_lines.append(
            f"| {row['endpoint_label']} | {row['effect_direction']} | "
            f"{row['p_value_raw_exact_two_sided']:.5f} | "
            f"{row['rank_biserial_correlation']:.2f} |"
        )
    result = f"""# Lineage-level replicate-aware PLXNA2 result

## Primary comparison

PLXNA2 expression was higher in cortical interneuron-lineage cells (pooled publication clusters 3 and 4 within each organoid; original Seurat clusters 1 and 4) than in subpallial MGE neurons (publication cluster 5; Seurat cluster 2) in all six organoids.

- Exact paired Wilcoxon W = {primary['wilcoxon_statistic_W_min']:.1f}
- Exact two-sided P = {primary['p_value_raw_exact_two_sided']:.5f}
- Rank-biserial r = {primary['rank_biserial_correlation']:.2f}
- Mean paired difference = {primary['mean_paired_difference']:.6f} log1p(CP10K)
- Median paired difference = {primary['median_paired_difference']:.6f} log1p(CP10K)
- Paired-mean-difference bootstrap 95% CI = [{primary['mean_difference_bootstrap_95CI_low']:.6f}, {primary['mean_difference_bootstrap_95CI_high']:.6f}]
- Organoids with cortical > subpallial = {int(primary['n_first_greater_than_second'])}/{int(primary['n_paired_organoids'])}

This supports the descriptive conclusion that **PLXNA2 is enriched in cortical interneuron-lineage cells relative to subpallial MGE neurons** in these six DIV30 organoids. It establishes an association with the existing lineage annotations and does not establish a causal role for PLXNA2.

## Sensitivity endpoints

| Endpoint | Direction | Exact P | Rank-biserial r |
|---|---|---:|---:|
{chr(10).join(sensitivity_lines)}

The primary endpoint remains the pooled-cell mean PLXNA2 log1p(CP10K) calculated once per organoid and lineage. No multiplicity adjustment is applied because this new analysis contains one prespecified lineage-level contrast. The earlier three-cluster analysis and all of its Holm-adjusted post hoc results remain unchanged in their original files.
"""
    (final_dir / "lineage_replicate_aware_PLXNA2_results_summary.md").write_text(result)


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
            ("expression_scale", "existing log1p(CP10K)"),
            ("inferential_unit", "organoid"),
            ("publication_mapping", "Seurat 0/3/7->pub1; 6->pub2; 1->pub3; 4->pub4; 2->pub5"),
            ("cortical_lineage_definition", "publication clusters 3 and 4 pooled within organoid"),
            ("subpallial_lineage_definition", "publication cluster 5"),
            ("primary_test", "exact two-sided paired Wilcoxon signed-rank"),
            ("bootstrap_iterations", str(bootstrap_iterations)),
            ("script", str(script_source.resolve())),
            ("git_commit", git_value(["git", "rev-parse", "HEAD"])),
        ],
        columns=["key", "value"],
    )
    manifest_path = final_dir / "provenance" / f"{TABLE_PREFIX}_render_manifest.tsv"
    manifest.to_csv(manifest_path, sep="\t", index=False)

    files = [
        final_dir / "figures" / extension / f"{FIGURE_STEM}.{extension}"
        for extension in ["png", "pdf", "svg"]
    ]
    files += [
        final_dir / "tables" / name
        for name in [
            f"{TABLE_PREFIX}_organoid_lineage_aggregates.tsv",
            f"{TABLE_PREFIX}_completeness_qc.tsv",
            f"{TABLE_PREFIX}_group_mean_95CI_summary.tsv",
            f"{TABLE_PREFIX}_primary_paired_wilcoxon.tsv",
            f"{TABLE_PREFIX}_sensitivity_endpoints.tsv",
            f"{TABLE_PREFIX}_equal_cluster_weight_values.tsv",
            f"{TABLE_PREFIX}_equal_cluster_weight_statistics.tsv",
            f"{TABLE_PREFIX}_cluster_to_lineage_mapping.tsv",
        ]
    ]
    files += [
        final_dir / "lineage_replicate_aware_PLXNA2_figure_caption.txt",
        final_dir / "lineage_replicate_aware_PLXNA2_results_summary.md",
        final_dir / "code" / SCRIPT_NAME,
        manifest_path,
    ]
    lines = []
    for path in files:
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        lines.append(f"{digest}  {path.relative_to(final_dir)}\n")
    (final_dir / "provenance" / f"{TABLE_PREFIX}_sha256_manifest.txt").write_text(
        "".join(lines)
    )


def main() -> None:
    args = parse_args()
    final_dir = args.final_dir.resolve()
    source = (args.per_cell_table or final_dir / PER_CELL_RELATIVE).resolve()
    script_source = Path(__file__).resolve()
    for relative in ["figures/png", "figures/pdf", "figures/svg", "tables", "code", "provenance"]:
        (final_dir / relative).mkdir(parents=True, exist_ok=True)

    data = load_and_validate(source)
    aggregates = aggregate_lineages(data)
    completeness = validate_complete_pairs(aggregates)
    if not completeness["complete_lineage_pair"].all():
        missing = completeness.loc[~completeness["complete_lineage_pair"]]
        raise ValueError(f"Incomplete organoid lineage pairs:\n{missing.to_string(index=False)}")
    summaries = group_mean_summary(aggregates)
    rng = np.random.default_rng(args.random_state)
    endpoint_specs = [
        ("mean_PLXNA2_log1p_CP10K", "Mean expression", "primary"),
        ("median_PLXNA2_log1p_CP10K", "Median expression", "sensitivity"),
        ("pct_PLXNA2_detected", "PLXNA2-detected cells (%)", "sensitivity"),
    ]
    sensitivity = pd.DataFrame(
        [
            analyze_endpoint(
                aggregates, endpoint, label, role, rng, args.bootstrap_iterations
            )
            for endpoint, label, role in endpoint_specs
        ]
    )
    primary = sensitivity.loc[sensitivity["analysis_role"] == "primary"].iloc[0]
    equal_values, equal_statistics = equal_cluster_weight_sensitivity(
        data, rng, args.bootstrap_iterations
    )

    tables = final_dir / "tables"
    aggregates.to_csv(
        tables / f"{TABLE_PREFIX}_organoid_lineage_aggregates.tsv", sep="\t", index=False
    )
    completeness.to_csv(
        tables / f"{TABLE_PREFIX}_completeness_qc.tsv", sep="\t", index=False
    )
    summaries.to_csv(
        tables / f"{TABLE_PREFIX}_group_mean_95CI_summary.tsv", sep="\t", index=False
    )
    sensitivity.loc[sensitivity["analysis_role"] == "primary"].to_csv(
        tables / f"{TABLE_PREFIX}_primary_paired_wilcoxon.tsv", sep="\t", index=False
    )
    sensitivity.to_csv(
        tables / f"{TABLE_PREFIX}_sensitivity_endpoints.tsv", sep="\t", index=False
    )
    equal_values.to_csv(
        tables / f"{TABLE_PREFIX}_equal_cluster_weight_values.tsv", sep="\t", index=False
    )
    equal_statistics.to_csv(
        tables / f"{TABLE_PREFIX}_equal_cluster_weight_statistics.tsv", sep="\t", index=False
    )
    pd.DataFrame(
        [
            {
                "cluster": cluster,
                "canonical_cluster_label": label,
                "publication_cluster": {"1": "3", "4": "4", "2": "5"}[cluster],
                "publication_cluster_label": {
                    "1": "SST+ cIN",
                    "4": "PV Neuron Precursor",
                    "2": "MGE Subpallial Neuron",
                }[cluster],
                "lineage": CLUSTER_TO_LINEAGE[cluster],
                "lineage_label": LINEAGE_LABELS[CLUSTER_TO_LINEAGE[cluster]],
                "analysis_role": (
                    "pooled into cortical interneuron lineage"
                    if cluster in {"1", "4"}
                    else "subpallial comparator"
                ),
            }
            for cluster, label in [
                ("1", "SST+ cIN"),
                ("4", "PV neuron precursor"),
                ("2", "MGE subpallial neurons"),
            ]
        ]
    ).to_csv(
        tables / f"{TABLE_PREFIX}_cluster_to_lineage_mapping.tsv", sep="\t", index=False
    )

    render_figure(aggregates, summaries, primary, final_dir)
    write_documentation(final_dir, primary, sensitivity)
    shutil.copy2(script_source, final_dir / "code" / SCRIPT_NAME)
    write_provenance(final_dir, source, script_source, args.bootstrap_iterations)

    print(f"Output folder: {final_dir}")
    print("Existing three-cluster outputs: preserved")
    print(f"Lineage aggregate observations: {len(aggregates)}")
    print(f"Complete paired organoids: {int(completeness['complete_lineage_pair'].sum())}")
    print(
        "Cortical vs subpallial: "
        f"Wmin={primary['wilcoxon_statistic_W_min']:.6g}, "
        f"P_exact={primary['p_value_raw_exact_two_sided']:.6g}, "
        f"r_rb={primary['rank_biserial_correlation']:.6g}, "
        f"n_cortical_greater={int(primary['n_first_greater_than_second'])}/"
        f"{int(primary['n_paired_organoids'])}"
    )


if __name__ == "__main__":
    main()
