#!/usr/bin/env python3
"""Summarize marker coexpression at multiple cutoffs from exported per-cell tables."""

from __future__ import annotations

import argparse
import gzip
from pathlib import Path
from textwrap import wrap

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


STATE_SPECS = [
    ("dcx_mki67", "DCX+ Ki67+", ["DCX", "MKI67"]),
    ("dcx_mki67_sox2", "DCX+ Ki67+ SOX2+", ["DCX", "MKI67", "SOX2"]),
    (
        "dcx_mki67_sox2_pcdh19",
        "DCX+ Ki67+ SOX2+ PCDH19+",
        ["DCX", "MKI67", "SOX2", "PCDH19"],
    ),
]

STUDY_SPECS = [
    {
        "study_id": "varela_div30",
        "study_label": "This Study, DIV 30",
        "prefix": "div30_marker",
        "sample_order": ["9583-MW-1", "9583-MW-2", "9583-MW-3", "9583-MW-4", "9583-MW-5", "9583-MW-6"],
    },
    {
        "study_id": "siebert_2026",
        "study_label": "Siebert 2026",
        "prefix": "siebert_marker",
        "sample_order": ["Young_1", "Young_2", "Old_1", "Old_2"],
    },
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--analysis-dir",
        default="/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/siebert_2026/analysis/siebert_marker_coexpression_umap_v3_sample_breakdown",
        help="Directory containing tables/ and plots/ from the marker coexpression analysis.",
    )
    parser.add_argument("--cutoffs", default="0.5,1.0", help="Comma-separated expression cutoffs.")
    return parser.parse_args()


def read_per_cell(path: Path) -> pd.DataFrame:
    usecols = ["cell_id", "sample", "DCX", "MKI67", "SOX2", "PCDH19"]
    with gzip.open(path, "rt") as handle:
      header = handle.readline().rstrip("\n").split("\t")
    cols = [col for col in usecols if col in header]
    missing = sorted(set(usecols) - set(cols))
    if missing:
        raise ValueError(f"{path} is missing required columns: {missing}")
    return pd.read_csv(path, sep="\t", usecols=cols)


def summarize_one(df: pd.DataFrame, study: dict[str, object], cutoff: float) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    sample_rows: list[dict[str, object]] = []
    overall_rows: list[dict[str, object]] = []
    sample_order = [s for s in study["sample_order"] if s in set(df["sample"].astype(str))]
    extras = sorted(set(df["sample"].astype(str)) - set(sample_order))
    samples = sample_order + extras

    for state_id, state_label, genes in STATE_SPECS:
        mask = np.logical_and.reduce([(df[gene].to_numpy() >= cutoff) for gene in genes])
        n_positive = int(mask.sum())
        overall_rows.append(
            {
                "study_id": study["study_id"],
                "study_label": study["study_label"],
                "threshold": cutoff,
                "state_id": state_id,
                "state_label": state_label,
                "genes": ",".join(genes),
                "n_positive": n_positive,
                "n_total": int(len(df)),
                "percent_cells": 100 * n_positive / len(df) if len(df) else np.nan,
            }
        )
        for sample in samples:
            sample_mask = df["sample"].astype(str).to_numpy() == sample
            n_total = int(sample_mask.sum())
            n_sample_positive = int((mask & sample_mask).sum())
            sample_rows.append(
                {
                    "study_id": study["study_id"],
                    "study_label": study["study_label"],
                    "sample": sample,
                    "threshold": cutoff,
                    "state_id": state_id,
                    "state_label": state_label,
                    "genes": ",".join(genes),
                    "n_positive": n_sample_positive,
                    "n_total": n_total,
                    "percent_cells": 100 * n_sample_positive / n_total if n_total else np.nan,
                }
            )
    return sample_rows, overall_rows


def write_tsv(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=False)


def label_wrap(text: str, width: int = 18) -> str:
    return "\n".join(wrap(text, width=width))


def plot_sample_heatmap(sample_df: pd.DataFrame, out_base: Path) -> None:
    states = [label for _, label, _ in STATE_SPECS]
    cutoffs = sorted(sample_df["threshold"].unique())
    row_records = []
    for study in STUDY_SPECS:
        study_df = sample_df[sample_df["study_id"] == study["study_id"]]
        order = [s for s in study["sample_order"] if s in set(study_df["sample"])]
        order += sorted(set(study_df["sample"]) - set(order))
        for sample in order:
            row_records.append((study["study_label"], sample))

    fig, axes = plt.subplots(1, len(cutoffs), figsize=(6.2 * len(cutoffs), 5.7), constrained_layout=True)
    if len(cutoffs) == 1:
        axes = [axes]
    vmax = max(1.0, float(sample_df["percent_cells"].max()))

    for ax, cutoff in zip(axes, cutoffs):
        matrix = np.full((len(row_records), len(states)), np.nan)
        for row_idx, (study_label, sample) in enumerate(row_records):
            sub = sample_df[
                (sample_df["study_label"] == study_label)
                & (sample_df["sample"] == sample)
                & (sample_df["threshold"] == cutoff)
            ]
            for col_idx, state in enumerate(states):
                vals = sub.loc[sub["state_label"] == state, "percent_cells"]
                if len(vals):
                    matrix[row_idx, col_idx] = vals.iloc[0]

        im = ax.imshow(matrix, aspect="auto", cmap="YlGnBu", vmin=0, vmax=vmax)
        ax.set_title(f"Cutoff >= {cutoff:g}", fontsize=12, fontweight="bold")
        ax.set_xticks(range(len(states)))
        ax.set_xticklabels([label_wrap(state, 16) for state in states], rotation=35, ha="right", fontsize=8)
        ax.set_yticks(range(len(row_records)))
        ax.set_yticklabels([f"{study}\n{sample}" for study, sample in row_records], fontsize=8)
        ax.tick_params(length=0)
        for row_idx in range(matrix.shape[0]):
            for col_idx in range(matrix.shape[1]):
                value = matrix[row_idx, col_idx]
                if np.isfinite(value):
                    color = "white" if value > 0.55 * vmax else "black"
                    ax.text(col_idx, row_idx, f"{value:.1f}%", ha="center", va="center", fontsize=7.5, color=color, fontweight="bold")
        for boundary in [idx - 0.5 for idx, (study, _) in enumerate(row_records) if idx > 0 and study != row_records[idx - 1][0]]:
            ax.axhline(boundary, color="black", linewidth=1.2)

    cbar = fig.colorbar(im, ax=axes, shrink=0.78, pad=0.02)
    cbar.set_label("% cells", fontsize=9)
    fig.suptitle("Sequential coexpression by study, sample, and cutoff", fontsize=14, fontweight="bold")

    for ext in ["png", "pdf", "svg"]:
        fig.savefig(out_base.with_suffix(f".{ext}"), dpi=450, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    args = parse_args()
    analysis_dir = Path(args.analysis_dir)
    table_dir = analysis_dir / "tables"
    plot_dir = analysis_dir / "plots"
    cutoffs = [float(x.strip()) for x in args.cutoffs.split(",") if x.strip()]

    all_sample_rows: list[dict[str, object]] = []
    all_overall_rows: list[dict[str, object]] = []
    for study in STUDY_SPECS:
        path = table_dir / f"{study['prefix']}_expression_and_coexpression_per_cell.tsv.gz"
        df = read_per_cell(path)
        for cutoff in cutoffs:
            sample_rows, overall_rows = summarize_one(df, study, cutoff)
            all_sample_rows.extend(sample_rows)
            all_overall_rows.extend(overall_rows)

    sample_df = pd.DataFrame(all_sample_rows)
    overall_df = pd.DataFrame(all_overall_rows)
    write_tsv(sample_df, table_dir / "combined_marker_coexpression_by_study_sample_cutoff.tsv")
    write_tsv(overall_df, table_dir / "combined_marker_coexpression_overall_by_study_cutoff.tsv")
    plot_sample_heatmap(sample_df, plot_dir / "combined_marker_coexpression_by_sample_cutoff_heatmap")


if __name__ == "__main__":
    main()
