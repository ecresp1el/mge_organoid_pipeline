#!/usr/bin/env python3
"""Plot winner-take-all MGE calls versus an MGE-score threshold."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


PROJECT_ROOT_DEFAULT = Path("/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder")
RUN_LABEL_DEFAULT = "cross_study_shi_seurat_label_transfer_v1"
RESULTS_DIRNAME = "cross_study_shi_seurat_label_transfer"

STUDIES = [
    ("varela_div30", "This study\nDIV30"),
    ("varela_div90", "This study\nDIV90"),
    ("siebert_2026", "Siebert\n2026"),
    ("walsh", "Walsh"),
    ("bershteyn_2025", "Bershteyn\n2025"),
    ("bershteyn_2023", "Bershteyn\n2023"),
    ("samarasinghe_2021", "Samarasinghe\n2021"),
]

WINNER_COLOR = "#2f6f73"
THRESHOLD_COLOR = "#d9893d"
LOST_COLOR = "#b8b8b8"


def output_paths(project_root: Path, run_label: str) -> tuple[Path, Path, Path]:
    run_dir = project_root / "results" / RESULTS_DIRNAME / run_label
    per_study_dir = run_dir / "tables" / "per_study"
    table_dir = run_dir / "tables"
    plot_dir = run_dir / "plots" / "summary"
    table_dir.mkdir(parents=True, exist_ok=True)
    plot_dir.mkdir(parents=True, exist_ok=True)
    return per_study_dir, table_dir, plot_dir


def load_tradeoff_summary(per_study_dir: Path, threshold: float) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for study_id, study_label in STUDIES:
        path = per_study_dir / f"{study_id}_shi_seurat_label_transfer_obs.tsv.gz"
        if not path.exists():
            raise FileNotFoundError(f"Missing per-study Shi obs table: {path}")
        obs = pd.read_csv(
            path,
            sep="\t",
            usecols=["shi_seurat_full_predicted_shi_label", "shi_seurat_full_prediction_score_MGE"],
        )
        total_cells = int(obs.shape[0])

        # Winner-take-all asks: "Did MGE have the highest support among all Shi
        # major labels?" It reproduces the default composition logic.
        top_mge = obs["shi_seurat_full_predicted_shi_label"].astype(str).eq("MGE")

        # Thresholded MGE asks a stricter question: "Is the MGE-specific support
        # at least this cutoff?" This is a QC/sensitivity subset, not the default
        # winner-take-all composition.
        mge_score = pd.to_numeric(obs["shi_seurat_full_prediction_score_MGE"], errors="coerce")
        threshold_mge = mge_score.ge(threshold)

        retained = top_mge & threshold_mge
        lost = top_mge & ~threshold_mge
        rows.append(
            {
                "study_id": study_id,
                "study_label": study_label.replace("\n", " "),
                "study_plot_label": study_label,
                "total_cells": total_cells,
                "winner_take_all_mge_cells": int(top_mge.sum()),
                "winner_take_all_mge_fraction": float(top_mge.mean()),
                f"mge_score_ge_{threshold:.2f}_cells": int(threshold_mge.sum()),
                f"mge_score_ge_{threshold:.2f}_fraction": float(threshold_mge.mean()),
                f"winner_mge_and_mge_score_ge_{threshold:.2f}_cells": int(retained.sum()),
                f"winner_mge_lost_by_{threshold:.2f}_cells": int(lost.sum()),
                f"winner_mge_lost_by_{threshold:.2f}_fraction_of_winner_mge": (
                    float(lost.sum() / top_mge.sum()) if int(top_mge.sum()) else np.nan
                ),
                f"winner_to_threshold_drop_{threshold:.2f}_percentage_points": float(
                    (top_mge.mean() - threshold_mge.mean()) * 100
                ),
            }
        )
    return pd.DataFrame(rows)


def format_axis_percent(ax: plt.Axes) -> None:
    ax.set_ylim(0, 1.02)
    ticks = np.linspace(0, 1, 6)
    ax.set_yticks(ticks)
    ax.set_yticklabels([f"{tick:.0%}" for tick in ticks])
    ax.grid(axis="y", color="#d8d8d8", linewidth=0.7)
    ax.set_axisbelow(True)


def save_figure(fig: plt.Figure, path: Path) -> None:
    fig.savefig(path.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(path.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def plot_fraction_comparison(summary: pd.DataFrame, threshold: float, plot_dir: Path) -> Path:
    labels = summary["study_plot_label"].tolist()
    x = np.arange(len(labels))
    width = 0.38
    threshold_col = f"mge_score_ge_{threshold:.2f}_fraction"

    fig, ax = plt.subplots(figsize=(10.5, 4.8))
    ax.bar(
        x - width / 2,
        summary["winner_take_all_mge_fraction"],
        width,
        color=WINNER_COLOR,
        label="Winner-take-all MGE",
    )
    ax.bar(
        x + width / 2,
        summary[threshold_col],
        width,
        color=THRESHOLD_COLOR,
        label=f"MGE score >= {threshold:.2f}",
    )
    format_axis_percent(ax)
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylabel("Fraction of all cells")
    ax.set_title(f"MGE calls: winner-take-all versus MGE score >= {threshold:.2f}")
    ax.legend(frameon=False, ncol=2, loc="upper right")
    fig.tight_layout()
    path = plot_dir / f"cross_study_shi_mge_winner_vs_score{int(threshold * 100):02d}_fraction_bar"
    save_figure(fig, path)
    return path.with_suffix(".png")


def plot_retention(summary: pd.DataFrame, threshold: float, plot_dir: Path) -> Path:
    labels = summary["study_plot_label"].tolist()
    x = np.arange(len(labels))
    retained_col = f"winner_mge_and_mge_score_ge_{threshold:.2f}_cells"
    lost_col = f"winner_mge_lost_by_{threshold:.2f}_cells"
    retained_fraction = summary[retained_col] / summary["total_cells"]
    lost_fraction = summary[lost_col] / summary["total_cells"]

    fig, ax = plt.subplots(figsize=(10.5, 4.8))
    ax.bar(x, retained_fraction, color=THRESHOLD_COLOR, label=f"Top MGE retained at >= {threshold:.2f}")
    ax.bar(x, lost_fraction, bottom=retained_fraction, color=LOST_COLOR, label="Top MGE below threshold")
    format_axis_percent(ax)
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylabel("Fraction of all cells")
    ax.set_title(f"What the {threshold:.2f} threshold removes from winner-take-all MGE")
    ax.legend(frameon=False, ncol=2, loc="upper right")
    fig.tight_layout()
    path = plot_dir / f"cross_study_shi_mge_score{int(threshold * 100):02d}_retention_stacked_bar"
    save_figure(fig, path)
    return path.with_suffix(".png")


def plot_drop_lollipop(summary: pd.DataFrame, threshold: float, plot_dir: Path) -> Path:
    labels = summary["study_plot_label"].tolist()
    x = np.arange(len(labels))
    threshold_col = f"mge_score_ge_{threshold:.2f}_fraction"

    winner = summary["winner_take_all_mge_fraction"].to_numpy()
    thresholded = summary[threshold_col].to_numpy()
    fig, ax = plt.subplots(figsize=(10.5, 4.8))
    for idx, low, high in zip(x, thresholded, winner, strict=True):
        ax.plot([idx, idx], [low, high], color="#777777", linewidth=2)
    ax.scatter(x, winner, color=WINNER_COLOR, s=58, label="Winner-take-all MGE", zorder=3)
    ax.scatter(x, thresholded, color=THRESHOLD_COLOR, s=58, label=f"MGE score >= {threshold:.2f}", zorder=3)
    for idx, low, high in zip(x, thresholded, winner, strict=True):
        ax.text(idx, high + 0.025, f"-{(high - low) * 100:.1f} pp", ha="center", va="bottom", fontsize=8)
    format_axis_percent(ax)
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylabel("Fraction of all cells")
    ax.set_title(f"Drop from winner-take-all MGE to MGE score >= {threshold:.2f}")
    ax.legend(frameon=False, ncol=2, loc="upper right")
    fig.tight_layout()
    path = plot_dir / f"cross_study_shi_mge_winner_to_score{int(threshold * 100):02d}_drop_lollipop"
    save_figure(fig, path)
    return path.with_suffix(".png")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, default=PROJECT_ROOT_DEFAULT)
    parser.add_argument("--run-label", default=RUN_LABEL_DEFAULT)
    parser.add_argument("--threshold", type=float, default=0.90)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    per_study_dir, table_dir, plot_dir = output_paths(args.project_root, args.run_label)
    summary = load_tradeoff_summary(per_study_dir, args.threshold)
    table_path = table_dir / f"cross_study_shi_mge_winner_vs_score{int(args.threshold * 100):02d}_summary.tsv"
    summary.drop(columns=["study_plot_label"]).to_csv(table_path, sep="\t", index=False)

    plot_paths = [
        plot_fraction_comparison(summary, args.threshold, plot_dir),
        plot_retention(summary, args.threshold, plot_dir),
        plot_drop_lollipop(summary, args.threshold, plot_dir),
    ]
    print(f"summary_table\t{table_path}", flush=True)
    for path in plot_paths:
        print(f"plot\t{path}", flush=True)


if __name__ == "__main__":
    main()
