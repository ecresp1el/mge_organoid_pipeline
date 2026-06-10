#!/usr/bin/env python3
"""Plot Shi label call modes and predicted-age sample composition.

Two call modes are used throughout this script:
  1. winner_take_all: the top transferred Shi label, with no score cutoff.
  2. score90_cutoff: an absolute support-score cutoff, defaulting to 0.90.

Three score scopes are summarized:
  1. all Shi major labels, one label at a time.
  2. MGE/LGE/CGE, using the summed support of those three labels.
  3. MGE only, using the MGE-specific support score.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import re


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

CONTROL_ONLY_STUDIES = {"samarasinghe_2021"}

SHI_LABELS = [
    "MGE",
    "LGE",
    "CGE",
    "progenitor",
    "Excitatory IPC",
    "Excitatory neuron",
    "Thalamic neurons",
    "Microglia",
    "OPC",
    "Endothelial",
]

LABEL_SETS = {
    "mge": ("mge", "MGE", ["MGE"]),
    "mge_lge_cge": ("mge_lge_cge", "MGE/LGE/CGE", ["MGE", "LGE", "CGE"]),
}

WINNER_COLOR = "#2f6f73"
THRESHOLD_COLOR = "#d9893d"
LOST_COLOR = "#b8b8b8"
GW_ORDER = ["GW09", "GW12", "GW13", "GW16", "GW18"]
GW_COLORS = {
    "GW09": "#4c6a9c",
    "GW12": "#56a0a6",
    "GW13": "#7dbd69",
    "GW16": "#e5b54a",
    "GW18": "#c45b45",
}


def label_token(label: str) -> str:
    return label.replace(" ", "_").replace("/", "_")


def score_col(label: str) -> str:
    return "shi_seurat_full_prediction_score_" + label_token(label)


def output_paths(project_root: Path, run_label: str) -> tuple[Path, Path, Path, Path]:
    run_dir = project_root / "results" / RESULTS_DIRNAME / run_label
    per_study_dir = run_dir / "tables" / "per_study"
    table_dir = run_dir / "tables"
    plot_dir = run_dir / "plots" / "summary"
    age_plot_dir = plot_dir / "predicted_age_sample_composition"
    table_dir.mkdir(parents=True, exist_ok=True)
    plot_dir.mkdir(parents=True, exist_ok=True)
    age_plot_dir.mkdir(parents=True, exist_ok=True)
    return per_study_dir, table_dir, plot_dir, age_plot_dir


def filter_analysis_obs(study_id: str, obs: pd.DataFrame) -> pd.DataFrame:
    """Apply study-specific analysis filters before calculating fractions."""
    if study_id not in CONTROL_ONLY_STUDIES:
        return obs
    sample_text = (
        obs.get("sample", pd.Series("", index=obs.index)).fillna("").astype(str)
        + " "
        + obs.get("sample_label", pd.Series("", index=obs.index)).fillna("").astype(str)
    )
    keep = sample_text.str.contains("Ctrl", case=False, na=False)
    if not bool(keep.any()):
        raise ValueError(f"{study_id}: control-only filter removed all cells.")
    return obs.loc[keep].copy()


def read_study_obs(per_study_dir: Path, study_id: str, labels: list[str]) -> pd.DataFrame:
    path = per_study_dir / f"{study_id}_shi_seurat_label_transfer_obs.tsv.gz"
    if not path.exists():
        raise FileNotFoundError(f"Missing per-study Shi obs table: {path}")
    cols = [
        "sample",
        "sample_label",
        "shi_seurat_full_predicted_shi_label",
        *[score_col(label) for label in labels],
    ]
    obs = pd.read_csv(path, sep="\t", usecols=cols)
    return filter_analysis_obs(study_id, obs)


def canonical_gw_label(value: object) -> str:
    """Collapse labels such as GW12_01 and GW12_02 to the plotted GW12 age bin."""
    match = re.search(r"GW\s*([0-9]+)", str(value).upper())
    if not match:
        return "unknown"
    return f"GW{int(match.group(1)):02d}"


def load_group_tradeoff_summary(
    per_study_dir: Path,
    threshold: float,
    group_slug: str,
    group_label: str,
    group_labels: list[str],
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for study_id, study_plot_label in STUDIES:
        obs = read_study_obs(per_study_dir, study_id, group_labels)
        total_cells = int(obs.shape[0])

        # Winner-take-all group membership asks whether the winning transferred
        # Shi major label is inside this group. This is the default composition
        # logic and does not impose any absolute score threshold.
        winner_group = obs["shi_seurat_full_predicted_shi_label"].astype(str).isin(group_labels)

        # Thresholded group membership asks whether the support assigned to this
        # group reaches the cutoff. For a multi-label group such as MGE/LGE/CGE,
        # the group support is the sum of the Seurat per-label supports for the
        # labels in that group; for MGE alone, this is simply the MGE score.
        group_score = sum(
            pd.to_numeric(obs[score_col(label)], errors="coerce").fillna(0.0)
            for label in group_labels
        )
        threshold_group = group_score.ge(threshold)

        retained = winner_group & threshold_group
        lost = winner_group & ~threshold_group
        rows.append(
            {
                "study_id": study_id,
                "study_label": study_plot_label.replace("\n", " "),
                "study_plot_label": study_plot_label,
                "label_set": group_label,
                "label_set_slug": group_slug,
                "labels_in_set": ",".join(group_labels),
                "call_mode_winner": "winner_take_all",
                "call_mode_threshold": f"score{int(threshold * 100):02d}_cutoff",
                "score_scope": group_label,
                "analysis_filter": "control_only" if study_id in CONTROL_ONLY_STUDIES else "all_cells",
                "total_cells": total_cells,
                "winner_take_all_cells": int(winner_group.sum()),
                "winner_take_all_fraction": float(winner_group.mean()),
                f"score_ge_{threshold:.2f}_cells": int(threshold_group.sum()),
                f"score_ge_{threshold:.2f}_fraction": float(threshold_group.mean()),
                f"winner_and_score_ge_{threshold:.2f}_cells": int(retained.sum()),
                f"winner_lost_by_{threshold:.2f}_cells": int(lost.sum()),
                f"winner_lost_by_{threshold:.2f}_fraction_of_winner": (
                    float(lost.sum() / winner_group.sum()) if int(winner_group.sum()) else np.nan
                ),
                f"winner_to_threshold_drop_{threshold:.2f}_percentage_points": float(
                    (winner_group.mean() - threshold_group.mean()) * 100
                ),
            }
        )
    return pd.DataFrame(rows)


def load_all_label_summary(per_study_dir: Path, threshold: float) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for study_id, study_plot_label in STUDIES:
        obs = read_study_obs(per_study_dir, study_id, SHI_LABELS)
        total_cells = int(obs.shape[0])
        predicted = obs["shi_seurat_full_predicted_shi_label"].astype(str)
        for label in SHI_LABELS:
            winner_label = predicted.eq(label)
            label_score = pd.to_numeric(obs[score_col(label)], errors="coerce")
            threshold_label = label_score.ge(threshold)
            rows.append(
                {
                    "study_id": study_id,
                    "study_label": study_plot_label.replace("\n", " "),
                    "study_plot_label": study_plot_label,
                    "shi_label": label,
                    "call_mode_winner": "winner_take_all",
                    "call_mode_threshold": f"score{int(threshold * 100):02d}_cutoff",
                    "score_scope": "all_shi_major_labels",
                    "analysis_filter": "control_only" if study_id in CONTROL_ONLY_STUDIES else "all_cells",
                    "total_cells": total_cells,
                    "winner_take_all_cells": int(winner_label.sum()),
                    "winner_take_all_fraction": float(winner_label.mean()),
                    f"score_ge_{threshold:.2f}_cells": int(threshold_label.sum()),
                    f"score_ge_{threshold:.2f}_fraction": float(threshold_label.mean()),
                    f"winner_lost_by_{threshold:.2f}_cells": int((winner_label & ~threshold_label).sum()),
                    f"winner_to_threshold_drop_{threshold:.2f}_percentage_points": float(
                        (winner_label.mean() - threshold_label.mean()) * 100
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
    group_label = str(summary["label_set"].iloc[0])
    group_slug = str(summary["label_set_slug"].iloc[0])
    x = np.arange(len(labels))
    width = 0.38
    threshold_col = f"score_ge_{threshold:.2f}_fraction"

    fig, ax = plt.subplots(figsize=(10.5, 4.8))
    ax.bar(
        x - width / 2,
        summary["winner_take_all_fraction"],
        width,
        color=WINNER_COLOR,
        label=f"Winner-take-all {group_label}",
    )
    ax.bar(
        x + width / 2,
        summary[threshold_col],
        width,
        color=THRESHOLD_COLOR,
        label=f"{group_label} score >= {threshold:.2f}",
    )
    format_axis_percent(ax)
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylabel("Fraction of all cells")
    ax.set_title(f"{group_label}: winner-take-all versus score >= {threshold:.2f}")
    ax.legend(frameon=False, ncol=2, loc="upper right")
    fig.tight_layout()
    path = plot_dir / f"cross_study_shi_{group_slug}_winner_vs_score{int(threshold * 100):02d}_fraction_bar"
    save_figure(fig, path)
    return path.with_suffix(".png")


def plot_retention(summary: pd.DataFrame, threshold: float, plot_dir: Path) -> Path:
    labels = summary["study_plot_label"].tolist()
    group_label = str(summary["label_set"].iloc[0])
    group_slug = str(summary["label_set_slug"].iloc[0])
    x = np.arange(len(labels))
    retained_col = f"winner_and_score_ge_{threshold:.2f}_cells"
    lost_col = f"winner_lost_by_{threshold:.2f}_cells"
    retained_fraction = summary[retained_col] / summary["total_cells"]
    lost_fraction = summary[lost_col] / summary["total_cells"]

    fig, ax = plt.subplots(figsize=(10.5, 4.8))
    ax.bar(x, retained_fraction, color=THRESHOLD_COLOR, label=f"Top {group_label} retained at >= {threshold:.2f}")
    ax.bar(x, lost_fraction, bottom=retained_fraction, color=LOST_COLOR, label=f"Top {group_label} below threshold")
    format_axis_percent(ax)
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylabel("Fraction of all cells")
    ax.set_title(f"What the {threshold:.2f} threshold removes from winner-take-all {group_label}")
    ax.legend(frameon=False, ncol=2, loc="upper right")
    fig.tight_layout()
    path = plot_dir / f"cross_study_shi_{group_slug}_score{int(threshold * 100):02d}_retention_stacked_bar"
    save_figure(fig, path)
    return path.with_suffix(".png")


def plot_drop_lollipop(summary: pd.DataFrame, threshold: float, plot_dir: Path) -> Path:
    labels = summary["study_plot_label"].tolist()
    group_label = str(summary["label_set"].iloc[0])
    group_slug = str(summary["label_set_slug"].iloc[0])
    x = np.arange(len(labels))
    threshold_col = f"score_ge_{threshold:.2f}_fraction"

    winner = summary["winner_take_all_fraction"].to_numpy()
    thresholded = summary[threshold_col].to_numpy()
    fig, ax = plt.subplots(figsize=(10.5, 4.8))
    for idx, low, high in zip(x, thresholded, winner, strict=True):
        ax.plot([idx, idx], [low, high], color="#777777", linewidth=2)
    ax.scatter(x, winner, color=WINNER_COLOR, s=58, label=f"Winner-take-all {group_label}", zorder=3)
    ax.scatter(x, thresholded, color=THRESHOLD_COLOR, s=58, label=f"{group_label} score >= {threshold:.2f}", zorder=3)
    for idx, low, high in zip(x, thresholded, winner, strict=True):
        ax.text(idx, min(1.0, high + 0.025), f"-{(high - low) * 100:.1f} pp", ha="center", va="bottom", fontsize=8)
    format_axis_percent(ax)
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylabel("Fraction of all cells")
    ax.set_title(f"Drop from winner-take-all {group_label} to score >= {threshold:.2f}")
    ax.legend(frameon=False, ncol=2, loc="upper right")
    fig.tight_layout()
    path = plot_dir / f"cross_study_shi_{group_slug}_winner_to_score{int(threshold * 100):02d}_drop_lollipop"
    save_figure(fig, path)
    return path.with_suffix(".png")


def plot_all_label_heatmap(
    all_summary: pd.DataFrame,
    value_col: str,
    title: str,
    path: Path,
    vmax: float,
    is_percentage_points: bool = False,
) -> Path:
    matrix = (
        all_summary.pivot(index="shi_label", columns="study_plot_label", values=value_col)
        .reindex(index=SHI_LABELS, columns=[label for _, label in STUDIES])
        .to_numpy(dtype=float)
    )
    fig, ax = plt.subplots(figsize=(10.8, 5.8))
    image = ax.imshow(matrix, aspect="auto", cmap="YlGnBu", vmin=0, vmax=vmax)
    ax.set_xticks(np.arange(len(STUDIES)))
    ax.set_xticklabels([label for _, label in STUDIES], rotation=0)
    ax.set_yticks(np.arange(len(SHI_LABELS)))
    ax.set_yticklabels(SHI_LABELS)
    for y in range(matrix.shape[0]):
        for x in range(matrix.shape[1]):
            value = matrix[y, x]
            text = f"{value:.1f}" if is_percentage_points else f"{value * 100:.1f}"
            ax.text(x, y, text, ha="center", va="center", fontsize=7, color="#1f1f1f")
    cbar = fig.colorbar(image, ax=ax, fraction=0.025, pad=0.02)
    cbar.set_label("Percentage points" if is_percentage_points else "Percent of cells")
    ax.set_title(title)
    fig.tight_layout()
    save_figure(fig, path)
    return path.with_suffix(".png")


def load_predicted_age_sample_composition(per_study_dir: Path) -> pd.DataFrame:
    rows: list[pd.DataFrame] = []
    for study_id, study_plot_label in STUDIES:
        path = per_study_dir / f"{study_id}_shi_seurat_label_transfer_obs.tsv.gz"
        obs = pd.read_csv(
            path,
            sep="\t",
            usecols=["sample", "sample_label", "shi_seurat_full_predicted_shi_week_label"],
        )
        obs = filter_analysis_obs(study_id, obs)
        obs["study_id"] = study_id
        obs["study_label"] = study_plot_label.replace("\n", " ")
        obs["study_plot_label"] = study_plot_label
        obs["sample"] = obs["sample"].astype(str)
        obs["sample_label"] = obs["sample_label"].fillna(obs["sample"]).astype(str)
        obs["predicted_age"] = obs["shi_seurat_full_predicted_shi_week_label"].map(canonical_gw_label)

        # This composition is also winner-take-all: each cell contributes to the
        # single Shi gestational-age label that won the week-level TransferData
        # call. It is not thresholded by age score.
        counts = (
            obs.groupby(["study_id", "study_label", "study_plot_label", "sample", "sample_label", "predicted_age"], observed=True)
            .size()
            .rename("n_cells")
            .reset_index()
        )
        sample_totals = counts.groupby(["study_id", "sample"], observed=True)["n_cells"].transform("sum")
        counts["fraction_of_sample"] = counts["n_cells"] / sample_totals
        counts["call_mode"] = "winner_take_all_predicted_age"
        counts["analysis_filter"] = "control_only" if study_id in CONTROL_ONLY_STUDIES else "all_cells"
        rows.append(counts)
    return pd.concat(rows, ignore_index=True)


def plot_age_composition_one_study(age_summary: pd.DataFrame, study_id: str, age_plot_dir: Path) -> Path:
    data = age_summary.loc[age_summary["study_id"] == study_id].copy()
    if data.empty:
        raise ValueError(f"No predicted-age composition rows for {study_id}")
    study_label = str(data["study_label"].iloc[0])
    sample_order = (
        data.groupby(["sample", "sample_label"], observed=True)["n_cells"]
        .sum()
        .sort_values(ascending=False)
        .reset_index()["sample_label"]
        .tolist()
    )
    pivot = (
        data.pivot_table(
            index="sample_label",
            columns="predicted_age",
            values="fraction_of_sample",
            aggfunc="sum",
            fill_value=0.0,
            observed=True,
        )
        .reindex(index=sample_order)
        .reindex(columns=GW_ORDER, fill_value=0.0)
    )

    height = max(4.2, 0.34 * len(pivot.index) + 1.8)
    fig, ax = plt.subplots(figsize=(9.2, height))
    left = np.zeros(len(pivot.index), dtype=float)
    y = np.arange(len(pivot.index))
    for gw in GW_ORDER:
        values = pivot[gw].to_numpy(dtype=float)
        ax.barh(y, values, left=left, color=GW_COLORS[gw], label=gw, height=0.78)
        left += values
    ax.set_yticks(y)
    ax.set_yticklabels(pivot.index, fontsize=8)
    ax.invert_yaxis()
    ax.set_xlim(0, 1.0)
    ax.set_xticks(np.linspace(0, 1, 6))
    ax.set_xticklabels([f"{tick:.0%}" for tick in np.linspace(0, 1, 6)])
    ax.set_xlabel("Percent of sample cells")
    ax.set_title(f"{study_label}: winner-take-all predicted Shi age by sample")
    ax.legend(title="Predicted age", frameon=False, ncol=min(len(GW_ORDER), 5), loc="lower center", bbox_to_anchor=(0.5, 1.02))
    ax.grid(axis="x", color="#d8d8d8", linewidth=0.7)
    ax.set_axisbelow(True)
    fig.tight_layout()
    path = age_plot_dir / f"{study_id}_shi_predicted_age_sample_composition_stacked_bar"
    save_figure(fig, path)
    return path.with_suffix(".png")


def write_predicted_age_sample_composition_outputs(
    per_study_dir: Path,
    table_dir: Path,
    age_plot_dir: Path,
) -> list[Path]:
    age_summary = load_predicted_age_sample_composition(per_study_dir)
    table_path = table_dir / "cross_study_shi_predicted_age_sample_composition.tsv"
    age_summary.drop(columns=["study_plot_label"]).to_csv(table_path, sep="\t", index=False)
    outputs = [table_path]
    for study_id, _ in STUDIES:
        outputs.append(plot_age_composition_one_study(age_summary, study_id, age_plot_dir))
    return outputs


def write_group_outputs(
    per_study_dir: Path,
    table_dir: Path,
    plot_dir: Path,
    threshold: float,
    group_slug: str,
    group_label: str,
    group_labels: list[str],
) -> list[Path]:
    summary = load_group_tradeoff_summary(per_study_dir, threshold, group_slug, group_label, group_labels)
    table_path = table_dir / f"cross_study_shi_{group_slug}_winner_vs_score{int(threshold * 100):02d}_summary.tsv"
    summary.drop(columns=["study_plot_label"]).to_csv(table_path, sep="\t", index=False)
    return [
        table_path,
        plot_fraction_comparison(summary, threshold, plot_dir),
        plot_retention(summary, threshold, plot_dir),
        plot_drop_lollipop(summary, threshold, plot_dir),
    ]


def write_all_label_outputs(
    per_study_dir: Path,
    table_dir: Path,
    plot_dir: Path,
    threshold: float,
) -> list[Path]:
    summary = load_all_label_summary(per_study_dir, threshold)
    table_path = table_dir / f"cross_study_shi_all_labels_winner_vs_score{int(threshold * 100):02d}_summary.tsv"
    summary.drop(columns=["study_plot_label"]).to_csv(table_path, sep="\t", index=False)
    threshold_col = f"score_ge_{threshold:.2f}_fraction"
    drop_col = f"winner_to_threshold_drop_{threshold:.2f}_percentage_points"
    return [
        table_path,
        plot_all_label_heatmap(
            summary,
            "winner_take_all_fraction",
            "Winner-take-all Shi major-label fractions",
            plot_dir / "cross_study_shi_all_labels_winner_take_all_fraction_heatmap",
            vmax=1.0,
        ),
        plot_all_label_heatmap(
            summary,
            threshold_col,
            f"Shi major-label score >= {threshold:.2f} fractions",
            plot_dir / f"cross_study_shi_all_labels_score{int(threshold * 100):02d}_fraction_heatmap",
            vmax=1.0,
        ),
        plot_all_label_heatmap(
            summary,
            drop_col,
            f"Drop from winner-take-all to score >= {threshold:.2f}",
            plot_dir / f"cross_study_shi_all_labels_winner_to_score{int(threshold * 100):02d}_drop_heatmap",
            vmax=100.0,
            is_percentage_points=True,
        ),
    ]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, default=PROJECT_ROOT_DEFAULT)
    parser.add_argument("--run-label", default=RUN_LABEL_DEFAULT)
    parser.add_argument("--threshold", type=float, default=0.90)
    parser.add_argument(
        "--which",
        choices=["all", "mge", "mge-lge-cge", "all-labels", "predicted-age"],
        default="all",
        help="Which tradeoff plot set to generate.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    per_study_dir, table_dir, plot_dir, age_plot_dir = output_paths(args.project_root, args.run_label)
    outputs: list[Path] = []
    if args.which in {"all", "mge"}:
        outputs.extend(write_group_outputs(per_study_dir, table_dir, plot_dir, args.threshold, *LABEL_SETS["mge"]))
    if args.which in {"all", "mge-lge-cge"}:
        outputs.extend(write_group_outputs(per_study_dir, table_dir, plot_dir, args.threshold, *LABEL_SETS["mge_lge_cge"]))
    if args.which in {"all", "all-labels"}:
        outputs.extend(write_all_label_outputs(per_study_dir, table_dir, plot_dir, args.threshold))
    if args.which in {"all", "predicted-age"}:
        outputs.extend(write_predicted_age_sample_composition_outputs(per_study_dir, table_dir, age_plot_dir))

    for path in outputs:
        kind = "summary_table" if path.suffix == ".tsv" else "plot"
        print(f"{kind}\t{path}", flush=True)


if __name__ == "__main__":
    main()
