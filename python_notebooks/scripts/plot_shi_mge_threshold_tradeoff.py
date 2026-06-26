#!/usr/bin/env python3
"""Plot Shi label-transfer score tradeoffs and predicted-age composition.

Upstream score source:
  Seurat TransferData is run twice in the R pipeline, once for Shi major labels
  and once for Shi gestational-week labels. For each cell, TransferData exports
  one score per reference class plus prediction.score.max. The R/Python schema
  checks require prediction.score.max to equal the row-wise maximum of the
  per-class scores, and uncertainty is defined as 1 - prediction.score.max.

Major-label tradeoff summaries:
  1. winner_take_all: cell membership is the top transferred Shi major label.
  2. score90_cutoff: cell membership is support score >= the cutoff, default
     0.90. For MGE this is the MGE score; for MGE/LGE/CGE this is the sum of
     those three major-label support scores; for all-label heatmaps this is each
     Shi major label evaluated separately.

Predicted-age sample-composition plots:
  These are winner-take-all week-label plots, not score90-cutoff plots. All-Shi
  age composition uses the whole-Shi week classifier. MGE and MGE/LGE/CGE age
  composition use the GE-only week classifier, which is produced by rerunning
  week TransferData against only Shi reference MGE/LGE/CGE cells.

Inputs:
  Reads per-study obs tables from
  results/cross_study_shi_seurat_label_transfer/<run_label>/tables/per_study/.

Outputs:
  Writes winner/threshold plots under
  plots/summary/threshold_tradeoff/ and predicted-age sample composition under
  plots/summary/predicted_age_sample_composition/{all_shi_major_labels,mge,mge_lge_cge}/.

Slurm entry point:
  slurm_templates/30_shi_threshold_tradeoff_plots.sbatch.template
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
RUN_LABEL_DEFAULT = "cross_study_shi_seurat_label_transfer_v2_ge_only_age"
RESULTS_DIRNAME = "cross_study_shi_seurat_label_transfer"

STUDIES = [
    ("varela_div30", "Varela et al.\nthis paper\nDIV30"),
    ("varela_div90", "Varela et al.\nthis paper\nDIV90"),
    ("siebert_2026", "Siebert\n2026"),
    ("walsh", "Walsh\n2025"),
    ("bershteyn_2023", "Bershteyn\n2023"),
    ("bershteyn_2025", "Bershteyn\n2025"),
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

AGE_CELL_SETS = [
    # These age-composition sets are denominator filters based on the winning
    # Shi major label. They do not apply the 0.90 support-score cutoff. The GE
    # subsets use the GE-only week classifier rather than whole-Shi week scores.
    (
        "all_shi_major_labels",
        "All Shi major-label cells",
        SHI_LABELS,
        "whole_shi_week_classifier",
        "Whole-Shi week classifier",
        "shi_seurat_full_predicted_shi_week_label",
    ),
    (
        "mge",
        "MGE winner-take-all cells",
        ["MGE"],
        "ge_only_week_classifier",
        "GE-only week classifier",
        "shi_seurat_ge_only_predicted_shi_week_label",
    ),
    (
        "mge_lge_cge",
        "MGE/LGE/CGE winner-take-all cells",
        ["MGE", "LGE", "CGE"],
        "ge_only_week_classifier",
        "GE-only week classifier",
        "shi_seurat_ge_only_predicted_shi_week_label",
    ),
]

WINNER_COLOR = "#2f6f73"
THRESHOLD_COLOR = "#d9893d"
LOST_COLOR = "#b8b8b8"
GW_ORDER = ["GW09", "GW12", "GW13", "GW16", "GW18"]
GW_COLORS = {
    "GW09": "#231611",
    "GW12": "#3F1C6A",
    "GW13": "#A02E6B",
    "GW16": "#EB5840",
    "GW18": "#FCC031",
}

BERSHTEYN_2023_SAMPLE_INFO = {
    # GSE208672 plus study design notes: D0 hESCs, D14 NPCs,
    # and week-6 End-of-Process cells (6 * 7 = DIV42). S/U denotes sorted
    # or unsorted, respectively.
    "D0": ("DIV0 hESC", 0.0, 0.0),
    "D14": ("DIV14 NPC", 14.0, 0.0),
    "MB460": ("DIV42 EOP L1 U", 42.0, 1.1),
    "MB461": ("DIV42 EOP L1 S", 42.0, 1.2),
    "MB279": ("DIV42 EOP L2 U", 42.0, 2.1),
    "MB280": ("DIV42 EOP L2 S", 42.0, 2.2),
    "MB527": ("DIV42 EOP L3 U", 42.0, 3.1),
    "MB528": ("DIV42 EOP L3 S", 42.0, 3.2),
    "MS35mock": ("DIV42 EOP B4 U (MS35mock)", 42.0, 4.0),
    "MS35r41": ("DIV42 EOP B5 U (MS35r41)", 42.0, 5.0),
    "mockv2dw": ("DIV42 EOP B6 U (mockv2dw)", 42.0, 6.0),
    "mockv2ym": ("DIV42 EOP B7 U (mockv2ym)", 42.0, 7.0),
    "r41v2dw": ("DIV42 EOP B8 U (r41v2dw)", 42.0, 8.0),
    "r41v2ym": ("DIV42 EOP B9 U (r41v2ym)", 42.0, 9.0),
}

BERSHTEYN_2025_SAMPLE_INFO = {
    # GSE283775 series matrix: all samples are sorted End-of-Process batches.
    "200319S": ("DIV unknown EOP batch 1 (200319S)", 999.0, 1.0),
    "010519S1": ("DIV unknown EOP batch 2 (010519S1)", 999.0, 2.0),
    "010519S2": ("DIV unknown EOP batch 3 (010519S2)", 999.0, 3.0),
    "111219S": ("DIV unknown EOP batch 4 (111219S)", 999.0, 4.0),
    "251219S": ("DIV unknown EOP batch 5 (251219S)", 999.0, 5.0),
    "150120S": ("DIV unknown EOP batch 6 (150120S)", 999.0, 6.0),
    "200520S1": ("DIV unknown EOP batch 7 (200520S1)", 999.0, 7.0),
    "200520S2": ("DIV unknown EOP batch 8 (200520S2)", 999.0, 8.0),
    "100620S": ("DIV unknown EOP batch 9 (100620S)", 999.0, 9.0),
    "010720S": ("DIV unknown EOP batch 10 (010720S)", 999.0, 10.0),
    "220720S1": ("DIV unknown EOP batch 11 (220720S1)", 999.0, 11.0),
    "220720S2": ("DIV unknown EOP batch 12 (220720S2)", 999.0, 12.0),
    "070120S": ("DIV unknown EOP batch 13 (070120S)", 999.0, 13.0),
    "280120S": ("DIV unknown EOP batch 14 (280120S)", 999.0, 14.0),
}

STUDY_ORDER = {study_id: idx for idx, (study_id, _) in enumerate(STUDIES)}


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


def first_number(pattern: str, text: str, default: float = 9999.0) -> float:
    match = re.search(pattern, text, flags=re.IGNORECASE)
    if not match:
        return default
    return float(match.group(1))


def sample_display_and_order(study_id: str, sample: object, sample_label: object) -> tuple[str, float, float]:
    sample_text = str(sample)
    label_text = str(sample_label) if pd.notna(sample_label) else sample_text
    display = label_text if label_text and label_text.lower() != "nan" else sample_text

    if study_id == "bershteyn_2023" and sample_text in BERSHTEYN_2023_SAMPLE_INFO:
        return BERSHTEYN_2023_SAMPLE_INFO[sample_text]
    if study_id == "bershteyn_2025" and sample_text in BERSHTEYN_2025_SAMPLE_INFO:
        return BERSHTEYN_2025_SAMPLE_INFO[sample_text]

    if study_id == "varela_div30":
        return f"DIV30 {sample_text}", 30.0, first_number(r"MW-(\d+)", sample_text)
    if study_id == "varela_div90":
        return f"DIV90 {sample_text}", 90.0, first_number(r"MW-(\d+)", sample_text)

    if study_id == "siebert_2026":
        age_order = 51.0 if sample_text.startswith("Young") else 164.0 if sample_text.startswith("Old") else 9999.0
        div_label = f"likely DIV{int(age_order)}" if age_order < 9999 else "DIV unknown"
        return f"{div_label} {sample_text}", age_order, first_number(r"_(\d+)", sample_text)

    if study_id == "walsh":
        within_order = 1.0 if "dFB" in label_text else 2.0 if "vFB" in label_text else first_number(r"GSM(\d+)", sample_text)
        div = first_number(r"d(\d+)", label_text, 75.0)
        return f"DIV{int(div)} {label_text.replace('_', ' ')}", div, within_order

    if study_id == "samarasinghe_2021":
        div = first_number(r"D(\d+)", display)
        return f"DIV{int(div)} {display}" if div < 9999 else f"DIV unknown {display}", div, first_number(r"_(\d+)_seurat", display)

    return display, first_number(r"D(?:IV)?\s*([0-9]+)", display), first_number(r"(\d+)", display)


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

        # Winner-take-all group membership uses only the transferred major-label
        # class with the largest Seurat support score. No absolute score cutoff
        # is imposed in this mode.
        winner_group = obs["shi_seurat_full_predicted_shi_label"].astype(str).isin(group_labels)

        # Thresholded group membership ignores which label won and asks whether
        # the Seurat support assigned to this group reaches the cutoff. For a
        # multi-label group such as MGE/LGE/CGE, group support is the sum of the
        # per-label support scores in that group; for MGE alone, it is the MGE
        # support score. This cutoff is on major-label support, not GW support.
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
            # In the all-label view, the 0.90 threshold is applied separately to
            # each Shi major label's own Seurat support score.
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


def load_predicted_age_sample_composition(
    per_study_dir: Path,
    cell_set_slug: str,
    cell_set_label: str,
    major_labels: list[str],
    age_classifier_slug: str,
    age_classifier_label: str,
    predicted_week_col: str,
) -> pd.DataFrame:
    rows: list[pd.DataFrame] = []
    for study_id, study_plot_label in STUDIES:
        path = per_study_dir / f"{study_id}_shi_seurat_label_transfer_obs.tsv.gz"
        obs = pd.read_csv(
            path,
            sep="\t",
            usecols=[
                "sample",
                "sample_label",
                "shi_seurat_full_predicted_shi_label",
                predicted_week_col,
            ],
        )
        obs = filter_analysis_obs(study_id, obs)
        # Age-composition subsets use the winner-take-all major-label call as
        # the denominator filter. They intentionally do not require the major
        # label support to be >= 0.90, and they do not threshold week scores.
        if set(major_labels) != set(SHI_LABELS):
            obs = obs.loc[obs["shi_seurat_full_predicted_shi_label"].astype(str).isin(major_labels)].copy()
        if obs.empty:
            continue
        obs["study_id"] = study_id
        obs["study_label"] = study_plot_label.replace("\n", " ")
        obs["study_plot_label"] = study_plot_label
        obs["age_cell_set_slug"] = cell_set_slug
        obs["age_cell_set_label"] = cell_set_label
        obs["age_classifier_slug"] = age_classifier_slug
        obs["age_classifier_label"] = age_classifier_label
        obs["sample"] = obs["sample"].astype(str)
        obs["sample_label"] = obs["sample_label"].fillna(obs["sample"]).astype(str)
        obs["predicted_age"] = obs[predicted_week_col].map(canonical_gw_label)

        sample_meta = obs[["sample", "sample_label"]].drop_duplicates().copy()
        sample_records = [
            sample_display_and_order(study_id, row.sample, row.sample_label)
            for row in sample_meta.itertuples(index=False)
        ]
        sample_meta["sample_display_label"] = [record[0] for record in sample_records]
        sample_meta["sample_age_order"] = [record[1] for record in sample_records]
        sample_meta["sample_order_within_study"] = [record[2] for record in sample_records]
        sample_meta["sample_plot_id"] = study_id + "::" + sample_meta["sample"]
        obs = obs.merge(sample_meta, on=["sample", "sample_label"], how="left", validate="many_to_one")

        # After the denominator filter above, each remaining cell contributes to
        # the single Shi gestational-age label that won the week-level
        # TransferData call. These stacked bars are not 0.90-cutoff age plots.
        counts = (
            obs.groupby(
                [
                    "study_id",
                    "study_label",
                    "study_plot_label",
                    "age_cell_set_slug",
                    "age_cell_set_label",
                    "age_classifier_slug",
                    "age_classifier_label",
                    "sample",
                    "sample_label",
                    "sample_display_label",
                    "sample_age_order",
                    "sample_order_within_study",
                    "sample_plot_id",
                    "predicted_age",
                ],
                observed=True,
            )
            .size()
            .rename("n_cells")
            .reset_index()
        )
        sample_totals = counts.groupby(["study_id", "sample"], observed=True)["n_cells"].transform("sum")
        counts["fraction_of_sample"] = counts["n_cells"] / sample_totals
        counts["call_mode"] = "winner_take_all_predicted_age"
        counts["predicted_week_source_column"] = predicted_week_col
        counts["major_label_filter"] = ",".join(major_labels)
        counts["analysis_filter"] = "control_only" if study_id in CONTROL_ONLY_STUDIES else "all_cells"
        rows.append(counts)
    if not rows:
        raise ValueError(f"No predicted-age composition rows for cell set {cell_set_slug}")
    return pd.concat(rows, ignore_index=True)


def plot_age_composition_one_study(age_summary: pd.DataFrame, study_id: str, age_plot_dir: Path) -> Path:
    data = age_summary.loc[age_summary["study_id"] == study_id].copy()
    if data.empty:
        raise ValueError(f"No predicted-age composition rows for {study_id}")
    cell_set_slug = str(data["age_cell_set_slug"].iloc[0])
    cell_set_label = str(data["age_cell_set_label"].iloc[0])
    age_classifier_label = str(data["age_classifier_label"].iloc[0])
    order = (
        data.groupby(
            ["sample_plot_id", "sample_display_label", "sample_age_order", "sample_order_within_study"],
            observed=True,
            as_index=False,
        )["n_cells"]
        .sum()
        .sort_values(["sample_age_order", "sample_order_within_study", "sample_display_label"], kind="stable")
    )
    pivot = (
        data.pivot_table(
            index="sample_plot_id",
            columns="predicted_age",
            values="fraction_of_sample",
            aggfunc="sum",
            fill_value=0.0,
            observed=True,
        )
        .reindex(index=order["sample_plot_id"].tolist())
        .reindex(columns=GW_ORDER, fill_value=0.0)
    )
    labels = order["sample_display_label"].tolist()

    width = max(7.5, 0.42 * len(pivot.index) + 2.6)
    fig, ax = plt.subplots(figsize=(width, 5.2))
    bottom = np.zeros(len(pivot.index), dtype=float)
    x = np.arange(len(pivot.index))
    for gw in GW_ORDER:
        values = pivot[gw].to_numpy(dtype=float) * 100.0
        ax.bar(x, values, bottom=bottom, color=GW_COLORS[gw], label=gw, width=0.82)
        bottom += values
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=90, ha="center", va="top", fontsize=8)
    ax.set_ylim(0, 100)
    ax.set_yticks(np.arange(0, 101, 20))
    ax.set_ylabel("% by predicted cell stage")
    ax.set_title(f"Sample composition\nby predicted cell stage\n{cell_set_label}\n{age_classifier_label}")
    ax.legend(title="Predicted stage", frameon=False, ncol=min(len(GW_ORDER), 5), loc="lower center", bbox_to_anchor=(0.5, 1.02))
    ax.grid(axis="y", color="#d8d8d8", linewidth=0.7)
    ax.set_axisbelow(True)
    fig.tight_layout()
    path = age_plot_dir / f"{study_id}_{cell_set_slug}_shi_predicted_age_sample_composition_stacked_bar"
    save_figure(fig, path)
    return path.with_suffix(".png")


def plot_age_composition_all_studies(age_summary: pd.DataFrame, age_plot_dir: Path) -> Path:
    data = age_summary.copy()
    cell_set_slug = str(data["age_cell_set_slug"].iloc[0])
    cell_set_label = str(data["age_cell_set_label"].iloc[0])
    age_classifier_label = str(data["age_classifier_label"].iloc[0])
    order = (
        data.groupby(
            [
                "study_id",
                "study_label",
                "study_plot_label",
                "sample_plot_id",
                "sample_display_label",
                "sample_age_order",
                "sample_order_within_study",
            ],
            observed=True,
            as_index=False,
        )["n_cells"]
        .sum()
    )
    order["study_plot_order"] = order["study_id"].map(STUDY_ORDER)
    order = order.sort_values(
        ["study_plot_order", "sample_age_order", "sample_order_within_study", "sample_display_label"],
        kind="stable",
    ).reset_index(drop=True)
    pivot = (
        data.pivot_table(
            index="sample_plot_id",
            columns="predicted_age",
            values="fraction_of_sample",
            aggfunc="sum",
            fill_value=0.0,
            observed=True,
        )
        .reindex(index=order["sample_plot_id"].tolist())
        .reindex(columns=GW_ORDER, fill_value=0.0)
    )

    study_gap = 1.6
    x_positions: list[float] = []
    current_x = 0.0
    previous_study: str | None = None
    for study_id in order["study_id"].tolist():
        if previous_study is not None and study_id != previous_study:
            current_x += study_gap
        x_positions.append(current_x)
        current_x += 1.0
        previous_study = study_id
    x = np.array(x_positions, dtype=float)

    width = max(13.0, 0.28 * float(x.max() + 1.0) + 3.2)
    fig, ax = plt.subplots(figsize=(width, 6.2))
    bottom = np.zeros(len(pivot.index), dtype=float)
    for gw in GW_ORDER:
        values = pivot[gw].to_numpy(dtype=float) * 100.0
        ax.bar(x, values, bottom=bottom, color=GW_COLORS[gw], label=gw, width=0.84)
        bottom += values

    boundaries = np.flatnonzero(order["study_id"].to_numpy()[1:] != order["study_id"].to_numpy()[:-1]) + 1
    for boundary in boundaries:
        ax.axvline(float(np.mean([x[boundary - 1], x[boundary]])), color="#4a4a4a", linewidth=0.8)

    for _, group in order.groupby("study_id", sort=False):
        center = float(np.mean(x[group.index.to_numpy()]))
        ax.text(
            center,
            1.02,
            str(group["study_label"].iloc[0]),
            transform=ax.get_xaxis_transform(),
            ha="center",
            va="bottom",
            fontsize=8.5,
        )

    ax.set_xticks(x)
    ax.set_xticklabels(order["sample_display_label"].tolist(), rotation=90, ha="center", va="top", fontsize=6.4)
    ax.set_xlim(x.min() - 1.0, x.max() + 1.0)
    ax.set_ylim(0, 100)
    ax.set_yticks(np.arange(0, 101, 20))
    ax.set_ylabel("% by predicted cell stage")
    handles, legend_labels = ax.get_legend_handles_labels()
    fig.suptitle(f"Sample composition\nby predicted cell stage\n{cell_set_label}\n{age_classifier_label}", y=0.985)
    fig.legend(
        handles,
        legend_labels,
        title="Predicted stage",
        frameon=False,
        ncol=min(len(GW_ORDER), 5),
        loc="upper center",
        bbox_to_anchor=(0.5, 0.855),
    )
    ax.grid(axis="y", color="#d8d8d8", linewidth=0.7)
    ax.set_axisbelow(True)
    fig.subplots_adjust(left=0.055, right=0.995, bottom=0.42, top=0.68)
    path = age_plot_dir / f"all_studies_{cell_set_slug}_shi_predicted_age_sample_composition_stacked_bar"
    save_figure(fig, path)
    return path.with_suffix(".png")


def write_predicted_age_sample_composition_outputs(
    per_study_dir: Path,
    table_dir: Path,
    age_plot_dir: Path,
) -> list[Path]:
    outputs: list[Path] = []
    for (
        cell_set_slug,
        cell_set_label,
        major_labels,
        age_classifier_slug,
        age_classifier_label,
        predicted_week_col,
    ) in AGE_CELL_SETS:
        age_summary = load_predicted_age_sample_composition(
            per_study_dir,
            cell_set_slug,
            cell_set_label,
            major_labels,
            age_classifier_slug,
            age_classifier_label,
            predicted_week_col,
        )
        table_path = table_dir / f"cross_study_shi_predicted_age_sample_composition_{cell_set_slug}.tsv"
        cell_set_plot_dir = age_plot_dir / cell_set_slug
        cell_set_plot_dir.mkdir(parents=True, exist_ok=True)
        age_summary.drop(columns=["study_plot_label"]).to_csv(table_path, sep="\t", index=False)
        outputs.append(table_path)
        for study_id, _ in STUDIES:
            if not age_summary["study_id"].eq(study_id).any():
                continue
            outputs.append(plot_age_composition_one_study(age_summary, study_id, cell_set_plot_dir))
        outputs.append(plot_age_composition_all_studies(age_summary, cell_set_plot_dir))
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
