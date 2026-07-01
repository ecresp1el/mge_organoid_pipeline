#!/usr/bin/env python3
"""Render a cross-study marker-expression and sample-positivity final panel."""

from __future__ import annotations

import argparse
from datetime import datetime
from pathlib import Path
import re
import shutil
import subprocess

import matplotlib

matplotlib.use("Agg")
matplotlib.rcParams.update(
    {
        "font.family": "Arial",
        "font.sans-serif": ["Arial", "Nimbus Sans", "Liberation Sans", "DejaVu Sans"],
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    }
)

import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.cm import ScalarMappable
import numpy as np
import pandas as pd

from mge_organoid_python.cross_study_marker_expression import (
    BACKGROUND_POINT_COLOR,
    PDF_EXPORT_DPI,
    PNG_EXPORT_DPI,
    SVG_EXPORT_DPI,
    _apply_internal_umap_plot_filters,
    _finite_quantile,
    add_marker_plot_coordinates,
    default_cross_study_marker_specs,
    load_marker_expression_tables,
    table_dir,
)


PROJECT_ROOT_DEFAULT = Path("/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder")
RUN_LABEL_DEFAULT = "cross_study_marker_expression_v12_pv_precursors_final_candidate"
FINAL_FOLDER_DEFAULT = "fig_cross_study_lhx6_expression_v12_candidate"
GENE_DEFAULT = "LHX6"
THRESHOLD_DEFAULT = 0.5
OUTPUT_STEM_TEMPLATE = "cross_study_{gene_token}_expression_umap_violin_sample_positive"
FINAL_FOLDER_TEMPLATE = "fig_cross_study_{gene_token}_expression_v12_candidate"

STUDY_COLORS = {
    "varela_div30": "#2B6CB0",
    "varela_div90": "#2F855A",
    "siebert_2026": "#805AD5",
    "walsh": "#DD6B20",
    "bershteyn_2025": "#B83280",
    "bershteyn_2023": "#319795",
    "samarasinghe_2021": "#4A5568",
}

BERSHTEYN_2023_SAMPLE_INFO = {
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


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, default=PROJECT_ROOT_DEFAULT)
    parser.add_argument("--run-label", default=RUN_LABEL_DEFAULT)
    parser.add_argument("--final-dir", type=Path, default=None)
    parser.add_argument("--gene", default=GENE_DEFAULT)
    parser.add_argument("--positive-threshold", type=float, default=THRESHOLD_DEFAULT)
    parser.add_argument("--max-umap-cells-per-study", type=int, default=0)
    parser.add_argument("--random-state", type=int, default=0)
    return parser.parse_args()


def safe_token(value: object) -> str:
    token = re.sub(r"[^A-Za-z0-9]+", "_", str(value).strip().lower())
    token = re.sub(r"_+", "_", token).strip("_")
    return token or "value"


def normalize_gene(value: str) -> str:
    text = str(value).strip()
    aliases = {
        "PV": "PVALB",
        "PVALB": "PVALB",
        "NKX2.1": "NKX2-1",
        "NKX2-1": "NKX2-1",
        "NKX2_1": "NKX2-1",
    }
    return aliases.get(text.upper(), text.upper())


def first_number(pattern: str, text: str, default: float = 9999.0) -> float:
    match = re.search(pattern, text, flags=re.IGNORECASE)
    if not match:
        return default
    return float(match.group(1))


def sample_display_and_order(study_id: str, sample: object) -> tuple[str, float, float]:
    sample_text = str(sample)
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
        within_order = 1.0 if "9671" in sample_text else 2.0 if "9672" in sample_text else first_number(r"GSM(\d+)", sample_text)
        label = "MEL1 dFB d75" if "9671" in sample_text else "MEL1 vFB d75" if "9672" in sample_text else sample_text
        return f"DIV75 {label}", 75.0, within_order
    if study_id == "samarasinghe_2021":
        div = first_number(r"D(\d+)", sample_text)
        compact = sample_text.replace("_seurat", "").replace("_", " ")
        return f"DIV{int(div)} {compact}" if div < 9999 else f"DIV unknown {compact}", div, first_number(r"_(\d+)_seurat", sample_text)
    return sample_text, first_number(r"D(?:IV)?\s*([0-9]+)", sample_text), first_number(r"(\d+)", sample_text)


def expression_floor_colormap(color_vmax: float, floor: float) -> LinearSegmentedColormap:
    floor_fraction = float(np.clip(floor / max(color_vmax, floor + 1e-6), 0.0, 1.0))
    return LinearSegmentedColormap.from_list(
        "grayFloorBlue",
        [
            (0.0, BACKGROUND_POINT_COLOR),
            (floor_fraction, BACKGROUND_POINT_COLOR),
            (1.0, "#0000ff"),
        ],
    )


def ordered_studies(data: pd.DataFrame, specs) -> pd.DataFrame:
    spec_order = {spec.study_id: idx for idx, spec in enumerate(specs)}
    studies = data[["study_id", "study_label"]].drop_duplicates().copy()
    studies["study_plot_order"] = studies["study_id"].map(spec_order)
    studies = studies.sort_values(["study_plot_order", "study_label"], kind="stable").reset_index(drop=True)
    return studies


def downsample_umap(data: pd.DataFrame, max_cells_per_study: int, random_state: int) -> pd.DataFrame:
    if max_cells_per_study <= 0:
        return data
    parts = []
    rng = np.random.default_rng(random_state)
    for _, group in data.groupby("study_id", sort=False):
        if group.shape[0] <= max_cells_per_study:
            parts.append(group)
            continue
        chosen = rng.choice(group.index.to_numpy(), size=max_cells_per_study, replace=False)
        parts.append(group.loc[np.sort(chosen)])
    return pd.concat(parts, ignore_index=True)


def add_sample_metadata(data: pd.DataFrame, studies: pd.DataFrame) -> pd.DataFrame:
    sample_meta = data[["study_id", "study_label", "sample"]].drop_duplicates().copy()
    display = [sample_display_and_order(row.study_id, row.sample) for row in sample_meta.itertuples(index=False)]
    sample_meta["sample_display_label"] = [value[0] for value in display]
    sample_meta["sample_age_order"] = [value[1] for value in display]
    sample_meta["sample_order_within_study"] = [value[2] for value in display]
    sample_meta["sample_plot_id"] = sample_meta["study_id"].astype(str) + "::" + sample_meta["sample"].astype(str)
    sample_meta = sample_meta.merge(studies[["study_id", "study_plot_order"]], on="study_id", how="left")
    return data.merge(sample_meta, on=["study_id", "study_label", "sample"], how="left", validate="many_to_one")


def summarize_by_study(data: pd.DataFrame, gene: str, threshold: float) -> pd.DataFrame:
    rows = []
    for (study_id, study_label), group in data.groupby(["study_id", "study_label"], sort=False):
        values = pd.to_numeric(group[gene], errors="coerce").to_numpy(dtype=float)
        finite = np.isfinite(values)
        finite_values = values[finite]
        positive = finite_values >= threshold
        rows.append(
            {
                "study_id": study_id,
                "study_label": study_label,
                "gene": gene,
                "positive_threshold_log_expression": threshold,
                "n_cells": int(group.shape[0]),
                "n_finite_expression_values": int(finite_values.size),
                "n_positive_cells": int(positive.sum()),
                "pct_positive_cells": float(positive.mean() * 100.0) if finite_values.size else np.nan,
                "expr_min": float(np.min(finite_values)) if finite_values.size else np.nan,
                "expr_q25": _finite_quantile(finite_values, 0.25),
                "expr_q50": _finite_quantile(finite_values, 0.50),
                "expr_q75": _finite_quantile(finite_values, 0.75),
                "expr_q90": _finite_quantile(finite_values, 0.90),
                "expr_q95": _finite_quantile(finite_values, 0.95),
                "expr_q99": _finite_quantile(finite_values, 0.99),
                "expr_max": float(np.max(finite_values)) if finite_values.size else np.nan,
            }
        )
    return pd.DataFrame(rows)


def summarize_by_sample(data: pd.DataFrame, gene: str, threshold: float) -> pd.DataFrame:
    rows = []
    group_cols = [
        "study_id",
        "study_label",
        "study_plot_order",
        "sample",
        "sample_display_label",
        "sample_age_order",
        "sample_order_within_study",
        "sample_plot_id",
    ]
    for keys, group in data.groupby(group_cols, sort=False, observed=True):
        values = pd.to_numeric(group[gene], errors="coerce").to_numpy(dtype=float)
        finite = np.isfinite(values)
        finite_values = values[finite]
        positive = finite_values >= threshold
        row = dict(zip(group_cols, keys, strict=True))
        row.update(
            {
                "gene": gene,
                "positive_threshold_log_expression": threshold,
                "n_cells": int(group.shape[0]),
                "n_finite_expression_values": int(finite_values.size),
                "n_positive_cells": int(positive.sum()),
                "fraction_positive_cells": float(positive.mean()) if finite_values.size else np.nan,
                "pct_positive_cells": float(positive.mean() * 100.0) if finite_values.size else np.nan,
                "expr_q50": _finite_quantile(finite_values, 0.50),
                "expr_q90": _finite_quantile(finite_values, 0.90),
                "expr_q95": _finite_quantile(finite_values, 0.95),
                "denominator": f"Visualization-filtered cells with finite {gene} expression in each sample",
            }
        )
        rows.append(row)
    out = pd.DataFrame(rows)
    return out.sort_values(
        ["study_plot_order", "sample_age_order", "sample_order_within_study", "sample_display_label"],
        kind="stable",
    ).reset_index(drop=True)


def save_figure(fig: plt.Figure, final_dir: Path, stem: str) -> list[Path]:
    outputs = []
    for ext, dpi in (("png", PNG_EXPORT_DPI), ("pdf", PDF_EXPORT_DPI), ("svg", SVG_EXPORT_DPI)):
        path = final_dir / "figures" / ext / f"{stem}.{ext}"
        path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(path, dpi=dpi, bbox_inches="tight")
        outputs.append(path)
    plt.close(fig)
    return outputs


def plot_lhx6_panel(
    data: pd.DataFrame,
    study_summary: pd.DataFrame,
    sample_summary: pd.DataFrame,
    studies: pd.DataFrame,
    *,
    gene: str,
    threshold: float,
    final_dir: Path,
    output_stem: str,
) -> list[Path]:
    n_studies = studies.shape[0]
    n_samples = sample_summary.shape[0]
    fig_width = max(2.05 * n_studies + 2.2, 0.31 * n_samples + 5.8, 14.0)
    fig_height = 8.9
    fig = plt.figure(figsize=(fig_width, fig_height), constrained_layout=False)
    gs = fig.add_gridspec(
        3,
        n_studies,
        height_ratios=[1.28, 0.92, 1.05],
        hspace=0.44,
        wspace=0.06,
    )
    umap_axes = [fig.add_subplot(gs[0, idx]) for idx in range(n_studies)]
    violin_ax = fig.add_subplot(gs[1, :])
    sample_ax = fig.add_subplot(gs[2, :])

    expr_values = pd.to_numeric(data[gene], errors="coerce").to_numpy(dtype=float)
    positive_expr = expr_values[np.isfinite(expr_values) & (expr_values > threshold)]
    vmax = max(float(np.quantile(positive_expr, 0.99)) if positive_expr.size else threshold, threshold + 1e-6)
    norm = Normalize(vmin=0.0, vmax=vmax, clip=True)
    cmap = expression_floor_colormap(vmax, threshold)

    for idx, study in enumerate(studies.itertuples(index=False)):
        ax = umap_axes[idx]
        subset = data.loc[data["study_id"].astype(str) == study.study_id].copy()
        x = pd.to_numeric(subset["UMAP1_plot"], errors="coerce").to_numpy(dtype=float)
        y = pd.to_numeric(subset["UMAP2_plot"], errors="coerce").to_numpy(dtype=float)
        expr = pd.to_numeric(subset[gene], errors="coerce").to_numpy(dtype=float)
        finite = np.isfinite(x) & np.isfinite(y)
        x = x[finite]
        y = y[finite]
        expr = expr[finite]
        colored = np.isfinite(expr) & (expr >= threshold)
        ax.scatter(x, y, s=0.12, c=BACKGROUND_POINT_COLOR, linewidths=0, rasterized=True)
        if colored.any():
            ax.scatter(
                x[colored],
                y[colored],
                s=0.38,
                c=expr[colored],
                cmap=cmap,
                norm=norm,
                linewidths=0,
                rasterized=True,
            )
        label = str(study.study_label).replace("This Study, ", "This Study\n")
        summary = study_summary.loc[study_summary["study_id"].astype(str).eq(study.study_id)].iloc[0]
        ax.set_title(f"{label}\n{summary.pct_positive_cells:.1f}% >= {threshold:g}", fontsize=8.2, pad=4)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_aspect("equal", adjustable="box")
        for spine in ax.spines.values():
            spine.set_visible(False)

    violin_values = []
    violin_positions = []
    for pos, study in enumerate(studies.itertuples(index=False), start=1):
        values = pd.to_numeric(data.loc[data["study_id"].astype(str).eq(study.study_id), gene], errors="coerce")
        values = values[np.isfinite(values)].to_numpy(dtype=float)
        violin_values.append(values)
        violin_positions.append(pos)
    violins = violin_ax.violinplot(
        violin_values,
        positions=violin_positions,
        widths=0.72,
        showmeans=False,
        showmedians=False,
        showextrema=False,
    )
    for body, study_id in zip(violins["bodies"], studies["study_id"].tolist(), strict=True):
        body.set_facecolor(STUDY_COLORS.get(str(study_id), "#777777"))
        body.set_edgecolor("#222222")
        body.set_linewidth(0.35)
        body.set_alpha(0.72)
    for pos, values in zip(violin_positions, violin_values, strict=True):
        if values.size == 0:
            continue
        q25, median, q75 = np.quantile(values, [0.25, 0.5, 0.75])
        violin_ax.vlines(pos, q25, q75, color="#202020", linewidth=1.2)
        violin_ax.scatter([pos], [median], color="#202020", s=8, zorder=3)
    violin_ax.axhline(threshold, color="#202020", linewidth=0.8, linestyle=(0, (3, 2)))
    violin_ax.set_xlim(0.4, n_studies + 0.6)
    violin_ax.set_ylabel(f"{gene} log1p(CP10K)", fontsize=8.5)
    violin_ax.set_xticks(violin_positions)
    violin_ax.set_xticklabels(studies["study_label"].str.replace("This Study, ", "", regex=False).tolist(), fontsize=7)
    violin_ax.tick_params(axis="y", labelsize=7)
    violin_ax.grid(axis="y", color="#d8d8d8", linewidth=0.6)
    violin_ax.set_axisbelow(True)
    violin_ax.spines[["top", "right"]].set_visible(False)

    sample_x = []
    current_x = 0.0
    previous_study = None
    for study_id in sample_summary["study_id"].astype(str).tolist():
        if previous_study is not None and study_id != previous_study:
            current_x += 1.4
        sample_x.append(current_x)
        current_x += 1.0
        previous_study = study_id
    sample_x = np.asarray(sample_x, dtype=float)
    sample_colors = [STUDY_COLORS.get(study_id, "#777777") for study_id in sample_summary["study_id"].astype(str)]
    sample_ax.bar(sample_x, sample_summary["pct_positive_cells"].to_numpy(dtype=float), width=0.82, color=sample_colors, edgecolor="#202020", linewidth=0.25)
    sample_ax.set_ylim(0, max(100.0, float(np.nanmax(sample_summary["pct_positive_cells"].to_numpy(dtype=float))) * 1.08))
    sample_ax.set_ylabel(f"% cells {gene} >= {threshold:g}", fontsize=8.5)
    sample_ax.set_xticks(sample_x)
    sample_ax.set_xticklabels(sample_summary["sample_display_label"].tolist(), rotation=90, ha="center", va="top", fontsize=5.7)
    sample_ax.tick_params(axis="y", labelsize=7)
    sample_ax.grid(axis="y", color="#d8d8d8", linewidth=0.6)
    sample_ax.set_axisbelow(True)
    sample_ax.spines[["top", "right"]].set_visible(False)
    boundaries = np.flatnonzero(sample_summary["study_id"].to_numpy()[1:] != sample_summary["study_id"].to_numpy()[:-1]) + 1
    for boundary in boundaries:
        sample_ax.axvline(float(np.mean([sample_x[boundary - 1], sample_x[boundary]])), color="#4a4a4a", linewidth=0.75)
    for _, group in sample_summary.groupby("study_id", sort=False):
        idx = group.index.to_numpy()
        center = float(np.mean(sample_x[idx]))
        label = str(group["study_label"].iloc[0])
        sample_ax.text(center, 1.02, label, transform=sample_ax.get_xaxis_transform(), ha="center", va="bottom", fontsize=7)

    fig.suptitle(f"{gene} expression and {gene}-positive sample composition across studies", fontsize=12, y=0.985)
    fig.subplots_adjust(left=0.055, right=0.965, top=0.90, bottom=0.30)
    cax = fig.add_axes([0.972, 0.705, 0.010, 0.17])
    sm = ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])
    cbar = fig.colorbar(sm, cax=cax)
    cbar.ax.tick_params(labelsize=6, length=2)
    cbar.set_ticks([0.0, threshold, vmax])
    cbar.set_ticklabels(["0", f"{threshold:g}", f"{vmax:.2g}"])
    cbar.set_label(f"{gene} log1p(CP10K)", fontsize=7)
    return save_figure(fig, final_dir, output_stem)


def write_readme(
    final_dir: Path,
    project_root: Path,
    run_label: str,
    gene: str,
    threshold: float,
    output_stem: str,
) -> None:
    text = f"""# Cross-study {gene} expression and sample positivity

Final-figure candidate focused on `{gene}` expression across studies.

Source prepared marker-expression run:
`{project_root}/results/cross_study_marker_expression/{run_label}`

Output figure:
- `figures/png/{output_stem}.png`
- `figures/pdf/{output_stem}.pdf`
- `figures/svg/{output_stem}.svg`

Panel definition:
- Top row: per-study UMAPs colored for `{gene}` expression.
- Middle row: per-study violin distributions of `{gene}` log1p(CP10K) expression.
- Bottom row: sample-level percentage of cells with `{gene}` expression >= `{threshold:g}`.

Denominator:
Visualization-filtered cells from the marker-expression figure workflow. Samarasinghe is controls-only, and DIV90 current clusters 6 and 7 are removed as stressed cells for plotting and summaries.

Reproducibility:
- `code/` contains the renderer used to make the package.
- `tables/` contains study summaries, sample-level positivity fractions, plot filters, and the output manifest.
- `provenance/` contains git state, render metadata, and checksums.
"""
    final_dir.mkdir(parents=True, exist_ok=True)
    (final_dir / "README.md").write_text(text)


def write_provenance(
    final_dir: Path,
    project_root: Path,
    run_label: str,
    gene: str,
    threshold: float,
    outputs: list[Path],
    output_stem: str,
) -> None:
    provenance = final_dir / "provenance"
    code_dir = final_dir / "code"
    provenance.mkdir(parents=True, exist_ok=True)
    code_dir.mkdir(parents=True, exist_ok=True)
    script_path = Path(__file__).resolve()
    shutil.copy2(script_path, code_dir / script_path.name)
    git_commit = subprocess.run(["git", "rev-parse", "HEAD"], check=False, capture_output=True, text=True)
    git_status = subprocess.run(["git", "status", "--short"], check=False, capture_output=True, text=True)
    (provenance / "git_commit.txt").write_text(git_commit.stdout if git_commit.returncode == 0 else git_commit.stderr)
    (provenance / "git_status_short.txt").write_text(git_status.stdout if git_status.returncode == 0 else git_status.stderr)
    rows = [
        ("rendered_at", datetime.now().astimezone().isoformat()),
        ("project_root", str(project_root)),
        ("run_label", run_label),
        ("gene", gene),
        ("positive_threshold_log_expression", f"{threshold:g}"),
        ("source_table_dir", str(table_dir(project_root, run_label))),
        ("output_stem", output_stem),
        ("png_dpi", str(PNG_EXPORT_DPI)),
        ("pdf_dpi", str(PDF_EXPORT_DPI)),
        ("svg_dpi", str(SVG_EXPORT_DPI)),
        ("div90_visualization_filter", "exclude current clusters 6 and 7 as stressed cells"),
        ("samarasinghe_visualization_filter", "controls only"),
        ("svg_text_policy", "editable text; matplotlib svg.fonttype none; Arial/Nimbus Sans fallback"),
    ]
    rows.extend((f"output_{idx + 1}", str(path)) for idx, path in enumerate(outputs))
    pd.DataFrame(rows, columns=["key", "value"]).to_csv(provenance / "render_manifest.tsv", sep="\t", index=False)
    files = [final_dir / "README.md", code_dir / script_path.name, provenance / "render_manifest.tsv", *outputs]
    files.extend(sorted((final_dir / "tables").glob("*.tsv")))
    sha_lines = []
    for path in files:
        if not path.exists():
            continue
        digest = subprocess.run(["sha256sum", str(path)], check=False, capture_output=True, text=True)
        if digest.returncode == 0:
            sha_lines.append(digest.stdout)
    (provenance / "sha256_manifest.txt").write_text("".join(sha_lines))


def main() -> None:
    args = parse_args()
    project_root = args.project_root
    gene = normalize_gene(args.gene)
    gene_token = safe_token(gene)
    output_stem = OUTPUT_STEM_TEMPLATE.format(gene_token=gene_token)
    final_folder = FINAL_FOLDER_DEFAULT if gene == GENE_DEFAULT else FINAL_FOLDER_TEMPLATE.format(gene_token=gene_token)
    final_dir = args.final_dir or (project_root / "final_figures" / final_folder)
    table_prefix = f"cross_study_{gene_token}_expression"
    for subdir in ("figures/png", "figures/pdf", "figures/svg", "tables", "logs", "provenance", "code"):
        (final_dir / subdir).mkdir(parents=True, exist_ok=True)

    specs = default_cross_study_marker_specs(project_root)
    raw = load_marker_expression_tables(specs, project_root, args.run_label, genes=[gene])
    filtered, filter_summary = _apply_internal_umap_plot_filters(raw)
    filtered = add_marker_plot_coordinates(filtered)
    finite_plot = np.isfinite(pd.to_numeric(filtered["UMAP1_plot"], errors="coerce")) & np.isfinite(
        pd.to_numeric(filtered["UMAP2_plot"], errors="coerce")
    )
    filtered = filtered.loc[finite_plot].copy()
    plot_data = downsample_umap(filtered, args.max_umap_cells_per_study, args.random_state)
    studies = ordered_studies(filtered, specs)
    annotated = add_sample_metadata(filtered, studies)

    study_summary = summarize_by_study(annotated, gene, args.positive_threshold)
    study_summary = study_summary.merge(studies[["study_id", "study_plot_order"]], on="study_id", how="left")
    study_summary = study_summary.sort_values("study_plot_order", kind="stable").reset_index(drop=True)
    sample_summary = summarize_by_sample(annotated, gene, args.positive_threshold)

    tables_dir = final_dir / "tables"
    filter_summary.to_csv(tables_dir / f"{table_prefix}_plot_filter_summary.tsv", sep="\t", index=False)
    study_summary.to_csv(tables_dir / f"{table_prefix}_study_summary.tsv", sep="\t", index=False)
    sample_summary.to_csv(tables_dir / f"{table_prefix}_sample_positive_fraction.tsv", sep="\t", index=False)

    plot_annotated = add_sample_metadata(plot_data, studies)
    outputs = plot_lhx6_panel(
        plot_annotated,
        study_summary,
        sample_summary,
        studies,
        gene=gene,
        threshold=args.positive_threshold,
        final_dir=final_dir,
        output_stem=output_stem,
    )
    manifest = pd.DataFrame(
        [
            {
                "output": str(path),
                "kind": path.suffix.removeprefix("."),
                "exists": path.exists(),
                "size_bytes": path.stat().st_size if path.exists() else 0,
                "gene": gene,
                "positive_threshold_log_expression": args.positive_threshold,
                "run_label": args.run_label,
                "source_table_dir": str(table_dir(project_root, args.run_label)),
                "max_umap_cells_per_study": args.max_umap_cells_per_study,
            }
            for path in outputs
        ]
    )
    manifest.to_csv(tables_dir / f"{table_prefix}_output_manifest.tsv", sep="\t", index=False)

    write_readme(final_dir, project_root, args.run_label, gene, args.positive_threshold, output_stem)
    write_provenance(final_dir, project_root, args.run_label, gene, args.positive_threshold, outputs, output_stem)
    print(manifest.to_string(index=False), flush=True)
    print(study_summary[["study_label", "n_cells", "n_positive_cells", "pct_positive_cells"]].to_string(index=False), flush=True)
    print(f"Final dir: {final_dir}", flush=True)


if __name__ == "__main__":
    main()
