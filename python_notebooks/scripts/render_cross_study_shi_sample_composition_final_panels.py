#!/usr/bin/env python3
"""Render final Shi sample-composition panels from the v1 combined table."""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from mge_organoid_python.cross_study_shi_prediction_plots import (
    FIGURE_EXPORT_DPI,
    LABEL_COLORS,
    SHI_LABEL_ORDER,
    STUDY_LABELS,
    STUDY_ORDER,
    add_plot_coordinates,
    apply_internal_umap_plot_filters,
    load_combined_table,
)


PROJECT_ROOT_DEFAULT = Path("/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder")
RUN_LABEL_DEFAULT = "cross_study_shi_seurat_label_transfer_v1"
FINAL_DIR_DEFAULT = (
    PROJECT_ROOT_DEFAULT
    / "final_figures"
    / "fig_cross_study_shi_label_transfer_v1_umap_score_grids_candidate"
)

GW_ORDER = ["GW09", "GW12", "GW13", "GW16", "GW18"]
GW_COLORS = {
    "GW09": "#231611",
    "GW12": "#3F1C6A",
    "GW13": "#A02E6B",
    "GW16": "#EB5840",
    "GW18": "#FCC031",
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

plt.rcParams.update(
    {
        "font.family": "Arial",
        "font.sans-serif": ["Arial", "Nimbus Sans", "Liberation Sans", "DejaVu Sans"],
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    }
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, default=PROJECT_ROOT_DEFAULT)
    parser.add_argument("--run-label", default=RUN_LABEL_DEFAULT)
    parser.add_argument("--final-dir", type=Path, default=FINAL_DIR_DEFAULT)
    return parser.parse_args()


def first_number(pattern: str, text: str, default: float = 9999.0) -> float:
    match = re.search(pattern, text, flags=re.IGNORECASE)
    if not match:
        return default
    return float(match.group(1))


def canonical_gw_label(value: object) -> str:
    match = re.search(r"GW\s*([0-9]+)", str(value).upper())
    if not match:
        return "unknown"
    return f"GW{int(match.group(1)):02d}"


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


def final_filtered_data(project_root: Path, run_label: str, final_dir: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    data = load_combined_table(project_root, run_label)
    filtered, filter_summary = apply_internal_umap_plot_filters(data)
    filtered = add_plot_coordinates(filtered)
    finite = np.isfinite(pd.to_numeric(filtered["UMAP1_plot"], errors="coerce")) & np.isfinite(
        pd.to_numeric(filtered["UMAP2_plot"], errors="coerce")
    )
    finite_summary = (
        filtered.assign(_finite_plot_coordinate=finite)
        .groupby("study_id", observed=True)
        .agg(
            n_cells_after_visualization_filter=("cell_id", "size"),
            n_cells_with_finite_plot_coordinates=("_finite_plot_coordinate", "sum"),
        )
        .reset_index()
    )
    filter_summary = filter_summary.merge(finite_summary, on="study_id", how="left")
    filter_summary["n_cells_final_denominator"] = filter_summary["n_cells_with_finite_plot_coordinates"].astype(int)
    filter_summary["n_cells_removed_for_missing_plot_coordinates"] = (
        filter_summary["n_cells_after_visualization_filter"] - filter_summary["n_cells_with_finite_plot_coordinates"]
    ).astype(int)
    tables_dir = final_dir / "tables"
    tables_dir.mkdir(parents=True, exist_ok=True)
    filter_summary.to_csv(tables_dir / "cross_study_shi_sample_composition_final_filter_summary.tsv", sep="\t", index=False)
    return filtered.loc[finite].copy(), filter_summary


def add_sample_metadata(data: pd.DataFrame) -> pd.DataFrame:
    out = data.copy()
    out["sample"] = out["sample"].astype(str)
    out["sample_label"] = out["sample_label"].fillna(out["sample"]).astype(str)
    sample_meta = out[["study_id", "sample", "sample_label"]].drop_duplicates().copy()
    sample_records = [
        sample_display_and_order(str(row.study_id), row.sample, row.sample_label)
        for row in sample_meta.itertuples(index=False)
    ]
    sample_meta["sample_display_label"] = [record[0] for record in sample_records]
    sample_meta["sample_age_order"] = [record[1] for record in sample_records]
    sample_meta["sample_order_within_study"] = [record[2] for record in sample_records]
    sample_meta["sample_plot_id"] = sample_meta["study_id"].astype(str) + "::" + sample_meta["sample"].astype(str)
    annotated = out.merge(sample_meta, on=["study_id", "sample", "sample_label"], how="left", validate="many_to_one")
    annotated["sample_label"] = annotated["sample_display_label"]
    return annotated


def sample_fraction_table(
    data: pd.DataFrame,
    *,
    value_col: str,
    value_name: str,
    category_order: list[str],
    table_kind: str,
) -> pd.DataFrame:
    major = data.loc[data["shi_seurat_full_predicted_shi_label"].astype(str).isin(SHI_LABEL_ORDER)].copy()
    major[value_name] = major[value_col].astype(str)
    if value_col == "shi_seurat_full_predicted_shi_week_label":
        major[value_name] = major[value_name].map(canonical_gw_label)
    major = major.loc[major[value_name].isin(category_order)].copy()
    counts = (
        major.groupby(
            [
                "study_id",
                "sample",
                "sample_label",
                "sample_display_label",
                "sample_age_order",
                "sample_order_within_study",
                "sample_plot_id",
                value_name,
            ],
            observed=True,
        )
        .size()
        .rename("n_cells")
        .reset_index()
    )
    sample_totals = counts.groupby(["study_id", "sample_plot_id"], observed=True)["n_cells"].transform("sum")
    counts["fraction_of_sample"] = counts["n_cells"] / sample_totals
    counts["study_label"] = counts["study_id"].map(STUDY_LABELS)
    counts["study_plot_order"] = counts["study_id"].map({study: idx for idx, study in enumerate(STUDY_ORDER)})
    counts["table_kind"] = table_kind
    counts["denominator"] = "All cells with a whole-Shi major-label winner after final visualization filters"
    counts["value_source_column"] = value_col
    return counts


def sample_order(table: pd.DataFrame) -> pd.DataFrame:
    order = (
        table.groupby(
            [
                "study_id",
                "study_label",
                "study_plot_order",
                "sample_plot_id",
                "sample_display_label",
                "sample_age_order",
                "sample_order_within_study",
            ],
            observed=True,
            as_index=False,
        )["n_cells"]
        .sum()
        .sort_values(["study_plot_order", "sample_age_order", "sample_order_within_study", "sample_display_label"], kind="stable")
        .reset_index(drop=True)
    )
    return order


def save_figure(fig: plt.Figure, stem: Path) -> list[Path]:
    outputs = []
    for ext in ("png", "pdf", "svg"):
        path = stem.parent / ext / f"{stem.name}.{ext}"
        path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(path, dpi=FIGURE_EXPORT_DPI, bbox_inches="tight")
        outputs.append(path)
    plt.close(fig)
    return outputs


def plot_all_studies_stacked_bar(
    table: pd.DataFrame,
    *,
    category_col: str,
    category_order: list[str],
    colors: dict[str, str],
    title: str,
    legend_title: str,
    output_stem: Path,
) -> list[Path]:
    order = sample_order(table)
    pivot = (
        table.pivot_table(
            index="sample_plot_id",
            columns=category_col,
            values="fraction_of_sample",
            aggfunc="sum",
            fill_value=0.0,
            observed=True,
        )
        .reindex(index=order["sample_plot_id"].tolist())
        .reindex(columns=category_order, fill_value=0.0)
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
    width = max(13.0, 0.30 * float(x.max() + 1.0) + 3.6)
    fig, ax = plt.subplots(figsize=(width, 6.4))
    bottom = np.zeros(len(pivot.index), dtype=float)
    for category in category_order:
        values = pivot[category].to_numpy(dtype=float) * 100.0
        ax.bar(x, values, bottom=bottom, color=colors.get(category, "#999999"), label=category, width=0.84)
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
    ax.set_ylabel("% of sample")
    handles, labels = ax.get_legend_handles_labels()
    fig.suptitle(title, y=0.985)
    fig.legend(
        handles,
        labels,
        title=legend_title,
        frameon=False,
        ncol=min(len(category_order), 5),
        loc="upper center",
        bbox_to_anchor=(0.5, 0.855),
    )
    ax.grid(axis="y", color="#d8d8d8", linewidth=0.7)
    ax.set_axisbelow(True)
    fig.subplots_adjust(left=0.055, right=0.995, bottom=0.42, top=0.68)
    return save_figure(fig, output_stem)


def main() -> None:
    args = parse_args()
    final_dir = args.final_dir
    figures_dir = final_dir / "figures"
    tables_dir = final_dir / "tables"
    data, filter_summary = final_filtered_data(args.project_root, args.run_label, final_dir)
    data = add_sample_metadata(data)

    age_table = sample_fraction_table(
        data,
        value_col="shi_seurat_full_predicted_shi_week_label",
        value_name="predicted_age",
        category_order=GW_ORDER,
        table_kind="whole_shi_predicted_age_sample_composition",
    )
    class_table = sample_fraction_table(
        data,
        value_col="shi_seurat_full_predicted_shi_label",
        value_name="predicted_shi_major_label",
        category_order=SHI_LABEL_ORDER,
        table_kind="whole_shi_major_class_sample_composition",
    )

    age_table_path = tables_dir / "cross_study_shi_final_all_shi_major_labels_predicted_age_sample_composition.tsv"
    class_table_path = tables_dir / "cross_study_shi_final_all_shi_major_labels_major_class_sample_composition.tsv"
    age_table.to_csv(age_table_path, sep="\t", index=False)
    class_table.to_csv(class_table_path, sep="\t", index=False)

    outputs = [age_table_path, class_table_path, tables_dir / "cross_study_shi_sample_composition_final_filter_summary.tsv"]
    outputs.extend(
        plot_all_studies_stacked_bar(
            age_table,
            category_col="predicted_age",
            category_order=GW_ORDER,
            colors=GW_COLORS,
            title="Sample composition by whole-Shi predicted stage\nAll Shi major-label cells",
            legend_title="Predicted stage",
            output_stem=figures_dir / "all_studies_all_shi_major_labels_shi_predicted_age_sample_composition_stacked_bar",
        )
    )
    outputs.extend(
        plot_all_studies_stacked_bar(
            class_table,
            category_col="predicted_shi_major_label",
            category_order=SHI_LABEL_ORDER,
            colors=LABEL_COLORS,
            title="Sample composition by whole-Shi predicted major class\nAll Shi major-label cells",
            legend_title="Predicted major class",
            output_stem=figures_dir / "all_studies_all_shi_major_labels_shi_major_class_sample_composition_stacked_bar",
        )
    )

    manifest = pd.DataFrame(
        [
            {
                "output": str(path),
                "kind": "table" if path.suffix == ".tsv" else path.suffix.removeprefix("."),
                "exists": path.exists(),
                "size_bytes": path.stat().st_size if path.exists() else 0,
            }
            for path in outputs
        ]
    )
    manifest.to_csv(tables_dir / "cross_study_shi_final_sample_composition_manifest.tsv", sep="\t", index=False)
    print(manifest.to_string(index=False), flush=True)
    print(filter_summary.to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
