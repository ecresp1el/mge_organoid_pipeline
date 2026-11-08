#!/usr/bin/env python3
"""Render cell-type and gestational-week prediction-score UMAP grids."""
from __future__ import annotations
import argparse
from pathlib import Path
import re
from typing import Sequence
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.colors import LinearSegmentedColormap, Normalize
import numpy as np
import pandas as pd

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Nimbus Sans", "Liberation Sans", "DejaVu Sans"],
    "svg.fonttype": "none", "pdf.fonttype": 42,
})
FIGURE_EXPORT_DPI = 150
EXPORT_FORMATS = ["png"]
WEEK_SCORE_PREFIX = "shi_seurat_full_week_prediction_score_"
WEEK_PREDICTION_SCORE_COL = "shi_seurat_full_week_prediction_score"

BACKGROUND_GREY = "#d0d0d0"

SCORE_BLUE = "#0000ff"

STUDY_ORDER = [
    "varela_div30",
    "varela_div90",
    "siebert_2026",
    "walsh",
    "bershteyn_2023",
    "bershteyn_2025",
    "samarasinghe_2021",
]

DEFAULT_EXCLUDED_PLOT_STUDY_IDS: set[str] = set()

STUDY_LABELS = {
    "varela_div30": "This Study, DIV 30",
    "varela_div90": "This Study, DIV 90",
    "siebert_2026": "Siebert et al. 2026",
    "walsh": "Walsh et al. 2025",
    "bershteyn_2025": "Bershteyn et al. 2025",
    "bershteyn_2023": "Bershteyn et al. 2023",
    "samarasinghe_2021": "Samarasinghe et al. 2021",
}

LABEL_SCORE_MAP = {
    "MGE": "shi_seurat_full_prediction_score_MGE",
    "LGE": "shi_seurat_full_prediction_score_LGE",
    "CGE": "shi_seurat_full_prediction_score_CGE",
    "progenitor": "shi_seurat_full_prediction_score_progenitor",
    "Excitatory IPC": "shi_seurat_full_prediction_score_Excitatory_IPC",
    "Excitatory neuron": "shi_seurat_full_prediction_score_Excitatory_neuron",
    "Thalamic neurons": "shi_seurat_full_prediction_score_Thalamic_neurons",
    "Microglia": "shi_seurat_full_prediction_score_Microglia",
    "OPC": "shi_seurat_full_prediction_score_OPC",
    "Endothelial": "shi_seurat_full_prediction_score_Endothelial",
}

GREY_BLUE_CMAP = LinearSegmentedColormap.from_list("greyBlue", [BACKGROUND_GREY, SCORE_BLUE])

def natural_sort_key(value: object) -> list[object]:
    """Return mixed text/numeric sort keys so label numbers follow numeric order."""
    parts = re.split(r"(\d+)", str(value))
    return [int(part) if part.isdigit() else part.lower() for part in parts]


def parse_gw_numeric(value: object) -> float:
    """Extract the gestational-week number from a column name or label."""
    match = re.search(r"GW\s*0*(\d+)", str(value), flags=re.IGNORECASE)
    if match:
        return float(int(match.group(1)))
    nums = re.findall(r"\d+", str(value))
    if not nums:
        return np.nan
    return float(int(nums[0]))


def canonical_gw_label(value: object) -> str:
    """Format a numeric gestational week as a zero-padded GW label."""
    numeric = parse_gw_numeric(value)
    if not np.isfinite(numeric):
        return str(value)
    return f"GW{int(numeric):02d}"


def label_score_columns(data: pd.DataFrame) -> list[str]:
    """Return cell-type score columns in biological display order, then other labels."""
    known = [col for col in LABEL_SCORE_MAP.values() if col in data.columns]
    other = [
        col
        for col in data.columns
        if col.startswith("shi_seurat_full_prediction_score_")
        and col not in known
        and col != "shi_seurat_full_prediction_score"
    ]
    return known + sorted(other, key=natural_sort_key)


def week_score_columns(data: pd.DataFrame) -> list[str]:
    """Return individual gestational-week score columns in ascending week order."""
    cols = prefixed_week_score_columns(data, WEEK_SCORE_PREFIX, WEEK_PREDICTION_SCORE_COL)
    return sorted(cols, key=lambda col: (parse_gw_numeric(col), natural_sort_key(col)))


def prefixed_week_score_columns(data: pd.DataFrame, prefix: str, max_score_col: str) -> list[str]:
    """Select individual week-score columns, excluding a summary maximum-score column."""
    return [col for col in data.columns if col.startswith(prefix) and col != max_score_col]


def label_from_score_col(col: str) -> str:
    """Return the biological display name associated with a cell-type score column."""
    reverse = {value: key for key, value in LABEL_SCORE_MAP.items()}
    if col in reverse:
        return reverse[col]
    return col.removeprefix("shi_seurat_full_prediction_score_").replace("_", " ")


def week_label_from_score_col(col: str) -> str:
    """Return the displayed GW label associated with a week-score column."""
    return canonical_gw_label(col.removeprefix(WEEK_SCORE_PREFIX))


def ordered_study_ids(data: pd.DataFrame) -> list[str]:
    """Order present studies according to the shared study display order."""
    present = set(data["study_id"].astype(str))
    return [study_id for study_id in STUDY_ORDER if study_id in present]


def apply_internal_umap_plot_filters(data: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Select cells for display and return per-study inclusion counts.
    
    DIV90 clusters 6 and 7 are omitted; Samarasinghe samples containing Ctrl are
    retained. Predictions and source UMAP columns are not recalculated."""
    rows = []
    parts = []
    for study_id, group in data.groupby("study_id", sort=False):
        group = group.copy()
        plot_filter = "none"
        plot_filter_label = ""
        if str(study_id) in DEFAULT_EXCLUDED_PLOT_STUDY_IDS:
            rows.append(
                {
                    "study_id": study_id,
                    "study_label": group["study_label"].iloc[0] if "study_label" in group.columns and not group.empty else STUDY_LABELS.get(str(study_id), str(study_id)),
                    "plot_filter": "excluded_study",
                    "plot_filter_label": "Excluded from plot",
                    "n_cells_before_plot_filter": int(group.shape[0]),
                    "n_cells_after_plot_filter": 0,
                    "n_cells_removed_by_plot_filter": int(group.shape[0]),
                    "sample_values_after_plot_filter": "",
                }
            )
            continue
        kept = group
        if str(study_id) == "varela_div90":
            cluster_values = kept["cluster"].astype(str)
            stressed_mask = cluster_values.isin({"6", "7", "6.0", "7.0"})
            kept = kept.loc[~stressed_mask].copy()
            if kept.empty:
                raise ValueError("DIV90 stressed-cell UMAP plot filter removed all cells.")
            plot_filter = "exclude_div90_stressed_clusters_6_7"
            plot_filter_label = "Stressed cells removed"
        if str(study_id) == "samarasinghe_2021":
            sample_values = kept["sample"].astype(str)
            control_mask = sample_values.str.contains("Ctrl", case=False, na=False)
            kept = kept.loc[control_mask].copy()
            if kept.empty:
                raise ValueError("Samarasinghe control-only UMAP plot filter removed all cells.")
            if not kept["sample"].astype(str).str.contains("Ctrl", case=False, na=False).all():
                raise ValueError("Samarasinghe UMAP plot filter retained non-control sample values.")
            plot_filter = "sample_contains_ctrl"
            plot_filter_label = "Controls only"
        kept["plot_filter"] = plot_filter
        kept["plot_filter_label"] = plot_filter_label
        rows.append(
            {
                "study_id": study_id,
                "study_label": kept["study_label"].iloc[0] if "study_label" in kept.columns and not kept.empty else STUDY_LABELS.get(str(study_id), str(study_id)),
                "plot_filter": plot_filter,
                "plot_filter_label": plot_filter_label,
                "n_cells_before_plot_filter": int(group.shape[0]),
                "n_cells_after_plot_filter": int(kept.shape[0]),
                "n_cells_removed_by_plot_filter": int(group.shape[0] - kept.shape[0]),
                "sample_values_after_plot_filter": ";".join(sorted(kept["sample"].astype(str).unique())) if "sample" in kept.columns else "",
            }
        )
        parts.append(kept)
    if not parts:
        raise ValueError("No prediction UMAP cells remain after internal plot filters.")
    return pd.concat(parts, ignore_index=True), pd.DataFrame(rows)


def _prep_axes(ax: plt.Axes) -> None:
    """Hide ticks and spines and use equal coordinate scaling for one UMAP panel."""
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_aspect("equal", adjustable="box")
    for spine in ax.spines.values():
        spine.set_visible(False)


def format_feature_title(label: object) -> str:
    """Wrap long biological labels onto compact column headers."""
    text = str(label)
    replacements = {
        "Excitatory IPC": "Excitatory\nIPC",
        "Excitatory neuron": "Excitatory\nneuron",
        "Thalamic neurons": "Thalamic\nneurons",
    }
    return replacements.get(text, text)


def add_plot_coordinates(data: pd.DataFrame) -> pd.DataFrame:
    """Copy the UMAP coordinates and reflect DIV90 UMAP2 for display.
    
    The source umap_1 and umap_2 columns remain unchanged."""
    out = data.copy()
    out["UMAP1_plot"] = pd.to_numeric(out["umap_1"], errors="coerce")
    out["UMAP2_plot"] = pd.to_numeric(out["umap_2"], errors="coerce")
    div90_mask = out["study_id"].astype(str) == "varela_div90"
    out.loc[div90_mask, "UMAP2_plot"] = -1.0 * out.loc[div90_mask, "UMAP2_plot"]
    return out


def _study_subsets(data: pd.DataFrame) -> tuple[list[str], dict[str, pd.DataFrame], dict[str, int]]:
    """Split cells by study and retain finite display coordinates.
    
    Returns study order, cell tables, and the plotted cell count for each study."""
    ids = ordered_study_ids(data)
    subsets = {}
    counts = {}
    for study_id in ids:
        subset = data.loc[data["study_id"].astype(str) == study_id].copy()
        if "UMAP1_plot" not in subset.columns or "UMAP2_plot" not in subset.columns:
            subset = add_plot_coordinates(subset)
        finite = np.isfinite(pd.to_numeric(subset["UMAP1_plot"], errors="coerce")) & np.isfinite(
            pd.to_numeric(subset["UMAP2_plot"], errors="coerce")
        )
        subset = subset.loc[finite].copy()
        subsets[study_id] = subset
        counts[study_id] = subset.shape[0]
    return ids, subsets, counts


def row_label_for_study(study_id: str, subset: pd.DataFrame, n_cells: int) -> str:
    """Compose the study label, applicable display note, and plotted cell count."""
    label = f"{STUDY_LABELS.get(study_id, study_id)}\n(n = {n_cells:,} cells)"
    if "plot_filter_label" in subset.columns:
        filter_labels = [value for value in subset["plot_filter_label"].dropna().astype(str).unique() if value]
        if filter_labels:
            label = f"{STUDY_LABELS.get(study_id, study_id)}\n{filter_labels[0]}\n(n = {n_cells:,} cells)"
    return label


def plot_continuous_umap_grid(
    data: pd.DataFrame,
    output_path: str | Path,
    features: Sequence[tuple[str, str, float, float, str]],
    title: str,
    *,
    point_size: float = 0.36,
    background_point_size: float = 0.15,
) -> pd.DataFrame:
    """Draw a study-by-feature grid of supplied prediction scores.
    
    Each panel shows a grey cell background and a continuous grey-to-blue score
    overlay on fixed 0-1 limits. Returns per-panel cell/finite-score counts;
    exports the assembled figure without modifying predictions or source coordinates."""
    study_ids, subsets, counts = _study_subsets(data)
    n_rows = len(study_ids)
    n_cols = len(features)
    fig_width = max(1.0 * n_cols + 1.85, 5.8)
    fig_height = max(1.05 * n_rows + 1.05, 4.0)
    fig, axes = plt.subplots(
        n_rows,
        n_cols,
        figsize=(fig_width, fig_height),
        squeeze=False,
        gridspec_kw={"wspace": 0.025, "hspace": 0.035},
    )
    rows = []
    for row_idx, study_id in enumerate(study_ids):
        subset = subsets[study_id]
        x = pd.to_numeric(subset["UMAP1_plot"], errors="coerce").to_numpy()
        y = pd.to_numeric(subset["UMAP2_plot"], errors="coerce").to_numpy()
        for col_idx, (col, label, vmin, vmax, cmap) in enumerate(features):
            ax = axes[row_idx, col_idx]
            values = pd.to_numeric(subset[col], errors="coerce").to_numpy(dtype=float)
            finite = np.isfinite(values)
            ax.scatter(x, y, s=background_point_size, c=BACKGROUND_GREY, linewidths=0, rasterized=True)
            if finite.any():
                ax.scatter(
                    x[finite],
                    y[finite],
                    s=point_size,
                    c=np.clip(values[finite], vmin, vmax),
                    cmap=GREY_BLUE_CMAP if str(cmap) in {"whiteBlue", "greyBlue", "magma"} else cmap,
                    norm=Normalize(vmin=vmin, vmax=vmax),
                    linewidths=0,
                    rasterized=True,
                )
            _prep_axes(ax)
            if row_idx == 0:
                ax.set_title(format_feature_title(label), fontsize=8, linespacing=0.92)
            rows.append(
                {
                    "study_id": study_id,
                    "feature": col,
                    "n_cells_plotted": counts[study_id],
                    "n_finite_values": int(finite.sum()),
                    "vmin": vmin,
                    "vmax": vmax,
                }
            )
        axes[row_idx, 0].set_ylabel(
            row_label_for_study(study_id, subset, counts[study_id]),
            fontsize=8,
            rotation=0,
            ha="right",
            va="center",
            )
    fig.suptitle(title, fontsize=12, y=0.995)
    fig.tight_layout(rect=(0.145, 0.15, 0.995, 0.965), w_pad=0.05, h_pad=0.08)
    fig.canvas.draw()
    cbar_y = max(0.040, axes[-1, 0].get_position().y0 - 0.072)
    for col_idx, (_, label, vmin, vmax, cmap) in enumerate(features):
        pos = axes[-1, col_idx].get_position()
        width = (pos.x1 - pos.x0) * 0.60
        x0 = pos.x0 + ((pos.x1 - pos.x0) - width) / 2
        cax = fig.add_axes([x0, cbar_y, width, 0.008])
        sm = ScalarMappable(norm=Normalize(vmin=vmin, vmax=vmax), cmap=GREY_BLUE_CMAP if str(cmap) in {"whiteBlue", "greyBlue", "magma"} else cmap)
        sm.set_array([])
        cbar = fig.colorbar(sm, cax=cax, orientation="horizontal")
        cbar.ax.tick_params(labelsize=4.8, length=1.2, pad=0.8)
        cbar.set_ticks([vmin, vmax])
        cbar.set_ticklabels([f"{vmin:g}", f"{vmax:g}"])
        cbar.outline.set_linewidth(0.4)
    save_figure(fig, output_path, also_pdf=True, also_svg=True)
    return pd.DataFrame(rows)

def save_figure(fig, output_path, *, also_pdf=True, also_svg=False):
    """Export a completed grid in the requested formats at the configured DPI.

    Input: matplotlib figure and output filename. PNG is written by default;
    optional PDF/SVG exports use the same layout and editable vector text.
    """
    path = Path(output_path)
    outputs = []
    for extension in EXPORT_FORMATS:
        target = path.with_suffix("." + extension)
        fig.savefig(target, dpi=FIGURE_EXPORT_DPI, bbox_inches="tight")
        outputs.append(target)
    plt.close(fig)
    return outputs


def load_scores(paths):
    """Read prepared cell-level score tables with matching UMAP coordinates.

    Inputs: TSV/TSV.gz files containing study_id, cell_id, sample, cluster,
    umap_1/umap_2, cell-type score columns, and individual GW score columns.
    Returns one combined table. Only columns used by these two grids are loaded.
    Duplicate study/cell keys and malformed finite scores raise an error.
    """
    required = {"study_id", "cell_id", "sample", "cluster", "umap_1", "umap_2"}
    parts = []
    for path in paths:
        columns = pd.read_csv(path, sep="\t", nrows=0).columns
        missing = required - set(columns)
        if missing:
            raise ValueError("Missing columns: " + ", ".join(sorted(missing)))
        score_columns = [c for c in columns if c.startswith("shi_seurat_full_prediction_score_")
                         or c.startswith(WEEK_SCORE_PREFIX)]
        usecols = list(required | set(score_columns) | ({"study_label"} & set(columns)))
        parts.append(pd.read_csv(path, sep="\t", usecols=usecols,
                                 dtype={"study_id": str, "cell_id": str, "sample": str, "cluster": str}))
    data = pd.concat(parts, ignore_index=True)
    if data.empty or data[["study_id", "cell_id"]].isna().any().any():
        raise ValueError("Provide nonempty tables with complete study and cell identifiers")
    if data.duplicated(["study_id", "cell_id"]).any():
        raise ValueError("Duplicate cell IDs within a study")
    unknown = set(data.study_id) - set(STUDY_ORDER)
    if unknown:
        raise ValueError("Add study display labels/order for: " + ", ".join(sorted(unknown)))
    for column in label_score_columns(data) + week_score_columns(data):
        values = pd.to_numeric(data[column], errors="raise")
        finite = values[np.isfinite(values)]
        if ((finite < -1e-8) | (finite > 1 + 1e-8)).any():
            raise ValueError("Prediction scores must be between 0 and 1: " + column)
        data[column] = values
    return data


def main():
    """Apply display inclusion rules and render the two prediction-score grids.

    CLI inputs are explicit file paths; output must be a new folder. Scores and
    source embeddings are read as supplied, with no label-transfer recomputation.
    """
    global FIGURE_EXPORT_DPI, EXPORT_FORMATS
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, nargs="+", required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--dpi", type=int, default=150)
    parser.add_argument("--formats", nargs="+", choices=["png", "pdf", "svg"], default=["png"])
    args = parser.parse_args()
    if args.outdir.exists():
        parser.error("--outdir must be a new folder")
    if args.dpi <= 0:
        parser.error("--dpi must be positive")
    FIGURE_EXPORT_DPI = args.dpi
    EXPORT_FORMATS = args.formats
    data = load_scores(args.input)
    data, filters = apply_internal_umap_plot_filters(data)
    data = add_plot_coordinates(data)
    cell_features = [(c, label_from_score_col(c), 0, 1, "greyBlue") for c in label_score_columns(data)]
    week_features = [(c, week_label_from_score_col(c), 0, 1, "greyBlue") for c in week_score_columns(data)]
    if not cell_features or not week_features:
        parser.error("Both cell-type and individual GW score columns are required")
    args.outdir.mkdir(parents=True)
    filters.to_csv(args.outdir / "display_cell_counts.tsv", sep="\t", index=False)
    for filename, features, title in [
        ("cell_type_prediction_scores.png", cell_features, "Cross-study all Shi label prediction scores"),
        ("gestational_week_prediction_scores.png", week_features, "Cross-study individual Shi GW prediction scores"),
    ]:
        manifest = plot_continuous_umap_grid(data, args.outdir / filename, features, title)
        manifest.to_csv(args.outdir / (Path(filename).stem + "_cell_counts.tsv"), sep="\t", index=False)


if __name__ == "__main__":
    main()
