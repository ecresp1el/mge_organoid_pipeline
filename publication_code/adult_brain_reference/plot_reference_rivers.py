#!/usr/bin/env python3
"""Adult brain reference assignments displayed as regional and fine-subtype rivers."""
from __future__ import annotations
import argparse
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams.update({"font.family": "DejaVu Sans", "svg.fonttype": "none", "pdf.fonttype": 42})
import matplotlib.pyplot as plt
from matplotlib.patches import PathPatch, Rectangle
from matplotlib.path import Path as MplPath
import numpy as np
import pandas as pd

ADULT_ORDER = [
    "Cortical SST+ LRP neurons",
    "Cortical SST+ nMt neurons",
    "Cortical SST+ Mt neurons",
    "Cortical PV+ basket neurons",
    "Cortical PV+ Chandelier neurons",
    "Subpallial SST+ LRP neurons",
    "Subpallial SST+ neurons",
    "Subpallial PV+ neurons",
    "Subpallial Cholinergic neurons",
]

# Regional colours are shared by coarse assignments, fine subtypes, and group strips.
REGION_COLORS = {"Pallial/cortical": "#8E6BB0", "Subpallial": "#D782A6"}


def reference_label_color(label: str) -> str:
    """Return the colour of a reference label's regional group.

    Cortical/pallial targets use purple and subpallial targets use pink.
    Unrecognized labels use neutral grey. Counts and assignments are unchanged.
    """
    return REGION_COLORS.get(label_group(label), "#777777")


PALLIAL_ORDER = ["Pallial/cortical", "Subpallial"]

MAJOR_SUBTYPE_ORDER = [
    "Pallial/cortical Pvalb",
    "Pallial/cortical Chandelier",
    "Pallial/cortical Sst",
    "Pallial/cortical Lamp5 Lhx6",
    "Pallial/cortical Vip",
    "Pallial/cortical Lamp5",
    "Pallial/cortical Sncg",
    "Pallial/cortical Pax6",
    "Pallial/cortical Medium spiny neuron",
    "Subpallial Pvalb",
    "Subpallial Chandelier",
    "Subpallial Sst",
    "Subpallial Lamp5 Lhx6",
    "Subpallial Vip",
    "Subpallial Lamp5",
    "Subpallial Medium spiny neuron",
    "Subpallial Eccentric medium spiny neuron",
    "Subpallial Cholinergic neurons",
    "Subpallial Cholinergic-GABA neurons",
    "Subpallial Pure cholinergic neurons",
]

FINE_SUBTYPE_ORDER = [
    "Cortical PV+ basket neurons",
    "Cortical PV+ Chandelier neurons",
    "Cortical SST+ LRP neurons",
    "Cortical SST+ Mt neurons",
    "Cortical SST+ nMt neurons",
    "Pallial/cortical Medium spiny neuron",
    "Pallial/cortical Vip",
    "Pallial/cortical Lamp5",
    "Pallial/cortical Sncg",
    "Pallial/cortical Pax6",
    "Pallial/cortical Lamp5 Lhx6",
    "Subpallial PV+ neurons",
    "Subpallial SST+ LRP neurons",
    "Subpallial SST+ neurons",
    "Subpallial Cholinergic neurons",
    "Subpallial Vip",
    "Subpallial Lamp5",
    "Subpallial Lamp5 Lhx6",
    "Subpallial Chandelier",
    "Subpallial Medium spiny neuron",
    "Subpallial Eccentric medium spiny neuron",
    "Subpallial Cholinergic-GABA neurons",
    "Subpallial Pure cholinergic neurons",
]

OTHER_SELECTED_REFERENCE = "Other selected reference"

def adult_label_order(labels: list[str]) -> list[str]:
    """Order known reference subtypes first, followed by other labels alphabetically.
    
    Input: subtype names. Returns: an ordered list used when no explicit target order is supplied."""
    known = [label for label in ADULT_ORDER if label in labels]
    extra = sorted(label for label in labels if label not in known)
    return known + extra


def draw_river(
    edges: pd.DataFrame,
    outpath: Path,
    right_order: list[str] | None = None,
    left_order: list[str] | None = None,
    right_title: str = "Adult brain reference subtype",
    plot_title: str = "DIV90 classes mapped to adult reference subtypes",
    right_group_getter=None,
    include_zero_right: bool = False,
    ax=None,
) -> None:
    """Draw count-proportional ribbons between query classes and adult reference labels.
    
    Inputs: div90_class/adult_subtype/n_cells edge table and optional node orders.
    Rectangles and ribbon widths use the supplied edge counts, with spacing between
    nodes. Target labels determine ribbon colours. Returns no value; draws on ax
    or saves a standalone plot when outpath is supplied."""
    if edges.empty:
        raise ValueError("No assigned cells available for river plot.")
    edges = edges.copy()
    edges["div90_class"] = edges["div90_class"].astype(str)
    edges["adult_subtype"] = edges["adult_subtype"].astype(str)
    left_sum = edges.groupby("div90_class")["n_cells"].sum()
    if left_order is None:
        left_order = left_sum.sort_values(ascending=False).index.astype(str).tolist()
    left_totals = left_sum.reindex([label for label in left_order if label in left_sum.index]).dropna()
    right_sum = edges.groupby("adult_subtype")["n_cells"].sum()
    if right_order is None:
        right_order = adult_label_order(right_sum.index.astype(str).tolist())
    if include_zero_right:
        right_totals = right_sum.reindex(right_order, fill_value=0)
    else:
        right_totals = right_sum.reindex([label for label in right_order if label in right_sum.index]).dropna()
    total = float(edges["n_cells"].sum())

    def positions(totals: pd.Series, group_getter=None) -> tuple[dict[str, tuple[float, float]], float]:
        """Allocate vertical node intervals from counts and between-group gaps.
        
        Returns: label-to-interval mapping and height per cell; left/right scales
        can differ because the numbers of gaps differ between the two sides."""
        labels = totals.index.astype(str).tolist()
        base_gap = 0.016
        group_gap = 0.052
        gaps = []
        for prev, cur in zip(labels, labels[1:]):
            if group_getter is not None and group_getter(prev) != group_getter(cur):
                gaps.append(group_gap)
            else:
                gaps.append(base_gap)
        usable = 1.0 - sum(gaps)
        if usable <= 0:
            raise ValueError("Too many river labels to fit with requested gaps.")
        scale = usable / total
        y = 1.0
        pos = {}
        for i, (label, value) in enumerate(totals.items()):
            h = scale * float(value)
            pos[str(label)] = (y - h, y)
            if i < len(gaps):
                y -= h + gaps[i]
        return pos, scale

    # Recompute geometry from displayed edges; hidden ribbons contribute
    # neither node height nor ribbon width in this panel.
    left_pos, left_scale = positions(left_totals)
    right_pos, right_scale = positions(right_totals, group_getter=right_group_getter)
    # Each ribbon consumes its share of the remaining source and target node.
    left_cursor = {k: v[1] for k, v in left_pos.items()}
    right_cursor = {k: v[1] for k, v in right_pos.items()}
    left_rank = {label: i for i, label in enumerate(left_totals.index.astype(str))}
    right_rank = {label: i for i, label in enumerate(right_totals.index.astype(str))}

    def text_positions(pos: dict[str, tuple[float, float]], min_sep: float = 0.055) -> dict[str, float]:
        """Separate label anchors vertically while preserving node order.
        
        Returns text y positions; leader lines connect displaced labels to their nodes."""
        labels = list(pos)
        desired = np.array([(pos[label][0] + pos[label][1]) / 2 for label in labels], dtype=float)
        if len(labels) <= 1:
            return {labels[0]: float(desired[0])} if labels else {}
        placed = desired.copy()
        placed[0] = min(0.98, placed[0])
        for i in range(1, len(placed)):
            placed[i] = min(placed[i], placed[i - 1] - min_sep)
        if placed[-1] < 0.02:
            placed += 0.02 - placed[-1]
            for i in range(len(placed) - 2, -1, -1):
                placed[i] = max(placed[i], placed[i + 1] + min_sep)
            if placed[0] > 0.98:
                placed -= placed[0] - 0.98
        return {label: float(y) for label, y in zip(labels, placed)}

    fig_height = max(7, 0.95 * max(len(left_totals), len(right_totals)))
    if ax is None:
        fig, ax = plt.subplots(figsize=(12, fig_height))
    else:
        fig = ax.figure
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    left_x0, left_x1 = 0.05, 0.12
    right_x0, right_x1 = 0.88, 0.95
    edges_sorted = edges.assign(
        _left_rank=edges["div90_class"].map(left_rank).fillna(10_000).astype(int),
        _right_rank=edges["adult_subtype"].map(right_rank).fillna(10_000).astype(int),
    ).sort_values(["_left_rank", "_right_rank", "div90_class", "adult_subtype"])
    for _, row in edges_sorted.iterrows():
        left = str(row["div90_class"])
        right = str(row["adult_subtype"])
        value = float(row["n_cells"])
        if value <= 0 or left not in left_pos or right not in right_pos:
            continue
        left_h = value * left_scale
        right_h = value * right_scale
        ly1 = left_cursor[left]
        ly0 = ly1 - left_h
        left_cursor[left] = ly0
        ry1 = right_cursor[right]
        ry0 = ry1 - right_h
        right_cursor[right] = ry0
        # Cubic Bezier boundaries connect equal cell counts at both ends.
        # Physical heights can differ because each side has different gaps.
        verts = [
            (left_x1, ly0),
            (0.38, ly0),
            (0.62, ry0),
            (right_x0, ry0),
            (right_x0, ry1),
            (0.62, ry1),
            (0.38, ly1),
            (left_x1, ly1),
            (left_x1, ly0),
        ]
        codes = [
            MplPath.MOVETO,
            MplPath.CURVE4,
            MplPath.CURVE4,
            MplPath.CURVE4,
            MplPath.LINETO,
            MplPath.CURVE4,
            MplPath.CURVE4,
            MplPath.CURVE4,
            MplPath.CLOSEPOLY,
        ]
        patch = PathPatch(
            MplPath(verts, codes),
            facecolor=reference_label_color(right),
            alpha=0.28,
            edgecolor="none",
        )
        ax.add_patch(patch)

    # Move text anchors only; node intervals and their counts remain fixed.
    left_text_y = text_positions(left_pos)
    right_text_y = text_positions(right_pos)
    for label, (y0, y1) in left_pos.items():
        center = (y0 + y1) / 2
        text_y = left_text_y[label]
        ax.add_patch(Rectangle((left_x0, y0), left_x1 - left_x0, y1 - y0, color="#C9CDD3", ec="white", lw=0.7))
        if abs(text_y - center) > 0.002:
            ax.plot([left_x0, left_x0 - 0.012], [center, text_y], color="#777777", lw=0.35, alpha=0.8)
        ax.text(left_x0 - 0.015, text_y, f"{label}\n{int(left_totals[label]):,}", ha="right", va="center", fontsize=8)
    for label, (y0, y1) in right_pos.items():
        center = (y0 + y1) / 2
        text_y = right_text_y[label]
        color = reference_label_color(label)
        ax.add_patch(Rectangle((right_x0, y0), right_x1 - right_x0, y1 - y0, color=color, ec="white", lw=0.7))
        if abs(text_y - center) > 0.002:
            ax.plot([right_x1, right_x1 + 0.012], [center, text_y], color="#777777", lw=0.35, alpha=0.8)
        ax.text(right_x1 + 0.015, text_y, f"{label}\n{int(right_totals[label]):,}", ha="left", va="center", fontsize=8)

    if right_group_getter is not None:
        group_spans: dict[str, list[float]] = {}
        for label, (y0, y1) in right_pos.items():
            group = right_group_getter(label)
            if group not in {"Pallial/cortical", "Subpallial"}:
                continue
            group_spans.setdefault(group, []).extend([y0, y1])
        group_colors = REGION_COLORS
        group_x0, group_x1 = 0.858, 0.875
        for group in ["Pallial/cortical", "Subpallial"]:
            if group not in group_spans:
                continue
            y0 = min(group_spans[group])
            y1 = max(group_spans[group])
            ax.add_patch(Rectangle((group_x0, y0), group_x1 - group_x0, y1 - y0, color=group_colors[group], ec="white", lw=0.6))
            ax.text(
                (group_x0 + group_x1) / 2,
                (y0 + y1) / 2,
                group.replace("/cortical", ""),
                ha="center",
                va="center",
                rotation=90,
                fontsize=8,
                color="white",
                fontweight="bold",
            )

    ax.text((left_x0 + left_x1) / 2, 1.03, "DIV90 class", ha="center", va="bottom", fontsize=11, fontweight="bold")
    ax.text((right_x0 + right_x1) / 2, 1.03, right_title, ha="center", va="bottom", fontsize=11, fontweight="bold")
    ax.set_title(plot_title, fontsize=13, pad=20)
    if outpath is not None:
        fig.savefig(outpath.with_suffix(".png"), dpi=300, bbox_inches="tight")
        fig.savefig(outpath.with_suffix(".pdf"), bbox_inches="tight")
        plt.close(fig)
def make_label_maps(ref: pd.DataFrame) -> tuple[dict[str, str], dict[str, str]]:
    """Build fine-to-major and fine-to-regional reference label lookups.
    
    Input: reference label hierarchy. Returns two dictionaries. Each fine label
    must have one parent at each level; conflicting hierarchies raise an error.
    The accounting label Other selected reference is excluded from these lookups."""
    required = ["unified_leaf_subtype", "unified_major_subtype_roi", "unified_pallial_subpallial_bin"]
    missing = [col for col in required if col not in ref.columns]
    if missing:
        raise ValueError(f"Reference metadata missing unified label columns: {missing}")
    maps = (
        ref.loc[~ref["unified_leaf_subtype"].eq(OTHER_SELECTED_REFERENCE), required]
        .drop_duplicates()
        .copy()
    )
    conflicts = maps.groupby("unified_leaf_subtype")[["unified_major_subtype_roi", "unified_pallial_subpallial_bin"]].nunique()
    conflicts = conflicts.loc[(conflicts > 1).any(axis=1)]
    if not conflicts.empty:
        raise ValueError(f"Unified leaf labels map to multiple parents: {conflicts.index.tolist()}")
    maps = maps.drop_duplicates("unified_leaf_subtype")
    leaf_to_major = dict(zip(maps["unified_leaf_subtype"].astype(str), maps["unified_major_subtype_roi"].astype(str)))
    leaf_to_bin = dict(zip(maps["unified_leaf_subtype"].astype(str), maps["unified_pallial_subpallial_bin"].astype(str)))
    return leaf_to_major, leaf_to_bin


def derive_obs_from_unified(unified_obs: pd.DataFrame, ref: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Derive regional and major labels from each cell's fine-label prediction.
    
    Inputs: cell-level predictions and the reference hierarchy. Returns regional,
    major, and fine prediction tables. Coarser labels inherit the winning fine
    label; there is no independent classifier or re-vote at the coarser levels."""
    leaf_to_major, leaf_to_bin = make_label_maps(ref)
    fine_obs = unified_obs.copy()
    subtype_obs = unified_obs.copy()
    pallial_obs = unified_obs.copy()
    # Aggregate the winning leaf annotation through the reference hierarchy.
    # No new neighbor search or confidence threshold is applied at these levels.
    subtype_obs["predicted.id"] = subtype_obs["predicted.id"].astype(str).map(leaf_to_major).fillna("Other")
    pallial_obs["predicted.id"] = pallial_obs["predicted.id"].astype(str).map(leaf_to_bin).fillna("Other")
    if subtype_obs["predicted.id"].eq("Other").any() or pallial_obs["predicted.id"].eq("Other").any():
        raise ValueError("Unified transfer produced labels that could not be mapped to major/bin parents.")
    return pallial_obs, subtype_obs, fine_obs


def make_edges(obs: pd.DataFrame, source_col: str = "div90_broad_class") -> pd.DataFrame:
    """Count cells for each query-class and predicted-reference-label pair.
    
    Input: a cell-level assignment table and source column. Returns a three-column
    edge table: div90_class, adult_subtype, and n_cells."""
    edges = (
        obs.groupby([source_col, "predicted.id"], dropna=False)
        .size()
        .reset_index(name="n_cells")
        .rename(columns={source_col: "div90_class", "predicted.id": "adult_subtype"})
    )
    return edges


def complete_order(order: list[str] | None, labels: pd.Series) -> list[str] | None:
    """Append unlisted target labels alphabetically to an explicit display order.
    
    Returns None when no order is supplied, otherwise the expanded label list."""
    if order is None:
        return None
    seen = set(order)
    extras = sorted(str(label) for label in labels.astype(str).unique() if str(label) not in seen)
    return list(order) + extras


def filter_edges_for_plot(edges: pd.DataFrame, min_cells: int) -> pd.DataFrame:
    """Keep ribbons supported by at least min_cells query cells.
    
    Returns a filtered copy of the edge table. This display threshold does not
    change cell predictions or the full edge counts."""
    if min_cells <= 1:
        return edges.copy()
    return edges.loc[edges["n_cells"].astype(float) >= float(min_cells)].copy()


def complete_left_order(edges: pd.DataFrame, left_order: list[str]) -> list[str]:
    """Append source classes absent from the common source ordering.
    
    Inputs: edge table and preferred order. Returns an order with extra classes
    ranked by descending displayed cell count."""
    seen = set(left_order)
    extras = (
        edges.loc[~edges["div90_class"].astype(str).isin(seen)]
        .groupby("div90_class")["n_cells"]
        .sum()
        .sort_values(ascending=False)
        .index.astype(str)
        .tolist()
    )
    return list(left_order) + extras


def left_order_by_pallial_fraction(pallial_edges: pd.DataFrame) -> list[str]:
    """Rank source classes by their fraction of Pallial/cortical assignments.
    
    Input: regional edges before display filtering. Returns descending fraction
    order, with total cell count and class name as deterministic tie breakers."""
    pivot = pallial_edges.pivot_table(
        index="div90_class",
        columns="adult_subtype",
        values="n_cells",
        aggfunc="sum",
        fill_value=0,
    )
    for col in PALLIAL_ORDER:
        if col not in pivot.columns:
            pivot[col] = 0
    pivot["total"] = pivot.sum(axis=1)
    pivot["pallial_fraction"] = np.where(
        pivot["total"].to_numpy() > 0,
        pivot["Pallial/cortical"].to_numpy(dtype=float) / pivot["total"].to_numpy(dtype=float),
        0.0,
    )
    ordered = (
        pivot.reset_index()
        .sort_values(["pallial_fraction", "total", "div90_class"], ascending=[False, False, True])
        ["div90_class"]
        .astype(str)
        .tolist()
    )
    return ordered


def label_group(label: str) -> str:
    """Assign a target-label prefix to Pallial/cortical, Subpallial, or Other.
    
    Returns a display group used for node spacing and the regional colour strip."""
    label = str(label)
    if label.startswith(("Pallial/cortical", "Cortical")):
        return "Pallial/cortical"
    if label.startswith("Subpallial"):
        return "Subpallial"
    if label in {"Medium spiny neuron", "Eccentric medium spiny neuron"}:
        return "Subpallial"
    return "Other"


def draw_combined_rivers(
    pallial_edges: pd.DataFrame,
    fine_edges: pd.DataFrame,
    out_prefix: Path,
    left_order: list[str],
) -> None:
    """Assemble regional and fine-subtype river panels.
    
    Inputs: displayed edge tables, output prefix, and one shared source order.
    Each panel recalculates node/ribbon geometry from its displayed edges.
    Outputs: one combined 300 dpi PNG and PDF/SVG versions."""
    panels = [
        (
            pallial_edges,
            PALLIAL_ORDER,
            "Adult brain reference region",
            "Pallial vs subpallial",
        ),
        (
            fine_edges,
            FINE_SUBTYPE_ORDER,
            "Adult brain reference subtype",
            "Fine cell types",
        ),
    ]
    fig, axes = plt.subplots(1, 2, figsize=(25, 13), constrained_layout=False)
    for ax, (edges, order, right_title, plot_title) in zip(axes, panels):
        order = complete_order(order, edges["adult_subtype"])
        plot_left_order = complete_left_order(edges, left_order)
        draw_river(
            edges,
            None,
            right_order=order,
            left_order=plot_left_order,
            right_title=right_title,
            plot_title=plot_title,
            right_group_getter=label_group,
            include_zero_right=False,
            ax=ax,
        )
    fig.subplots_adjust(left=0.04, right=0.985, top=0.9, bottom=0.04, wspace=0.55)
    fig.savefig(out_prefix.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(out_prefix.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(out_prefix.with_suffix(".svg"), bbox_inches="tight")
    plt.close(fig)

def main() -> None:
    """Read cell assignments and label hierarchy; export one combined river figure.

    Arguments supply the files directly, independent of their storage layout.
    Full and displayed edge tables retain the counts used to size each ribbon.
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--predictions", type=Path, required=True, help="Cell-level fine-label predictions, TSV or TSV.gz")
    parser.add_argument("--label-map", type=Path, required=True, help="Fine, major, and regional reference-label hierarchy")
    parser.add_argument("--outdir", type=Path, required=True, help="New output folder")
    parser.add_argument("--source-class-col", default="div90_broad_class")
    parser.add_argument("--min-edge-cells", type=int, default=10)
    args = parser.parse_args()
    if args.outdir.exists():
        parser.error("--outdir must be a new folder")
    if args.min_edge_cells < 1:
        parser.error("--min-edge-cells must be positive")
    obs = pd.read_csv(args.predictions, sep="\t")
    required = [args.source_class_col, "predicted.id"]
    missing = set(required) - set(obs)
    if missing:
        parser.error("Missing prediction columns: " + ", ".join(sorted(missing)))
    if obs.empty or obs[required].isna().any().any():
        parser.error("Predictions must contain cells with nonmissing source and target labels")
    ref = pd.read_csv(args.label_map, sep="\t")
    label_columns = ["unified_leaf_subtype", "unified_major_subtype_roi", "unified_pallial_subpallial_bin"]
    if set(label_columns) - set(ref) or ref[label_columns].isna().any().any():
        parser.error("Label map needs complete fine, major, and regional label columns")
    pallial, major, fine = derive_obs_from_unified(obs, ref)
    levels = {"regional": pallial, "major_subtype": major, "fine_subtype": fine}
    # Compute a single source-cluster order from all regional assignments before
    # ribbon filtering; use that same order in both panels.
    order = left_order_by_pallial_fraction(make_edges(pallial, args.source_class_col))
    edges = {}
    totals = []
    tables = args.outdir / "tables"
    tables.mkdir(parents=True)
    for level, cells in levels.items():
        full = make_edges(cells, args.source_class_col)
        shown = filter_edges_for_plot(full, args.min_edge_cells)
        if shown.empty:
            raise ValueError(f"No {level} ribbons meet the minimum cell count")
        full.to_csv(tables / f"{level}_edges.tsv", sep="\t", index=False)
        shown.to_csv(tables / f"{level}_displayed_edges.tsv", sep="\t", index=False)
        totals.append({"level": level, "all_cells": int(full.n_cells.sum()),
                       "displayed_cells": int(shown.n_cells.sum()),
                       "minimum_edge_cells": args.min_edge_cells})
        edges[level] = shown
    pd.DataFrame(totals).to_csv(tables / "displayed_cell_counts.tsv", sep="\t", index=False)
    pd.DataFrame({"source_class": order}).to_csv(tables / "source_class_order.tsv", sep="\t", index=False)
    figures = args.outdir / "figures"
    figures.mkdir()
    draw_combined_rivers(edges["regional"], edges["fine_subtype"],
                         figures / "adult_brain_reference_assignments", order)


if __name__ == "__main__":
    main()
