#!/usr/bin/env python3
"""Render Step 02 overlap figures from cached mapping tables only."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.gridspec import GridSpec
import numpy as np
import pandas as pd
from PIL import Image
import seaborn as sns


STUDIES = [
    "varela_div30",
    "varela_div90",
    "walsh",
    "bershteyn_2025",
    "bershteyn_2023",
    "siebert_2026",
]
LABELS = {
    "varela_div30": "Varela DIV30",
    "varela_div90": "Varela DIV90",
    "walsh": "Walsh",
    "bershteyn_2025": "Bershteyn 2025",
    "bershteyn_2023": "Bershteyn 2023",
    "siebert_2026": "Siebert 2026",
}
COLORS = {
    "mapped": "#367B9B",
    "ambiguous": "#E5A84B",
    "unresolved": "#C85A5A",
    "strict": "#3A9D78",
    "duplicate_loss": "#9A6FB0",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument("--dpi", type=int, default=300)
    parser.add_argument("--svg", choices=["true", "false"], default="false")
    parser.add_argument("--max-patterns", type=int, default=30)
    return parser.parse_args()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(16 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def set_style() -> None:
    sns.set_theme(style="whitegrid", context="notebook")
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "axes.titleweight": "bold",
            "axes.labelcolor": "#222222",
            "text.color": "#222222",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def save_figure(fig: plt.Figure, run_dir: Path, stem: str, dpi: int, svg: bool) -> list[dict[str, object]]:
    rows = []
    formats = ["png", "pdf"] + (["svg"] if svg else [])
    for file_format in formats:
        path = run_dir / "figures" / file_format / f"{stem}.{file_format}"
        path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(path, dpi=dpi, bbox_inches="tight", facecolor="white")
        rows.append(
            {
                "figure_id": stem,
                "format": file_format,
                "relative_path": str(path.relative_to(run_dir)),
                "dpi_requested": dpi,
                "size_bytes": path.stat().st_size,
                "sha256": sha256_file(path),
            }
        )
    plt.close(fig)
    return rows


def pattern_table(presence: pd.DataFrame, suffix: str) -> pd.DataFrame:
    columns = [f"{study}_{suffix}" for study in STUDIES]
    patterns = (
        presence[columns]
        .astype(bool)
        .groupby(columns, dropna=False)
        .size()
        .rename("n_genes")
        .reset_index()
    )
    patterns["n_studies"] = patterns[columns].sum(axis=1).astype(int)
    patterns["pattern"] = patterns[columns].apply(
        lambda row: "|".join(study for study, value in zip(STUDIES, row) if value), axis=1
    )
    # An UpSet universe consists only of identities present in at least one set.
    # Strict filtering can make an identity ineligible in all six; those rows
    # remain available in the source presence table but are not an intersection.
    patterns = patterns.loc[patterns["n_studies"].gt(0)]
    return patterns.sort_values(["n_genes", "n_studies"], ascending=[False, False]).reset_index(drop=True)


def select_patterns(patterns: pd.DataFrame, max_patterns: int) -> pd.DataFrame:
    suffix_columns = [column for column in patterns if column.endswith("_present") or column.endswith("_strict_one_to_one")]
    required = patterns.loc[patterns["n_studies"].isin([1, 5, 6])]
    top = patterns.head(max_patterns)
    selected = pd.concat([required, top]).drop_duplicates("pattern")
    selected = selected.sort_values(["n_genes", "n_studies"], ascending=[False, False]).head(max_patterns)
    return selected.reset_index(drop=True), suffix_columns


def render_upset(
    patterns: pd.DataFrame,
    presence: pd.DataFrame,
    suffix: str,
    title: str,
    subtitle: str,
    max_patterns: int,
) -> tuple[plt.Figure, pd.DataFrame]:
    selected, columns = select_patterns(patterns, max_patterns)
    fig = plt.figure(figsize=(17, 10.5), constrained_layout=True)
    grid = GridSpec(2, 2, figure=fig, width_ratios=[1.25, 4.8], height_ratios=[3.1, 2.0])
    set_ax = fig.add_subplot(grid[:, 0])
    bar_ax = fig.add_subplot(grid[0, 1])
    matrix_ax = fig.add_subplot(grid[1, 1], sharex=bar_ax)

    set_sizes = [int(presence[column].astype(bool).sum()) for column in columns]
    y = np.arange(len(STUDIES))
    set_ax.barh(y, set_sizes, color="#4C7899")
    set_ax.set_yticks(y, [LABELS[study] for study in STUDIES])
    set_ax.invert_yaxis()
    set_ax.set_xlabel("Mapped gene identities")
    set_ax.set_title("Set sizes")
    for index, value in enumerate(set_sizes):
        set_ax.text(value, index, f" {value:,}", va="center", fontsize=9)
    set_ax.spines[["top", "right", "left"]].set_visible(False)

    x = np.arange(selected.shape[0])
    bar_ax.bar(x, selected["n_genes"], color="#2F6F8F", width=0.76)
    bar_ax.set_ylabel("Genes in exact intersection pattern")
    bar_ax.set_title(f"{title}\n{subtitle}", loc="left", pad=12)
    bar_ax.tick_params(axis="x", labelbottom=False)
    bar_ax.spines[["top", "right"]].set_visible(False)
    for index, value in enumerate(selected["n_genes"]):
        bar_ax.text(index, value, f"{int(value):,}", ha="center", va="bottom", rotation=90, fontsize=7)

    matrix_ax.set_ylim(-0.75, len(STUDIES) - 0.25)
    matrix_ax.set_yticks(y, [LABELS[study] for study in STUDIES])
    matrix_ax.invert_yaxis()
    matrix_ax.set_xlabel(f"Top {selected.shape[0]} exact presence/absence patterns (forced to include 1-, 5-, and 6-study patterns)")
    matrix_ax.set_xticks(x, [str(i + 1) for i in x], fontsize=8)
    matrix_ax.grid(False)
    for col_index, (_, row) in enumerate(selected.iterrows()):
        active = [bool(row[column]) for column in columns]
        active_y = [i for i, value in enumerate(active) if value]
        matrix_ax.scatter([col_index] * len(STUDIES), y, s=28, color="#D5D5D5", zorder=2)
        if active_y:
            matrix_ax.plot([col_index, col_index], [min(active_y), max(active_y)], color="#222222", lw=1.4, zorder=3)
            matrix_ax.scatter([col_index] * len(active_y), active_y, s=42, color="#222222", zorder=4)
    matrix_ax.spines[["top", "right", "left"]].set_visible(False)
    return fig, selected


def pairwise_matrices(pairwise: pd.DataFrame, representation: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    subset = pairwise.loc[pairwise["representation"].eq(representation)]
    counts = pd.DataFrame(index=STUDIES, columns=STUDIES, dtype=float)
    jaccard = pd.DataFrame(index=STUDIES, columns=STUDIES, dtype=float)
    for study in STUDIES:
        counts.loc[study, study] = np.nan
        jaccard.loc[study, study] = 1.0
    for row in subset.itertuples(index=False):
        counts.loc[row.left_study, row.right_study] = row.n_intersection
        counts.loc[row.right_study, row.left_study] = row.n_intersection
        jaccard.loc[row.left_study, row.right_study] = row.jaccard
        jaccard.loc[row.right_study, row.left_study] = row.jaccard
    labels = [LABELS[study] for study in STUDIES]
    counts.index = counts.columns = labels
    jaccard.index = jaccard.columns = labels
    return counts, jaccard


def render_pairwise(pairwise: pd.DataFrame, outcome: pd.DataFrame) -> plt.Figure:
    mapped_counts, mapped_jaccard = pairwise_matrices(pairwise, "mapped_common_gene_id")
    strict_counts, strict_jaccard = pairwise_matrices(pairwise, "strict_one_to_one_common_gene_id")
    mapped_sizes = outcome.set_index("study_id")["mapped_unique_common_gene_ids"].astype(int)
    strict_sizes = outcome.set_index("study_id")["strict_one_to_one_common_gene_ids"].astype(int)
    for study in STUDIES:
        label = LABELS[study]
        mapped_counts.loc[label, label] = mapped_sizes[study]
        strict_counts.loc[label, label] = strict_sizes[study]

    fig, axes = plt.subplots(2, 2, figsize=(17, 14), constrained_layout=True)
    cmap_count = LinearSegmentedColormap.from_list("count", ["#F4F7F8", "#307A9A"])
    cmap_jaccard = LinearSegmentedColormap.from_list("jac", ["#F7F2E8", "#D9983D", "#337A5B"])
    for ax, matrix, title, cmap, fmt, vmax in [
        (axes[0, 0], mapped_counts, "Mapped common identities: pairwise intersection", cmap_count, ".0f", None),
        (axes[0, 1], mapped_jaccard, "Mapped common identities: Jaccard overlap", cmap_jaccard, ".3f", 1),
        (axes[1, 0], strict_counts, "Strict one-to-one identities: pairwise intersection", cmap_count, ".0f", None),
        (axes[1, 1], strict_jaccard, "Strict one-to-one identities: Jaccard overlap", cmap_jaccard, ".3f", 1),
    ]:
        sns.heatmap(matrix, annot=True, fmt=fmt, cmap=cmap, vmin=0, vmax=vmax, ax=ax, square=True,
                    linewidths=0.6, linecolor="white", cbar_kws={"shrink": 0.72})
        ax.set_title(title, pad=12)
        ax.set_xlabel("")
        ax.set_ylabel("")
        ax.tick_params(axis="x", rotation=40)
        ax.tick_params(axis="y", rotation=0)
    fig.suptitle(
        "Pairwise gene overlap across the six canonical studies\n"
        "Diagonal cells are study set sizes; off-diagonal cells show shared identities",
        fontsize=18,
        fontweight="bold",
    )
    return fig


def render_coverage(outcome: pd.DataFrame, presence: pd.DataFrame) -> plt.Figure:
    outcome = outcome.set_index("study_id").loc[STUDIES]
    y = np.arange(len(STUDIES))
    labels = [LABELS[study] for study in STUDIES]
    fig, axes = plt.subplots(1, 3, figsize=(19, 7.5), constrained_layout=True)

    left = np.zeros(len(STUDIES))
    for column, label, color in [
        ("mapped_features", "Mapped", COLORS["mapped"]),
        ("ambiguous_features", "Ambiguous", COLORS["ambiguous"]),
        ("unresolved_features", "Unresolved", COLORS["unresolved"]),
    ]:
        values = outcome[column].astype(int).to_numpy()
        axes[0].barh(y, values, left=left, color=color, label=label)
        left += values
    axes[0].set_yticks(y, labels)
    axes[0].invert_yaxis()
    axes[0].set_xlabel("Input features")
    axes[0].set_title("Mapping status by study")
    axes[0].legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.12), ncol=3)

    mapped = outcome["mapped_unique_common_gene_ids"].astype(int).to_numpy()
    strict = outcome["strict_one_to_one_common_gene_ids"].astype(int).to_numpy()
    axes[1].barh(y - 0.18, mapped, height=0.34, color=COLORS["mapped"], label="Mapped identity set")
    axes[1].barh(y + 0.18, strict, height=0.34, color=COLORS["strict"], label="Strict one-to-one set")
    axes[1].set_yticks(y, labels)
    axes[1].invert_yaxis()
    axes[1].set_xlabel("Unique common gene identities")
    axes[1].set_title("Duplicate-sensitive set sizes")
    axes[1].legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.12), ncol=2)
    for index, (full, clean) in enumerate(zip(mapped, strict)):
        axes[1].text(full, index - 0.18, f" {full:,}", va="center", fontsize=8)
        axes[1].text(clean, index + 0.18, f" {clean:,}", va="center", fontsize=8)

    degree = presence["n_studies_present"].value_counts().reindex(range(1, 7), fill_value=0)
    axes[2].bar(degree.index, degree.values, color=["#C85A5A", "#D67B55", "#DEA953", "#84A85D", "#4D9875", "#306F8E"])
    axes[2].set_xticks(range(1, 7))
    axes[2].set_xlabel("Number of studies containing mapped identity")
    axes[2].set_ylabel("Gene identities")
    axes[2].set_title("How broadly identities overlap")
    for x, value in degree.items():
        axes[2].text(x, value, f"{value:,}", ha="center", va="bottom", fontsize=9)

    fig.suptitle(
        "Gene-mapping coverage and non-overlap diagnostics\n"
        "Report-only summaries; no expression matrix was changed",
        fontsize=18,
        fontweight="bold",
    )
    for ax in axes:
        ax.spines[["top", "right"]].set_visible(False)
    return fig


def update_readme(run_dir: Path, dpi: int, svg: bool) -> None:
    readme_path = run_dir / "README.md"
    text = readme_path.read_text(encoding="utf-8")
    start = "<!-- GENE_OVERLAP_FIGURES_START -->"
    end = "<!-- GENE_OVERLAP_FIGURES_END -->"
    section = f"""{start}

## Plot-only overlap visualization extension

The completed mapping tables were rendered into overlap figures without
reopening canonical objects or reading expression matrices. PNG and PDF are
mandatory and were written at {dpi} dpi. SVG generation was `{str(svg).lower()}`.

- `gene_overlap_upset_identity_level`: exact multi-study presence/absence
  combinations for mapped common identities.
- `gene_overlap_upset_strict_one_to_one`: the same view after excluding
  duplicate-to-one identities.
- `gene_overlap_pairwise_heatmaps`: pairwise intersection counts and Jaccard
  overlap for both mapped and strict representations.
- `gene_mapping_coverage_and_nonoverlap`: per-study mapped/ambiguous/unresolved
  coverage, duplicate-sensitive set sizes, and the number of studies sharing
  each mapped identity.

The exact plotted pattern counts are in `tables/`. This extension does not
change the Step 02 review stop or choose an intersection policy.

{end}"""
    if start in text and end in text:
        before = text.split(start, 1)[0].rstrip()
        after = text.split(end, 1)[1].lstrip()
        text = before + "\n\n" + section + "\n\n" + after
    else:
        text = text.rstrip() + "\n\n" + section + "\n"
    readme_path.write_text(text, encoding="utf-8")


def write_package_checksums(run_dir: Path) -> None:
    rows = []
    excluded = {Path("provenance/package_checksums.tsv"), Path("SUCCESS.txt"), Path("FIGURES_SUCCESS.txt")}
    for path in sorted(run_dir.rglob("*")):
        relative = path.relative_to(run_dir)
        if path.is_file() and relative.parts[0] != "logs" and relative not in excluded:
            rows.append(
                {"relative_path": str(relative), "size_bytes": path.stat().st_size, "sha256": sha256_file(path)}
            )
    pd.DataFrame(rows).to_csv(run_dir / "provenance/package_checksums.tsv", sep="\t", index=False)


def main() -> None:
    args = parse_args()
    run_dir = args.run_dir.resolve()
    if not (run_dir / "SUCCESS.txt").is_file():
        raise FileNotFoundError("Completed Step 02 SUCCESS.txt is required")
    tables = run_dir / "tables"
    presence = pd.read_csv(tables / "mapped_gene_presence_matrix.tsv.gz", sep="\t")
    pairwise = pd.read_csv(tables / "pairwise_gene_overlap.tsv", sep="\t")
    outcome = pd.read_csv(tables / "gene_mapping_outcome_by_dataset.tsv", sep="\t")
    if outcome["study_id"].tolist() != STUDIES:
        raise ValueError("Unexpected study order/content")

    set_style()
    manifest_rows: list[dict[str, object]] = []
    identity_patterns = pattern_table(presence, "present")
    strict_patterns = pattern_table(presence, "strict_one_to_one")
    identity_patterns.to_csv(tables / "intersection_pattern_counts_identity_level.tsv", sep="\t", index=False)
    strict_patterns.to_csv(tables / "intersection_pattern_counts_strict_one_to_one.tsv", sep="\t", index=False)

    fig, plotted = render_upset(
        identity_patterns,
        presence,
        "present",
        "Mapped common gene identities: exact multi-study intersections",
        "A dot marks study membership; missing dots show which studies do not contain that identity pattern",
        args.max_patterns,
    )
    plotted.to_csv(tables / "plotted_intersection_patterns_identity_level.tsv", sep="\t", index=False)
    manifest_rows.extend(save_figure(fig, run_dir, "gene_overlap_upset_identity_level", args.dpi, args.svg == "true"))

    fig, plotted = render_upset(
        strict_patterns,
        presence,
        "strict_one_to_one",
        "Strict one-to-one gene identities: exact multi-study intersections",
        "Duplicate-to-one identities are excluded before counting exact overlap patterns",
        args.max_patterns,
    )
    plotted.to_csv(tables / "plotted_intersection_patterns_strict_one_to_one.tsv", sep="\t", index=False)
    manifest_rows.extend(
        save_figure(fig, run_dir, "gene_overlap_upset_strict_one_to_one", args.dpi, args.svg == "true")
    )

    manifest_rows.extend(
        save_figure(
            render_pairwise(pairwise, outcome),
            run_dir,
            "gene_overlap_pairwise_heatmaps",
            args.dpi,
            args.svg == "true",
        )
    )
    manifest_rows.extend(
        save_figure(
            render_coverage(outcome, presence),
            run_dir,
            "gene_mapping_coverage_and_nonoverlap",
            args.dpi,
            args.svg == "true",
        )
    )

    figure_manifest = pd.DataFrame(manifest_rows)
    figure_manifest.to_csv(tables / "gene_overlap_figure_manifest.tsv", sep="\t", index=False)
    for row in figure_manifest.loc[figure_manifest["format"].eq("png")].itertuples(index=False):
        with Image.open(run_dir / row.relative_path) as image:
            observed = image.info.get("dpi", (None, None))
            if observed[0] is None or abs(float(observed[0]) - args.dpi) > 0.1:
                raise ValueError(f"PNG DPI validation failed for {row.relative_path}: {observed}")

    update_readme(run_dir, args.dpi, args.svg == "true")
    completed = datetime.now(timezone.utc).astimezone().isoformat()
    (run_dir / "FIGURES_SUCCESS.txt").write_text(
        "PASS\n"
        f"completed={completed}\n"
        "source=cached_step02_mapping_tables_only\n"
        f"png_pdf_dpi={args.dpi}\n"
        f"svg={args.svg}\n"
        f"figures={figure_manifest['figure_id'].nunique()}\n"
        "review_stop=YES\n",
        encoding="utf-8",
    )
    write_package_checksums(run_dir)
    print(f"Rendered {figure_manifest['figure_id'].nunique()} gene-overlap figures", flush=True)


if __name__ == "__main__":
    main()
