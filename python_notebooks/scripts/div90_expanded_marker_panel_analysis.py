#!/usr/bin/env python3
"""Expanded DIV90 marker panels with UMAPs, DotPlots, heatmaps, and summaries."""

from __future__ import annotations

import argparse
import os
import sys
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.cm import ScalarMappable
import numpy as np
import pandas as pd


PROJECT_ROOT_DEFAULT = "/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder"


@dataclass(frozen=True)
class MarkerPanel:
    panel_id: str
    label: str
    goal: str
    genes: tuple[str, ...]


PANELS = [
    MarkerPanel(
        panel_id="panel1_developmental_migration",
        label="Developmental state / migration markers",
        goal=(
            "Determine whether clusters 1, 3, and 11 contain distinct migratory "
            "versus post-migratory interneuron populations."
        ),
        genes=("DCX", "SOX11", "STMN2", "DCLK1", "CXCR4", "ROBO1", "ROBO2", "CXCL12R"),
    ),
    MarkerPanel(
        panel_id="panel2_mge_interneuron_maturation",
        label="MGE interneuron maturation markers",
        goal="Determine whether the cortical-fated interneuron compartment contains separable maturation programs.",
        genes=("LHX6", "SOX6", "SATB1", "MAF", "MAFB", "ERBB4", "GAD1", "GAD2"),
    ),
    MarkerPanel(
        panel_id="panel3_pv_lineage_maturation",
        label="PV-lineage maturation markers",
        goal="Determine whether putative PV precursors separate into distinct maturation states despite low PVALB expression.",
        genes=("PVALB", "PPARGC1A", "KCNC1", "KCNC2", "GPR149", "ERBB4", "SOX6"),
    ),
    MarkerPanel(
        panel_id="panel4_striatal_pallidal_specification",
        label="Striatal versus pallidal specification markers",
        goal=(
            "Test whether the large MGE striatal/GP-fated compartment contains unresolved "
            "striatal versus pallidal programs."
        ),
        genes=("NKX2-1", "LHX8", "GBX2", "ISL1", "EBF1", "FOXP1", "BCL11B", "TAC1"),
    ),
]


GENE_ALIASES = {
    "CXCL12R": ("ACKR3", "CXCR7", "CXCL12R"),
    "NKX2-1": ("NKX2-1", "NKX2.1"),
}


def repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def add_src_to_path() -> None:
    src = repo_root() / "python_notebooks" / "src"
    if str(src) not in sys.path:
        sys.path.insert(0, str(src))


def unique_preserve_order(values: list[str] | tuple[str, ...]) -> list[str]:
    seen = set()
    out = []
    for value in values:
        if value in seen:
            continue
        seen.add(value)
        out.append(value)
    return out


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--project-root",
        default=os.environ.get("PROJECT_ROOT", PROJECT_ROOT_DEFAULT),
        help="Project/data root containing results/python_anndata/varela_div90.h5ad.",
    )
    parser.add_argument(
        "--h5ad",
        default=None,
        help="Optional explicit DIV90 H5AD path.",
    )
    parser.add_argument(
        "--outdir",
        default=None,
        help="Output directory. Default: PROJECT_ROOT/results/div90_expanded_marker_panel_analysis/div90_expanded_marker_panel_analysis_v1",
    )
    parser.add_argument("--ncols", type=int, default=4, help="UMAP gene panels per row.")
    parser.add_argument("--max-cells", type=int, default=None, help="Optional UMAP plotting downsample.")
    parser.add_argument("--vmax-quantile", type=float, default=0.99, help="Positive-expression color scale quantile.")
    parser.add_argument("--point-size", type=float, default=0.32, help="UMAP expressing-cell point size.")
    parser.add_argument("--background-point-size", type=float, default=0.12, help="UMAP background point size.")
    parser.add_argument(
        "--enrichment-z",
        type=float,
        default=1.0,
        help="Cluster-restricted enrichment z-score threshold for marker calls.",
    )
    return parser.parse_args()


def positive_vmax(values: pd.Series, quantile: float) -> float:
    arr = pd.to_numeric(values, errors="coerce").to_numpy(dtype=float)
    positive = arr[np.isfinite(arr) & (arr > 0)]
    if positive.size == 0:
        return 1.0
    vmax = float(np.quantile(positive, quantile))
    if not np.isfinite(vmax) or vmax <= 0:
        vmax = float(np.max(positive))
    return max(vmax, 1e-6)


def zscore_by_gene(summary: pd.DataFrame) -> pd.DataFrame:
    out = summary.copy()
    out["mean_expr_z_by_gene"] = np.nan
    for gene, idx in out.groupby("gene").groups.items():
        values = out.loc[idx, "mean_expr"].to_numpy(dtype=float)
        sd = float(np.nanstd(values))
        if not np.isfinite(sd) or sd == 0:
            out.loc[idx, "mean_expr_z_by_gene"] = 0.0
        else:
            out.loc[idx, "mean_expr_z_by_gene"] = (values - float(np.nanmean(values))) / sd
    return out


def cluster_sort_key(value: str) -> tuple[int, str]:
    text = str(value)
    try:
        return (0, f"{int(float(text)):05d}")
    except ValueError:
        return (1, text)


def ordered_clusters(values: pd.Series) -> list[str]:
    return sorted(values.astype(str).unique().tolist(), key=cluster_sort_key)


def panel_display_gene(gene: str, gene_match: pd.DataFrame) -> str:
    row = gene_match.loc[gene_match["gene"] == gene]
    if row.empty:
        return gene
    matched = str(row.iloc[0].get("matched_feature", ""))
    if matched and matched != gene:
        return f"{gene}\n({matched})"
    return gene


def expression_colormap() -> LinearSegmentedColormap:
    return LinearSegmentedColormap.from_list("grey_to_blue", ["#d7d7d7", "#3157ff", "#07008f"])


def plot_umap_grid(
    data: pd.DataFrame,
    genes: list[str],
    gene_match: pd.DataFrame,
    output_path: Path,
    *,
    title: str,
    ncols: int,
    max_cells: int | None,
    vmax_quantile: float,
    point_size: float,
    background_point_size: float,
) -> pd.DataFrame:
    plot_data = data.copy()
    if max_cells and plot_data.shape[0] > max_cells:
        plot_data = plot_data.sample(n=max_cells, random_state=0).sort_index()
    ncols = max(1, min(ncols, len(genes)))
    nrows = int(np.ceil(len(genes) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(3.25 * ncols, 3.15 * nrows + 0.5), squeeze=False)
    x = pd.to_numeric(plot_data["umap_1"], errors="coerce").to_numpy(dtype=float)
    y = pd.to_numeric(plot_data["umap_2"], errors="coerce").to_numpy(dtype=float)
    finite_umap = np.isfinite(x) & np.isfinite(y)
    x = x[finite_umap]
    y = y[finite_umap]
    plot_data = plot_data.loc[finite_umap].copy()
    x_pad = (np.nanmax(x) - np.nanmin(x)) * 0.03
    y_pad = (np.nanmax(y) - np.nanmin(y)) * 0.03
    xlim = (np.nanmin(x) - x_pad, np.nanmax(x) + x_pad)
    ylim = (np.nanmin(y) - y_pad, np.nanmax(y) + y_pad)
    cmap = expression_colormap()
    rows = []

    for i, gene in enumerate(genes):
        ax = axes[i // ncols, i % ncols]
        expr = pd.to_numeric(plot_data[gene], errors="coerce").to_numpy(dtype=float)
        positive = np.isfinite(expr) & (expr > 0)
        vmax = positive_vmax(plot_data[gene], vmax_quantile)
        norm = Normalize(vmin=0.0, vmax=vmax, clip=True)
        ax.scatter(x, y, s=background_point_size, c="#d0d0d0", linewidths=0, rasterized=True)
        if positive.any():
            ax.scatter(x[positive], y[positive], s=point_size, c=expr[positive], cmap=cmap, norm=norm, linewidths=0, rasterized=True)
        ax.set_title(panel_display_gene(gene, gene_match), fontsize=10)
        ax.set_xlim(*xlim)
        ax.set_ylim(*ylim)
        ax.set_aspect("equal", adjustable="box")
        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_visible(False)
        cbar = fig.colorbar(ScalarMappable(norm=norm, cmap=cmap), ax=ax, orientation="horizontal", fraction=0.045, pad=0.025)
        cbar.set_ticks([0.0, vmax])
        cbar.set_ticklabels(["0", f"{vmax:.2g}"])
        cbar.ax.tick_params(labelsize=7, length=1.5, pad=1)
        cbar.outline.set_linewidth(0.3)
        rows.append(
            {
                "gene": gene,
                "n_cells_plotted": int(plot_data.shape[0]),
                "n_positive_cells": int(positive.sum()),
                "pct_positive_cells": float(positive.mean() * 100),
                "color_scale_max": vmax,
                "color_scale_max_rule": f"q{vmax_quantile:g}_positive_expression",
            }
        )

    for i in range(len(genes), nrows * ncols):
        axes[i // ncols, i % ncols].axis("off")

    fig.suptitle(title, fontsize=13)
    fig.text(0.5, 0.015, "Expression scale: log1p(CP10K); each gene panel has its own scale", ha="center", fontsize=8)
    fig.tight_layout(rect=(0.01, 0.04, 0.99, 0.93), w_pad=0.4, h_pad=0.7)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=300)
    fig.savefig(output_path.with_suffix(".pdf"))
    plt.close(fig)
    return pd.DataFrame(rows)


def summarize_by_cluster(data: pd.DataFrame, genes: list[str]) -> pd.DataFrame:
    rows = []
    for cluster, group in data.groupby("cluster_id", sort=False):
        cluster_name = group["cluster_number_name"].iloc[0] if "cluster_number_name" in group else str(cluster)
        for gene in genes:
            values = pd.to_numeric(group[gene], errors="coerce").to_numpy(dtype=float)
            finite = values[np.isfinite(values)]
            positive = finite[finite > 0]
            rows.append(
                {
                    "cluster_id": str(cluster),
                    "cluster_number_name": cluster_name,
                    "gene": gene,
                    "n_cells": int(group.shape[0]),
                    "n_positive_cells": int(positive.size),
                    "pct_expressed": float(positive.size / finite.size * 100) if finite.size else np.nan,
                    "mean_expr": float(np.nanmean(finite)) if finite.size else np.nan,
                    "median_expr": float(np.nanmedian(finite)) if finite.size else np.nan,
                    "positive_mean_expr": float(np.nanmean(positive)) if positive.size else np.nan,
                }
            )
    return zscore_by_gene(pd.DataFrame(rows))


def plot_dotplot(summary: pd.DataFrame, genes: list[str], output_path: Path, title: str) -> None:
    clusters = ordered_clusters(summary["cluster_id"])
    x_lookup = {gene: i for i, gene in enumerate(genes)}
    y_lookup = {cluster: i for i, cluster in enumerate(clusters)}
    plot_df = summary.copy()
    fig, ax = plt.subplots(figsize=(max(7.5, 0.62 * len(genes) + 2.0), max(5.5, 0.36 * len(clusters) + 1.6)))
    sizes = np.clip(plot_df["pct_expressed"].to_numpy(dtype=float), 0, 100)
    sc = ax.scatter(
        [x_lookup[g] for g in plot_df["gene"]],
        [y_lookup[str(c)] for c in plot_df["cluster_id"]],
        s=np.maximum(8, sizes * 4.2),
        c=plot_df["mean_expr"].to_numpy(dtype=float),
        cmap=expression_colormap(),
        linewidths=0.25,
        edgecolors="#222222",
    )
    ax.set_xticks(range(len(genes)))
    ax.set_xticklabels(genes, rotation=45, ha="right")
    ax.set_yticks(range(len(clusters)))
    ax.set_yticklabels(clusters)
    ax.invert_yaxis()
    ax.set_xlabel("Gene")
    ax.set_ylabel("DIV90 cluster_id")
    ax.set_title(title)
    cbar = fig.colorbar(sc, ax=ax, fraction=0.028, pad=0.025)
    cbar.set_label("Mean expression")
    for pct, xpos in zip([10, 40, 80], [0.72, 0.80, 0.88]):
        ax.scatter([], [], s=pct * 4.2, c="white", edgecolors="#222222", label=f"{pct}%")
    ax.legend(title="% expressed", loc="upper left", bbox_to_anchor=(1.03, 0.98), frameon=False)
    fig.tight_layout()
    fig.savefig(output_path, dpi=300)
    fig.savefig(output_path.with_suffix(".pdf"))
    plt.close(fig)


def plot_heatmap(summary: pd.DataFrame, genes: list[str], output_path: Path, title: str) -> None:
    clusters = ordered_clusters(summary["cluster_id"])
    matrix = summary.pivot(index="cluster_id", columns="gene", values="mean_expr_z_by_gene").reindex(index=clusters, columns=genes)
    fig, ax = plt.subplots(figsize=(max(7.5, 0.62 * len(genes) + 2.0), max(5.5, 0.36 * len(clusters) + 1.5)))
    im = ax.imshow(matrix.to_numpy(dtype=float), aspect="auto", cmap="RdBu_r", vmin=-2.0, vmax=2.0)
    ax.set_xticks(range(len(genes)))
    ax.set_xticklabels(genes, rotation=45, ha="right")
    ax.set_yticks(range(len(clusters)))
    ax.set_yticklabels(clusters)
    ax.set_xlabel("Gene")
    ax.set_ylabel("DIV90 cluster_id")
    ax.set_title(title)
    cbar = fig.colorbar(im, ax=ax, fraction=0.028, pad=0.025)
    cbar.set_label("Cluster mean z-score per gene")
    fig.tight_layout()
    fig.savefig(output_path, dpi=300)
    fig.savefig(output_path.with_suffix(".pdf"))
    plt.close(fig)


def write_report(
    outdir: Path,
    panel_rows: list[dict[str, str]],
    enriched: pd.DataFrame,
    target_enriched: pd.DataFrame,
    gene_match: pd.DataFrame,
) -> None:
    lines = [
        "# DIV90 Expanded Marker Panel Analysis",
        "",
        "This run used the cached DIV90 AnnData and exported Seurat UMAP coordinates. No Slurm jobs or Seurat RDS loads were used.",
        "",
        "## Panels",
        "",
    ]
    for row in panel_rows:
        lines.extend(
            [
                f"### {row['panel_label']}",
                "",
                f"- Panel ID: `{row['panel_id']}`",
                f"- Genes: `{row['genes']}`",
                f"- Goal: {row['goal']}",
                "",
            ]
        )
    missing = gene_match.loc[~gene_match["matched"].astype(bool), "gene"].tolist()
    lines.extend(["## Gene Matching", ""])
    if missing:
        lines.append("- Missing genes: `" + "`, `".join(missing) + "`")
    else:
        lines.append("- All requested panel genes matched the DIV90 H5AD feature set.")
    cx = gene_match.loc[gene_match["gene"] == "CXCL12R"]
    if not cx.empty:
        lines.append(f"- `CXCL12R` was plotted using matched feature `{cx.iloc[0]['matched_feature']}`.")
    lines.extend(["", "## Cluster-Restricted Marker Signals", ""])
    if enriched.empty:
        lines.append("- No marker passed the enrichment threshold.")
    else:
        for panel_id, group in enriched.groupby("panel_id", sort=False):
            lines.append(f"### {panel_id}")
            top = group.sort_values(["cluster_id", "mean_expr_z_by_gene"], ascending=[True, False])
            for cluster, cluster_group in top.groupby("cluster_id", sort=False):
                genes = ", ".join(
                    f"{r.gene} (z={r.mean_expr_z_by_gene:.2f}, {r.pct_expressed:.1f}%)"
                    for r in cluster_group.itertuples()
                )
                lines.append(f"- Cluster `{cluster}`: {genes}")
            lines.append("")
    lines.extend(["## Candidate Subcluster Flags For Clusters 0, 1, 3, 11", ""])
    if target_enriched.empty:
        lines.append("- No cluster-restricted markers passed the threshold within clusters `0`, `1`, `3`, or `11`.")
    else:
        for cluster, group in target_enriched.groupby("cluster_id", sort=False):
            pieces = []
            for panel_id, panel_group in group.groupby("panel_id", sort=False):
                genes = ", ".join(panel_group.sort_values("mean_expr_z_by_gene", ascending=False)["gene"].tolist())
                pieces.append(f"{panel_id}: {genes}")
            lines.append(f"- Cluster `{cluster}` candidate domains: " + "; ".join(pieces))
    (outdir / "div90_expanded_marker_panel_report.md").write_text("\n".join(lines) + "\n")


def main() -> int:
    add_src_to_path()
    from mge_organoid_python.cross_study_marker_expression import (
        CrossStudyMarkerSpec,
        extract_marker_expression_from_h5ad,
    )

    args = parse_args()
    project_root = Path(args.project_root).expanduser().resolve()
    h5ad_path = Path(args.h5ad).expanduser().resolve() if args.h5ad else project_root / "results/python_anndata/varela_div90.h5ad"
    outdir = (
        Path(args.outdir).expanduser().resolve()
        if args.outdir
        else project_root / "results/div90_expanded_marker_panel_analysis/div90_expanded_marker_panel_analysis_v1"
    )
    plot_dir = outdir / "plots"
    table_dir = outdir / "tables"
    plot_dir.mkdir(parents=True, exist_ok=True)
    table_dir.mkdir(parents=True, exist_ok=True)

    all_genes = unique_preserve_order([gene for panel in PANELS for gene in panel.genes])
    extraction_genes = []
    for gene in all_genes:
        extraction_genes.extend(GENE_ALIASES.get(gene, (gene,)))
    extraction_genes = unique_preserve_order(tuple(extraction_genes))

    spec = CrossStudyMarkerSpec(
        study_id="varela_div90",
        study_label="This Study, DIV 90",
        seurat_path="/nfs/turbo/umms-parent/Manny_test/ventral_sosrs_output/umap_props_output/clustered_day90_with_cluster_names_2.rds",
        h5ad_path=str(h5ad_path),
        sample_col="orig.ident",
        cluster_col="cluster_id",
    )
    expression_raw, match_raw = extract_marker_expression_from_h5ad(
        spec=spec,
        output_path=table_dir / "div90_expanded_marker_expression_raw_alias_columns.tsv.gz",
        project_root=project_root,
        genes=extraction_genes,
        obsm_keys=("X_umap_seurat", "X_umap"),
    )

    gene_match_rows = []
    expression = expression_raw[["cell_id", "study_id", "study_label", "sample", "cluster", "umap_1", "umap_2"]].copy()
    expression = expression.rename(columns={"cluster": "cluster_id"})
    for gene in all_genes:
        candidates = list(GENE_ALIASES.get(gene, (gene,)))
        matched_feature = ""
        matched_col = None
        match_type = "missing"
        for candidate in candidates:
            row = match_raw.loc[match_raw["gene"] == candidate]
            if not row.empty and bool(row.iloc[0]["matched"]):
                matched_feature = str(row.iloc[0]["matched_feature"])
                matched_col = candidate
                match_type = "exact" if candidate == gene else "alias"
                break
        expression[gene] = np.nan if matched_col is None else expression_raw[matched_col]
        gene_match_rows.append(
            {
                "gene": gene,
                "matched_feature": matched_feature,
                "matched": matched_col is not None,
                "match_type": match_type,
                "aliases_considered": ",".join(candidates),
            }
        )
    gene_match = pd.DataFrame(gene_match_rows)

    cluster_names = (
        expression_raw[["cell_id"]]
        .merge(read_obs_columns(h5ad_path, ["cluster_id", "cluster_number_name"]), on="cell_id", how="left")
    )
    expression = expression.merge(cluster_names, on="cell_id", how="left")
    expression["cluster_id"] = expression["cluster_id_y"].fillna(expression["cluster_id_x"]).astype(str)
    expression = expression.drop(columns=[col for col in ["cluster_id_x", "cluster_id_y"] if col in expression.columns])
    expression["cluster_number_name"] = expression["cluster_number_name"].fillna(expression["cluster_id"])

    gene_match.to_csv(table_dir / "div90_expanded_marker_gene_match_table.tsv", sep="\t", index=False)
    expression.to_csv(table_dir / "div90_expanded_marker_expression_table.tsv.gz", sep="\t", index=False)

    panel_rows = []
    all_summaries = []
    all_umap_manifest = []
    for panel in PANELS:
        panel_dir = plot_dir / panel.panel_id
        panel_table_dir = table_dir / panel.panel_id
        panel_dir.mkdir(parents=True, exist_ok=True)
        panel_table_dir.mkdir(parents=True, exist_ok=True)
        genes = list(panel.genes)
        panel_rows.append(
            {
                "panel_id": panel.panel_id,
                "panel_label": panel.label,
                "goal": panel.goal,
                "genes": ",".join(genes),
            }
        )
        umap_manifest = plot_umap_grid(
            expression,
            genes,
            gene_match,
            panel_dir / f"{panel.panel_id}_umap_expression_grid.png",
            title=panel.label,
            ncols=args.ncols,
            max_cells=args.max_cells,
            vmax_quantile=args.vmax_quantile,
            point_size=args.point_size,
            background_point_size=args.background_point_size,
        )
        umap_manifest["panel_id"] = panel.panel_id
        all_umap_manifest.append(umap_manifest)
        summary = summarize_by_cluster(expression, genes)
        summary["panel_id"] = panel.panel_id
        summary["panel_label"] = panel.label
        summary.to_csv(panel_table_dir / f"{panel.panel_id}_cluster_expression_summary.tsv", sep="\t", index=False)
        all_summaries.append(summary)
        plot_dotplot(summary, genes, panel_dir / f"{panel.panel_id}_cluster_dotplot.png", f"{panel.label}: cluster DotPlot")
        plot_heatmap(summary, genes, panel_dir / f"{panel.panel_id}_cluster_average_heatmap.png", f"{panel.label}: cluster-average heatmap")

    panel_table = pd.DataFrame(panel_rows)
    panel_table.to_csv(table_dir / "div90_expanded_marker_panels.tsv", sep="\t", index=False)
    combined_summary = pd.concat(all_summaries, ignore_index=True)
    combined_summary.to_csv(table_dir / "div90_expanded_marker_cluster_expression_summary.tsv", sep="\t", index=False)
    pd.concat(all_umap_manifest, ignore_index=True).to_csv(table_dir / "div90_expanded_marker_umap_manifest.tsv", sep="\t", index=False)

    enriched = combined_summary.loc[
        (combined_summary["mean_expr_z_by_gene"] >= args.enrichment_z) & (combined_summary["pct_expressed"] >= 5)
    ].copy()
    enriched = enriched.sort_values(["panel_id", "cluster_id", "mean_expr_z_by_gene"], ascending=[True, True, False])
    enriched.to_csv(table_dir / "div90_expanded_marker_cluster_restricted_enrichment.tsv", sep="\t", index=False)
    target_enriched = enriched.loc[enriched["cluster_id"].isin(["0", "1", "3", "11"])].copy()
    target_enriched.to_csv(table_dir / "div90_expanded_marker_candidate_subcluster_flags_clusters_0_1_3_11.tsv", sep="\t", index=False)

    write_report(outdir, panel_rows, enriched, target_enriched, gene_match)
    print(f"Wrote analysis to: {outdir}")
    print(f"Report: {outdir / 'div90_expanded_marker_panel_report.md'}")
    return 0


def read_obs_columns(h5ad_path: Path, columns: list[str]) -> pd.DataFrame:
    import anndata as ad

    backed = ad.read_h5ad(h5ad_path, backed="r")
    try:
        obs = backed.obs.copy()
        if "cell_id" not in obs.columns:
            obs = obs.reset_index(names="cell_id")
        else:
            obs = obs.reset_index(drop=True)
        keep = ["cell_id"] + [col for col in columns if col in obs.columns]
        return obs[keep].copy()
    finally:
        backed.file.close()


if __name__ == "__main__":
    raise SystemExit(main())
