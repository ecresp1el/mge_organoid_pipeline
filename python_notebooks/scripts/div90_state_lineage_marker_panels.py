#!/usr/bin/env python3
"""DIV90 marker panels for developmental state and lineage diversification."""

from __future__ import annotations

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
H5AD_RELATIVE = "results/python_anndata/varela_div90.h5ad"


@dataclass(frozen=True)
class Panel:
    panel_id: str
    label: str
    genes: tuple[str, ...]
    goal: str


PANELS = [
    Panel(
        panel_id="panel1_developmental_migration",
        label="Panel 1: developmental state / migration",
        genes=("DCX", "SOX11", "STMN2", "DCLK1", "CXCR4", "ROBO1", "ROBO2", "ACKR3"),
        goal=(
            "Determine whether clusters 1, 3, and 11 contain distinct migratory versus "
            "post-migratory interneuron populations."
        ),
    ),
    Panel(
        panel_id="panel2_mge_interneuron_maturation",
        label="Panel 2: MGE interneuron maturation",
        genes=("LHX6", "SOX6", "SATB1", "MAF", "MAFB", "ERBB4", "GAD1", "GAD2"),
        goal="Determine whether cortical-fated interneurons contain separable maturation programs.",
    ),
    Panel(
        panel_id="panel3_pv_lineage_maturation",
        label="Panel 3: PV-lineage maturation",
        genes=("PVALB", "PPARGC1A", "KCNC1", "KCNC2", "GPR149", "ERBB4", "SOX6"),
        goal="Determine whether putative PV precursors separate into distinct maturation states despite low PVALB.",
    ),
    Panel(
        panel_id="panel4_striatal_pallidal_specification",
        label="Panel 4: striatal versus pallidal specification",
        genes=("NKX2-1", "LHX8", "GBX2", "ISL1", "EBF1", "FOXP1", "BCL11B", "TAC1"),
        goal="Test whether the MGE striatal/GP-fated compartment contains unresolved specification programs.",
    ),
]

TARGET_CLUSTERS = ("0", "1", "3", "11")
DISPLAY_LABELS = {"ACKR3": "ACKR3/CXCR7"}


def repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def add_src_to_path() -> None:
    src = repo_root() / "python_notebooks" / "src"
    if str(src) not in sys.path:
        sys.path.insert(0, str(src))


def unique(values: list[str] | tuple[str, ...]) -> list[str]:
    seen: set[str] = set()
    out: list[str] = []
    for value in values:
        if value in seen:
            continue
        seen.add(value)
        out.append(value)
    return out


def natural_cluster_order(values: pd.Series) -> list[str]:
    numeric = pd.to_numeric(values.astype(str), errors="coerce")
    ordered = pd.DataFrame({"cluster": values.astype(str), "numeric": numeric}).drop_duplicates()
    ordered = ordered.sort_values(["numeric", "cluster"], na_position="last")
    return ordered["cluster"].tolist()


def positive_vmax(values: pd.Series, quantile: float = 0.99) -> float:
    arr = pd.to_numeric(values, errors="coerce").to_numpy(dtype=float)
    positive = arr[np.isfinite(arr) & (arr > 0)]
    if positive.size == 0:
        return 1.0
    vmax = float(np.quantile(positive, quantile))
    if not np.isfinite(vmax) or vmax <= 0:
        vmax = float(np.max(positive))
    return vmax if vmax > 0 else 1.0


def plot_umap_panel(data: pd.DataFrame, genes: list[str], title: str, path: Path) -> None:
    ncols = min(4, len(genes))
    nrows = int(np.ceil(len(genes) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(3.1 * ncols, 3.05 * nrows), squeeze=False)
    cmap = LinearSegmentedColormap.from_list("grey_blue", ["#d9d9d9", "#2e49ff", "#08008f"])
    x = data["umap_1"].to_numpy(dtype=float)
    y = data["umap_2"].to_numpy(dtype=float)
    finite = np.isfinite(x) & np.isfinite(y)
    data = data.loc[finite].copy()
    x = x[finite]
    y = y[finite]
    xpad = (x.max() - x.min()) * 0.035
    ypad = (y.max() - y.min()) * 0.035
    for idx, gene in enumerate(genes):
        ax = axes[idx // ncols][idx % ncols]
        expr = pd.to_numeric(data[gene], errors="coerce").to_numpy(dtype=float)
        positive = np.isfinite(expr) & (expr > 0)
        vmax = positive_vmax(data[gene])
        norm = Normalize(vmin=0, vmax=vmax, clip=True)
        ax.scatter(x, y, s=0.12, c="#cfcfcf", linewidths=0, rasterized=True)
        if positive.any():
            ax.scatter(x[positive], y[positive], s=0.42, c=expr[positive], cmap=cmap, norm=norm, linewidths=0, rasterized=True)
        ax.set_title(DISPLAY_LABELS.get(gene, gene), fontsize=10)
        ax.set_xlim(x.min() - xpad, x.max() + xpad)
        ax.set_ylim(y.min() - ypad, y.max() + ypad)
        ax.set_aspect("equal", adjustable="box")
        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_visible(False)
        cbar = fig.colorbar(ScalarMappable(norm=norm, cmap=cmap), ax=ax, orientation="horizontal", fraction=0.045, pad=0.02)
        cbar.set_ticks([0, vmax])
        cbar.set_ticklabels(["0", f"{vmax:.2g}"])
        cbar.ax.tick_params(labelsize=6, length=1)
        cbar.outline.set_linewidth(0.25)
    for idx in range(len(genes), nrows * ncols):
        axes[idx // ncols][idx % ncols].axis("off")
    fig.suptitle(title, fontsize=13)
    fig.text(0.5, 0.012, "log1p(CP10K); each marker uses its own positive-expression q99 color scale", ha="center", fontsize=8)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout(rect=(0.01, 0.035, 0.99, 0.93), w_pad=0.35, h_pad=0.55)
    fig.savefig(path, dpi=300)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def cluster_summary(data: pd.DataFrame, genes: list[str], cluster_order: list[str]) -> pd.DataFrame:
    rows = []
    for cluster in cluster_order:
        sub = data.loc[data["cluster"].astype(str) == cluster]
        for gene in genes:
            values = pd.to_numeric(sub[gene], errors="coerce").to_numpy(dtype=float)
            finite = values[np.isfinite(values)]
            positive = finite[finite > 0]
            rows.append(
                {
                    "cluster": cluster,
                    "gene": gene,
                    "gene_label": DISPLAY_LABELS.get(gene, gene),
                    "n_cells": int(sub.shape[0]),
                    "n_positive": int(positive.size),
                    "pct_expressed": float(positive.size / finite.size * 100) if finite.size else np.nan,
                    "mean_expr": float(np.mean(finite)) if finite.size else np.nan,
                    "median_expr": float(np.median(finite)) if finite.size else np.nan,
                    "positive_mean_expr": float(np.mean(positive)) if positive.size else np.nan,
                    "q90_expr": float(np.quantile(finite, 0.9)) if finite.size else np.nan,
                }
            )
    return pd.DataFrame(rows)


def plot_dotplot(summary: pd.DataFrame, genes: list[str], cluster_order: list[str], title: str, path: Path) -> None:
    matrix = summary.pivot(index="cluster", columns="gene", values="mean_expr").reindex(index=cluster_order, columns=genes)
    pct = summary.pivot(index="cluster", columns="gene", values="pct_expressed").reindex(index=cluster_order, columns=genes)
    fig, ax = plt.subplots(figsize=(max(6.5, 0.55 * len(genes) + 2.0), max(4.8, 0.28 * len(cluster_order) + 1.2)))
    xs, ys, sizes, colors = [], [], [], []
    for yidx, cluster in enumerate(cluster_order):
        for xidx, gene in enumerate(genes):
            xs.append(xidx)
            ys.append(yidx)
            sizes.append(max(float(pct.loc[cluster, gene]), 0.0) * 3.2)
            colors.append(float(matrix.loc[cluster, gene]))
    sc = ax.scatter(xs, ys, s=sizes, c=colors, cmap="viridis", edgecolors="#303030", linewidths=0.15)
    ax.set_xticks(range(len(genes)))
    ax.set_xticklabels([DISPLAY_LABELS.get(g, g) for g in genes], rotation=45, ha="right")
    ax.set_yticks(range(len(cluster_order)))
    ax.set_yticklabels(cluster_order)
    ax.invert_yaxis()
    ax.set_xlabel("Marker")
    ax.set_ylabel("DIV90 cluster_id")
    ax.set_title(title)
    cbar = fig.colorbar(sc, ax=ax, fraction=0.035, pad=0.02)
    cbar.set_label("Mean expression")
    for size in (10, 30, 60):
        ax.scatter([], [], s=size * 3.2, c="#d0d0d0", edgecolors="#303030", linewidths=0.15, label=f"{size}%")
    ax.legend(title="% expressed", loc="center left", bbox_to_anchor=(1.12, 0.5), frameon=False)
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=300)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def plot_heatmap(summary: pd.DataFrame, genes: list[str], cluster_order: list[str], title: str, path: Path) -> None:
    means = summary.pivot(index="cluster", columns="gene", values="mean_expr").reindex(index=cluster_order, columns=genes)
    z = means.copy()
    for gene in genes:
        vals = means[gene].to_numpy(dtype=float)
        sd = np.nanstd(vals)
        z[gene] = 0 if sd == 0 or not np.isfinite(sd) else (vals - np.nanmean(vals)) / sd
    fig, ax = plt.subplots(figsize=(max(6.5, 0.55 * len(genes) + 2.0), max(4.8, 0.28 * len(cluster_order) + 1.0)))
    im = ax.imshow(z.to_numpy(dtype=float), aspect="auto", cmap="RdBu_r", vmin=-2, vmax=2)
    ax.set_xticks(range(len(genes)))
    ax.set_xticklabels([DISPLAY_LABELS.get(g, g) for g in genes], rotation=45, ha="right")
    ax.set_yticks(range(len(cluster_order)))
    ax.set_yticklabels(cluster_order)
    ax.set_xlabel("Marker")
    ax.set_ylabel("DIV90 cluster_id")
    ax.set_title(title)
    cbar = fig.colorbar(im, ax=ax, fraction=0.035, pad=0.02)
    cbar.set_label("Cluster mean z-score per marker")
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=300)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def enrichment_table(summary: pd.DataFrame, genes: list[str]) -> pd.DataFrame:
    rows = []
    for gene in genes:
        sub = summary.loc[summary["gene"] == gene].sort_values("mean_expr", ascending=False).reset_index(drop=True)
        vals = sub["mean_expr"].to_numpy(dtype=float)
        sd = np.nanstd(vals)
        top = sub.iloc[0]
        second = sub.iloc[1] if sub.shape[0] > 1 else top
        top_z = np.nan if sd == 0 else (float(top["mean_expr"]) - float(np.nanmean(vals))) / sd
        rows.append(
            {
                "gene": gene,
                "gene_label": DISPLAY_LABELS.get(gene, gene),
                "top_cluster": str(top["cluster"]),
                "top_cluster_mean_expr": float(top["mean_expr"]),
                "top_cluster_pct_expressed": float(top["pct_expressed"]),
                "second_cluster": str(second["cluster"]),
                "second_cluster_mean_expr": float(second["mean_expr"]),
                "top_minus_second_mean_expr": float(top["mean_expr"] - second["mean_expr"]),
                "top_to_second_mean_ratio": float((top["mean_expr"] + 1e-9) / (second["mean_expr"] + 1e-9)),
                "top_cluster_z": top_z,
                "cluster_restricted_enrichment": bool(top_z >= 1.25 and (top["mean_expr"] - second["mean_expr"]) > 0.05),
            }
        )
    return pd.DataFrame(rows)


def candidate_domains(enrichment: pd.DataFrame, summary: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for cluster in TARGET_CLUSTERS:
        enriched = enrichment.loc[
            (enrichment["top_cluster"].astype(str) == cluster)
            & (enrichment["cluster_restricted_enrichment"])
        ]
        for _, row in enriched.iterrows():
            s = summary.loc[(summary["cluster"].astype(str) == cluster) & (summary["gene"] == row["gene"])].iloc[0]
            rows.append(
                {
                    "candidate_cluster": cluster,
                    "supporting_gene": row["gene"],
                    "supporting_gene_label": row["gene_label"],
                    "mean_expr": s["mean_expr"],
                    "pct_expressed": s["pct_expressed"],
                    "top_minus_second_mean_expr": row["top_minus_second_mean_expr"],
                    "interpretation": "cluster-level marker-restricted enrichment; inspect UMAP panel for within-cluster domain shape",
                }
            )
    return pd.DataFrame(rows)


def write_report(outdir: Path, panel_rows: list[dict[str, object]], candidates: pd.DataFrame, gene_matches: pd.DataFrame) -> None:
    lines = [
        "# DIV90 State And Lineage Marker Panel Analysis",
        "",
        "This local run used the cached DIV90 H5AD and exported Seurat UMAP coordinates. No Slurm jobs were submitted.",
        "",
        "## Gene Matching",
        "",
    ]
    missing = gene_matches.loc[~gene_matches["matched"].astype(bool), "gene"].tolist()
    lines.append("- Missing requested markers: " + (", ".join(missing) if missing else "none"))
    ackr3 = gene_matches.loc[gene_matches["gene"] == "ACKR3"]
    if not ackr3.empty:
        lines.append("- `CXCL12R/CXCR7` was represented as `ACKR3`; aliases considered were `ACKR3,CXCR7,CXCL12R`.")
    lines.extend(["", "## Panel Outputs", ""])
    for row in panel_rows:
        lines.extend(
            [
                f"### {row['label']}",
                "",
                f"- Goal: {row['goal']}",
                f"- UMAP: `{row['umap_png']}`",
                f"- DotPlot: `{row['dotplot_png']}`",
                f"- Heatmap: `{row['heatmap_png']}`",
                f"- Cluster summary: `{row['summary_tsv']}`",
                f"- Enrichment summary: `{row['enrichment_tsv']}`",
                "",
            ]
        )
    lines.extend(["## Candidate Cluster-Level Domains", ""])
    if candidates.empty:
        lines.append("No requested marker met the cluster-restricted enrichment rule within target clusters 0, 1, 3, or 11.")
    else:
        for cluster, sub in candidates.groupby("candidate_cluster", sort=False):
            genes = ", ".join(sub["supporting_gene_label"].astype(str).tolist())
            lines.append(f"- Cluster `{cluster}`: candidate marker-restricted domain supported by {genes}.")
    lines.extend(["", "Rule used: top cluster z-score >= 1.25 and top-minus-second mean expression > 0.05. Treat this as a screening flag, not a replacement for subclustering."])
    (outdir / "DIV90_state_lineage_marker_panel_report.md").write_text("\n".join(lines) + "\n")


def main() -> int:
    add_src_to_path()
    from mge_organoid_python import cross_study_marker_expression as csm

    csm.GENE_ALIASES = dict(csm.GENE_ALIASES)
    csm.GENE_ALIASES["ACKR3"] = ["ACKR3", "CXCR7", "CXCL12R"]

    project_root = Path(os.environ.get("PROJECT_ROOT", PROJECT_ROOT_DEFAULT))
    h5ad = project_root / H5AD_RELATIVE
    outdir = project_root / "results/div90_state_lineage_marker_panels/div90_state_lineage_marker_panels_v1"
    table_dir = outdir / "tables"
    plot_dir = outdir / "plots"
    table_dir.mkdir(parents=True, exist_ok=True)
    plot_dir.mkdir(parents=True, exist_ok=True)

    all_genes = unique([gene for panel in PANELS for gene in panel.genes])
    spec = csm.CrossStudyMarkerSpec(
        study_id="varela_div90",
        study_label="This Study, DIV 90",
        seurat_path="/nfs/turbo/umms-parent/Manny_test/ventral_sosrs_output/umap_props_output/clustered_day90_with_cluster_names_2.rds",
        h5ad_path=str(h5ad),
        sample_col="orig.ident",
        cluster_col="cluster_id",
    )
    expression_path = table_dir / "div90_state_lineage_all_panel_expression.tsv.gz"
    expression, gene_matches = csm.extract_marker_expression_from_h5ad(
        spec=spec,
        output_path=expression_path,
        project_root=project_root,
        genes=all_genes,
        obsm_keys=("X_umap_seurat", "X_umap"),
    )
    gene_matches.to_csv(table_dir / "div90_state_lineage_gene_match_table.tsv", sep="\t", index=False)
    expression["cluster"] = expression["cluster"].astype(str)
    cluster_order = natural_cluster_order(expression["cluster"])

    all_summaries = []
    all_enrichment = []
    all_candidates = []
    panel_rows = []
    for panel in PANELS:
        genes = list(panel.genes)
        panel_plot_dir = plot_dir / panel.panel_id
        panel_table_dir = table_dir / panel.panel_id
        panel_plot_dir.mkdir(parents=True, exist_ok=True)
        panel_table_dir.mkdir(parents=True, exist_ok=True)
        umap_path = panel_plot_dir / f"{panel.panel_id}_umap_grid.png"
        dotplot_path = panel_plot_dir / f"{panel.panel_id}_cluster_dotplot.png"
        heatmap_path = panel_plot_dir / f"{panel.panel_id}_cluster_average_heatmap.png"

        plot_umap_panel(expression, genes, panel.label, umap_path)
        summary = cluster_summary(expression, genes, cluster_order)
        summary.insert(0, "panel_id", panel.panel_id)
        summary.insert(1, "panel_label", panel.label)
        summary_path = panel_table_dir / f"{panel.panel_id}_cluster_expression_summary.tsv"
        summary.to_csv(summary_path, sep="\t", index=False)

        enrich = enrichment_table(summary, genes)
        enrich.insert(0, "panel_id", panel.panel_id)
        enrich.insert(1, "panel_label", panel.label)
        enrich_path = panel_table_dir / f"{panel.panel_id}_cluster_restricted_enrichment.tsv"
        enrich.to_csv(enrich_path, sep="\t", index=False)

        plot_dotplot(summary, genes, cluster_order, panel.label + " DotPlot", dotplot_path)
        plot_heatmap(summary, genes, cluster_order, panel.label + " cluster-average heatmap", heatmap_path)

        candidates = candidate_domains(enrich, summary)
        if not candidates.empty:
            candidates.insert(0, "panel_id", panel.panel_id)
            candidates.insert(1, "panel_label", panel.label)
            all_candidates.append(candidates)
        all_summaries.append(summary)
        all_enrichment.append(enrich)
        panel_rows.append(
            {
                "panel_id": panel.panel_id,
                "label": panel.label,
                "goal": panel.goal,
                "genes": ",".join(genes),
                "umap_png": str(umap_path),
                "dotplot_png": str(dotplot_path),
                "heatmap_png": str(heatmap_path),
                "summary_tsv": str(summary_path),
                "enrichment_tsv": str(enrich_path),
            }
        )

    combined_summary = pd.concat(all_summaries, ignore_index=True)
    combined_enrichment = pd.concat(all_enrichment, ignore_index=True)
    combined_candidates = pd.concat(all_candidates, ignore_index=True) if all_candidates else pd.DataFrame()
    combined_summary.to_csv(table_dir / "div90_state_lineage_combined_cluster_expression_summary.tsv", sep="\t", index=False)
    combined_enrichment.to_csv(table_dir / "div90_state_lineage_combined_cluster_restricted_enrichment.tsv", sep="\t", index=False)
    combined_candidates.to_csv(table_dir / "div90_state_lineage_candidate_domains_clusters_0_1_3_11.tsv", sep="\t", index=False)
    pd.DataFrame(panel_rows).to_csv(table_dir / "div90_state_lineage_panel_manifest.tsv", sep="\t", index=False)
    write_report(outdir, panel_rows, combined_candidates, gene_matches)

    print(f"Wrote report: {outdir / 'DIV90_state_lineage_marker_panel_report.md'}")
    print(f"Wrote panel manifest: {table_dir / 'div90_state_lineage_panel_manifest.tsv'}")
    print(f"Wrote candidate table: {table_dir / 'div90_state_lineage_candidate_domains_clusters_0_1_3_11.tsv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
