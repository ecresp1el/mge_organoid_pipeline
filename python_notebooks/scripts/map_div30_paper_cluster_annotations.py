#!/usr/bin/env python3
"""Map DIV30 Seurat clusters to paper cell-state annotations and plot UMAP.

This is a provenance-preserving sidecar step: it does not modify the source
DIV30 Seurat object or the Jia score table. It takes the numeric cluster IDs
already present in the Jia scoring export (`seurat_clusters` and
`RNA_snn_res.0.2`), verifies that those two cluster columns agree cell-by-cell,
and attaches the manual/paper labels that were used for the DIV30 object:

* Radial glia: Seurat clusters 0, 3, 7
* Inhibitory progenitors: Seurat cluster 6
* SST+ cIN: Seurat cluster 1
* PV neuron precursor: Seurat cluster 4
* MGE subpallial neurons: Seurat cluster 2

The UMAP plotted here comes from the converted Seurat reduction table under
`results/python_anndata/varela_div30_2f0j5mwk/umap.tsv`, which traces back to
`results/varela_this_paper/varela_this_paper_seurat.rds`.
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd


DEFAULT_PROJECT_ROOT = Path("/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder")
DEFAULT_RUN_LABEL = "div30_paper_cluster_annotations_v1"

PAPER_CLUSTER_ANNOTATION = {
    "0": "Radial glia",
    "3": "Radial glia",
    "7": "Radial glia",
    "6": "Inhibitory progenitors",
    "1": "SST+ cIN",
    "4": "PV neuron precursor",
    "2": "MGE subpallial neurons",
}

PAPER_DEVELOPMENTAL_STAGE = {
    "Radial glia": "root_progenitor",
    "Inhibitory progenitors": "intermediate_progenitor",
    "PV neuron precursor": "neuronal_precursor",
    "SST+ cIN": "interneuron",
    "MGE subpallial neurons": "subpallial_neuron",
}

PAPER_DEVELOPMENTAL_ORDER = {
    "Radial glia": 0,
    "Inhibitory progenitors": 1,
    "PV neuron precursor": 2,
    "SST+ cIN": 3,
    "MGE subpallial neurons": 3,
}

ANNOTATION_ORDER = [
    "Radial glia",
    "Inhibitory progenitors",
    "PV neuron precursor",
    "SST+ cIN",
    "MGE subpallial neurons",
]

ANNOTATION_COLORS = {
    "Radial glia": "#4c78a8",
    "Inhibitory progenitors": "#f58518",
    "PV neuron precursor": "#54a24b",
    "SST+ cIN": "#b279a2",
    "MGE subpallial neurons": "#e45756",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--project-root",
        type=Path,
        default=Path(os.environ.get("PROJECT_ROOT", DEFAULT_PROJECT_ROOT)),
    )
    parser.add_argument("--run-label", default=os.environ.get("RUN_LABEL", DEFAULT_RUN_LABEL))
    parser.add_argument("--jia-scores", type=Path, default=None)
    parser.add_argument("--umap", type=Path, default=None)
    parser.add_argument("--outdir", type=Path, default=None)
    return parser.parse_args()


def default_paths(project_root: Path, run_label: str) -> dict[str, Path]:
    return {
        "jia_scores": project_root
        / "results"
        / "jia_program_div30_scoring"
        / "jia_program_div30_scoring_v1"
        / "tables"
        / "div30_jia_program_scores_obs.tsv",
        "umap": project_root / "results" / "python_anndata" / "varela_div30_2f0j5mwk" / "umap.tsv",
        "outdir": project_root / "results" / "div30_paper_cluster_annotations" / run_label,
    }


def natural_cluster_key(value: object) -> tuple[int, str]:
    text = str(value)
    return (int(text), text) if text.isdigit() else (10**9, text)


def add_paper_annotations(scores: pd.DataFrame) -> pd.DataFrame:
    required = {"cell_id", "seurat_clusters", "RNA_snn_res.0.2", "orig.ident"}
    missing = sorted(required.difference(scores.columns))
    if missing:
        raise ValueError(f"Jia score table is missing required columns: {missing}")

    out = scores.copy()
    out["seurat_clusters"] = out["seurat_clusters"].astype(str)
    out["RNA_snn_res.0.2"] = out["RNA_snn_res.0.2"].astype(str)
    cluster_match = out["seurat_clusters"] == out["RNA_snn_res.0.2"]
    if not bool(cluster_match.all()):
        mismatch_count = int((~cluster_match).sum())
        raise ValueError(
            "seurat_clusters and RNA_snn_res.0.2 do not match for "
            f"{mismatch_count} cells; refusing to map ambiguous cluster labels."
        )

    out["paper_cluster_annotation"] = out["seurat_clusters"].map(PAPER_CLUSTER_ANNOTATION)
    unmapped = sorted(out.loc[out["paper_cluster_annotation"].isna(), "seurat_clusters"].unique())
    if unmapped:
        raise ValueError(f"Unmapped DIV30 Seurat clusters: {unmapped}")

    out["paper_developmental_stage"] = out["paper_cluster_annotation"].map(PAPER_DEVELOPMENTAL_STAGE)
    out["paper_developmental_order"] = out["paper_cluster_annotation"].map(PAPER_DEVELOPMENTAL_ORDER)
    out["paper_cluster_annotation"] = pd.Categorical(
        out["paper_cluster_annotation"],
        categories=ANNOTATION_ORDER,
        ordered=True,
    )
    return out


def write_mapping_tables(annotated: pd.DataFrame, table_dir: Path) -> list[Path]:
    outputs: list[Path] = []
    mapping = pd.DataFrame(
        [
            {
                "seurat_clusters": cluster,
                "RNA_snn_res.0.2": cluster,
                "paper_cluster_annotation": annotation,
                "paper_developmental_stage": PAPER_DEVELOPMENTAL_STAGE[annotation],
                "paper_developmental_order": PAPER_DEVELOPMENTAL_ORDER[annotation],
            }
            for cluster, annotation in sorted(PAPER_CLUSTER_ANNOTATION.items(), key=lambda item: natural_cluster_key(item[0]))
        ]
    )
    path = table_dir / "div30_paper_cluster_annotation_mapping.tsv"
    mapping.to_csv(path, sep="\t", index=False)
    outputs.append(path)

    count_df = (
        annotated.groupby(
            ["seurat_clusters", "RNA_snn_res.0.2", "paper_cluster_annotation"],
            observed=True,
        )
        .size()
        .reset_index(name="n_cells")
    )
    count_df["paper_developmental_stage"] = count_df["paper_cluster_annotation"].map(PAPER_DEVELOPMENTAL_STAGE)
    count_df["paper_developmental_order"] = count_df["paper_cluster_annotation"].map(PAPER_DEVELOPMENTAL_ORDER)
    count_df["_cluster_sort"] = count_df["seurat_clusters"].map(lambda value: natural_cluster_key(value)[0])
    count_df = count_df.sort_values(["paper_developmental_order", "_cluster_sort"]).drop(columns="_cluster_sort")
    count_path = table_dir / "div30_paper_cluster_annotation_counts.tsv"
    count_df.to_csv(count_path, sep="\t", index=False)
    outputs.append(count_path)

    sample_counts = pd.crosstab(
        annotated["orig.ident"],
        annotated["paper_cluster_annotation"],
    ).reindex(columns=ANNOTATION_ORDER, fill_value=0)
    sample_counts.index.name = "orig.ident"
    sample_path = table_dir / "div30_paper_cluster_annotation_counts_by_sample.tsv"
    sample_counts.to_csv(sample_path, sep="\t")
    outputs.append(sample_path)

    cell_cols = [
        "obs_name",
        "cell_id",
        "orig.ident",
        "seurat_clusters",
        "RNA_snn_res.0.2",
        "paper_cluster_annotation",
        "paper_developmental_stage",
        "paper_developmental_order",
        "jia_score_IPC",
        "jia_score_RGC1",
        "jia_score_RGC2",
    ]
    available = [col for col in cell_cols if col in annotated.columns]
    cell_path = table_dir / "div30_jia_scores_with_paper_cluster_annotations.tsv.gz"
    annotated[available].to_csv(cell_path, sep="\t", index=False, compression="gzip")
    outputs.append(cell_path)

    return outputs


def plot_umap(annotated: pd.DataFrame, umap_path: Path, plot_dir: Path) -> Path:
    if not umap_path.exists():
        raise FileNotFoundError(f"Missing UMAP table: {umap_path}")
    umap = pd.read_csv(umap_path, sep="\t")
    required = {"cell_id", "UMAP_1", "UMAP_2"}
    missing = sorted(required.difference(umap.columns))
    if missing:
        raise ValueError(f"UMAP table is missing required columns: {missing}")

    plot_df = umap.merge(
        annotated[["cell_id", "seurat_clusters", "paper_cluster_annotation"]],
        on="cell_id",
        how="inner",
        validate="one_to_one",
    )
    if len(plot_df) != len(annotated):
        raise ValueError(
            f"UMAP join retained {len(plot_df)} cells, expected {len(annotated)}. "
            "Check cell IDs between Jia score and UMAP tables."
        )

    fig, axes = plt.subplots(1, 2, figsize=(12, 5.4), constrained_layout=True)
    for ax in axes:
        ax.set_aspect("equal", adjustable="box")
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_xlabel("")
        ax.set_ylabel("")
        for spine in ax.spines.values():
            spine.set_visible(False)

    for label in ANNOTATION_ORDER:
        sub = plot_df[plot_df["paper_cluster_annotation"].astype(str) == label]
        axes[0].scatter(
            sub["UMAP_1"],
            sub["UMAP_2"],
            s=0.6,
            linewidths=0,
            alpha=0.85,
            color=ANNOTATION_COLORS[label],
            label=f"{label} ({len(sub):,})",
            rasterized=True,
        )
    axes[0].set_title("DIV30 paper cluster annotations", fontsize=12)
    axes[0].legend(loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=False, markerscale=5, fontsize=8)

    cluster_values = sorted(plot_df["seurat_clusters"].astype(str).unique(), key=natural_cluster_key)
    palette = plt.get_cmap("tab10")
    for idx, cluster in enumerate(cluster_values):
        sub = plot_df[plot_df["seurat_clusters"].astype(str) == cluster]
        axes[1].scatter(
            sub["UMAP_1"],
            sub["UMAP_2"],
            s=0.6,
            linewidths=0,
            alpha=0.85,
            color=palette(idx % 10),
            label=f"{cluster} ({len(sub):,})",
            rasterized=True,
        )
    axes[1].set_title("DIV30 Seurat clusters", fontsize=12)
    axes[1].legend(loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=False, markerscale=5, fontsize=8)

    out = plot_dir / "div30_umap_paper_cluster_annotations_and_seurat_clusters.png"
    fig.savefig(out, dpi=300)
    plt.close(fig)
    return out


def main() -> None:
    args = parse_args()
    project_root = args.project_root.expanduser().resolve()
    paths = default_paths(project_root, args.run_label)
    jia_scores_path = (args.jia_scores or paths["jia_scores"]).expanduser().resolve()
    umap_path = (args.umap or paths["umap"]).expanduser().resolve()
    outdir = (args.outdir or paths["outdir"]).expanduser().resolve()
    table_dir = outdir / "tables"
    plot_dir = outdir / "plots"
    table_dir.mkdir(parents=True, exist_ok=True)
    plot_dir.mkdir(parents=True, exist_ok=True)

    print(f"[Div30PaperClusterAnnotations] project_root={project_root}", flush=True)
    print(f"[Div30PaperClusterAnnotations] jia_scores={jia_scores_path}", flush=True)
    print(f"[Div30PaperClusterAnnotations] umap={umap_path}", flush=True)
    print(f"[Div30PaperClusterAnnotations] outdir={outdir}", flush=True)

    if not jia_scores_path.exists():
        raise FileNotFoundError(f"Missing Jia score table: {jia_scores_path}")
    scores = pd.read_csv(jia_scores_path, sep="\t")
    annotated = add_paper_annotations(scores)

    outputs = write_mapping_tables(annotated, table_dir)
    plot_path = plot_umap(annotated, umap_path, plot_dir)
    outputs.append(plot_path)

    manifest = pd.DataFrame(
        [
            {
                "path": str(path),
                "kind": "plot" if path.suffix.lower() == ".png" else "table",
                "bytes": path.stat().st_size,
            }
            for path in outputs
        ]
    )
    manifest_path = table_dir / "div30_paper_cluster_annotation_manifest.tsv"
    manifest.to_csv(manifest_path, sep="\t", index=False)

    print("[Div30PaperClusterAnnotations] counts by paper annotation", flush=True)
    counts = annotated["paper_cluster_annotation"].value_counts().reindex(ANNOTATION_ORDER, fill_value=0)
    print(counts.to_string(), flush=True)
    print(f"[Div30PaperClusterAnnotations] complete: {outdir}", flush=True)


if __name__ == "__main__":
    main()
