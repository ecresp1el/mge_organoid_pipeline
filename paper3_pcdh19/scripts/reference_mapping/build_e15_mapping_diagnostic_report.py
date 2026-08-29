#!/usr/bin/env python3
"""Join both diagnostic mappings, plot unchanged UMAPs, and write the report."""

from __future__ import annotations

import argparse
import gzip
import json
import math
import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd


class ReportContractError(RuntimeError):
    pass


def natural_sample_key(value):
    match = re.search(r"-(\d+)$", str(value))
    return int(match.group(1)) if match else str(value)


class ExistingUmapPlotter:
    def __init__(self, data: pd.DataFrame, figure_dir: Path):
        self.data = data
        self.figure_dir = figure_dir
        self.figure_dir.mkdir(parents=True, exist_ok=True)
        self.samples = sorted(data["sample_id"].unique(), key=natural_sample_key)

    def categorical(self, column: str, stem: str, title: str):
        labels = sorted(self.data[column].fillna("NA").astype(str).unique())
        palette = {label: plt.get_cmap("tab20")(i % 20) for i, label in enumerate(labels)}
        if "Unassigned" in palette:
            palette["Unassigned"] = (0.7, 0.7, 0.7, 1.0)
        fig, axes = plt.subplots(3, 4, figsize=(19, 14), constrained_layout=True)
        for axis, sample in zip(axes.flat, self.samples):
            frame = self.data[self.data["sample_id"] == sample]
            values = frame[column].fillna("NA").astype(str)
            axis.scatter(
                frame["vendor_umap_1"], frame["vendor_umap_2"],
                c=[palette[x] for x in values], s=0.45, linewidths=0,
                alpha=0.75, rasterized=True,
            )
            axis.set_title(f"{sample} (n={len(frame):,})", fontsize=10)
            axis.set_xticks([]); axis.set_yticks([])
            axis.set_aspect("equal", adjustable="datalim")
        handles = [Line2D([0], [0], marker="o", color="none", markerfacecolor=palette[x],
                          markeredgecolor="none", markersize=6, label=x) for x in labels]
        fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, -0.02),
                   ncol=min(6, len(labels)), fontsize=8)
        fig.suptitle(title + "\nOriginal Cell Ranger UMAP; each panel has an independent coordinate system", fontsize=15)
        for suffix in ("png", "pdf"):
            fig.savefig(self.figure_dir / f"{stem}.{suffix}", dpi=300 if suffix == "png" else 200)
        plt.close(fig)

    def continuous(self, column: str, stem: str, title: str):
        values = pd.to_numeric(self.data[column], errors="coerce")
        fig, axes = plt.subplots(3, 4, figsize=(18, 13), constrained_layout=True)
        normalization = Normalize(vmin=0, vmax=1)
        artist = None
        for axis, sample in zip(axes.flat, self.samples):
            mask = self.data["sample_id"] == sample
            artist = axis.scatter(
                self.data.loc[mask, "vendor_umap_1"], self.data.loc[mask, "vendor_umap_2"],
                c=values.loc[mask], cmap="viridis", norm=normalization,
                s=0.45, linewidths=0, alpha=0.8, rasterized=True,
            )
            axis.set_title(f"{sample} (n={int(mask.sum()):,})", fontsize=10)
            axis.set_xticks([]); axis.set_yticks([])
            axis.set_aspect("equal", adjustable="datalim")
        fig.colorbar(artist, ax=axes, shrink=0.75, label=column)
        fig.suptitle(title + "\nOriginal Cell Ranger UMAP; each panel has an independent coordinate system", fontsize=15)
        for suffix in ("png", "pdf"):
            fig.savefig(self.figure_dir / f"{stem}.{suffix}", dpi=300 if suffix == "png" else 200)
        plt.close(fig)


class MappingDiagnosticReport:
    def __init__(self, args):
        self.args = args
        self.output_dir = args.output_dir
        self.table_dir = self.output_dir / "tables"
        self.figure_dir = self.output_dir / "figures"
        self.table_dir.mkdir(parents=True, exist_ok=True)

    @staticmethod
    def find_column(columns, candidates):
        lower = {name.lower(): name for name in columns}
        for candidate in candidates:
            if candidate.lower() in lower:
                return lower[candidate.lower()]
        for name in columns:
            if any(name.lower().endswith(candidate.lower()) for candidate in candidates):
                return name
        return None

    def load(self):
        mind = pd.read_csv(self.args.mind, sep="\t", compression="gzip", low_memory=False)
        mmc = pd.read_csv(self.args.mapmycells, sep="\t", compression="gzip", low_memory=False)
        for frame, label in ((mind, "MIND"), (mmc, "MapMyCells")):
            if frame["cell_id"].duplicated().any():
                raise ReportContractError(f"Duplicate {label} cell IDs")
        if set(mind.cell_id) != set(mmc.cell_id):
            raise ReportContractError("MIND and MapMyCells cell identities differ")
        mmc_only = [name for name in mmc.columns if name.startswith("mmc_")]
        base = [
            "cell_id", "sample_id", "submitted_sample_name", "genotype", "sex", "design_group",
            "existing_cluster", "existing_cluster_numeric", "vendor_umap_1", "vendor_umap_2",
            "MIND_class", "MIND_class_raw", "MIND_class_confidence",
            "MIND_cluster", "MIND_cluster_raw", "MIND_cluster_confidence",
        ]
        score_columns = [name for name in mind.columns if name.startswith("mind_") and name not in base]
        self.data = mind[base + score_columns].merge(
            mmc[["cell_id"] + mmc_only], on="cell_id", validate="one_to_one"
        )
        self.mmc_class = self.find_column(
            self.data.columns, ["mmc_class_name", "class_name"]
        )
        self.mmc_subclass = self.find_column(
            self.data.columns, ["mmc_subclass_name", "subclass_name"]
        )
        self.mmc_cluster = self.find_column(
            self.data.columns, ["mmc_cluster_name", "cluster_name"]
        )
        self.mmc_class_conf = self.find_column(
            self.data.columns,
            ["mmc_class_bootstrapping_probability", "class_bootstrapping_probability"]
        )
        if not all((self.mmc_class, self.mmc_subclass, self.mmc_cluster, self.mmc_class_conf)):
            raise ReportContractError(
                "Could not identify class/subclass/cluster and class confidence in MapMyCells verbose output"
            )
        self.data[self.mmc_class_conf] = pd.to_numeric(
            self.data[self.mmc_class_conf], errors="coerce"
        )

    def tables(self):
        concordance_path = self.table_dir / "per_cell_mapping_concordance.tsv.gz"
        self.data.to_csv(concordance_path, sep="\t", index=False, compression="gzip")
        group_columns = [
            "sample_id", "existing_cluster", "MIND_class", "MIND_cluster",
            self.mmc_class, self.mmc_subclass, self.mmc_cluster,
        ]
        grouped = self.data.groupby(group_columns, dropna=False, observed=True).agg(
            cells=("cell_id", "size"),
            mean_MIND_class_confidence=("MIND_class_confidence", "mean"),
            mean_MIND_cluster_confidence=("MIND_cluster_confidence", "mean"),
            mean_MapMyCells_class_confidence=(self.mmc_class_conf, "mean"),
        ).reset_index().sort_values(["sample_id", "existing_cluster", "cells"], ascending=[True, True, False])
        grouped.to_csv(self.table_dir / "cluster_label_concordance.tsv.gz", sep="\t", index=False,
                       compression="gzip")

        summary = self.data.groupby("sample_id", observed=True).agg(
            cells=("cell_id", "size"),
            MIND_class_assigned=("MIND_class", lambda x: int((x != "Unassigned").sum())),
            MIND_cluster_assigned=("MIND_cluster", lambda x: int((x != "Unassigned").sum())),
            mean_MIND_class_confidence=("MIND_class_confidence", "mean"),
            mean_MIND_cluster_confidence=("MIND_cluster_confidence", "mean"),
            mean_MapMyCells_class_confidence=(self.mmc_class_conf, "mean"),
        ).reset_index()
        summary["MIND_class_assigned_fraction"] = summary.MIND_class_assigned / summary.cells
        summary["MIND_cluster_assigned_fraction"] = summary.MIND_cluster_assigned / summary.cells
        summary.to_csv(self.table_dir / "sample_mapping_summary.tsv", sep="\t", index=False)

        non_neural_pattern = re.compile(
            r"microgl|macroph|immune|vascular|endothel|pericy|oligodend|opc|astro|ependym|fibro",
            re.IGNORECASE,
        )
        label_text = self.data[[self.mmc_class, self.mmc_subclass]].fillna("").astype(str).agg(" | ".join, axis=1)
        non_neural = self.data.loc[label_text.str.contains(non_neural_pattern),
                                   ["sample_id", self.mmc_class, self.mmc_subclass]].copy()
        non_neural.groupby(["sample_id", self.mmc_class, self.mmc_subclass], dropna=False).size().rename(
            "cells"
        ).reset_index().to_csv(self.table_dir / "mapmycells_named_non_neural_populations.tsv",
                               sep="\t", index=False)
        return summary, non_neural

    def plots(self):
        plotter = ExistingUmapPlotter(self.data, self.figure_dir)
        plotter.categorical("MIND_class", "01_mind_class_existing_umap",
                            "Bandler/MIND broad-class transfer (fixed abstention applied)")
        plotter.categorical("MIND_cluster", "02_mind_12_state_existing_umap",
                            "Bandler/MIND 12-state transfer (fixed abstention applied)")
        plotter.continuous("MIND_cluster_confidence", "03_mind_cluster_confidence_existing_umap",
                           "Bandler/MIND winning 12-state prediction score")
        plotter.categorical(self.mmc_class, "04_mapmycells_class_existing_umap",
                            "MapMyCells Allen WMB class comparator")
        plotter.categorical(self.mmc_subclass, "05_mapmycells_subclass_existing_umap",
                            "MapMyCells Allen WMB subclass comparator")
        plotter.continuous(self.mmc_class_conf, "06_mapmycells_class_confidence_existing_umap",
                           "MapMyCells class bootstrapping probability")
        self.mind_hierarchy_heatmap()

    def mind_hierarchy_heatmap(self):
        composition = pd.crosstab(
            self.data["MIND_class"].fillna("NA"),
            self.data["MIND_cluster"].fillna("NA"),
            dropna=False,
        )
        composition.to_csv(self.table_dir / "mind_class_by_state_composition.tsv", sep="\t")
        width = max(14, 0.9 * composition.shape[1] + 5)
        height = max(4, 0.8 * composition.shape[0] + 2)
        fig, axis = plt.subplots(figsize=(width, height), constrained_layout=True)
        image = axis.imshow(np.log1p(composition.values), cmap="Blues", aspect="auto")
        for row in range(composition.shape[0]):
            for column in range(composition.shape[1]):
                count = int(composition.iat[row, column])
                if count:
                    axis.text(column, row, f"{count:,}", ha="center", va="center", fontsize=7,
                              color="white" if np.log1p(count) > 0.65 * np.log1p(composition.values.max()) else "black")
        axis.set_xticks(range(composition.shape[1]), labels=composition.columns, rotation=45, ha="right")
        axis.set_yticks(range(composition.shape[0]), labels=composition.index)
        axis.set_xlabel("Bandler/MIND 12-state call")
        axis.set_ylabel("Bandler/MIND broad-class call")
        axis.set_title("Independently transferred MIND hierarchy: broad class × 12-state composition")
        fig.colorbar(image, ax=axis, label="log(1 + cells)", shrink=0.8)
        for suffix in ("png", "pdf"):
            fig.savefig(self.figure_dir / f"07_mind_hierarchy_composition.{suffix}",
                        dpi=300 if suffix == "png" else 200)
        plt.close(fig)

    def report(self, summary, non_neural):
        threshold = float(self.args.threshold)
        cells = len(self.data)
        mind_class_assigned = int((self.data.MIND_class != "Unassigned").sum())
        mind_cluster_assigned = int((self.data.MIND_cluster != "Unassigned").sum())
        mind_states = sorted(self.data.loc[self.data.MIND_cluster != "Unassigned", "MIND_cluster"].unique())
        broad = sorted(self.data.loc[self.data.MIND_class != "Unassigned", "MIND_class"].unique())
        broad_counts = self.data["MIND_class"].value_counts()
        state_counts = self.data["MIND_cluster"].value_counts()
        dominant_state = state_counts.drop(labels=["Unassigned"], errors="ignore").index[0]
        dominant_state_cells = int(state_counts[dominant_state])
        mitotic_cells = int(broad_counts.get("Mitotic", 0))
        broad_minimum = float(pd.to_numeric(self.data.MIND_class_confidence).min())
        mmc_classes = int(self.data[self.mmc_class].nunique(dropna=True))
        mmc_subclasses = int(self.data[self.mmc_subclass].nunique(dropna=True))
        mmc_clusters = int(self.data[self.mmc_cluster].nunique(dropna=True))
        mmc_class_counts = self.data[self.mmc_class].value_counts()
        top_mmc_class = str(mmc_class_counts.index[0])
        top_mmc_class_cells = int(mmc_class_counts.iloc[0])
        named_non_neural = non_neural.groupby(self.mmc_subclass, dropna=False).size().sort_values(ascending=False)
        named_non_neural_text = ", ".join(
            f"{name}: {int(count):,}" for name, count in named_non_neural.head(12).items()
        )
        min_state_fraction = float(summary.MIND_cluster_assigned_fraction.min())
        max_state_fraction = float(summary.MIND_cluster_assigned_fraction.max())
        report = f"""# E15 MGE diagnostic reference mapping report

## Scope and invariant

This is the first fast diagnostic mapping pass for all **{cells:,}** Paper 3 cells. It did not recluster, integrate, manually annotate, or recompute an embedding. Every plot uses the original Cell Ranger UMAP coordinates and is faceted by sample because each sample's coordinate system was computed independently. `existing_cluster` is likewise sample-qualified (`sample_id:graph_cluster`).

## Diagnostic answer

**Bandler/MIND recovers developmental-MGE similarity, but this pass is not a clean or balanced resolution of the reference hierarchy.** It produces high 12-state assignment coverage in every sample ({min_state_fraction:.1%}–{max_state_fraction:.1%}), but only 11 of the 12 reference states are selected after abstention and the output is strongly concentrated: **{mitotic_cells:,}/{cells:,} ({mitotic_cells/cells:.1%})** are called `Mitotic`, while **{dominant_state_cells:,}/{cells:,} ({dominant_state_cells/cells:.1%})** are assigned `{dominant_state}`. This concentration must be checked against marker expression, sample QC, and technical compatibility before treating the transferred states as annotations.

**MapMyCells provides the intended independent sanity check.** It assigns **{top_mmc_class_cells:,}/{cells:,} ({top_mmc_class_cells/cells:.1%})** cells to adult-taxonomy `{top_mmc_class}`, which supports a broad immature GABAergic program but is not an MGE-specific developmental identity. Separately, it flags **{len(non_neural):,}** explicitly named non-neural cells. The largest such subclasses are {named_non_neural_text}.

## Bandler/MIND transfer

The reference is the **4,481 definitively barcode-resolved CA301 WT E15.5 MGE cells** with the later interactive MIND atlas hierarchy. These labels are not represented as the original paper's unrecovered 21-label taxonomy. The transfer preserves the raw winner, every Seurat prediction-score column, and the maximum score.

A fixed, predeclared threshold of **{threshold:.2f}** maps lower-scoring winners to `Unassigned`; it was not optimized against these results.

- Broad-class assigned: **{mind_class_assigned:,}/{cells:,} ({mind_class_assigned/cells:.1%})**
- Twelve-state assigned: **{mind_cluster_assigned:,}/{cells:,} ({mind_cluster_assigned/cells:.1%})**
- Broad labels observed after abstention: {', '.join(broad) if broad else 'none'}
- State labels observed after abstention ({len(mind_states)}): {', '.join(mind_states) if mind_states else 'none'}

The transfer is a developmental-MGE similarity diagnostic. It cannot discover a genuinely non-MGE class that is absent from CA301; low confidence and `Unassigned` are therefore biologically important outputs.

Because the broad transfer has only two normalized candidate labels, its winning score is mathematically at least 0.50. The observed minimum was **{broad_minimum:.6f}**, so the predeclared 0.50 rule cannot abstain at the broad level and all broad calls remain assigned. This is a limitation of this diagnostic cutoff, not evidence that every cell is truly represented by CA301. The same fixed cutoff remains informative for the 12-state transfer, where 23,199 cells are `Unassigned`.

## MapMyCells comparator

The independent comparator used pinned `cell_type_mapper` 1.7.4 with the official Allen WMB-10X CCN20230722 assets, raw-count normalization, hierarchical (not flat) mapping, and the software defaults of 100 bootstrap iterations and 0.5 marker downsampling. The official workflow does not directly map `CCN20230722_SUPT` during hierarchical traversal; MapMyCells retains that level and backfills it from the lower assignment. Every hierarchy level returned by the software—including that inferred/backfilled supertype—plus identifiers/names, bootstrapping probabilities, correlations, aggregate probabilities, direct-assignment flags, and runner-up fields remains in the per-cell table and extended JSON.

- Allen labels represented: **{mmc_classes} classes, {mmc_subclasses} subclasses, {mmc_clusters} clusters**
- Cells matching an explicit non-neural name pattern at class/subclass: **{len(non_neural):,}**

This adult whole-mouse-brain taxonomy is a sanity-check comparator, not developmental E15 ground truth. Named microglial, vascular, oligodendroglial, astroglial, or related calls are reported as flags for review, not accepted annotations.

## Primary deliverables

- `tables/per_cell_mapping_concordance.tsv.gz`: requested cell-level cluster × MIND hierarchy/scores × full MapMyCells hierarchy/metrics.
- `tables/cluster_label_concordance.tsv.gz`: aggregated label combinations with cell counts and mean confidences.
- `tables/sample_mapping_summary.tsv`: assignment and confidence summary per sample.
- `tables/mapmycells_named_non_neural_populations.tsv`: explicit broad out-of-reference sanity-check calls.
- `figures/01` through `06`: both mappings and confidence on unchanged, sample-faceted UMAPs.
- `figures/07_mind_hierarchy_composition`: broad-class × 12-state composition, retaining discordance between the two independently transferred levels rather than forcing it away.
- `Bandler_MIND/paper3_query_with_mind_and_mapmycells_labels.rds`: one derived Seurat query containing original counts, existing metadata/UMAP, all MIND predictions/scores, and all flattened MapMyCells hierarchy/confidence fields. The temporary MIND-only RDS is removed only after this combined object validates.
- `MapMyCells/mapmycells_extended.json.gz`: untouched semantic content of the official extended output, gzip-compressed after successful parsing.

## Interpretation boundary

This report answers whether expected MGE developmental structure and obvious broad comparator populations are visible. It does not select a final reference, set biological identities, reconstruct the old 21 Bandler labels, or tune mapping parameters.
"""
        (self.output_dir / "E15_MGE_MAPPING_DIAGNOSTIC_REPORT.md").write_text(report, encoding="utf-8")
        (self.output_dir / "SUCCESS.txt").write_text(
            "PASS\n"
            f"cells={cells}\n"
            "reclustered=false\nintegrated=false\nexisting_umap_recomputed=false\n",
            encoding="utf-8",
        )

    def run(self):
        self.load()
        summary, non_neural = self.tables()
        self.plots()
        self.report(summary, non_neural)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--mind", required=True, type=Path)
    parser.add_argument("--mapmycells", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--threshold", required=True, type=float)
    args = parser.parse_args()
    MappingDiagnosticReport(args).run()


if __name__ == "__main__":
    main()
