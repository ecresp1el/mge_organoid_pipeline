#!/usr/bin/env python3
"""Publish Bandler sample, annotation, artifact-scope, and visual evidence."""

from __future__ import annotations

import argparse
import csv
import os
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Mapping, Sequence

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def atomic_tsv(path: Path, columns: Sequence[str], rows: Iterable[Mapping[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(prefix=f".{path.name}.", dir=path.parent, text=True)
    temporary = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(columns), delimiter="\t", lineterminator="\n")
            writer.writeheader()
            for row in rows:
                writer.writerow({column: row.get(column, "") for column in columns})
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    finally:
        if temporary.exists():
            temporary.unlink()


def atomic_text(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(prefix=f".{path.name}.", dir=path.parent, text=True)
    temporary = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
            handle.write(content)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    finally:
        if temporary.exists():
            temporary.unlink()


def clean(value: object) -> str:
    if pd.isna(value):
        return "NA"
    if isinstance(value, float) and value.is_integer():
        return str(int(value))
    return str(value).strip()


def read_marker_table(path: Path) -> pd.DataFrame:
    frame = pd.read_excel(path, header=1)
    frame.columns = [str(column).strip() for column in frame.columns]
    required = {"cluster", "gene"}
    if not required.issubset(frame.columns):
        raise RuntimeError(f"Marker table lacks {sorted(required)}: {path}")
    frame = frame[frame["cluster"].notna() & frame["gene"].notna()].copy()
    unnamed = [column for column in frame.columns if column.startswith("Unnamed:")]
    return frame.drop(columns=unnamed, errors="ignore")


@dataclass(frozen=True)
class PublishedArtifactPaths:
    sample_ledger: Path
    broad_markers: Path
    refined_markers: Path
    embryonic_markers: Path

    @classmethod
    def from_source_root(cls, source_root: Path) -> "PublishedArtifactPaths":
        root = source_root / "Bandler2022" / "metadata_sources"
        return cls(
            root / "41586_2021_4237_MOESM2_ESM.xlsx",
            root / "41586_2021_4237_MOESM3_ESM.xlsx",
            root / "41586_2021_4237_MOESM4_ESM.xlsx",
            root / "41586_2021_4237_MOESM5_ESM.xlsx",
        )

    def validate(self) -> None:
        for path in (self.sample_ledger, self.broad_markers, self.refined_markers, self.embryonic_markers):
            if not path.is_file():
                raise FileNotFoundError(path)


class PublishedSampleLedger:
    """Extract the seven datasets used in the published embryonic integration."""

    GEO = {
        "CA298": "GSM5684879", "CA299": "GSM5684878", "CA300": "GSM5684877",
        "CA301": "GSM5684876", "CA302": "GSM5684875", "CA303": "GSM5684874",
        "MUC28072": "GSM5684900",
    }

    def __init__(self, workbook: Path):
        self.workbook = workbook

    def rows(self) -> list[dict[str, object]]:
        transcriptomes = pd.read_excel(self.workbook, sheet_name="Transcriptome Datasets", header=0)
        processing = pd.read_excel(self.workbook, sheet_name="Datasets processing", header=2)
        transcriptomes.columns = [str(column).strip() for column in transcriptomes.columns]
        processing.columns = [str(column).strip() for column in processing.columns]
        sample_column = "Sample ID Seq"
        if sample_column not in transcriptomes or "Dataset" not in processing:
            raise RuntimeError("Unexpected Supplementary Data 1 schema")
        result = []
        for _, row in transcriptomes.iterrows():
            sample_seq = clean(row.get(sample_column))
            sample = next((token for token in self.GEO if sample_seq.startswith(token)), None)
            if sample is None:
                continue
            matched = processing[processing["Dataset"].astype(str).str.startswith(sample_seq, na=False)]
            if matched.empty:
                matched = processing[processing["Dataset"].astype(str).str.startswith(sample, na=False)]
            processed = matched.iloc[0] if not matched.empty else pd.Series(dtype=object)
            result.append({
                "sample_id": sample,
                "sample_id_seq": sample_seq,
                "geo_accession": self.GEO[sample],
                "experiment": clean(row.get("Experiment")),
                "stage_injected": clean(row.get("Stage Injected")),
                "age_of_collection": clean(row.get("Age of collection")),
                "anatomical_region": clean(row.get("Anatomical region collected")),
                "pooled_brains": clean(row.get("Pooled brains")),
                "facs_selection": clean(row.get("FACSorted on")),
                "tenx_version": clean(row.get("10x Version")),
                "sequencing_info": clean(row.get("Sequencing info")),
                "library_type": clean(row.get("Barcode or transcriptome library")),
                "cells_before_filter": clean(processed.get("Estimated Number of Cells before filter")),
                "doublet_finder": clean(processed.get("DoubletFinder")),
                "nfeature_min": clean(processed.get("nFeature_RNA min")),
                "nfeature_max": clean(processed.get("nFeature_RNA max")),
                "ncount_max": clean(processed.get("nCount_RNA")),
                "percent_mito_max": clean(processed.get("percent.mito")),
                "cells_after_filter_published": clean(processed.get("Estimated Number of Cells after filter")),
                "source": "Bandler Supplementary Data 1",
            })
        order = {sample: index for index, sample in enumerate(self.GEO)}
        return sorted(result, key=lambda row: order[str(row["sample_id"])])


class PublishedAnnotationTables:
    """Preserve publisher marker rows and publish one inventory per hierarchy level."""

    def __init__(self, paths: PublishedArtifactPaths, metadata_dir: Path):
        self.paths = paths
        self.metadata_dir = metadata_dir

    @staticmethod
    def marker_rows(frame: pd.DataFrame) -> list[dict[str, object]]:
        return [{column: clean(value) for column, value in row.items()} for row in frame.to_dict("records")]

    @staticmethod
    def inventory(frame: pd.DataFrame, level: str) -> list[dict[str, object]]:
        rows = []
        for cluster, group in frame.groupby("cluster", sort=False):
            genes = [clean(value) for value in group["gene"].tolist()]
            rows.append({
                "annotation_level": level,
                "published_label": clean(cluster),
                "marker_rows": len(group),
                "marker_genes": "|".join(genes),
                "cell_counts_available": "no",
            })
        return rows

    def publish(self) -> dict[str, list[dict[str, object]]]:
        frames = {
            "postnatal_broad_class": read_marker_table(self.paths.broad_markers),
            "postnatal_refined_cluster": read_marker_table(self.paths.refined_markers),
            "embryonic_cluster": read_marker_table(self.paths.embryonic_markers),
        }
        for stem, frame in frames.items():
            marker_rows = self.marker_rows(frame)
            atomic_tsv(
                self.metadata_dir / f"published_{stem}_markers.tsv",
                tuple(frame.columns), marker_rows,
            )
            inventory = self.inventory(frame, stem)
            atomic_tsv(
                self.metadata_dir / f"published_{stem}_inventory.tsv",
                ("annotation_level", "published_label", "marker_rows", "marker_genes", "cell_counts_available"),
                inventory,
            )
        return {name: self.inventory(frame, name) for name, frame in frames.items()}


class BandlerHierarchyFigure:
    """Contrast observed STICR cell composition with published embryonic labels."""

    def __init__(self, class_cluster_path: Path, embryonic_inventory: Sequence[Mapping[str, object]], output_stem: Path):
        self.class_cluster_path = class_cluster_path
        self.embryonic_inventory = embryonic_inventory
        self.output_stem = output_stem

    def publish(self) -> None:
        cross = pd.read_csv(self.class_cluster_path, sep="\t")
        totals = cross.groupby("broad_class", sort=False)["cells"].sum().sort_values(ascending=True)
        fig = plt.figure(figsize=(16, 10), constrained_layout=True)
        grid = fig.add_gridspec(1, 2, width_ratios=(1.15, 1))
        left = fig.add_subplot(grid[0, 0])
        palette = plt.get_cmap("tab20")
        for row_index, broad_class in enumerate(totals.index):
            subset = cross[cross["broad_class"] == broad_class].sort_values("cells", ascending=False)
            left_edge = 0
            for cluster_index, (_, row) in enumerate(subset.iterrows()):
                left.barh(row_index, row["cells"], left=left_edge, color=palette(cluster_index % 20), height=0.72, linewidth=0)
                left_edge += row["cells"]
            left.text(left_edge + totals.max() * 0.012, row_index, f"{int(left_edge):,}", va="center", fontsize=9)
        left.set_yticks(np.arange(len(totals)), totals.index)
        left.set_xlabel("Cells in recovered author Seurat object")
        left.set_title("A. Recovered STICR object\nobserved broad class → refined cluster composition", loc="left", weight="bold")
        left.spines[["top", "right"]].set_visible(False)

        right = fig.add_subplot(grid[0, 1])
        right.axis("off")
        mitotic = [row for row in self.embryonic_inventory if str(row["published_label"]).startswith("m_")]
        immature = [row for row in self.embryonic_inventory if str(row["published_label"]).startswith("i_")]
        y = 0.97
        right.text(0, y, "B. Published embryonic annotation vocabulary", transform=right.transAxes, weight="bold", fontsize=13, va="top")
        y -= 0.065
        right.text(0, y, "Mitotic/progenitor-associated states", transform=right.transAxes, weight="bold", color="#7b3294", va="top")
        y -= 0.042
        for row in mitotic:
            genes = str(row["marker_genes"]).split("|")[:4]
            right.text(0.02, y, f"{row['published_label']}  ·  {', '.join(genes)}", transform=right.transAxes, fontsize=9.5, va="top")
            y -= 0.039
        y -= 0.02
        right.text(0, y, "Postmitotic/immature states", transform=right.transAxes, weight="bold", color="#008837", va="top")
        y -= 0.042
        for row in immature:
            genes = str(row["marker_genes"]).split("|")[:4]
            right.text(0.02, y, f"{row['published_label']}  ·  {', '.join(genes)}", transform=right.transAxes, fontsize=9, va="top")
            y -= 0.034
        right.text(
            0, 0.005,
            "The embryonic labels come from Supplementary Data 4 marker rows.\nNo per-cell CA301 assignments or embryonic UMAP were present in that workbook.",
            transform=right.transAxes, fontsize=9, color="#555555", va="bottom",
        )
        fig.suptitle("Bandler reference artifacts: observed cells versus published label definitions", fontsize=16, weight="bold")
        self.output_stem.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(self.output_stem.with_suffix(".png"), dpi=220, bbox_inches="tight")
        fig.savefig(self.output_stem.with_suffix(".pdf"), bbox_inches="tight")
        plt.close(fig)


class BandlerEvidencePublisher:
    """Coordinate the evidence package without treating aggregate markers as cell labels."""

    def __init__(self, run_dir: Path, source_root: Path):
        self.run_dir = run_dir.resolve()
        self.source_root = source_root.resolve()
        self.audit_dir = self.run_dir / "Bandler2022" / "author_object_audit"
        self.metadata_dir = self.run_dir / "Bandler2022" / "metadata"
        self.figure_dir = self.run_dir / "Bandler2022" / "figures"
        self.paths = PublishedArtifactPaths.from_source_root(self.source_root)

    def _artifact_scope(self) -> list[dict[str, object]]:
        return [
            {"artifact": "CA301/GSM5684876 deposited counts", "cells": 4516, "cell_labels": "no", "embedding": "no", "contains_CA301": "yes", "role": "Exact WT E15.5 MGE expression matrix; stable CA301-prefixed barcodes."},
            {"artifact": "Recovered STICR.seuratobject.RDS", "cells": 65700, "cell_labels": "yes", "embedding": "yes", "contains_CA301": "no", "role": "Postnatal forebrain reference with broad classes and refined clusters."},
            {"artifact": "Supplementary Data 4", "cells": "NA", "cell_labels": "definitions only", "embedding": "no", "contains_CA301": "not mappable", "role": "Twenty-one embryonic cluster definitions and marker genes; no cell IDs."},
            {"artifact": "2025 Mayer-lab interactive GE atlas", "cells": "not exported", "cell_labels": "yes in live app", "embedding": "yes in live app", "contains_CA301": "not yet proven", "role": "Later Mayer-lab reanalysis; underlying EXCIT_INHIBIT_cleaned_sub.rds is not publicly linked."},
        ]

    def _report(self, samples: Sequence[Mapping[str, object]], inventories: Mapping[str, Sequence[Mapping[str, object]]]) -> str:
        ca301 = next(row for row in samples if row["sample_id"] == "CA301")
        broad_count = len(inventories["postnatal_broad_class"])
        refined_count = len(inventories["postnatal_refined_cluster"])
        embryonic_count = len(inventories["embryonic_cluster"])
        sample_lines = [
            f"| `{row['sample_id']}` | `{row['geo_accession']}` | {row['age_of_collection']} | {row['anatomical_region']} | {row['pooled_brains']} | {row['tenx_version']} | {row['cells_after_filter_published']} |"
            for row in samples
        ]
        return "\n".join([
            "# Bandler author-object and embryonic-annotation recovery",
            "",
            "## Result",
            "",
            "The recovered 1.02-GB author artifact is a genuine Seurat object with 21,051 features, 65,700 cells, 11 observed broad classes, 51 observed refined clusters, and PCA, Harmony, and UMAP reductions. It is the **postnatal STICR reference**, not the integrated embryonic object: its 18 samples exclude CA298–CA303, and exact CA301 cell-ID overlap is zero.",
            "",
            f"The publisher supplements independently define {broad_count} postnatal broad classes, {refined_count} postnatal refined clusters, and {embryonic_count} embryonic clusters. These are marker tables, not cell-level metadata; they cannot be joined to CA301 barcodes without an additional artifact.",
            "",
            "## Embryonic integration sample ledger",
            "",
            "| Sample | GEO | Collection age | Region | Pooled brains | 10x | Published cells after filtering |",
            "| --- | --- | --- | --- | ---: | --- | ---: |",
            *sample_lines,
            "",
            f"For CA301, Supplementary Data 1 reports {ca301['cells_after_filter_published']} cells after filtering, whereas the deposited count matrix contains 4,516 columns. The four-cell difference is retained as an unresolved version/filtering discrepancy.",
            "",
            "## What is and is not recovered",
            "",
            "- Recovered: exact author STICR Seurat object, postnatal broad/fine annotations, author UMAP, official sample design, all four marker supplements, and TrackerSeq lineage metadata.",
            "- Proven: CA301 is WT MGE collected at E15.5, made from six pooled brains using 10x v2 and NovaSeq/Broad sequencing.",
            "- Proven: the paper used CA298–CA303 plus MUC28072 in its integrated embryonic analysis and published 21 embryonic cluster names.",
            "- Not recovered: the original integrated embryonic Seurat object or a barcode-to-cluster/UMAP table for CA301.",
            "- Later lead: the Mayer lab's 2025 interactive atlas uses a local `EXCIT_INHIBIT_cleaned_sub.rds` containing Bandler embryonic cells, stages, UMAP, classes, clusters, study, and cell IDs. The GitHub deployment code names this object but does not publish it or its generated TSV/HDF5 files.",
            "",
            "Therefore CA301 remains the strongest anatomical/age match, but its 4,516 deposited cells must not be assigned the 21 published embryonic labels until a barcode-preserving author or later-lab artifact is obtained.",
            "",
            "![Bandler recovered and published annotation structure](figures/01_bandler_recovered_and_published_annotation_structure.png)",
            "",
            "## Machine-readable outputs",
            "",
            "- `metadata/published_embryonic_sample_design.tsv`",
            "- `metadata/published_embryonic_cluster_inventory.tsv` and `published_embryonic_cluster_markers.tsv`",
            "- `metadata/published_postnatal_broad_class_inventory.tsv` and marker table",
            "- `metadata/published_postnatal_refined_cluster_inventory.tsv` and marker table",
            "- `author_object_audit/author_artifact_scope.tsv`",
            "- `author_object_audit/author_seurat_class_by_cluster.tsv`",
            "- `author_object_audit/author_seurat_sample_inventory.tsv`",
            "",
        ])

    def publish(self) -> None:
        self.paths.validate()
        required = self.audit_dir / "author_seurat_class_by_cluster.tsv"
        if not required.is_file():
            raise FileNotFoundError(required)
        samples = PublishedSampleLedger(self.paths.sample_ledger).rows()
        atomic_tsv(
            self.metadata_dir / "published_embryonic_sample_design.tsv",
            tuple(samples[0]), samples,
        )
        inventories = PublishedAnnotationTables(self.paths, self.metadata_dir).publish()
        scope = self._artifact_scope()
        atomic_tsv(
            self.audit_dir / "author_artifact_scope.tsv",
            ("artifact", "cells", "cell_labels", "embedding", "contains_CA301", "role"), scope,
        )
        figure_stem = self.figure_dir / "01_bandler_recovered_and_published_annotation_structure"
        BandlerHierarchyFigure(required, inventories["embryonic_cluster"], figure_stem).publish()
        atomic_text(
            self.run_dir / "Bandler2022" / "BANDLER_AUTHOR_OBJECT_RECOVERY_REPORT.md",
            self._report(samples, inventories),
        )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", required=True, type=Path)
    parser.add_argument("--source-root", required=True, type=Path)
    args = parser.parse_args()
    BandlerEvidencePublisher(args.run_dir, args.source_root).publish()


if __name__ == "__main__":
    main()
