#!/usr/bin/env python3
"""Build the human-readable report from a completed reference-curation run."""

from __future__ import annotations

import argparse
import csv
import hashlib
import os
import platform
import shlex
import sys
import tempfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Sequence


STUDIES = ("LaManno2021", "Bandler2022", "Mayer2018")
ANNOTATION_DICTIONARY_COLUMNS = (
    "annotation_column", "published_label", "parent_label_if_known",
    "annotation_level", "n_cells_global", "n_cells_MGE",
    "n_cells_age_matched_MGE", "paper_figure_or_source", "notes",
)
OBSERVED_SUMMARY_COLUMNS = (
    "paper", "P0_file", "actual_object_type", "dimensions",
    "published_samples", "P0_linked_samples", "P0_linked_sample_ids",
    "age_match", "MGE_match", "author_embedding", "embedded_annotations",
    "immediate_cell_type_visibility", "selection_bias", "next_required_action",
)
LEVEL_BY_COLUMN = {
    "Age": "developmental_age",
    "Class": "broad_class",
    "Subclass": "subclass",
    "CellType": "fine_cell_type",
    "cluster_id": "cluster",
    "dissection": "anatomical_region",
    "development_stage_ontology_term_id": "developmental_age_ontology",
    "tissue_ontology_term_id": "tissue_ontology",
    "suspension_type": "technical_metadata",
    "tissue_type": "technical_metadata",
}


def read_tsv(path: Path) -> List[Dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return [dict(row) for row in csv.DictReader(handle, delimiter="\t")]


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def atomic_write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(prefix=f".{path.name}.", dir=str(path.parent), text=True)
    temporary = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
            handle.write(text)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    finally:
        if temporary.exists():
            temporary.unlink()


def atomic_write_tsv(path: Path, columns: Sequence[str], rows: Iterable[Mapping[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(prefix=f".{path.name}.", dir=str(path.parent), text=True)
    temporary = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, delimiter="\t", fieldnames=list(columns), lineterminator="\n")
            writer.writeheader()
            for row in rows:
                writer.writerow({column: row.get(column, "") for column in columns})
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    finally:
        if temporary.exists():
            temporary.unlink()


class AnnotationDictionaryPublisher:
    """Publish only labels that were observed in an inspected cell object."""

    def __init__(self, run_dir: Path):
        self.run_dir = run_dir

    def publish(self) -> Dict[str, List[Dict[str, str]]]:
        dictionaries: Dict[str, List[Dict[str, str]]] = {}
        for study in STUDIES:
            rows: List[Dict[str, str]] = []
            counts_dir = self.run_dir / study / "audit" / "obs_value_counts"
            if counts_dir.is_dir():
                for path in sorted(counts_dir.glob("*.tsv")):
                    column = path.stem
                    for observed in read_tsv(path):
                        rows.append({
                            "annotation_column": column,
                            "published_label": observed.get("published_label", ""),
                            "parent_label_if_known": "NA",
                            "annotation_level": LEVEL_BY_COLUMN.get(column, "author_metadata"),
                            "n_cells_global": observed.get("n_cells", ""),
                            "n_cells_MGE": "NA",
                            "n_cells_age_matched_MGE": "NA",
                            "paper_figure_or_source": f"{study}/audit/obs_value_counts/{path.name}",
                            "notes": "Directly observed in the P0 object; parent and MGE-specific counts were not inferred.",
                        })
            atomic_write_tsv(
                self.run_dir / study / "metadata" / "annotation_dictionary.tsv",
                ANNOTATION_DICTIONARY_COLUMNS, rows,
            )
            dictionaries[study] = rows
        return dictionaries


class ReferenceCurationReportBuilder:
    """Summarize observed object contents without upgrading paper claims to cell evidence."""

    def __init__(self, run_dir: Path):
        self.run_dir = run_dir.resolve()
        success = self.run_dir / "SUCCESS.txt"
        if not success.is_file():
            raise RuntimeError(f"Completed SUCCESS marker is required: {success}")
        self.checkpoints = {row["paper"]: row for row in read_tsv(self.run_dir / "tables" / "early_processed_object_checkpoint.tsv")}
        self.summaries = {row["paper"]: row for row in read_tsv(self.run_dir / "tables" / "study_sample_summary.tsv")}

    @staticmethod
    def _count_lookup(dictionary: Sequence[Mapping[str, str]], column: str, label: str) -> str:
        for row in dictionary:
            if row.get("annotation_column") == column and row.get("published_label") == label:
                return row.get("n_cells_global", "")
        return "NOT_OBSERVED"

    def _observed_summary(self) -> List[Dict[str, str]]:
        rows = []
        for study in STUDIES:
            checkpoint = self.checkpoints[study]
            summary = self.summaries[study]
            if study == "LaManno2021":
                age_match = "Exact e15.0 and e15.5 labels are present; MGE-by-age count is not yet proven"
                mge_match = "Ventral forebrain dissections exist, but a direct MGE definition is not yet proven"
                visibility = "YES: Class, Subclass, CellType, and cluster_id are embedded"
                selection = "No selection field recovered at this checkpoint"
            elif study == "Bandler2022":
                age_match = "Exact E15 sample"
                mge_match = "Exact WT MGE sample CA301/GSM5684876"
                visibility = "NO: deposited P0 is a counts-only sparse matrix"
                selection = "GEO metadata includes GFP/FACS-selected samples; CA301 is registered WT MGE"
            else:
                age_match = "E13.5, not E15"
                mge_match = "MGE sample is directly identifiable"
                visibility = "NO: deposited P0 is a counts-only CSV"
                selection = "MGE_E13.5_Lhx6pos; strong positive-selection bias"
            rows.append({
                "paper": study,
                "P0_file": checkpoint.get("P0_file", ""),
                "actual_object_type": checkpoint.get("actual_object_type", ""),
                "dimensions": checkpoint.get("dimensions", ""),
                "published_samples": summary.get("published_sample_count", ""),
                "P0_linked_samples": summary.get("processed_object_sample_count", ""),
                "P0_linked_sample_ids": summary.get("processed_object_sample_ids", ""),
                "age_match": age_match,
                "MGE_match": mge_match,
                "author_embedding": checkpoint.get("author_embedding_present", ""),
                "embedded_annotations": checkpoint.get("annotation_columns_present", ""),
                "immediate_cell_type_visibility": visibility,
                "selection_bias": selection,
                "next_required_action": checkpoint.get("next_minimal_action", ""),
            })
        return rows

    def _report_text(self, dictionaries: Mapping[str, Sequence[Mapping[str, str]]]) -> str:
        lamanno = dictionaries["LaManno2021"]
        count = lambda column, label: self._count_lookup(lamanno, column, label)
        lines = [
            "# Developing-mouse MGE reference curation: observed-object report",
            "",
            f"Run: `{self.run_dir.name}`",
            "",
            "This report describes the files that were actually downloaded and inspected. It does not treat labels mentioned in a paper as cell-level annotations unless those labels were present in the P0 object or joined by a stable identifier. No Paper 3 cells were loaded, mapped, integrated, or annotated.",
            "",
            "## Main result",
            "",
            "| Candidate | Actual P0 object | Cells | Cell annotations in P0? | Author embedding in P0? | Immediate implication |",
            "| --- | --- | ---: | --- | --- | --- |",
            "| La Manno 2021 | Backed AnnData H5AD, 31,053 genes | 292,495 | Yes: `Class`, `Subclass`, `CellType`, `cluster_id`, age and dissection | Yes: `X_UMAP`, `X_tSNE` | Only candidate immediately viewable with author cell labels; MGE subset still needs an author-supported definition. |",
            "| Bandler 2022 | `dgCMatrix`, 19,808 genes | 4,516 | No | No | Exact WT E15 MGE sample, but published embryonic labels require a barcode join. |",
            "| Mayer 2018 | Gzip CSV, 19,272 data rows | 42,418 | No | No | Samples are identifiable, including 6,515 enriched E13.5 MGE cells, but published labels/embedding require another artifact. |",
            "",
            "The checkpoint therefore does **not** yet justify choosing a final reference. La Manno is most immediately inspectable; Bandler is the closest biological match; Mayer is traceable but younger and Lhx6-positive enriched.",
            "",
            "## La Manno 2021",
            "",
            "The H5AD contains 292,495 cells, 31,053 genes, 93 object-linked sample IDs, 20 exact age labels, 24 broad `Class` labels, 133 `Subclass` labels, 748 `CellType` categories including the deposited `nan` category, and 798 `cluster_id` values. Exact e15.0 and e15.5 cells are present globally (18,443 and 8,796 cells, respectively).",
            "",
            "Relevant exact author labels observed in the object:",
            "",
            "| Requested biology | Exact observed author field/label | Global cells | Current limitation |",
            "| --- | --- | ---: | --- |",
            f"| Radial glia/apical progenitor | `Class = Radial glia` | {count('Class', 'Radial glia')} | MGE-specific subset not yet counted |",
            f"| Intermediate progenitor | `Subclass = Neuronal intermediate progenitor` | {count('Subclass', 'Neuronal intermediate progenitor')} | Parent hierarchy not inferred |",
            f"| Neuroblast | `Class = Neuroblast` | {count('Class', 'Neuroblast')} | MGE-specific subset not yet counted |",
            f"| Neuron | `Class = Neuron` | {count('Class', 'Neuron')} | Includes many regions/lineages |",
            f"| Forebrain GABAergic | `Subclass = Forebrain GABAergic` | {count('Subclass', 'Forebrain GABAergic')} | MGE versus other forebrain sources not yet separated |",
            f"| Glioblast | `Class = Glioblast` | {count('Class', 'Glioblast')} | Requires subtype review |",
            f"| Forebrain astrocyte | `Subclass = Forebrain astrocyte` | {count('Subclass', 'Forebrain astrocyte')} | MGE-specific subset not yet counted |",
            f"| Immune | `Class = Immune` | {count('Class', 'Immune')} | Broad class |",
            f"| Non-cycling microglia | `Subclass = Non-cycling microglia` | {count('Subclass', 'Non-cycling microglia')} | Global count |",
            f"| Cycling microglia | `Subclass = Cycling microglia` | {count('Subclass', 'Cycling microglia')} | Not a neural progenitor label |",
            f"| Oligodendrocyte lineage | `Class = Oligodendrocyte` | {count('Class', 'Oligodendrocyte')} | Global count |",
            f"| OPC | `Subclass = Oligodendrocyte precursor cell` | {count('Subclass', 'Oligodendrocyte precursor cell')} | Global count |",
            f"| PreOPC | `Subclass = PreOPC` | {count('Subclass', 'PreOPC')} | Global count |",
            f"| Vascular | `Class = Vascular` | {count('Class', 'Vascular')} | Broad class |",
            f"| Endothelial | `Subclass = Endothelial` | {count('Subclass', 'Endothelial')} | Global count |",
            f"| Pericyte | `Subclass = Pericyte` | {count('Subclass', 'Pericyte')} | Global count |",
            "",
            "No explicit author label containing `MGE`, `ganglionic`, `NKX2-1`, `LHX6`, `NPY`, `NXPH1`, `SIX3`, or `GUCY1A3` was found by exact string search in the exported annotation vocabularies. `dissection` includes `ForebrainVentral` (22,642 cells) and `ForebrainVentroLateral` (17,277), but those must not be relabeled MGE without the author taxonomy/code definition. An explicit cycling **neural** progenitor label was also not recovered by string match; `Neuronal intermediate progenitor` is present, while cycling labels currently observed refer to immune populations.",
            "",
            "## Bandler 2022",
            "",
            "The deposited CA301 file has two gzip layers and resolves to a sparse `dgCMatrix` with 19,808 gene rows and 4,516 uniquely named cell columns. Cell IDs retain the `CA301_` prefix followed by a 10x barcode, providing a stable basis for a future join. The object contains no metadata columns, Seurat reductions, UMAP/tSNE coordinates, or cell-type labels. Therefore labels such as published embryonic states cannot yet be assigned to these 4,516 cells from this file alone.",
            "",
            "Biological strength: CA301/GSM5684876 is the exact WT E15 MGE sample. Immediate limitation: find the smallest author/supplement/intermediate table mapping `CA301_<barcode>` to the published embryonic annotation and embedding.",
            "",
            "## Mayer 2018",
            "",
            "The deposited 10x CSV contains 42,418 cell columns across six directly recoverable sample tokens:",
            "",
            "| Sample token | Cells | Relevance |",
            "| --- | ---: | --- |",
            "| `CGE_E13.5_Lhx6neg` | 7,106 | CGE comparator |",
            "| `MGE_E13.5_Lhx6pos` | 6,515 | Relevant MGE population, but Lhx6-positive selected and younger than E15 |",
            "| `Cortex_E18.5` | 8,625 | Later cortex |",
            "| `Subcortex_E18.5` | 8,237 | Later subcortex |",
            "| `Cortex_P10` | 6,346 | Postnatal cortex |",
            "| `Subcortex_P10` | 5,589 | Postnatal subcortex |",
            "",
            "The CSV contains expression only: no cell metadata columns and no saved embedding. It proves which sample each cell column belongs to, but not the published progenitor/neuronal label for each cell. The next action is to locate the smallest author metadata or serialized analysis object that maps these exact column IDs to the published annotations and tSNE/UMAP coordinates.",
            "",
            "## What can and cannot be concluded now",
            "",
            "- La Manno clearly contains the broad dissection classes requested, including radial glia, intermediate progenitors, neuroblasts, neurons, astroglial, immune/microglial, oligodendrocyte/OPC, endothelial, vascular and pericyte labels.",
            "- La Manno is the only candidate whose author annotations and embedding coordinates are immediately available in the inspected P0 object.",
            "- Bandler is the best exact age/anatomy match, but its deposited P0 is not annotation-ready.",
            "- Mayer's relevant MGE component is E13.5 and Lhx6-positive selected, so it is not a census of an unbiased whole E15 MGE dissection.",
            "- Exact MGE-specific and age-matched label counts, author-embedding plots, and a final readiness ranking remain the next reviewed stage; they were not inferred during this checkpoint.",
            "",
            "## Reproducibility and outputs",
            "",
            "The exact submitted Python, R, configuration, source registry, handoff, shell submitter and SLURM files are under `code/` and `config/`. The report generator executed from `code/generate_reference_curation_report.py`. Detailed evidence is in:",
            "",
            "- `tables/early_processed_object_checkpoint.tsv`",
            "- `tables/all_candidate_reference_samples.tsv`",
            "- `tables/study_sample_summary.tsv`",
            "- `tables/reference_curation_requirements_ledger.tsv`",
            "- `tables/reference_annotation_availability.tsv`",
            "- `<study>/metadata/annotation_dictionary.tsv`",
            "- `<study>/audit/` structure, identifiers, embedding and value-count files.",
            "",
            "The Turbo filesystem reports the cached source files as writable despite the workflow's read-only chmod request. Their exact byte sizes and SHA-256 values are recorded in each `<study>/audit/source_file.tsv`; the workflow did not modify their contents after download.",
        ]
        return "\n".join(lines) + "\n"

    def run(self) -> None:
        dictionaries = AnnotationDictionaryPublisher(self.run_dir).publish()
        summary_path = self.run_dir / "tables" / "reference_annotation_availability.tsv"
        report_path = self.run_dir / "REFERENCE_CURATION_REPORT.md"
        atomic_write_tsv(
            summary_path,
            OBSERVED_SUMMARY_COLUMNS, self._observed_summary(),
        )
        atomic_write_text(report_path, self._report_text(dictionaries))
        generated_paths = [
            report_path, summary_path,
            *(self.run_dir / study / "metadata" / "annotation_dictionary.tsv" for study in STUDIES),
        ]
        atomic_write_tsv(
            self.run_dir / "provenance" / "report_output_manifest.tsv",
            ("relative_path", "bytes", "sha256"),
            ({
                "relative_path": str(path.relative_to(self.run_dir)),
                "bytes": path.stat().st_size,
                "sha256": sha256_file(path),
            } for path in generated_paths),
        )
        script_path = Path(__file__).resolve()
        atomic_write_tsv(
            self.run_dir / "provenance" / "report_generation.tsv",
            ("generated_utc", "python", "command", "executed_script", "script_sha256"),
            ({
                "generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
                "python": platform.python_version(),
                "command": shlex.join(sys.argv),
                "executed_script": str(script_path),
                "script_sha256": sha256_file(script_path),
            },),
        )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", required=True)
    args = parser.parse_args()
    ReferenceCurationReportBuilder(Path(args.run_dir)).run()


if __name__ == "__main__":
    main()
