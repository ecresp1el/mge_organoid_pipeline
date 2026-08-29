#!/usr/bin/env python3
"""Plot observed annotation hierarchy and E15-MGE evidence for the candidates."""

from __future__ import annotations

import argparse
import csv
import hashlib
import os
import platform
import shlex
import sys
import tempfile
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Iterable, Mapping, Sequence

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


TARGET_CLASSES = (
    "Radial glia", "Neuroblast", "Neuron", "Glioblast", "Immune",
    "Oligodendrocyte", "Vascular",
)
TARGET_COLOR = "#2878B5"
OTHER_COLOR = "#B8C2CC"
PALETTE = ("#2878B5", "#E07A5F", "#3D9970", "#8E6C8A", "#E9C46A", "#6C8EAD")


def decode(value: object) -> str:
    return value.decode("utf-8") if isinstance(value, bytes) else str(value)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


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


def atomic_save_figure(figure, path: Path, **kwargs) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    suffix = path.suffix
    descriptor, temporary_name = tempfile.mkstemp(prefix=f".{path.stem}.", suffix=suffix, dir=str(path.parent))
    os.close(descriptor)
    temporary = Path(temporary_name)
    try:
        figure.savefig(temporary, **kwargs)
        os.replace(temporary, path)
    finally:
        if temporary.exists():
            temporary.unlink()


class H5adCategoricalReader:
    """Read selected H5AD observation categoricals without loading expression."""

    def __init__(self, path: Path):
        self.path = path

    @staticmethod
    def _categorical(group) -> np.ndarray:
        categories = np.asarray([decode(value) for value in group["categories"][:]], dtype=object)
        codes = group["codes"][:]
        values = np.full(codes.shape, "NA", dtype=object)
        valid = codes >= 0
        values[valid] = categories[codes[valid]]
        return values

    def read(self, columns: Sequence[str]) -> Dict[str, np.ndarray]:
        with h5py.File(self.path, "r") as handle:
            obs = handle["obs"]
            missing = [column for column in columns if column not in obs]
            if missing:
                raise RuntimeError(f"Required H5AD obs columns are absent: {missing}")
            result = {}
            for column in columns:
                node = obs[column]
                if not isinstance(node, h5py.Group) or "categories" not in node or "codes" not in node:
                    raise RuntimeError(f"Expected categorical H5AD encoding for obs/{column}")
                result[column] = self._categorical(node)
            lengths = {len(values) for values in result.values()}
            if len(lengths) != 1:
                raise RuntimeError(f"Observation columns have inconsistent lengths: {lengths}")
            return result


class HierarchyComposition:
    """Compute exact cell-level co-occurrence across author annotation levels."""

    def __init__(self, obs: Mapping[str, np.ndarray]):
        self.obs = obs
        self.class_counts = Counter(obs["Class"])
        self.subclass_by_class: Dict[str, Counter] = defaultdict(Counter)
        self.celltype_by_subclass: Dict[str, Counter] = defaultdict(Counter)
        self.age_by_dissection: Dict[str, Counter] = defaultdict(Counter)
        for broad, subclass, cell_type, age, dissection in zip(
            obs["Class"], obs["Subclass"], obs["CellType"], obs["Age"], obs["dissection"],
        ):
            self.subclass_by_class[str(broad)][str(subclass)] += 1
            self.celltype_by_subclass[str(subclass)][str(cell_type)] += 1
            self.age_by_dissection[str(age)][str(dissection)] += 1

    def publish_tables(self, metadata_dir: Path) -> list[Path]:
        class_path = metadata_dir / "class_by_subclass_composition.tsv"
        age_path = metadata_dir / "age_by_dissection_composition.tsv"
        fine_path = metadata_dir / "forebrain_gabaergic_celltype_composition.tsv"
        atomic_write_tsv(class_path, ("Class", "Subclass", "n_cells", "percent_within_Class"), (
            {
                "Class": broad, "Subclass": subclass, "n_cells": count,
                "percent_within_Class": f"{100 * count / self.class_counts[broad]:.6f}",
            }
            for broad in sorted(self.subclass_by_class)
            for subclass, count in self.subclass_by_class[broad].most_common()
        ))
        atomic_write_tsv(age_path, ("Age", "dissection", "n_cells"), (
            {"Age": age, "dissection": dissection, "n_cells": count}
            for age in sorted(self.age_by_dissection)
            for dissection, count in self.age_by_dissection[age].most_common()
        ))
        atomic_write_tsv(fine_path, ("Subclass", "CellType", "n_cells", "percent_within_Subclass"), (
            {
                "Subclass": "Forebrain GABAergic", "CellType": cell_type,
                "n_cells": count,
                "percent_within_Subclass": f"{100 * count / sum(self.celltype_by_subclass['Forebrain GABAergic'].values()):.6f}",
            }
            for cell_type, count in self.celltype_by_subclass["Forebrain GABAergic"].most_common()
        ))
        return [class_path, age_path, fine_path]


class HierarchyPlotter:
    """Render broad, subclass, and fine-label composition without overplotting."""

    def __init__(self, composition: HierarchyComposition):
        self.composition = composition

    def build(self):
        figure = plt.figure(figsize=(17, 12), constrained_layout=True)
        grid = figure.add_gridspec(2, 2, width_ratios=(0.9, 1.35), height_ratios=(1.0, 1.0))
        broad_ax = figure.add_subplot(grid[:, 0])
        nested_ax = figure.add_subplot(grid[0, 1])
        fine_ax = figure.add_subplot(grid[1, 1])

        broad_rows = self.composition.class_counts.most_common()
        labels = [row[0] for row in broad_rows][::-1]
        values = [row[1] for row in broad_rows][::-1]
        colors = [TARGET_COLOR if label in TARGET_CLASSES else OTHER_COLOR for label in labels]
        broad_ax.barh(labels, values, color=colors, edgecolor="white", linewidth=0.4)
        broad_ax.set_title("A. Author Class composition", loc="left", fontweight="bold")
        broad_ax.set_xlabel("Cells")
        broad_ax.grid(axis="x", color="#E5E7EB", linewidth=0.7)
        broad_ax.set_axisbelow(True)
        for index, value in enumerate(values):
            broad_ax.text(value + max(values) * 0.008, index, f"{value:,}", va="center", fontsize=7)

        y_positions = np.arange(len(TARGET_CLASSES))
        for y, broad in zip(y_positions, TARGET_CLASSES):
            counter = self.composition.subclass_by_class[broad]
            total = sum(counter.values())
            top = counter.most_common(4)
            other = total - sum(value for _, value in top)
            pieces = top + ([('Other subclasses', other)] if other else [])
            left = 0.0
            for index, (subclass, count) in enumerate(pieces):
                percent = 100 * count / total
                color = PALETTE[index % len(PALETTE)] if subclass != "Other subclasses" else "#D7DCE2"
                nested_ax.barh(y, percent, left=left, color=color, height=0.68, edgecolor="white", linewidth=0.7)
                if percent >= 9:
                    label = subclass if len(subclass) <= 23 else subclass[:21] + "…"
                    nested_ax.text(left + percent / 2, y, f"{label}\n{percent:.1f}%", ha="center", va="center", fontsize=6.5)
                left += percent
            nested_ax.text(101.2, y, f"n={total:,}", va="center", fontsize=7)
        nested_ax.set_yticks(y_positions, TARGET_CLASSES)
        nested_ax.invert_yaxis()
        nested_ax.set_xlim(0, 113)
        nested_ax.set_xlabel("Percent within Class")
        nested_ax.set_title("B. Observed Class → Subclass composition\n(top four subclasses plus remainder)", loc="left", fontweight="bold")
        nested_ax.spines[["top", "right", "left"]].set_visible(False)

        fine_rows = self.composition.celltype_by_subclass["Forebrain GABAergic"].most_common(15)
        fine_labels = [row[0] for row in fine_rows][::-1]
        fine_values = [row[1] for row in fine_rows][::-1]
        fine_ax.barh(fine_labels, fine_values, color="#3D9970")
        fine_ax.set_title("C. Forebrain GABAergic → CellType\n(top 15 exact author labels)", loc="left", fontweight="bold")
        fine_ax.set_xlabel("Cells")
        fine_ax.grid(axis="x", color="#E5E7EB", linewidth=0.7)
        fine_ax.set_axisbelow(True)
        for index, value in enumerate(fine_values):
            fine_ax.text(value + max(fine_values) * 0.01, index, f"{value:,}", va="center", fontsize=7)

        figure.suptitle(
            "La Manno 2021: author annotation hierarchy and composition",
            fontsize=17, fontweight="bold",
        )
        figure.text(
            0.5, 0.005,
            "Class/Subclass/CellType nesting is computed from the same P0 cells. Counts are atlas-wide, not MGE-specific.",
            ha="center", fontsize=9, color="#4B5563",
        )
        return figure


class CandidateEvidencePlotter:
    """Show E15-MGE and annotation availability as yes/no/unproven evidence."""

    ROWS = ("La Manno 2021", "Bandler 2022", "Mayer 2018")
    COLUMNS = ("Exact E15\npresent", "MGE directly\nidentified", "Exact E15 MGE\nproven", "Cell labels\nin P0", "Author embedding\nin P0")
    VALUES = np.asarray((
        (2, 1, 1, 2, 2),
        (2, 2, 2, 0, 0),
        (0, 2, 0, 0, 0),
    ))
    TEXT = (
        ("Yes, global", "Not explicit\n(ventral forebrain)", "Unproven", "Yes", "Yes"),
        ("Yes", "Yes", "Yes", "No", "No"),
        ("No (E13.5)", "Yes", "No", "No", "No"),
    )

    def build(self):
        from matplotlib.colors import ListedColormap

        figure, ax = plt.subplots(figsize=(12, 6.3))
        figure.subplots_adjust(left=0.17, right=0.98, bottom=0.23, top=0.72)
        cmap = ListedColormap(("#D1D5DB", "#F2C14E", "#4C956C"))
        ax.imshow(self.VALUES, cmap=cmap, vmin=0, vmax=2, aspect="auto")
        ax.set_xticks(np.arange(len(self.COLUMNS)), self.COLUMNS)
        ax.set_yticks(np.arange(len(self.ROWS)), self.ROWS)
        ax.tick_params(top=True, bottom=False, labeltop=True, labelbottom=False, length=0)
        for row in range(self.VALUES.shape[0]):
            for column in range(self.VALUES.shape[1]):
                ax.text(column, row, self.TEXT[row][column], ha="center", va="center", fontsize=10,
                        color="white" if self.VALUES[row, column] == 2 else "#1F2937", fontweight="bold")
        ax.set_xticks(np.arange(-0.5, len(self.COLUMNS), 1), minor=True)
        ax.set_yticks(np.arange(-0.5, len(self.ROWS), 1), minor=True)
        ax.grid(which="minor", color="white", linewidth=2)
        ax.tick_params(which="minor", bottom=False, left=False)
        ax.set_title("Candidate reference evidence: E15 MGE and immediate annotation readiness", fontsize=15, fontweight="bold", pad=55)
        figure.text(
            0.5, 0.045,
            "Green = directly observed/proven; amber = partial or unproven; gray = absent from the inspected P0.\n"
            "Bandler: 4,516 exact WT E15 MGE cells. Mayer: 6,515 MGE E13.5 Lhx6-positive cells. La Manno: exact E15 MGE count not yet established.",
            ha="center", fontsize=9, color="#374151",
        )
        return figure


class VisualAuditWorkflow:
    """Coordinate evidence extraction, plotting, manifests, and provenance."""

    def __init__(self, run_dir: Path):
        self.run_dir = run_dir.resolve()
        source_rows = []
        with (self.run_dir / "LaManno2021" / "audit" / "source_file.tsv").open("r", encoding="utf-8", newline="") as handle:
            source_rows = list(csv.DictReader(handle, delimiter="\t"))
        if len(source_rows) != 1:
            raise RuntimeError("Expected one La Manno source record")
        self.h5ad = Path(source_rows[0]["resolved_path"])

    def run(self) -> None:
        obs = H5adCategoricalReader(self.h5ad).read(("Class", "Subclass", "CellType", "Age", "dissection"))
        composition = HierarchyComposition(obs)
        metadata_paths = composition.publish_tables(self.run_dir / "LaManno2021" / "metadata")
        hierarchy = HierarchyPlotter(composition).build()
        evidence = CandidateEvidencePlotter().build()
        figure_paths = []
        for stem, figure, directory in (
            ("01_author_annotation_hierarchy_composition", hierarchy, self.run_dir / "LaManno2021" / "figures"),
            ("02_candidate_e15_mge_annotation_evidence", evidence, self.run_dir / "figures"),
        ):
            for suffix in (".png", ".pdf"):
                path = directory / f"{stem}{suffix}"
                atomic_save_figure(figure, path, dpi=300, bbox_inches="tight")
                figure_paths.append(path)
            plt.close(figure)
        status_path = self.run_dir / "tables" / "candidate_e15_mge_status.tsv"
        atomic_write_tsv(status_path, (
            "paper", "exact_E15_present", "MGE_directly_identified", "exact_E15_MGE_proven",
            "relevant_cells", "evidence", "interpretation",
        ), (
            {"paper": "LaManno2021", "exact_E15_present": "yes_global", "MGE_directly_identified": "no_explicit_label", "exact_E15_MGE_proven": "no", "relevant_cells": "not_defined", "evidence": "Age plus dissection fields", "interpretation": "E15 and ventral forebrain coexist, but MGE is not an explicit author label"},
            {"paper": "Bandler2022", "exact_E15_present": "yes", "MGE_directly_identified": "yes", "exact_E15_MGE_proven": "yes", "relevant_cells": "4516", "evidence": "CA301/GSM5684876 WT E15 MGE sample", "interpretation": "Exact match; counts-only P0 requires label join"},
            {"paper": "Mayer2018", "exact_E15_present": "no", "MGE_directly_identified": "yes", "exact_E15_MGE_proven": "no", "relevant_cells": "6515", "evidence": "MGE_E13.5_Lhx6pos sample token", "interpretation": "MGE is E13.5 and experimentally enriched"},
        ))
        generated = [*metadata_paths, *figure_paths, status_path]
        atomic_write_tsv(self.run_dir / "provenance" / "hierarchy_plot_output_manifest.tsv", (
            "relative_path", "bytes", "sha256",
        ), ({
            "relative_path": str(path.relative_to(self.run_dir)),
            "bytes": path.stat().st_size,
            "sha256": sha256_file(path),
        } for path in generated))
        script_path = Path(__file__).resolve()
        atomic_write_tsv(self.run_dir / "provenance" / "hierarchy_plot_generation.tsv", (
            "generated_utc", "python", "command", "executed_script", "script_sha256",
        ), ({
            "generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
            "python": platform.python_version(),
            "command": shlex.join(sys.argv),
            "executed_script": str(script_path),
            "script_sha256": sha256_file(script_path),
        },))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", required=True)
    args = parser.parse_args()
    VisualAuditWorkflow(Path(args.run_dir)).run()


if __name__ == "__main__":
    main()
