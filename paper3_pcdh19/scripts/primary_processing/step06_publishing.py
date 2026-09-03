"""Atomic publication, provenance, reporting, and approval-ledger updates."""

from __future__ import annotations

import os
import shutil
import tempfile
from datetime import datetime, timezone
from importlib import metadata
from pathlib import Path

import pandas as pd

from .step06_validation import sha256


def utc_now() -> str:
    """Return a timezone-aware UTC timestamp."""

    return datetime.now(timezone.utc).isoformat()


class AtomicStep06Publisher:
    """Keep every Step 06 output private until validation succeeds."""

    def __init__(self, run_dir: Path):
        """Create run-scoped staging and publication directories."""

        self.run_dir = run_dir
        self.stage_dir = Path(tempfile.mkdtemp(prefix=".step06_staging_", dir=run_dir))
        self.objects_dir = self.stage_dir / "objects"
        self.tables_dir = self.stage_dir / "tables"
        self.figures_dir = self.stage_dir / "figures"
        for directory in (self.objects_dir, self.tables_dir, self.figures_dir):
            directory.mkdir()

    def publish(self) -> None:
        """Atomically expose all validated output groups and status files."""

        targets = {
            self.objects_dir: self.run_dir / "objects",
            self.tables_dir: self.run_dir / "tables",
            self.figures_dir: self.run_dir / "figures",
            self.stage_dir / "STEP06_DIAGNOSTIC_REPORT.md": self.run_dir / "STEP06_DIAGNOSTIC_REPORT.md",
            self.stage_dir / "STEP_STATUS.tsv": self.run_dir / "STEP_STATUS.tsv",
        }
        for target in targets.values():
            if target.exists():
                raise FileExistsError(f"Refusing to replace published Step 06 asset: {target}")
        for source, target in targets.items():
            os.replace(source, target)
        self.stage_dir.rmdir()

    def discard(self) -> None:
        """Remove only unpublished staging outputs after an error."""

        if self.stage_dir.exists():
            shutil.rmtree(self.stage_dir)


class Step06ProvenancePublisher:
    """Create software, status, manifest, and human-readable report assets."""

    packages = (
        "anndata", "h5py", "igraph", "leidenalg", "matplotlib", "numpy",
        "pandas", "scanpy", "scikit-learn", "scikit-misc", "scipy",
        "seaborn", "umap-learn",
    )

    def software_versions(self) -> pd.DataFrame:
        """Return exact interpreter and package versions."""

        rows = [{"component": "python", "version": os.sys.version.replace("\n", " ")}]
        rows.extend({"component": name, "version": metadata.version(name)} for name in self.packages)
        return pd.DataFrame(rows)

    def output_manifest(self, stage_dir: Path) -> pd.DataFrame:
        """Hash every staged file except the manifest itself."""

        rows = []
        for path in sorted(item for item in stage_dir.rglob("*") if item.is_file()):
            if path.name == "output_manifest.tsv":
                continue
            rows.append({
                "relative_path": str(path.relative_to(stage_dir)),
                "bytes": path.stat().st_size,
                "sha256": sha256(path),
            })
        return pd.DataFrame(rows)

    def status_frame(self, run_id: str, input_run_id: str, cells: int, genes: int) -> pd.DataFrame:
        """Return the run-local IN_REVIEW status record."""

        return pd.DataFrame([{
            "step": "06_technical_sample_batch_diagnostics",
            "run_id": run_id,
            "status": "IN_REVIEW",
            "input_checkpoint": input_run_id,
            "output_checkpoint": "objects/pcdh19_step06_unintegrated_diagnostics.h5ad",
            "completed_utc": utc_now(),
            "cells_before": cells,
            "cells_after": cells,
            "genes_before": genes,
            "genes_after": genes,
            "approval_decision": "",
            "approved_run_id": "",
        }])

    def report(self, run_id: str, input_run_id: str, results, checks: pd.DataFrame, settings) -> str:
        """Render the Step 06 methods, findings, and review boundary."""

        failures = int((checks["status"] != "PASS").sum())
        lines = [
            "# Step 06 unintegrated technical/sample/batch diagnostic report",
            "",
            f"- Run ID: `{run_id}`.",
            "- Status: **IN_REVIEW**.",
            f"- Exact approved input: `{input_run_id}`.",
            "- Lineage: Steps 03–05 bypassed; no scDblFinder or ambient/contaminant output consumed.",
            "- Runtime: Python-only AnnData/Scanpy ecosystem; no R runtime or R-based package invoked.",
            f"- Cells/genes retained: **{settings.expected_cells:,} x {settings.expected_genes:,}; zero removed**.",
            f"- Normalization: total-count normalization to {settings.target_sum:g}, followed by log1p.",
            f"- HVGs: {settings.hvg_flavor}, {settings.n_top_genes:,}, selected with technical sample as batch key.",
            f"- PCA: {settings.pca_components} components; neighbors: k={settings.n_neighbors}, metric={settings.neighbor_metric}, n_jobs={settings.n_jobs}.",
            f"- UMAP: min_dist={settings.umap_min_dist}, spread={settings.umap_spread}, seed={settings.random_seed}.",
            f"- Leiden: resolution={settings.leiden_resolution}, igraph flavor, seed={settings.random_seed}.",
            f"- Rendering: deterministic sample-balanced subset capped at {settings.render_max_cells:,}; all numerical analyses used all cells.",
            f"- Validation: **{len(checks)-failures} PASS, {failures} FAIL**.",
            "",
            "## Provisional diagnostic outcome",
            "",
            f"**Outcome {results.outcome_code}: {results.outcome_label}.**",
            "",
            results.outcome_rationale,
            "",
            "This outcome is descriptive and does not authorize integration or batch correction.",
            "",
            "## Primary review asset",
            "",
            "`figures/step06_primary_diagnostic_report.pdf` contains the coherent A–L report.",
            "Machine-readable companion tables and the rendering-cell ledger are under `tables/`.",
            "",
            "## Review boundary",
            "",
            "The raw Step 02 checkpoint remains unchanged. Step 06 adds diagnostic metadata,",
            "coordinates, and graph state to a new checkpoint, removes no cells or genes, and",
            "stops before integration, correction, or final annotation.",
            "",
        ]
        return "\n".join(lines)


class Step06ApprovalLedger:
    """Append a completed Step 06 run without implying approval."""

    def update(self, path: Path, run_id: str, input_run_id: str, settings, results, code_version: str) -> None:
        """Insert one IN_REVIEW row atomically after verifying lineage."""

        ledger = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
        approved = ledger.loc[(ledger["run_id"] == input_run_id) & (ledger["status"] == "APPROVED")]
        if len(approved) != 1:
            raise ValueError("Approval ledger lacks exactly one approved Step 02 input")
        ledger = ledger.loc[ledger["run_id"] != run_id].copy()
        row = {column: "" for column in ledger.columns}
        row.update({
            "step": "06_technical_sample_batch_diagnostics",
            "run_id": run_id,
            "status": "IN_REVIEW",
            "input_checkpoint": input_run_id,
            "output_checkpoint": "objects/pcdh19_step06_unintegrated_diagnostics.h5ad",
            "date_run_timestamp": utc_now(),
            "code_config_version": code_version,
            "summary": "Unintegrated A-L technical/sample/batch diagnostic report; zero cells removed; awaiting review.",
            "cells_before": str(settings.expected_cells),
            "cells_after": str(settings.expected_cells),
            "genes_before": str(settings.expected_genes),
            "genes_after": str(settings.expected_genes),
            "important_parameters": f"normalize_total={settings.target_sum:g}; log1p; {settings.n_top_genes} {settings.hvg_flavor} HVGs; PCA={settings.pca_components}; k={settings.n_neighbors}; unintegrated UMAP; Leiden={settings.leiden_resolution}",
            "important_findings": f"Provisional outcome {results.outcome_code}: {results.outcome_label}.",
            "outstanding_questions": "Human review of sample structure and whether integration is justified.",
        })
        ledger = pd.concat([ledger, pd.DataFrame([row])], ignore_index=True)
        temporary = path.with_suffix(path.suffix + ".tmp")
        ledger.to_csv(temporary, sep="\t", index=False)
        os.replace(temporary, path)
