"""Atomic publication, manifests, reports, and approval-ledger management."""

from __future__ import annotations

import hashlib
import os
import shutil
import tempfile
from datetime import datetime, timezone
from importlib import metadata
from pathlib import Path

import pandas as pd


def _utc_now() -> str:
    """Return an ISO-8601 UTC timestamp for machine-readable provenance."""

    return datetime.now(timezone.utc).isoformat()


def _sha256(path: Path) -> str:
    """Calculate a streaming SHA-256 digest for one published file."""

    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


class AtomicRunPublisher:
    """Stage scientific outputs and expose them only after complete validation."""

    def __init__(self, run_dir: Path):
        """Create a private staging directory inside the run filesystem."""

        self.run_dir = run_dir
        self.stage_dir = Path(tempfile.mkdtemp(prefix=".step00_staging_", dir=run_dir))
        self.objects_dir = self.stage_dir / "objects"
        self.tables_dir = self.stage_dir / "tables"
        self.objects_dir.mkdir()
        self.tables_dir.mkdir()

    def publish(self) -> None:
        """Atomically rename staged objects, tables, and report into the run."""

        targets = {
            self.objects_dir: self.run_dir / "objects",
            self.tables_dir: self.run_dir / "tables",
            self.stage_dir / "STEP00_VALIDATION_REPORT.md": self.run_dir / "STEP00_VALIDATION_REPORT.md",
            self.stage_dir / "STEP_STATUS.tsv": self.run_dir / "STEP_STATUS.tsv",
        }
        for target in targets.values():
            if target.exists():
                raise FileExistsError(f"Refusing to replace published Step 00 asset: {target}")
        for source, target in targets.items():
            os.replace(source, target)
        self.stage_dir.rmdir()

    def discard(self) -> None:
        """Remove unpublished staging files after a failed computation."""

        if self.stage_dir.exists():
            shutil.rmtree(self.stage_dir)


class ProvenancePublisher:
    """Write environment, manifest, report, and approval-state artifacts."""

    PACKAGES = ("anndata", "h5py", "numpy", "pandas", "scanpy", "scipy")

    def software_versions(self) -> pd.DataFrame:
        """Return exact Python-package versions used by the frozen run."""

        rows = [{"component": "python", "version": os.sys.version.replace("\n", " ")}]
        rows.extend({"component": package, "version": metadata.version(package)} for package in self.PACKAGES)
        return pd.DataFrame(rows)

    def output_manifest(self, stage_dir: Path) -> pd.DataFrame:
        """Manifest every staged output except the manifest itself."""

        rows = []
        for path in sorted(item for item in stage_dir.rglob("*") if item.is_file()):
            if path.name == "output_manifest.tsv":
                continue
            rows.append({"relative_path": str(path.relative_to(stage_dir)), "bytes": path.stat().st_size, "sha256": _sha256(path)})
        return pd.DataFrame(rows)

    def status_frame(self, run_id: str, checkpoint: str, n_cells: int, n_genes: int) -> pd.DataFrame:
        """Create the run-local IN_REVIEW status record required by the handoff."""

        return pd.DataFrame([{"step": "00_input_validation_and_canonical_anndata", "run_id": run_id, "status": "IN_REVIEW", "output_checkpoint": checkpoint, "completed_utc": _utc_now(), "cells_after": n_cells, "genes_after": n_genes, "approval_decision": "", "approved_run_id": ""}])

    def report(
        self,
        run_id: str,
        checkpoint: str,
        n_cells: int,
        n_genes: int,
        total_nnz: int,
        sample_summary: pd.DataFrame,
        checks: pd.DataFrame,
    ) -> str:
        """Render the concise human review report for Step 00."""

        failures = int((checks["status"] != "PASS").sum())
        lines = [
            "# Step 00 validation report",
            "",
            f"- Run: `{run_id}`",
            "- Status: **IN_REVIEW**",
            f"- Canonical checkpoint: `{checkpoint}`",
            f"- Dimensions: **{n_cells:,} cells × {n_genes:,} genes**",
            f"- Nonzero count entries: **{total_nnz:,}**",
            f"- Validation checks: **{len(checks) - failures} PASS, {failures} FAIL**",
            "- Matrix state: sparse, unnormalized integer Cell Ranger counts in `.X`.",
            "- Reductions/graphs/annotations: none created.",
            "- Raw matrices: registered for possible Step 04 use; uncalled droplets excluded.",
            "",
            "## Sample composition",
            "",
            "| Sample | Genotype | Sex | Design group | Cells | Genes | Nonzero entries |",
            "| --- | --- | --- | --- | ---: | ---: | ---: |",
        ]
        for row in sample_summary.itertuples(index=False):
            lines.append(f"| `{row.technical_sample_id}` | {row.genotype} | {row.sex} | `{row.design_group}` | {int(row.n_cells):,} | {int(row.n_genes):,} | {int(row.nnz):,} |")
        lines.extend([
            "",
            "## Review boundary",
            "",
            "This computation is not approved merely because it succeeded. Review the",
            "tables and checkpoint summary, request any additions or reruns, and provide",
            "an explicit YES before Step 00 can become the input to Step 01.",
            "",
        ])
        return "\n".join(lines)


class ApprovalLedger:
    """Maintain the cross-run primary-processing approval ledger atomically."""

    COLUMNS = (
        "step",
        "run_id",
        "status",
        "input_checkpoint",
        "output_checkpoint",
        "date_run_timestamp",
        "code_config_version",
        "summary",
        "cells_before",
        "cells_after",
        "genes_before",
        "genes_after",
        "important_parameters",
        "important_findings",
        "additional_analyses_requested",
        "reruns_modifications",
        "outstanding_questions",
        "approval_decision",
        "approved_run_id",
    )

    def update_in_review(
        self,
        path: Path,
        run_id: str,
        checkpoint: str,
        n_cells: int,
        n_genes: int,
        code_version: str,
    ) -> None:
        """Insert or replace one Step 00 run as IN_REVIEW without approving it."""

        if path.exists():
            ledger = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
        else:
            ledger = pd.DataFrame(columns=self.COLUMNS)
        ledger = ledger.loc[ledger.get("run_id", pd.Series(dtype=str)) != run_id].copy()
        row = {column: "" for column in self.COLUMNS}
        row.update({
            "step": "00_input_validation_and_canonical_anndata",
            "run_id": run_id,
            "status": "IN_REVIEW",
            "input_checkpoint": "12 per-sample Cell Ranger filtered feature-barcode matrices",
            "output_checkpoint": checkpoint,
            "date_run_timestamp": _utc_now(),
            "code_config_version": code_version,
            "summary": "Canonical raw-count AnnData constructed and structurally validated; awaiting user review.",
            "cells_before": str(n_cells),
            "cells_after": str(n_cells),
            "genes_before": str(n_genes),
            "genes_after": str(n_genes),
            "important_parameters": "Gene Expression only; .X raw integer counts; no filtering or transformation",
            "important_findings": "See run STEP00_VALIDATION_REPORT.md and validation_checks.tsv",
            "outstanding_questions": "User approval of exact Step 00 run",
        })
        ledger = pd.concat([ledger, pd.DataFrame([row])], ignore_index=True)
        path.parent.mkdir(parents=True, exist_ok=True)
        temporary = path.with_suffix(path.suffix + ".tmp")
        ledger.loc[:, list(self.COLUMNS)].to_csv(temporary, sep="\t", index=False)
        os.replace(temporary, path)
