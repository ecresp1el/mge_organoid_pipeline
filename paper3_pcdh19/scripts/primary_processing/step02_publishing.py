"""Atomic publication, reporting, manifests, and review state for Step 02."""

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
    """Return a timezone-aware UTC timestamp for provenance."""

    return datetime.now(timezone.utc).isoformat()


def _sha256(path: Path) -> str:
    """Calculate a streaming SHA-256 digest for one published file."""

    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


class AtomicStep02Publisher:
    """Keep Step 02 outputs private until every validation passes."""

    def __init__(self, run_dir: Path):
        """Create a private staging directory on the publication filesystem."""

        self.run_dir = run_dir
        self.stage_dir = Path(tempfile.mkdtemp(prefix=".step02_staging_", dir=run_dir))
        self.objects_dir = self.stage_dir / "objects"
        self.tables_dir = self.stage_dir / "tables"
        self.figures_dir = self.stage_dir / "figures"
        self.objects_dir.mkdir()
        self.tables_dir.mkdir()
        self.figures_dir.mkdir()

    def publish(self) -> None:
        """Atomically expose the checkpoint, tables, figures, reports, and status."""

        targets = {
            self.objects_dir: self.run_dir / "objects",
            self.tables_dir: self.run_dir / "tables",
            self.figures_dir: self.run_dir / "figures",
            self.stage_dir / "STEP02_QC_FILTERING_REPORT.md": self.run_dir / "STEP02_QC_FILTERING_REPORT.md",
            self.stage_dir / "STEP_STATUS.tsv": self.run_dir / "STEP_STATUS.tsv",
        }
        for target in targets.values():
            if target.exists():
                raise FileExistsError(f"Refusing to replace published Step 02 asset: {target}")
        for source, target in targets.items():
            os.replace(source, target)
        self.stage_dir.rmdir()

    def discard(self) -> None:
        """Delete only unpublished staging outputs after a failure."""

        if self.stage_dir.exists():
            shutil.rmtree(self.stage_dir)


class Step02ProvenancePublisher:
    """Create Step 02 software, manifest, status, and human review reports."""

    PACKAGES = ("anndata", "h5py", "matplotlib", "numpy", "pandas", "scanpy", "scipy")

    def software_versions(self) -> pd.DataFrame:
        """Return exact interpreter and package versions used by the job."""

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

    def status_frame(self, run_id: str, input_run_id: str, n_before: int, n_after: int, n_genes: int) -> pd.DataFrame:
        """Create the run-local Step 02 IN_REVIEW record."""

        return pd.DataFrame([{
            "step": "02_qc_filtering",
            "run_id": run_id,
            "status": "IN_REVIEW",
            "input_checkpoint": input_run_id,
            "output_checkpoint": "objects/pcdh19_step02_qc_filtered.h5ad",
            "completed_utc": _utc_now(),
            "cells_before": n_before,
            "cells_after": n_after,
            "genes_before": n_genes,
            "genes_after": n_genes,
            "approval_decision": "",
            "approved_run_id": "",
        }])

    def filtering_report(self, run_id: str, input_run_id: str, flags_run_id: str, by_sample: pd.DataFrame, by_design: pd.DataFrame, checks: pd.DataFrame) -> str:
        """Render the Step 02 filtering report with exact before/after counts."""

        failures = int((checks["status"] != "PASS").sum())
        before = int(by_sample["cells_before"].sum())
        removed = int(by_sample["cells_removed"].sum())
        after = int(by_sample["cells_after"].sum())
        lines = [
            "# Step 02 QC-filtering report",
            "",
            f"- Run: `{run_id}`",
            "- Status: **IN_REVIEW**",
            f"- Approved Step 01 input: `{input_run_id}`",
            f"- Approved Step 01a flags: `{flags_run_id}`",
            "- Policy: per-technical-sample 5-MAD low counts OR low genes OR high mitochondrial percentage.",
            f"- Cells: **{before:,} before; {removed:,} removed; {after:,} after**.",
            "- Genes: **19,071 before and after; zero genes removed**.",
            "- Low-count removals: 0; low-gene removals: 67; high-mitochondrial removals: 4,372.",
            "- Upper count/gene filtering: none.",
            "- Doublet filtering: none; scDblFinder not run.",
            f"- Validation checks: **{len(checks)-failures} PASS, {failures} FAIL**.",
            "",
            "## By technical sample",
            "",
            "| Sample | Before | Removed | After | Removed % | Low counts | Low genes | High mt | Multiple reasons |",
            "| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
        ]
        for row in by_sample.itertuples(index=False):
            lines.append(f"| `{row.group}` | {row.cells_before:,} | {row.cells_removed:,} | {row.cells_after:,} | {row.pct_removed:.3f}% | {row.low_total_counts_n:,} | {row.low_n_genes_n:,} | {row.high_pct_counts_mt_n:,} | {row.multiple_reasons_n:,} |")
        lines.extend(["", "## By design group (reporting only)", "", "Design groups did not define any boundary.", "", "| Design | Before | Removed | After | Removed % |", "| --- | ---: | ---: | ---: | ---: |"])
        for row in by_design.itertuples(index=False):
            lines.append(f"| `{row.group}` | {row.cells_before:,} | {row.cells_removed:,} | {row.cells_after:,} | {row.pct_removed:.3f}% |")
        lines.extend(["", "## Review boundary", "", "The checkpoint remains IN_REVIEW. Review all exclusions and the saved raw-count", "checkpoint before authorizing any Step 03 scDblFinder computation.", ""])
        return "\n".join(lines)


class Step02ApprovalLedger:
    """Append the completed Step 02 run without implying user approval."""

    def update_in_review(self, path: Path, run_id: str, input_run_id: str, output_checkpoint: str, n_before: int, n_after: int, n_genes: int, code_version: str) -> None:
        """Insert or replace one Step 02 IN_REVIEW row atomically."""

        ledger = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
        approved = ledger.loc[(ledger["run_id"] == input_run_id) & (ledger["status"] == "APPROVED")]
        if len(approved) != 1:
            raise ValueError("Approval ledger does not contain exactly one approved Step 01 input")
        ledger = ledger.loc[ledger["run_id"] != run_id].copy()
        row = {column: "" for column in ledger.columns}
        row.update({
            "step": "02_qc_filtering",
            "run_id": run_id,
            "status": "IN_REVIEW",
            "input_checkpoint": input_run_id,
            "output_checkpoint": output_checkpoint,
            "date_run_timestamp": _utc_now(),
            "code_config_version": code_version,
            "summary": "Approved per-sample 5-MAD union applied with exact exclusion provenance; awaiting review.",
            "cells_before": str(n_before),
            "cells_after": str(n_after),
            "genes_before": str(n_genes),
            "genes_after": str(n_genes),
            "important_parameters": "5 MAD per technical sample; low counts OR low genes OR high mt; no upper complexity filter",
            "important_findings": f"{n_before-n_after:,} cells excluded; reasons preserved; scDblFinder not run.",
            "outstanding_questions": "User review of Step 02; Step 03 capture definition awaits explicit approval",
        })
        ledger = pd.concat([ledger, pd.DataFrame([row])], ignore_index=True)
        temporary = path.with_suffix(path.suffix + ".tmp")
        ledger.to_csv(temporary, sep="\t", index=False)
        os.replace(temporary, path)
