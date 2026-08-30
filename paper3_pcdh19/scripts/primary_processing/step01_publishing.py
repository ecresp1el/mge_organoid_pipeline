"""Atomic publication, reporting, manifests, and review-ledger updates for Step 01."""

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
    """Return a timezone-aware UTC timestamp for Step 01 provenance."""

    return datetime.now(timezone.utc).isoformat()


def _sha256(path: Path) -> str:
    """Compute a streaming SHA-256 digest for one staged output."""

    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


class AtomicStep01Publisher:
    """Keep Step 01 scientific outputs private until every validation passes."""

    def __init__(self, run_dir: Path):
        """Create run-local staging directories on the publication filesystem."""

        self.run_dir = run_dir
        self.stage_dir = Path(tempfile.mkdtemp(prefix=".step01_staging_", dir=run_dir))
        self.objects_dir = self.stage_dir / "objects"
        self.tables_dir = self.stage_dir / "tables"
        self.figures_dir = self.stage_dir / "figures"
        self.objects_dir.mkdir()
        self.tables_dir.mkdir()
        self.figures_dir.mkdir()

    def publish(self) -> None:
        """Atomically move validated objects, tables, figures, report, and status."""

        targets = {
            self.objects_dir: self.run_dir / "objects",
            self.tables_dir: self.run_dir / "tables",
            self.figures_dir: self.run_dir / "figures",
            self.stage_dir / "STEP01_QC_METRICS_REPORT.md": self.run_dir / "STEP01_QC_METRICS_REPORT.md",
            self.stage_dir / "STEP_STATUS.tsv": self.run_dir / "STEP_STATUS.tsv",
        }
        for target in targets.values():
            if target.exists():
                raise FileExistsError(f"Refusing to replace published Step 01 asset: {target}")
        for source, target in targets.items():
            os.replace(source, target)
        self.stage_dir.rmdir()

    def discard(self) -> None:
        """Delete only unpublished staging data after a failed computation."""

        if self.stage_dir.exists():
            shutil.rmtree(self.stage_dir)


class Step01ProvenancePublisher:
    """Create Step 01 environment, report, status, and output-manifest assets."""

    PACKAGES = ("anndata", "h5py", "matplotlib", "numpy", "pandas", "scanpy", "scipy")

    def software_versions(self) -> pd.DataFrame:
        """Return exact interpreter and dependency versions used by the job."""

        rows = [{"component": "python", "version": os.sys.version.replace("\n", " ")}]
        rows.extend({"component": package, "version": metadata.version(package)} for package in self.PACKAGES)
        return pd.DataFrame(rows)

    def output_manifest(self, stage_dir: Path) -> pd.DataFrame:
        """Manifest every staged Step 01 output except the manifest itself."""

        rows = []
        for path in sorted(item for item in stage_dir.rglob("*") if item.is_file()):
            if path.name == "output_manifest.tsv":
                continue
            rows.append({"relative_path": str(path.relative_to(stage_dir)), "bytes": path.stat().st_size, "sha256": _sha256(path)})
        return pd.DataFrame(rows)

    def status_frame(self, run_id: str, input_checkpoint: str, output_checkpoint: str, n_cells: int, n_genes: int) -> pd.DataFrame:
        """Create the run-local IN_REVIEW record without implying approval."""

        return pd.DataFrame([{
            "step": "01_qc_metrics",
            "run_id": run_id,
            "status": "IN_REVIEW",
            "input_checkpoint": input_checkpoint,
            "output_checkpoint": output_checkpoint,
            "completed_utc": _utc_now(),
            "cells_before": n_cells,
            "cells_after": n_cells,
            "genes_before": n_genes,
            "genes_after": n_genes,
            "approval_decision": "",
            "approved_run_id": "",
        }])

    def report(
        self,
        run_id: str,
        input_run_id: str,
        output_checkpoint: str,
        n_genes: int,
        pooled: pd.DataFrame,
        sample_summary: pd.DataFrame,
        checks: pd.DataFrame,
        mitochondrial_genes: int,
        ribosomal_genes: int,
        plot_count: int,
    ) -> str:
        """Render the human-facing Step 01 review report from computed summaries."""

        row = pooled.iloc[0]
        failures = int((checks["status"] != "PASS").sum())
        lines = [
            "# Step 01 QC-metrics report",
            "",
            f"- Run: `{run_id}`",
            "- Status: **IN_REVIEW**",
            f"- Approved Step 00 input: `{input_run_id}`",
            f"- Output checkpoint: `{output_checkpoint}`",
            f"- Dimensions: **{int(row['n_cells']):,} cells × {n_genes:,} genes**",
            f"- Validation checks: **{len(checks) - failures} PASS, {failures} FAIL**",
            f"- Diagnostic figure files: **{plot_count}**",
            "- Calculation: `scanpy.pp.calculate_qc_metrics()` with no percent-top calculation.",
            "- Filtering/thresholds/exclusions: **none**.",
            "",
            "## Pooled descriptive values",
            "",
            f"- Total counts: median {row['total_counts_p50']:,.0f}; 5th–95th percentile {row['total_counts_p05']:,.0f}–{row['total_counts_p95']:,.0f}.",
            f"- Detected genes: median {row['n_genes_by_counts_p50']:,.0f}; 5th–95th percentile {row['n_genes_by_counts_p05']:,.0f}–{row['n_genes_by_counts_p95']:,.0f}.",
            f"- Mitochondrial percentage: median {row['pct_counts_mt_p50']:.3f}%; 5th–95th percentile {row['pct_counts_mt_p05']:.3f}–{row['pct_counts_mt_p95']:.3f}%.",
            "",
            "## QC gene sets",
            "",
            f"- Mitochondrial: {mitochondrial_genes} delivered genes with symbols beginning `mt-`.",
            f"- Cytosolic ribosomal: {ribosomal_genes} delivered genes beginning `Rpl`/`Rps`.",
        ]
        if ribosomal_genes == 0:
            lines.append("- Ribosomal count/fraction metrics were not fabricated because this Flex panel contains no Rpl/Rps features.")
        lines.extend([
            "",
            "## Per-sample cell counts",
            "",
            "| Sample | Cells | Median counts | Median genes | Median mt % |",
            "| --- | ---: | ---: | ---: | ---: |",
        ])
        for sample in sample_summary.itertuples(index=False):
            lines.append(f"| `{sample.group}` | {int(sample.n_cells):,} | {sample.total_counts_p50:,.0f} | {sample.n_genes_by_counts_p50:,.0f} | {sample.pct_counts_mt_p50:.3f} |")
        lines.extend([
            "",
            "## Review boundary",
            "",
            "These distributions describe the complete approved cell set. They do not",
            "define high/low quality, pass/fail status, or a filtering threshold. Request",
            "additional Step 01 diagnostics or explicitly approve this run before any",
            "Step 02 threshold proposal.",
            "",
        ])
        return "\n".join(lines)


class Step01ApprovalLedger:
    """Append the computed Step 01 run while preserving approved Step 00 state."""

    def update_in_review(
        self,
        path: Path,
        run_id: str,
        input_run_id: str,
        output_checkpoint: str,
        n_cells: int,
        n_genes: int,
        code_version: str,
    ) -> None:
        """Insert or replace one Step 01 IN_REVIEW row atomically."""

        ledger = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
        required_step00 = ledger.loc[
            (ledger["run_id"] == input_run_id) & (ledger["status"] == "APPROVED")
        ]
        if len(required_step00) != 1:
            raise ValueError("Approval ledger does not contain exactly one approved Step 00 input")
        ledger = ledger.loc[ledger["run_id"] != run_id].copy()
        row = {column: "" for column in ledger.columns}
        row.update({
            "step": "01_qc_metrics",
            "run_id": run_id,
            "status": "IN_REVIEW",
            "input_checkpoint": input_run_id,
            "output_checkpoint": output_checkpoint,
            "date_run_timestamp": _utc_now(),
            "code_config_version": code_version,
            "summary": "Scanpy QC metrics calculated without filtering; awaiting user review.",
            "cells_before": str(n_cells),
            "cells_after": str(n_cells),
            "genes_before": str(n_genes),
            "genes_after": str(n_genes),
            "important_parameters": "scanpy.pp.calculate_qc_metrics; mt prefix mt-; percent_top=None; no thresholds",
            "important_findings": "See STEP01_QC_METRICS_REPORT.md, QC summaries, and diagnostic figures.",
            "outstanding_questions": "User review and explicit approval of exact Step 01 run",
        })
        ledger = pd.concat([ledger, pd.DataFrame([row])], ignore_index=True)
        temporary = path.with_suffix(path.suffix + ".tmp")
        ledger.to_csv(temporary, sep="\t", index=False)
        os.replace(temporary, path)
