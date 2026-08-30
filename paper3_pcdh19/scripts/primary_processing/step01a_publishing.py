"""Atomic publication, reporting, manifests, and ledger updates for Step 01a."""

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


class AtomicStep01aPublisher:
    """Keep Step 01a tables and figures private until validation passes."""

    def __init__(self, run_dir: Path):
        """Create a private staging directory on the run filesystem."""

        self.run_dir = run_dir
        self.stage_dir = Path(tempfile.mkdtemp(prefix=".step01a_staging_", dir=run_dir))
        self.tables_dir = self.stage_dir / "tables"
        self.figures_dir = self.stage_dir / "figures"
        self.tables_dir.mkdir()
        self.figures_dir.mkdir()

    def publish(self) -> None:
        """Atomically expose validated tables, figures, report, and status."""

        targets = {
            self.tables_dir: self.run_dir / "tables",
            self.figures_dir: self.run_dir / "figures",
            self.stage_dir / "STEP01A_MAD_SENSITIVITY_REPORT.md": self.run_dir / "STEP01A_MAD_SENSITIVITY_REPORT.md",
            self.stage_dir / "STEP_STATUS.tsv": self.run_dir / "STEP_STATUS.tsv",
        }
        for target in targets.values():
            if target.exists():
                raise FileExistsError(f"Refusing to replace published Step 01a asset: {target}")
        for source, target in targets.items():
            os.replace(source, target)
        self.stage_dir.rmdir()

    def discard(self) -> None:
        """Delete only unpublished staging outputs after a failure."""

        if self.stage_dir.exists():
            shutil.rmtree(self.stage_dir)


class Step01aProvenancePublisher:
    """Create environment, manifest, report, status, and workflow-ledger assets."""

    PACKAGES = ("anndata", "h5py", "matplotlib", "numpy", "pandas", "scanpy", "scipy")

    def software_versions(self) -> pd.DataFrame:
        """Return exact interpreter and dependency versions used by the job."""

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

    def status_frame(self, run_id: str, input_run_id: str, n_cells: int, n_genes: int) -> pd.DataFrame:
        """Create the run-local IN_REVIEW status without implying filtering approval."""

        return pd.DataFrame([{
            "step": "01a_qc_mad_sensitivity",
            "run_id": run_id,
            "status": "IN_REVIEW",
            "input_checkpoint": input_run_id,
            "output_checkpoint": "tables/per_cell_mad_candidate_flags.tsv.gz",
            "completed_utc": _utc_now(),
            "cells_before": n_cells,
            "cells_after": n_cells,
            "genes_before": n_genes,
            "genes_after": n_genes,
            "approval_decision": "",
            "approved_run_id": "",
        }])

    def report(self, run_id: str, input_run_id: str, summary: pd.DataFrame, checks: pd.DataFrame, plot_count: int) -> str:
        """Render the human-facing sensitivity report from computed summaries."""

        failures = int((checks["status"] != "PASS").sum())
        lines = [
            "# Step 01a per-sample MAD sensitivity report",
            "",
            f"- Run: `{run_id}`",
            "- Status: **IN_REVIEW**",
            f"- Exact Step 01 input: `{input_run_id}`",
            "- Scope: technical-outlier sensitivity only; **no cells removed**.",
            "- Boundary groups: 12 technical samples independently; no design/genotype grouping.",
            "- Stringencies: 3, 4, and 5 scaled MAD.",
            "- Scaled MAD: `1.4826 × median(abs(x - median(x)))`.",
            "- Low counts/genes use `log1p` values and are back-transformed to original units.",
            "- High mitochondrial percentage uses the observed percentage scale.",
            "- Upper counts/genes: not calculated; high-complexity assessment is reserved for Step 03 scDblFinder.",
            f"- Validation checks: **{len(checks) - failures} PASS, {failures} FAIL**.",
            f"- Diagnostic figure files: **{plot_count}**.",
            "",
            "## Candidate percentages by technical sample",
            "",
            "| Sample | MAD | Low counts | Low genes | High mt | Any | All three |",
            "| --- | ---: | ---: | ---: | ---: | ---: | ---: |",
        ]
        for row in summary.itertuples(index=False):
            lines.append(
                f"| `{row.technical_sample_id}` | {row.stringency_mad:g} | "
                f"{row.low_total_counts_n:,} ({row.low_total_counts_pct:.2f}%) | "
                f"{row.low_n_genes_n:,} ({row.low_n_genes_pct:.2f}%) | "
                f"{row.high_pct_counts_mt_n:,} ({row.high_pct_counts_mt_pct:.2f}%) | "
                f"{row.any_candidate_n:,} ({row.any_candidate_pct:.2f}%) | "
                f"{row.all_three_candidate_n:,} ({row.all_three_candidate_pct:.2f}%) |"
            )
        lines.extend([
            "",
            "## Review boundary",
            "",
            "These are candidate sensitivity flags, not QC calls. The input H5AD was",
            "opened read-only and no replacement H5AD was created. Review the boundary",
            "plots and tables, request further Step 01 diagnostics if needed, and stop",
            "before Step 02 until the combined Step 01/01a evidence is explicitly approved.",
            "",
        ])
        return "\n".join(lines)


class Step01aApprovalLedger:
    """Register Step 01a as an in-review amendment without approving Step 01."""

    def update_in_review(self, path: Path, run_id: str, input_run_id: str, n_cells: int, n_genes: int, code_version: str) -> None:
        """Insert or replace one Step 01a row while preserving all prior states."""

        ledger = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
        upstream = ledger.loc[(ledger["run_id"] == input_run_id) & (ledger["status"] == "IN_REVIEW")]
        if len(upstream) != 1:
            raise ValueError("Approval ledger does not contain exactly one in-review Step 01 input")
        ledger = ledger.loc[ledger["run_id"] != run_id].copy()
        row = {column: "" for column in ledger.columns}
        row.update({
            "step": "01a_qc_mad_sensitivity",
            "run_id": run_id,
            "status": "IN_REVIEW",
            "input_checkpoint": input_run_id,
            "output_checkpoint": "tables/per_cell_mad_candidate_flags.tsv.gz",
            "date_run_timestamp": _utc_now(),
            "code_config_version": code_version,
            "summary": "Per-sample 3/4/5 scaled-MAD candidate sensitivity analysis; no filtering.",
            "cells_before": str(n_cells),
            "cells_after": str(n_cells),
            "genes_before": str(n_genes),
            "genes_after": str(n_genes),
            "important_parameters": "per technical sample; low log1p counts; low log1p genes; high mt percentage; scaled MAD 1.4826; 3/4/5 MAD",
            "important_findings": "See STEP01A_MAD_SENSITIVITY_REPORT.md, candidate summaries, and boundary figures.",
            "outstanding_questions": "User review and explicit approval of combined Step 01/01a evidence",
        })
        ledger = pd.concat([ledger, pd.DataFrame([row])], ignore_index=True)
        temporary = path.with_suffix(path.suffix + ".tmp")
        ledger.to_csv(temporary, sep="\t", index=False)
        os.replace(temporary, path)
