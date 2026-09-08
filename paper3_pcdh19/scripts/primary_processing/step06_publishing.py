"""Stage and publish Step 06 outputs with provenance and review status.

``AtomicStep06Publisher`` owns temporary objects/tables/figures directories.
``Step06ProvenancePublisher`` builds versions, hashes, status, and report text.
``Step06ApprovalLedger`` updates the shared TSV only after outputs publish.
All successful computations are IN_REVIEW, with approval fields left blank.

Each rename is atomic on the filesystem; the full collection and ledger are
not a single transaction. Staging cleanup cannot roll back assets already
renamed. These classes manage files, not scientific tuning.
"""

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
    """Return a timezone-aware UTC timestamp.

    Returns
    -------
    str
        Timezone-aware UTC ISO-8601 timestamp.

    Notes
    -----
    Timestamp only; no file access.
    """

    return datetime.now(timezone.utc).isoformat()


class AtomicStep06Publisher:
    """Stage outputs and publish validated groups with individual atomic renames.

    Notes
    -----
    Own a unique staging directory beneath the new run. The caller writes and
    validates assets there before publish() renames each group. discard() cleans
    only remaining staging; it cannot undo already-published groups after a partial
    failure.
    """

    def __init__(self, run_dir: Path):
        """Create run-scoped staging and publication directories.

        Parameters
        ----------
        run_dir : pathlib.Path
            Run-scoped directory; never the approved upstream checkpoint directory.

        Notes
        -----
        Create a unique .step06_staging_* child of existing run_dir with objects, tables and
        figures directories. Filesystem errors propagate.
        """

        self.run_dir = run_dir
        self.stage_dir = Path(tempfile.mkdtemp(prefix=".step06_staging_", dir=run_dir))
        self.objects_dir = self.stage_dir / "objects"
        self.tables_dir = self.stage_dir / "tables"
        self.figures_dir = self.stage_dir / "figures"
        for directory in (self.objects_dir, self.tables_dir, self.figures_dir):
            directory.mkdir()

    def publish(self) -> None:
        """Atomically expose all validated output groups and status files.

        Returns
        -------
        None
            Effects are described below.

        Raises
        ------
        FileExistsError or RuntimeError
            A target already exists or staging contains unexpected top-level entries.

        Notes
        -----
        Quarantine .DS_Store, reject existing destination groups and unexpected
        staging entries, then rename objects/tables/figures/report/status one by one.
        Each os.replace is atomic; the collection is NOT a single transaction.
        Refuses overwrite of published assets. Does not itself validate scientific
        results: the workflow must pass its checks before calling this method.
        """

        self._quarantine_incidental_files()
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
        expected_sources = set(targets)
        unexpected = sorted(
            str(path.relative_to(self.stage_dir))
            for path in self.stage_dir.iterdir()
            if path not in expected_sources
        )
        if unexpected:
            raise RuntimeError(
                f"Refusing partial publication with unexpected staging entries: {unexpected}"
            )
        # Atomic per group, not one transaction. A later rename failure can leave earlier groups published.
        for source, target in targets.items():
            os.replace(source, target)
        self.stage_dir.rmdir()

    def _quarantine_incidental_files(self) -> None:
        """Preserve browser metadata outside staging before publication.

        Returns
        -------
        None
            Effects are described below.

        Notes
        -----
        Move .DS_Store files from staging into run provenance/incidental_files with encoded
        relative names. They are excluded from the scientific output manifest.
        """

        destination = self.run_dir / "provenance" / "incidental_files"
        for index, path in enumerate(sorted(self.stage_dir.rglob(".DS_Store")), start=1):
            destination.mkdir(parents=True, exist_ok=True)
            relative = str(path.relative_to(self.stage_dir)).replace("/", "__")
            target = destination / f"{index:03d}__{relative}"
            os.replace(path, target)

    def discard(self) -> None:
        """Remove only unpublished staging outputs after an error.

        Returns
        -------
        None
            Effects are described below.

        Notes
        -----
        Recursively delete only this publisher's remaining staging directory. It cannot undo
        already-renamed output groups or a ledger update.
        """

        if self.stage_dir.exists():
            shutil.rmtree(self.stage_dir)


class Step06ProvenancePublisher:
    """Create software, status, manifest, and human-readable report assets.

    Notes
    -----
    Build metadata and report values without recomputing scientific results. The
    workflow writes returned tables/text. Manifest hashes cover staged scientific
    files, while code/configuration/live event provenance are retained separately in
    the run.
    """

    packages = (
        "anndata", "h5py", "igraph", "leidenalg", "matplotlib", "numpy",
        "pandas", "scanpy", "scikit-learn", "scikit-misc", "scipy",
        "seaborn", "umap-learn",
    )

    def software_versions(self) -> pd.DataFrame:
        """Return exact interpreter and package versions.

        Returns
        -------
        pandas.DataFrame
            component/version rows for Python and the declared package list.

        Notes
        -----
        Read installed distribution metadata; missing required packages raise
        PackageNotFoundError. No package installation or environment mutation occurs.
        """

        rows = [{"component": "python", "version": os.sys.version.replace("\n", " ")}]
        rows.extend({"component": name, "version": metadata.version(name)} for name in self.packages)
        return pd.DataFrame(rows)

    def output_manifest(self, stage_dir: Path) -> pd.DataFrame:
        """Hash every staged file except the manifest itself.

        Parameters
        ----------
        stage_dir : pathlib.Path
            Unpublished output tree used for relative paths and content hashes.

        Returns
        -------
        pandas.DataFrame
            Relative path, bytes and SHA-256 for each staged scientific file.

        Notes
        -----
        Read/hash staged files; exclude output_manifest.tsv itself and .DS_Store. Frozen code,
        config and live progress outside staging are not included in this manifest.
        """

        rows = []
        for path in sorted(item for item in stage_dir.rglob("*") if item.is_file()):
            if path.name == "output_manifest.tsv":
                continue
            if path.name == ".DS_Store":
                continue
            rows.append({
                "relative_path": str(path.relative_to(stage_dir)),
                "bytes": path.stat().st_size,
                "sha256": sha256(path),
            })
        return pd.DataFrame(rows)

    def status_frame(self, run_id: str, input_run_id: str, cells: int, genes: int) -> pd.DataFrame:
        """Return the run-local IN_REVIEW status record.

        Parameters
        ----------
        run_id : str
            Identity of the Step 06 output run.
        input_run_id : str
            Identity of the approved Step 02 input run.
        cells : int
            Retained cell count, recorded identically before and after diagnostics.
        genes : int
            Retained gene count, recorded identically before and after diagnostics.

        Returns
        -------
        pandas.DataFrame
            One-row run-local status with IN_REVIEW and empty approval fields.

        Notes
        -----
        Build an in-memory record only; caller writes STEP_STATUS.tsv. Successful computation
        does not imply approval.
        """

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
        """Render the Step 06 methods, findings, and review boundary.

        Parameters
        ----------
        run_id : str
            Identity of the Step 06 output run.
        input_run_id : str
            Identity of the approved Step 02 input run.
        results : Step06Results
            Companion tables plus the automatic, unreviewed diagnostic outcome.
        checks : pandas.DataFrame
            Validation ledger with a status column containing PASS or FAIL.
        settings : Step06Settings
            Resolved scientific/rendering controls; see the settings class and tuning
            guide.

        Returns
        -------
        str
            Markdown methods, validation count, automatic outcome and review boundary.

        Notes
        -----
        Format existing results without recomputation. The headline is copied from _outcome
        and remains heuristic. This concise report is not a complete parameter inventory;
        consult frozen code/config and the tuning guide.
        """

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
    """Append a completed Step 06 run without implying approval.

    Notes
    -----
    Own the final shared-ledger update after scientific assets are published.
    Preserve the approved upstream lineage and record this run as IN_REVIEW. The
    update uses a temporary-file rename but does not provide concurrent-writer
    locking.
    """

    def update(self, path: Path, run_id: str, input_run_id: str, settings, results, code_version: str) -> None:
        """Insert one IN_REVIEW row atomically after verifying lineage.

        Parameters
        ----------
        path : pathlib.Path
            File to read or update, as described in Notes.
        run_id : str
            Identity of the Step 06 output run.
        input_run_id : str
            Identity of the approved Step 02 input run.
        settings : Step06Settings
            Resolved scientific/rendering controls; see the settings class and tuning
            guide.
        results : Step06Results
            Companion tables plus the automatic, unreviewed diagnostic outcome.
        code_version : str
            Frozen repository identity to store verbatim in the approval ledger.

        Returns
        -------
        None
            Effects are described below.

        Raises
        ------
        ValueError
            The input run is not uniquely approved in the live ledger.

        Notes
        -----
        Read the shared ledger at path; require exactly one approved input-run row.
        Replace any row for run_id with a completed IN_REVIEW record, write a .tmp
        sibling, then atomically replace the ledger. This is not a concurrent-writer
        locking protocol. It neither grants approval nor starts a downstream step.
        """

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
