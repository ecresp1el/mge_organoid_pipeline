"""Atomic publication, manifests, reports, and review state for Step 03."""

from __future__ import annotations

import hashlib
import os
import shutil
import tempfile
from datetime import datetime, timezone
from importlib import metadata
from pathlib import Path

import pandas as pd


def utc_now() -> str:
    """Return one timezone-aware UTC timestamp."""

    return datetime.now(timezone.utc).isoformat()


def sha256_file(path: Path) -> str:
    """Calculate one published file's streaming SHA-256 digest."""

    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


class AtomicStep03Publisher:
    """Keep all scientific Step 03 outputs private until validation passes."""

    def __init__(self, run_dir: Path):
        """Create one private staging directory on the destination filesystem."""

        self.run_dir = run_dir
        self.stage_dir = Path(tempfile.mkdtemp(prefix=".step03_staging_", dir=run_dir))
        self.objects_dir = self.stage_dir / "objects"
        self.tables_dir = self.stage_dir / "tables"
        self.figures_dir = self.stage_dir / "figures"
        self.objects_dir.mkdir()
        self.tables_dir.mkdir()
        self.figures_dir.mkdir()

    def publish(self) -> None:
        """Atomically expose the complete validated Step 03 review package."""

        targets = {
            self.objects_dir: self.run_dir / "objects",
            self.tables_dir: self.run_dir / "tables",
            self.figures_dir: self.run_dir / "figures",
            self.stage_dir / "STEP03_SCDBLFINDER_REPORT.md": self.run_dir / "STEP03_SCDBLFINDER_REPORT.md",
            self.stage_dir / "STEP_STATUS.tsv": self.run_dir / "STEP_STATUS.tsv",
        }
        for target in targets.values():
            if target.exists():
                raise FileExistsError(f"Refusing to replace published Step 03 asset: {target}")
        for source, target in targets.items():
            os.replace(source, target)
        self.stage_dir.rmdir()

    def discard(self) -> None:
        """Remove only unpublished staged assets after failure."""

        if self.stage_dir.exists():
            shutil.rmtree(self.stage_dir)


class Step03ReportBuilder:
    """Build exact summaries, object-state documentation, and the review report."""

    @staticmethod
    def summaries(frame: pd.DataFrame, sample_field: str, design_field: str) -> dict[str, pd.DataFrame]:
        """Create overall, sample, design, cluster, score, and concordance tables."""

        def grouped(field: str) -> pd.DataFrame:
            """Summarize primary calls and scores for one metadata field."""

            table = frame.groupby(field, observed=True).agg(
                cells=("cell_id", "size"),
                called_doublets=("primary_class", lambda x: int((x == "doublet").sum())),
                median_score=("primary_score", "median"),
                mean_score=("primary_score", "mean"),
            ).reset_index()
            table["called_doublet_pct"] = 100 * table["called_doublets"] / table["cells"]
            return table

        primary_doublets = int((frame["primary_class"] == "doublet").sum())
        replicate_doublets = int((frame["replicate_class"] == "doublet").sum())
        agreement = int((frame["primary_class"] == frame["replicate_class"]).sum())
        overall = pd.DataFrame([{
            "cells": len(frame),
            "primary_called_doublets": primary_doublets,
            "primary_called_doublet_pct": 100 * primary_doublets / len(frame),
            "replicate_called_doublets": replicate_doublets,
            "replicate_called_doublet_pct": 100 * replicate_doublets / len(frame),
            "class_agreement_n": agreement,
            "class_agreement_pct": 100 * agreement / len(frame),
            "score_pearson": float(np_corr(frame["primary_score"], frame["replicate_score"])),
            "empirical_lowest_doublet_score": float(frame.loc[frame["primary_class"] == "doublet", "primary_score"].min()),
            "cells_removed": 0,
        }])
        cluster = grouped("primary_cluster").sort_values("primary_cluster")
        quantiles = frame.groupby(sample_field, observed=True)["primary_score"].quantile([0, 0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99, 1]).unstack().reset_index()
        quantiles.columns = [sample_field, "q000", "q001", "q005", "q025", "q050", "q075", "q095", "q099", "q100"]
        confusion = pd.crosstab(frame["primary_class"], frame["replicate_class"]).rename_axis("primary_class").reset_index()
        reproducibility_sample = Step03ReportBuilder._reproducibility_by_group(frame, sample_field)
        reproducibility_design = Step03ReportBuilder._reproducibility_by_group(frame, design_field)
        return {
            "overall": overall,
            "sample": grouped(sample_field),
            "design": grouped(design_field),
            "cluster": cluster,
            "score_quantiles": quantiles,
            "confusion": confusion,
            "reproducibility_sample": reproducibility_sample,
            "reproducibility_design": reproducibility_design,
        }

    @staticmethod
    def _reproducibility_by_group(frame: pd.DataFrame, field: str) -> pd.DataFrame:
        """Summarize concordant and discordant calls for one metadata grouping."""

        rows = []
        for group, subset in frame.groupby(field, observed=True):
            primary_doublet = subset["primary_class"].eq("doublet")
            replicate_doublet = subset["replicate_class"].eq("doublet")
            agreement = primary_doublet.eq(replicate_doublet)
            rows.append({
                field: group,
                "cells": len(subset),
                "call_agreement_n": int(agreement.sum()),
                "call_agreement_pct": 100 * float(agreement.mean()),
                "both_doublet": int((primary_doublet & replicate_doublet).sum()),
                "primary_only_doublet": int((primary_doublet & ~replicate_doublet).sum()),
                "replicate_only_doublet": int((~primary_doublet & replicate_doublet).sum()),
                "both_singlet": int((~primary_doublet & ~replicate_doublet).sum()),
                "score_pearson": np_corr(subset["primary_score"], subset["replicate_score"]),
            })
        return pd.DataFrame(rows)

    @staticmethod
    def data_dictionary() -> pd.DataFrame:
        """Describe every canonical Step 03 checkpoint addition and object state."""

        rows = [
            ("pcdh19_step03_scdblfinder.h5ad", "446,349 cells x 19,071 genes; raw sparse integer counts remain in .X; no cells or genes removed", "IN_REVIEW checkpoint"),
            ("capture_id", "Constant GEX_1; the independently processed capture passed to scDblFinder samples", "approved"),
            ("scDblFinder_score", "Primary-seed scDblFinder score", "review-only; no filtering"),
            ("scDblFinder_class", "Primary-seed package call: singlet or doublet", "review-only; no filtering"),
            ("scDblFinder_cluster", "Fast cluster generated internally because clusters=TRUE", "diagnostic"),
            ("scDblFinder_replicate_score", "Second-seed score with otherwise identical model arguments", "reproducibility diagnostic"),
            ("scDblFinder_replicate_class", "Second-seed package call", "reproducibility diagnostic"),
            ("scDblFinder_call_reproduced", "Whether primary and second-seed package calls agree", "reproducibility diagnostic"),
            ("obsm['X_scdblfinder_pca']", "Real-cell coordinates from the primary run's internal normalized PCA; artificial doublets excluded", "diagnostic expression space, not an integrated UMAP"),
            ("layers/.raw/obsp", "Empty; no normalized expression layer, raw alias, or neighbor graph is saved", "object state"),
            ("scdblfinder_per_cell_details.tsv.gz", "Complete real-cell native-R output fields plus primary and repeat results", "audit evidence"),
            ("scdblfinder_primary_stats.rds", "Package-generated primary threshold/model statistics only; no expression matrix", "native-R audit object"),
            ("scdblfinder_internal_pca.rds", "Exact native-R real-cell internal PCA matrix", "native-R diagnostic object; no counts"),
        ]
        return pd.DataFrame(rows, columns=["asset_or_field", "definition", "decision_status"])

    @staticmethod
    def report(run_id: str, input_run_id: str, summaries: dict[str, pd.DataFrame], checks: pd.DataFrame, pca_dimensions: int) -> str:
        """Render the human-readable Step 03 review report and mandatory stop."""

        overall = summaries["overall"].iloc[0]
        failures = int((checks["status"] != "PASS").sum())
        lines = [
            "# Step 03 scDblFinder report",
            "",
            f"- Run ID: `{run_id}`.",
            f"- Approved Step 02 input: `{input_run_id}`.",
            "- Capture definition: one independently processed capture, `GEX_1`; the 12 Probe Barcode samples remain reporting metadata.",
            "- Scientific arguments: `samples='capture_id'`, `clusters=TRUE`, `dbr.sd=1`, no supplied `dbr`; other model parameters were package defaults.",
            "- Output-only return settings: primary `returnType='full'` to retain the exact internal PCA/model diagnostics; second seed `returnType='scores'` for reproducibility. These do not change classifier training or thresholding.",
            f"- Cells/genes retained: **{int(overall.cells):,} cells x 19,071 genes; zero removed**.",
            f"- Primary calls: **{int(overall.primary_called_doublets):,} doublets ({overall.primary_called_doublet_pct:.3f}%)**.",
            f"- Second-seed calls: **{int(overall.replicate_called_doublets):,} doublets ({overall.replicate_called_doublet_pct:.3f}%)**.",
            f"- Reproducibility: **{overall.class_agreement_pct:.3f}% call agreement**; score Pearson correlation **{overall.score_pearson:.5f}**.",
            f"- Preserved expression-space diagnostic: {pca_dimensions} internal PCA coordinates for every real cell; this is not a final normalized, integrated, or clustered analysis object.",
            f"- Validation: **{len(checks)-failures} PASS, {failures} FAIL**.",
            "",
            "## Primary calls by technical sample",
            "",
            "The detector was fit to one capture. These sample summaries are post hoc diagnostics; they did not define separate doublet-rate priors or thresholds.",
            "",
            "| Technical sample | Cells | Called doublets | Called % | Median score |",
            "| --- | ---: | ---: | ---: | ---: |",
        ]
        for row in summaries["sample"].itertuples(index=False):
            lines.append(f"| `{getattr(row, summaries['sample'].columns[0])}` | {row.cells:,} | {row.called_doublets:,} | {row.called_doublet_pct:.3f}% | {row.median_score:.5f} |")
        lines.extend([
            "",
            "## Review boundary",
            "",
            "This run remains **IN_REVIEW**. No predicted doublet was removed. Review the score distribution, called fraction, second-seed reproducibility, sample/design composition, generated-cluster composition, and internal-PCA localization before deciding whether any call should become an exclusion in a later approved step.",
            "",
            "The withdrawn Howitt preprint is not used or cited. Methods provenance is limited to Germain et al. and the current scDblFinder documentation.",
            "",
        ])
        return "\n".join(lines)


def np_corr(left: pd.Series, right: pd.Series) -> float:
    """Calculate a finite Pearson correlation without another dependency."""

    return float(left.astype(float).corr(right.astype(float), method="pearson"))


class Step03ProvenancePublisher:
    """Create Python versions, output manifests, and Step 03 review status."""

    PACKAGES = ("anndata", "h5py", "matplotlib", "numpy", "pandas", "scanpy", "scipy")

    def python_versions(self) -> pd.DataFrame:
        """Return exact Python interpreter and package versions used to finalize."""

        rows = [{"component": "python", "version": os.sys.version.replace("\n", " ")}]
        rows.extend({"component": package, "version": metadata.version(package)} for package in self.PACKAGES)
        return pd.DataFrame(rows)

    def output_manifest(self, stage_dir: Path) -> pd.DataFrame:
        """Manifest every staged output except the manifest itself."""

        rows = []
        for path in sorted(item for item in stage_dir.rglob("*") if item.is_file()):
            if path.name == "output_manifest.tsv":
                continue
            rows.append({"relative_path": str(path.relative_to(stage_dir)), "bytes": path.stat().st_size, "sha256": sha256_file(path)})
        return pd.DataFrame(rows)

    def status_frame(self, run_id: str, input_run_id: str, n_cells: int, n_genes: int) -> pd.DataFrame:
        """Create the run-local Step 03 IN_REVIEW record without exclusions."""

        return pd.DataFrame([{
            "step": "03_scdblfinder",
            "run_id": run_id,
            "status": "IN_REVIEW",
            "input_checkpoint": input_run_id,
            "output_checkpoint": "objects/pcdh19_step03_scdblfinder.h5ad",
            "completed_utc": utc_now(),
            "cells_before": n_cells,
            "cells_after": n_cells,
            "genes_before": n_genes,
            "genes_after": n_genes,
            "approval_decision": "",
            "approved_run_id": "",
        }])


class Step03ApprovalLedger:
    """Append one completed Step 03 run without implying user approval."""

    def update_in_review(self, path: Path, run_id: str, input_run_id: str, calls: int, call_pct: float, code_version: str) -> None:
        """Insert one atomic IN_REVIEW ledger row after successful publication."""

        ledger = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
        approved = ledger.loc[(ledger["run_id"] == input_run_id) & (ledger["status"] == "APPROVED") & (ledger["approved_run_id"] == input_run_id)]
        if len(approved) != 1:
            raise ValueError("Approval ledger does not contain exactly one approved Step 02 input")
        ledger = ledger.loc[ledger["run_id"] != run_id].copy()
        row = {column: "" for column in ledger.columns}
        row.update({
            "step": "03_scdblfinder",
            "run_id": run_id,
            "status": "IN_REVIEW",
            "input_checkpoint": input_run_id,
            "output_checkpoint": "objects/pcdh19_step03_scdblfinder.h5ad",
            "date_run_timestamp": utc_now(),
            "code_config_version": code_version,
            "summary": "scDblFinder detection on one approved GEX_1 capture; scores/calls retained; no cells removed.",
            "cells_before": "446349",
            "cells_after": "446349",
            "genes_before": "19071",
            "genes_after": "19071",
            "important_parameters": "samples=GEX_1 capture; clusters=TRUE; dbr.sd=1; dbr omitted; otherwise model defaults",
            "important_findings": f"{calls:,} primary called doublets ({call_pct:.3f}%); no exclusions; see reproducibility and PCA diagnostics.",
            "outstanding_questions": "User review of score separation, called fraction, reproducibility, sample composition, generated clusters, and internal-PCA localization",
        })
        temporary = path.with_suffix(path.suffix + ".tmp")
        pd.concat([ledger, pd.DataFrame([row])], ignore_index=True).to_csv(temporary, sep="\t", index=False)
        os.replace(temporary, path)
