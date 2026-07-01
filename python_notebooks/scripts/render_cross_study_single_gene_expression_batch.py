#!/usr/bin/env python3
"""Render multiple single-gene cross-study expression panels into one package."""

from __future__ import annotations

import argparse
from datetime import datetime
from pathlib import Path
import shutil
import subprocess
import sys

import pandas as pd

from render_cross_study_lhx6_expression_final_panel import (
    OUTPUT_STEM_TEMPLATE,
    normalize_gene,
    safe_token,
)


PROJECT_ROOT_DEFAULT = Path("/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder")
RUN_LABEL_DEFAULT = "cross_study_marker_expression_v12_pv_precursors_final_candidate"
FINAL_FOLDER_DEFAULT = "fig_cross_study_single_gene_expression_v12_candidate"
DEFAULT_GENES = ["NKX2.1", "LHX6", "LHX8", "MAFB", "MEF2C", "CRABP1", "TAC1", "VIPR2", "PV"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, default=PROJECT_ROOT_DEFAULT)
    parser.add_argument("--run-label", default=RUN_LABEL_DEFAULT)
    parser.add_argument("--final-dir", type=Path, default=None)
    parser.add_argument("--genes", nargs="+", default=DEFAULT_GENES)
    parser.add_argument("--positive-threshold", type=float, default=0.5)
    parser.add_argument("--fresh", action="store_true", help="Delete the consolidated final-dir before rendering.")
    parser.add_argument("--skip-existing", action="store_true", help="Reuse existing per-gene outputs/tables when present.")
    return parser.parse_args()


def write_batch_readme(final_dir: Path, project_root: Path, run_label: str, genes: list[str], threshold: float) -> None:
    gene_labels = ", ".join(genes)
    text = f"""# Cross-study single-gene marker expression batch

Consolidated final-figure candidate package for repeated single-gene panels:
{gene_labels}.

Source prepared marker-expression run:
`{project_root}/results/cross_study_marker_expression/{run_label}`

Each gene panel uses the same three-row layout:
- Top row: per-study UMAPs colored for marker expression.
- Middle row: per-study violin distributions of marker log1p(CP10K) expression.
- Bottom row: sample-level percentage of cells with marker expression >= `{threshold:g}`.

Shared denominator:
Visualization-filtered cells from the marker-expression workflow. Samarasinghe
is controls-only, and DIV90 current clusters 6 and 7 are removed as stressed
cells for plotting and summaries.

Aliases:
- `NKX2.1` uses the marker table column `NKX2-1`.
- `PV` uses the marker table column `PVALB`.

Package contents:
- `figures/png|pdf|svg/`: one output per gene.
- `tables/`: per-gene study summaries, sample-positive fractions, filter
  summaries, output manifests, plus a combined study summary.
- `code/`: render scripts used for this batch.
- `provenance/`: batch manifest, git state, handoff snapshots, and checksums.
"""
    (final_dir / "README.md").write_text(text)


def main() -> None:
    args = parse_args()
    final_dir = args.final_dir or (args.project_root / "final_figures" / FINAL_FOLDER_DEFAULT)
    if args.fresh and final_dir.exists():
        shutil.rmtree(final_dir)
    for subdir in ("figures/png", "figures/pdf", "figures/svg", "tables", "logs", "provenance", "code"):
        (final_dir / subdir).mkdir(parents=True, exist_ok=True)

    script_dir = Path(__file__).resolve().parent
    single_script = script_dir / "render_cross_study_lhx6_expression_final_panel.py"
    rendered = []
    for gene_arg in args.genes:
        gene = normalize_gene(gene_arg)
        token = safe_token(gene)
        expected_outputs = [
            final_dir / "figures" / ext / f"{OUTPUT_STEM_TEMPLATE.format(gene_token=token)}.{ext}"
            for ext in ("png", "pdf", "svg")
        ]
        expected_table = final_dir / "tables" / f"cross_study_{token}_expression_study_summary.tsv"
        if args.skip_existing and expected_table.exists() and all(path.exists() and path.stat().st_size > 0 for path in expected_outputs):
            print(f"Skipping existing complete gene output: {gene}", flush=True)
            rendered.append(
                {
                    "requested_gene": gene_arg,
                    "canonical_gene": gene,
                    "gene_token": token,
                    "output_stem": OUTPUT_STEM_TEMPLATE.format(gene_token=token),
                }
            )
            continue
        cmd = [
            sys.executable,
            str(single_script),
            "--project-root",
            str(args.project_root),
            "--run-label",
            args.run_label,
            "--final-dir",
            str(final_dir),
            "--gene",
            gene_arg,
            "--positive-threshold",
            str(args.positive_threshold),
        ]
        print("Running:", " ".join(cmd), flush=True)
        subprocess.run(cmd, check=True)
        rendered.append(
            {
                "requested_gene": gene_arg,
                "canonical_gene": gene,
                "gene_token": token,
                "output_stem": OUTPUT_STEM_TEMPLATE.format(gene_token=token),
            }
        )

    write_batch_readme(final_dir, args.project_root, args.run_label, [row["canonical_gene"] for row in rendered], args.positive_threshold)
    shutil.copy2(single_script, final_dir / "code" / single_script.name)
    shutil.copy2(Path(__file__).resolve(), final_dir / "code" / Path(__file__).name)

    combined_study = []
    combined_manifest = []
    for row in rendered:
        token = row["gene_token"]
        study_path = final_dir / "tables" / f"cross_study_{token}_expression_study_summary.tsv"
        manifest_path = final_dir / "tables" / f"cross_study_{token}_expression_output_manifest.tsv"
        if study_path.exists():
            combined_study.append(pd.read_csv(study_path, sep="\t"))
        if manifest_path.exists():
            combined_manifest.append(pd.read_csv(manifest_path, sep="\t"))

    if combined_study:
        pd.concat(combined_study, ignore_index=True).to_csv(
            final_dir / "tables" / "cross_study_single_gene_expression_combined_study_summary.tsv",
            sep="\t",
            index=False,
        )
    if combined_manifest:
        pd.concat(combined_manifest, ignore_index=True).to_csv(
            final_dir / "tables" / "cross_study_single_gene_expression_combined_output_manifest.tsv",
            sep="\t",
            index=False,
        )

    pd.DataFrame(rendered).to_csv(
        final_dir / "provenance" / "batch_gene_manifest.tsv",
        sep="\t",
        index=False,
    )
    pd.DataFrame(
        [
            ("rendered_at", datetime.now().astimezone().isoformat()),
            ("project_root", str(args.project_root)),
            ("run_label", args.run_label),
            ("final_dir", str(final_dir)),
            ("positive_threshold_log_expression", f"{args.positive_threshold:g}"),
            ("fresh_render", str(bool(args.fresh))),
            ("genes", ",".join(args.genes)),
        ],
        columns=["key", "value"],
    ).to_csv(final_dir / "provenance" / "batch_render_manifest.tsv", sep="\t", index=False)

    stale_single_manifest = final_dir / "provenance" / "render_manifest.tsv"
    if stale_single_manifest.exists():
        stale_single_manifest.unlink()

    for handoff in (
        Path("python_notebooks/HANDOFF_transition_to_final_figs.md"),
        Path("python_notebooks/HANDOFF_cross_study_marker_synthesis_concept.md"),
    ):
        if handoff.exists():
            shutil.copy2(handoff, final_dir / "provenance" / handoff.name)

    git_status = subprocess.run(["git", "status", "--short"], check=False, capture_output=True, text=True)
    (final_dir / "provenance" / "git_status_short.txt").write_text(git_status.stdout if git_status.returncode == 0 else git_status.stderr)
    git_commit = subprocess.run(["git", "rev-parse", "HEAD"], check=False, capture_output=True, text=True)
    (final_dir / "provenance" / "git_commit.txt").write_text(git_commit.stdout if git_commit.returncode == 0 else git_commit.stderr)

    checksum_inputs = [
        final_dir / "README.md",
        *sorted((final_dir / "code").glob("*")),
        *sorted((final_dir / "figures" / "png").glob("*.png")),
        *sorted((final_dir / "figures" / "pdf").glob("*.pdf")),
        *sorted((final_dir / "figures" / "svg").glob("*.svg")),
        *sorted((final_dir / "tables").glob("*.tsv")),
        *sorted((final_dir / "provenance").glob("*.tsv")),
        *sorted((final_dir / "provenance").glob("HANDOFF_*.md")),
        final_dir / "provenance" / "git_status_short.txt",
        final_dir / "provenance" / "git_commit.txt",
    ]
    sha_lines = []
    for path in checksum_inputs:
        if path.exists() and path.is_file():
            digest = subprocess.run(["sha256sum", str(path)], check=False, capture_output=True, text=True)
            if digest.returncode == 0:
                sha_lines.append(digest.stdout)
    (final_dir / "provenance" / "sha256_manifest.txt").write_text("".join(sha_lines))

    print(f"Consolidated package: {final_dir}", flush=True)
    print(pd.DataFrame(rendered).to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
