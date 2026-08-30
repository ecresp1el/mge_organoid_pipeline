"""Command-line entry point for approved Step 02 QC filtering."""

from __future__ import annotations

import argparse
from pathlib import Path

from .step02_models import Step02Paths, Step02Settings
from .step02_workflow import Step02Workflow


def _parser() -> argparse.ArgumentParser:
    """Construct the explicit Step 02 input, policy, and evidence contract."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-h5ad", type=Path, required=True)
    parser.add_argument("--input-step01-status", type=Path, required=True)
    parser.add_argument("--input-step01-manifest", type=Path, required=True)
    parser.add_argument("--input-flags", type=Path, required=True)
    parser.add_argument("--input-step01a-status", type=Path, required=True)
    parser.add_argument("--input-step01a-manifest", type=Path, required=True)
    parser.add_argument("--approval-ledger-snapshot", type=Path, required=True)
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument("--workflow-root", type=Path, required=True)
    parser.add_argument("--frozen-code-dir", type=Path, required=True)
    parser.add_argument("--expected-step01-run-id", required=True)
    parser.add_argument("--expected-step01-bytes", type=int, required=True)
    parser.add_argument("--expected-step01-sha256", required=True)
    parser.add_argument("--expected-step01a-run-id", required=True)
    parser.add_argument("--expected-flags-bytes", type=int, required=True)
    parser.add_argument("--expected-flags-sha256", required=True)
    parser.add_argument("--expected-cells-before", type=int, default=450_788)
    parser.add_argument("--expected-cells-after", type=int, default=446_349)
    parser.add_argument("--expected-removed", type=int, default=4_439)
    parser.add_argument("--expected-genes", type=int, default=19_071)
    parser.add_argument("--selected-stringency", type=int, default=5)
    parser.add_argument("--sample-field", default="technical_sample_id")
    parser.add_argument("--design-field", default="design_group")
    parser.add_argument("--plot-dpi", type=int, default=180)
    parser.add_argument("--compression", default="lzf")
    return parser


def main() -> None:
    """Apply only approved filtering, publish review assets, and stop."""

    arguments = _parser().parse_args()
    settings = Step02Settings(
        expected_cells_before=arguments.expected_cells_before,
        expected_cells_after=arguments.expected_cells_after,
        expected_removed=arguments.expected_removed,
        expected_genes=arguments.expected_genes,
        selected_stringency=arguments.selected_stringency,
        sample_field=arguments.sample_field,
        design_field=arguments.design_field,
        plot_dpi=arguments.plot_dpi,
        compression=arguments.compression,
    )
    paths = Step02Paths(
        input_h5ad=arguments.input_h5ad,
        input_step01_status=arguments.input_step01_status,
        input_step01_manifest=arguments.input_step01_manifest,
        input_flags=arguments.input_flags,
        input_step01a_status=arguments.input_step01a_status,
        input_step01a_manifest=arguments.input_step01a_manifest,
        approval_ledger_snapshot=arguments.approval_ledger_snapshot,
        run_dir=arguments.run_dir,
        workflow_root=arguments.workflow_root,
        frozen_code_dir=arguments.frozen_code_dir,
        expected_step01_run_id=arguments.expected_step01_run_id,
        expected_step01_bytes=arguments.expected_step01_bytes,
        expected_step01_sha256=arguments.expected_step01_sha256,
        expected_step01a_run_id=arguments.expected_step01a_run_id,
        expected_flags_bytes=arguments.expected_flags_bytes,
        expected_flags_sha256=arguments.expected_flags_sha256,
    )
    output = Step02Workflow(paths, settings).run()
    print(f"Published Step 02 filtered checkpoint: {output}")
    print("Status: IN_REVIEW; scDblFinder not run")


if __name__ == "__main__":
    main()
