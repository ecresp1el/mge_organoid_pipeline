"""Command-line entry point for Step 01a technical-outlier sensitivity."""

from __future__ import annotations

import argparse
from pathlib import Path

from .step01a_models import Step01aPaths, Step01aSettings
from .step01a_workflow import Step01aWorkflow


def _parse_stringencies(value: str) -> tuple[float, ...]:
    """Parse an ordered comma-delimited list of positive MAD multipliers."""

    parsed = tuple(float(item) for item in value.split(","))
    if not parsed or any(item <= 0 for item in parsed) or tuple(sorted(parsed)) != parsed:
        raise argparse.ArgumentTypeError("stringencies must be ordered positive values")
    return parsed


def _parser() -> argparse.ArgumentParser:
    """Construct the explicit read-only Step 01a command-line contract."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-h5ad", type=Path, required=True)
    parser.add_argument("--input-step-status", type=Path, required=True)
    parser.add_argument("--input-output-manifest", type=Path, required=True)
    parser.add_argument("--approval-ledger-snapshot", type=Path, required=True)
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument("--workflow-root", type=Path, required=True)
    parser.add_argument("--frozen-code-dir", type=Path, required=True)
    parser.add_argument("--expected-input-run-id", required=True)
    parser.add_argument("--expected-input-bytes", type=int, required=True)
    parser.add_argument("--expected-input-sha256", required=True)
    parser.add_argument("--expected-cells", type=int, default=450_788)
    parser.add_argument("--expected-genes", type=int, default=19_071)
    parser.add_argument("--expected-samples", type=int, default=12)
    parser.add_argument("--sample-field", default="technical_sample_id")
    parser.add_argument("--stringencies", type=_parse_stringencies, default=(3.0, 4.0, 5.0))
    parser.add_argument("--mad-scale", type=float, default=1.4826)
    parser.add_argument("--histogram-bins", type=int, default=80)
    parser.add_argument("--plot-dpi", type=int, default=180)
    return parser


def main() -> None:
    """Run Step 01a only and print the published per-cell candidate table."""

    arguments = _parser().parse_args()
    settings = Step01aSettings(
        expected_cells=arguments.expected_cells,
        expected_genes=arguments.expected_genes,
        expected_samples=arguments.expected_samples,
        sample_field=arguments.sample_field,
        stringencies=arguments.stringencies,
        mad_scale=arguments.mad_scale,
        histogram_bins=arguments.histogram_bins,
        plot_dpi=arguments.plot_dpi,
    )
    paths = Step01aPaths(
        input_h5ad=arguments.input_h5ad,
        input_step_status=arguments.input_step_status,
        input_output_manifest=arguments.input_output_manifest,
        approval_ledger_snapshot=arguments.approval_ledger_snapshot,
        run_dir=arguments.run_dir,
        workflow_root=arguments.workflow_root,
        frozen_code_dir=arguments.frozen_code_dir,
        expected_input_run_id=arguments.expected_input_run_id,
        expected_input_bytes=arguments.expected_input_bytes,
        expected_input_sha256=arguments.expected_input_sha256,
    )
    output = Step01aWorkflow(paths, settings).run()
    print(f"Published Step 01a candidate table: {output}")
    print("Status: IN_REVIEW; no cells removed")


if __name__ == "__main__":
    main()

