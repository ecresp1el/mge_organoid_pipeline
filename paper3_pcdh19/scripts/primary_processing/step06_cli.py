"""Command-line entry point for Step 06 unintegrated diagnostics."""

from __future__ import annotations

import argparse
from pathlib import Path

from .step06_models import Step06Paths, Step06Settings
from .step06_workflow import Step06Workflow


def parser() -> argparse.ArgumentParser:
    """Construct the explicit approved-input and diagnostic parameter contract."""

    command = argparse.ArgumentParser(description=__doc__)
    command.add_argument("--input-h5ad", type=Path, required=True)
    command.add_argument("--input-status", type=Path, required=True)
    command.add_argument("--input-manifest", type=Path, required=True)
    command.add_argument("--approval-ledger-snapshot", type=Path, required=True)
    command.add_argument("--bypass-decision", type=Path, required=True)
    command.add_argument("--run-dir", type=Path, required=True)
    command.add_argument("--workflow-root", type=Path, required=True)
    command.add_argument("--frozen-code-dir", type=Path, required=True)
    command.add_argument("--expected-step02-run-id", required=True)
    command.add_argument("--expected-input-bytes", type=int, required=True)
    command.add_argument("--expected-input-sha256", required=True)
    for name in ("expected-cells", "expected-genes", "n-top-genes", "pca-components", "centroid-components", "n-neighbors", "n-jobs", "random-seed", "render-max-cells", "plot-dpi"):
        command.add_argument(f"--{name}", type=int, required=True)
    for name in ("target-sum", "umap-min-dist", "umap-spread", "leiden-resolution"):
        command.add_argument(f"--{name}", type=float, required=True)
    for name in ("sample-field", "genotype-field", "sex-field", "design-field", "hvg-flavor", "neighbor-metric", "compression"):
        command.add_argument(f"--{name}", required=True)
    return command


def main() -> None:
    """Run the frozen Step 06 workflow and report its review boundary."""

    arguments = vars(parser().parse_args())
    path_names = {
        "input_h5ad", "input_status", "input_manifest", "approval_ledger_snapshot",
        "bypass_decision", "run_dir", "workflow_root", "frozen_code_dir",
        "expected_step02_run_id", "expected_input_bytes", "expected_input_sha256",
    }
    path_values = {name: arguments.pop(name) for name in path_names}
    paths = Step06Paths(**path_values)
    settings = Step06Settings(**arguments)
    output = Step06Workflow(paths, settings).run()
    print(f"Published Step 06 unintegrated diagnostic checkpoint: {output}")
    print("Status: IN_REVIEW; no integration, correction, or cell removal")


if __name__ == "__main__":
    main()
