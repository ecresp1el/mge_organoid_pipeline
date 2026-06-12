#!/usr/bin/env python
"""Run the Seurat object inventory helper for one configured study."""

import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
SRC_ROOT = REPO_ROOT / "python_notebooks" / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from mge_organoid_python import cached_h5ad_path, default_studies, resolve_project_root


def parse_args():
    parser = argparse.ArgumentParser(description="Run Seurat inventory for one study.")
    parser.add_argument("--study-id", required=True, help="Study ID from default_studies().")
    parser.add_argument(
        "--project-root",
        default=os.environ.get("PROJECT_ROOT"),
        help="Runtime PROJECT_ROOT. Defaults to env PROJECT_ROOT or the documented Great Lakes path.",
    )
    parser.add_argument(
        "--results-dirname",
        default=os.environ.get("SEURAT_INVENTORY_RESULTS_DIRNAME", "seurat_anndata_umap_inventory"),
        help="Results dirname under PROJECT_ROOT/results.",
    )
    parser.add_argument(
        "--run-label",
        default=os.environ.get("SEURAT_INVENTORY_RUN_LABEL", "seurat_anndata_umap_inventory_v1"),
        help="Run label under the inventory results directory.",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Rerun the R inventory even if the completion marker already exists.",
    )
    parser.add_argument(
        "--rscript",
        default=os.environ.get("RSCRIPT") or shutil.which("Rscript"),
        help="Rscript executable. Defaults to env RSCRIPT or PATH lookup.",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    if not args.rscript:
        raise SystemExit("Rscript was not found")

    project_root = resolve_project_root(args.project_root)
    studies_by_id = {study.study_id: study for study in default_studies()}
    if args.study_id not in studies_by_id:
        raise SystemExit(
            "Unknown study ID '{}'. Valid IDs: {}".format(
                args.study_id,
                ", ".join(sorted(studies_by_id)),
            )
        )
    study = studies_by_id[args.study_id]

    source = Path(study.seurat_path).expanduser()
    if not source.exists():
        raise FileNotFoundError("Missing Seurat source for {}: {}".format(study.study_id, source))
    h5ad = cached_h5ad_path(study, project_root=project_root)
    if not h5ad.exists():
        raise FileNotFoundError("Missing cached H5AD for {}: {}".format(study.study_id, h5ad))

    run_dir = project_root / "results" / args.results_dirname / args.run_label
    table_dir = run_dir / "tables" / "seurat_object_inventory" / study.study_id
    log_dir = run_dir / "logs"
    table_dir.mkdir(parents=True, exist_ok=True)
    log_dir.mkdir(parents=True, exist_ok=True)

    marker = table_dir / "seurat_inventory_complete.tsv"
    log_path = log_dir / "inspect_seurat_object.{}.log".format(study.study_id)
    if marker.exists() and not args.overwrite:
        print("Using existing Seurat inventory marker: {}".format(marker), flush=True)
        return 0

    r_helper = REPO_ROOT / "python_notebooks" / "scripts" / "inspect_seurat_object.R"
    if not r_helper.exists():
        raise FileNotFoundError("Missing R helper: {}".format(r_helper))

    cmd = [
        args.rscript,
        str(r_helper),
        "--study_id",
        study.study_id,
        "--label",
        study.label,
        "--seurat",
        str(source),
        "--h5ad",
        str(h5ad),
        "--outdir",
        str(table_dir),
    ]

    print("=== Seurat inventory one-study run ===", flush=True)
    print("study_id:", study.study_id, flush=True)
    print("source:", source, flush=True)
    print("h5ad:", h5ad, flush=True)
    print("outdir:", table_dir, flush=True)
    print("log:", log_path, flush=True)
    print("command:", " ".join(cmd), flush=True)

    with log_path.open("w") as log:
        proc = subprocess.Popen(
            cmd,
            cwd=str(REPO_ROOT),
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,
        )
        assert proc.stdout is not None
        for line in proc.stdout:
            print(line, end="", flush=True)
            log.write(line)
        returncode = proc.wait()

    if returncode != 0:
        raise RuntimeError(
            "Seurat inventory failed for {} with exit code {}. See {}".format(
                study.study_id,
                returncode,
                log_path,
            )
        )
    if not marker.exists():
        raise RuntimeError("Seurat inventory completed without marker: {}".format(marker))

    print("Completed Seurat inventory marker: {}".format(marker), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
