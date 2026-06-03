#!/usr/bin/env python
"""Run a Seurat cluster-resolution sweep for one configured study."""

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

from mge_organoid_python import default_studies, resolve_project_root


def parse_args():
    parser = argparse.ArgumentParser(description="Run Seurat cluster-resolution sweep for one study.")
    parser.add_argument("--study-id", required=True, help="Study ID from default_studies().")
    parser.add_argument(
        "--project-root",
        default=os.environ.get("PROJECT_ROOT"),
        help="Runtime PROJECT_ROOT. Defaults to env PROJECT_ROOT or the documented Great Lakes path.",
    )
    parser.add_argument(
        "--results-dirname",
        default=os.environ.get("SEURAT_RESOLUTION_SWEEP_RESULTS_DIRNAME", "seurat_cluster_resolution_sweep"),
        help="Results dirname under PROJECT_ROOT/results.",
    )
    parser.add_argument(
        "--run-label",
        default=os.environ.get("SEURAT_RESOLUTION_SWEEP_RUN_LABEL", "seurat_cluster_resolution_sweep_v1"),
        help="Run label under the sweep results directory.",
    )
    parser.add_argument(
        "--resolutions",
        default=os.environ.get("SEURAT_RESOLUTION_SWEEP_RESOLUTIONS", "0.1,0.2,0.3,0.4,0.5,0.6,0.8,1.0,1.2"),
        help="Comma/colon/space-separated resolution grid.",
    )
    parser.add_argument("--graph", default=os.environ.get("SEURAT_RESOLUTION_SWEEP_GRAPH", "RNA_snn"))
    parser.add_argument("--reduction", default=os.environ.get("SEURAT_RESOLUTION_SWEEP_REDUCTION", "umap"))
    parser.add_argument("--algorithm", default=os.environ.get("SEURAT_RESOLUTION_SWEEP_ALGORITHM", "1"))
    parser.add_argument("--random-seed", default=os.environ.get("SEURAT_RESOLUTION_SWEEP_RANDOM_SEED", "0"))
    parser.add_argument("--overwrite", action="store_true", help="Rerun even if completion marker exists.")
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

    run_dir = project_root / "results" / args.results_dirname / args.run_label
    study_dir = run_dir / study.study_id
    log_dir = run_dir / "logs"
    study_dir.mkdir(parents=True, exist_ok=True)
    log_dir.mkdir(parents=True, exist_ok=True)

    marker = study_dir / "resolution_sweep_complete.tsv"
    log_path = log_dir / "seurat_cluster_resolution_sweep.{}.log".format(study.study_id)
    if marker.exists() and not args.overwrite:
        print("Using existing resolution sweep marker: {}".format(marker), flush=True)
        return 0

    r_helper = REPO_ROOT / "python_notebooks" / "scripts" / "seurat_cluster_resolution_sweep.R"
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
        "--outdir",
        str(study_dir),
        "--graph",
        args.graph,
        "--reduction",
        args.reduction,
        "--resolutions",
        args.resolutions,
        "--algorithm",
        str(args.algorithm),
        "--random_seed",
        str(args.random_seed),
        "--overwrite",
        "true" if args.overwrite else "false",
    ]

    print("=== Seurat cluster-resolution sweep ===", flush=True)
    print("study_id:", study.study_id, flush=True)
    print("source:", source, flush=True)
    print("outdir:", study_dir, flush=True)
    print("log:", log_path, flush=True)
    print("resolutions:", args.resolutions, flush=True)
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
            "Resolution sweep failed for {} with exit code {}. See {}".format(
                study.study_id,
                returncode,
                log_path,
            )
        )
    if not marker.exists():
        raise RuntimeError("Resolution sweep completed without marker: {}".format(marker))

    print("Completed resolution sweep marker: {}".format(marker), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
