#!/usr/bin/env python
"""Batch entry point for cached Seurat RDS to AnnData H5AD conversion."""

import argparse
import os
import platform
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
SRC_ROOT = REPO_ROOT / "python_notebooks" / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from mge_organoid_python import SeuratToAnnDataConverter, default_studies, resolve_project_root


def parse_args():
    parser = argparse.ArgumentParser(description="Convert configured Seurat objects to cached H5AD files.")
    parser.add_argument(
        "--project-root",
        default=os.environ.get("PROJECT_ROOT"),
        help="Runtime PROJECT_ROOT. Defaults to env PROJECT_ROOT or the documented Great Lakes path.",
    )
    parser.add_argument(
        "--studies",
        nargs="+",
        default=["all"],
        help="Study IDs to convert. Use all, or one or more of: shi_2019_paper_qc varela_div30 varela_div90.",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Rebuild H5AD files even if cached outputs already exist.",
    )
    parser.add_argument(
        "--allow-login-node",
        action="store_true",
        help="Allow conversion on gl-login nodes. Not recommended.",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    project_root = resolve_project_root(args.project_root)
    repo_root = REPO_ROOT

    if args.allow_login_node:
        os.environ["MGE_ALLOW_LOGIN_NODE_CONVERSION"] = "true"

    all_studies = default_studies()
    studies_by_id = {study.study_id: study for study in all_studies}
    requested = args.studies
    if requested == ["all"]:
        studies = all_studies
    else:
        missing = [study_id for study_id in requested if study_id not in studies_by_id]
        if missing:
            raise SystemExit("Unknown study ID(s): {}".format(", ".join(missing)))
        studies = [studies_by_id[study_id] for study_id in requested]

    print("=== Seurat to AnnData batch conversion ===", flush=True)
    print("hostname:", platform.node(), flush=True)
    print("PROJECT_ROOT:", project_root, flush=True)
    print("repo_root:", repo_root, flush=True)
    print("studies:", ", ".join(study.study_id for study in studies), flush=True)
    print("overwrite:", args.overwrite, flush=True)

    converter = SeuratToAnnDataConverter(
        project_root=project_root,
        repo_root=repo_root,
        overwrite=args.overwrite,
        verbose=True,
    )
    paths = converter.convert_many_files(studies, overwrite=args.overwrite)

    print("=== Conversion outputs ===", flush=True)
    for study_id, path in paths.items():
        exists = Path(path).exists()
        size = Path(path).stat().st_size if exists else 0
        print("{}\t{}\texists={}\tsize_bytes={}".format(study_id, path, exists, size), flush=True)


if __name__ == "__main__":
    main()
