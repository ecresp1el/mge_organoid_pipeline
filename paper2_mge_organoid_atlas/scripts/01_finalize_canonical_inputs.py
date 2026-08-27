#!/usr/bin/env python3
"""Validate six canonical input pairs, write frozen manifests, and publish atomically."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import os
from pathlib import Path
import stat

import anndata as ad
import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--staging-dir", type=Path, required=True)
    parser.add_argument("--canonical-dir", type=Path, required=True)
    return parser.parse_args()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(16 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_key_values(path: Path) -> dict[str, str]:
    table = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    key_column = "key" if "key" in table.columns else "check"
    if key_column not in table.columns or "value" not in table.columns:
        raise ValueError(f"Expected key/check and value columns in {path}")
    return dict(zip(table[key_column], table["value"], strict=True))


def compare_summaries(study_dir: Path) -> None:
    rds = pd.read_csv(study_dir / "matrix_summary_rds.tsv", sep="\t")
    h5 = pd.read_csv(study_dir / "matrix_summary_h5ad.tsv", sep="\t")
    if list(rds["matrix"]) != list(h5["matrix"]):
        raise ValueError(f"Matrix layer lists differ in {study_dir}")
    for matrix_name in rds["matrix"]:
        left = rds.loc[rds["matrix"] == matrix_name].iloc[0]
        right = h5.loc[h5["matrix"] == matrix_name].iloc[0]
        for column in ["n_cells", "n_features", "nnz"]:
            if int(left[column]) != int(right[column]):
                raise ValueError(f"{study_dir.name} {matrix_name} {column} mismatch")
        for column in ["value_sum", "value_sum_squares", "value_min_nonzero", "value_max"]:
            if not pd.isna(left[column]) or not pd.isna(right[column]):
                if not abs(float(left[column]) - float(right[column])) <= 1e-8 * max(1.0, abs(float(left[column]))):
                    raise ValueError(f"{study_dir.name} {matrix_name} {column} mismatch")


def attempt_remove_write_bits(root: Path) -> bool:
    """Best-effort POSIX hardening; some Turbo NFS exports restore rw modes."""
    for path in sorted(root.rglob("*"), reverse=True):
        mode = stat.S_IMODE(path.stat().st_mode)
        path.chmod(mode & ~0o222)
    mode = stat.S_IMODE(root.stat().st_mode)
    root.chmod(mode & ~0o222)
    return not any(
        stat.S_IMODE(path.stat().st_mode) & 0o222
        for path in [root, *root.rglob("*")]
    )


def main() -> None:
    args = parse_args()
    staging = args.staging_dir.resolve()
    canonical = args.canonical_dir.resolve()
    if staging.parent != canonical.parent:
        raise ValueError("Staging and canonical directories must share the same parent")
    if not staging.name.startswith(".canonical_build_"):
        raise ValueError(f"Unsafe staging directory name: {staging.name}")
    if canonical.name != "canonical":
        raise ValueError(f"Canonical destination must be named canonical: {canonical}")
    if canonical.exists():
        raise FileExistsError(f"Frozen canonical directory already exists: {canonical}")

    manifest = pd.read_csv(args.manifest, sep="\t", dtype=str, keep_default_na=False)
    if manifest.shape[0] != 6 or manifest["study_id"].nunique() != 6:
        raise ValueError("Canonical manifest must contain exactly six unique studies")

    dataset_rows: list[dict[str, object]] = []
    file_rows: list[dict[str, object]] = []
    for record in manifest.to_dict(orient="records"):
        study_id = record["study_id"]
        study_dir = staging / study_id
        rds = study_dir / f"{study_id}_minimal.rds"
        h5ad = study_dir / f"{study_id}_minimal.h5ad"
        required = [
            rds,
            h5ad,
            study_dir / "canonical_manifest.tsv",
            study_dir / "matrix_summary_rds.tsv",
            study_dir / "matrix_summary_h5ad.tsv",
            study_dir / "rds_h5ad_equivalence.tsv",
            study_dir / "r_sessionInfo.txt",
            study_dir / "python_environment.tsv",
            study_dir / "SUCCESS.txt",
        ]
        missing = [str(path) for path in required if not path.is_file() or path.stat().st_size == 0]
        if missing:
            raise FileNotFoundError(f"Missing/empty canonical outputs for {study_id}: {missing}")
        if (study_dir / ".bridge").exists():
            raise ValueError(f"Temporary bridge was not removed for {study_id}")
        rds_manifest = read_key_values(study_dir / "canonical_manifest.tsv")
        expected_rds_values = {
            "object_class": "Seurat",
            "assays_present": "RNA",
            "reductions_present": "false",
            "graphs_present": "false",
            "neighbors_present": "false",
            "commands_present": "false",
            "tools_present": "false",
            "images_present": "false",
            "scale_data_present": "false",
            "integrated_assays_present": "false",
        }
        bad_rds_checks = [
            key for key, expected in expected_rds_values.items() if rds_manifest.get(key) != expected
        ]
        if bad_rds_checks:
            raise ValueError(f"Minimal Seurat validation failed for {study_id}: {bad_rds_checks}")
        checks = read_key_values(study_dir / "rds_h5ad_equivalence.tsv")
        required_passes = [
            "rds_exists",
            "h5ad_reopened",
            "cell_ids_and_order",
            "feature_ids_and_order",
            "cell_metadata_values",
            "feature_metadata_values",
            "counts_exact",
            "analysis_embeddings_absent",
            "analysis_pairwise_arrays_absent",
        ]
        failures = [name for name in required_passes if checks.get(name) != "PASS"]
        if failures:
            raise ValueError(f"Equivalence validation failed for {study_id}: {failures}")
        compare_summaries(study_dir)

        adata = ad.read_h5ad(h5ad, backed="r")
        try:
            if str(adata.uns.get("study_id", "")) != study_id:
                raise ValueError(f"H5AD study_id does not match directory for {study_id}")
            if not all(str(cell).startswith(f"{study_id}::") for cell in adata.obs_names):
                raise ValueError(f"H5AD canonical cell IDs have the wrong study prefix for {study_id}")
            if adata.obsm or adata.varm or adata.obsp or adata.varp:
                raise ValueError(f"Analysis-specific arrays found in {h5ad}")
            layers = ["counts (X)"] + list(adata.layers.keys())
            dataset_rows.append(
                {
                    "study_id": study_id,
                    "study_label": record["display_name"],
                    "n_cells": adata.n_obs,
                    "n_features": adata.n_vars,
                    "expression_layers": ";".join(layers),
                    "source_rds": record["source_rds"],
                    "source_sha256": record["source_sha256"],
                    "R_version": rds_manifest["R_version"],
                    "Seurat_version": rds_manifest["Seurat_version"],
                    "SeuratObject_version": rds_manifest["SeuratObject_version"],
                    "Python_version": str(adata.uns["software_versions"]["Python"]),
                    "Scanpy_version": str(adata.uns["software_versions"]["Scanpy"]),
                    "AnnData_version": str(adata.uns["software_versions"]["AnnData"]),
                    "rds_h5ad_equivalence": "PASS",
                    "cell_or_gene_subset_performed": "false",
                    "normalization_performed": "false",
                    "gene_harmonization_performed": "false",
                    "celltype_harmonization_performed": "false",
                }
            )
        finally:
            adata.file.close()
        for format_name, path in [("rds", rds), ("h5ad", h5ad)]:
            file_rows.append(
                {
                    "study_id": study_id,
                    "format": format_name,
                    "relative_path": str(path.relative_to(staging)),
                    "size_bytes": path.stat().st_size,
                    "sha256": sha256_file(path),
                }
            )

    provenance = staging / "provenance"
    provenance.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(dataset_rows).to_csv(provenance / "canonical_dataset_manifest.tsv", sep="\t", index=False)
    pd.DataFrame(file_rows).to_csv(provenance / "canonical_file_checksums.tsv", sep="\t", index=False)
    manifest.to_csv(provenance / "source_object_manifest.tsv", sep="\t", index=False)

    completed = datetime.now(timezone.utc).astimezone().isoformat()
    readme = f"""# Frozen Paper 2 canonical inputs

This directory is the only allowed input layer for downstream Paper 2 analyses.
It contains six minimal, non-integrated organoid datasets. Each study has a
current stable Seurat RDS and a current stable AnnData H5AD representing the
same cells, genes, raw counts, optional source RNA normalized data, and selected
provenance metadata.

The H5AD `X` matrix is raw counts. When the source RNA assay had a normalized
`data` layer, it is preserved as `layers["lognorm"]`; no normalization was
calculated during this build. The Seurat RDS stores the same matrices in the
RNA `counts` and optional `data` layers.

No cells or genes were subset. No genes, cell types, or metadata values were
harmonized. No UMAP, PCA, reductions, graphs, neighbors, SCT models, integrated
assays, commands, or plotting artifacts are present.

Every RDS was reloaded before its H5AD bridge was written. Each H5AD was then
reopened and checked for exact sparse-matrix equality, identical ordered cell
and feature IDs, and equivalent selected metadata. File and source checksums
are under `provenance/`.

Frozen at: `{completed}`

This directory is logically frozen. The build launcher refuses to overwrite
an existing `inputs/canonical/`, and checksums detect content changes. POSIX
write-bit removal is attempted as additional hardening, but the Turbo NFS
export may preserve project-group write modes. Historical source objects
remain provenance records and must not be read directly by downstream Paper 2
workflows.
"""
    (staging / "README.md").write_text(readme, encoding="utf-8")
    (staging / "FROZEN.txt").write_text(
        f"status=FROZEN\nfrozen_at={completed}\ndatasets=6\nequivalence=PASS\n",
        encoding="utf-8",
    )
    (staging / "SUCCESS.txt").write_text(
        f"PASS\ncompleted={completed}\ndatasets=6\ncanonical_destination={canonical}\n",
        encoding="utf-8",
    )

    # Checksum every small provenance/config/code artifact after the large RDS
    # and H5AD checksums have already been recorded in canonical_file_checksums.
    package_rows = []
    for path in sorted(staging.rglob("*")):
        relative = path.relative_to(staging)
        if (
            path.is_file()
            and path.suffix not in {".rds", ".h5ad"}
            and relative.parts[0] != "logs"
            and relative != Path("provenance/supporting_file_checksums.tsv")
        ):
            package_rows.append(
                {
                    "relative_path": str(path.relative_to(staging)),
                    "size_bytes": path.stat().st_size,
                    "sha256": sha256_file(path),
                }
            )
    pd.DataFrame(package_rows).to_csv(provenance / "supporting_file_checksums.tsv", sep="\t", index=False)

    os.replace(staging, canonical)
    permissions_hardened = attempt_remove_write_bits(canonical)
    if not permissions_hardened:
        print(
            "WARNING: Turbo retained write bits; freeze is enforced by the "
            "launcher overwrite guard, FROZEN marker, and checksums.",
            flush=True,
        )
    print(f"Frozen canonical input directory published: {canonical}", flush=True)


if __name__ == "__main__":
    main()
