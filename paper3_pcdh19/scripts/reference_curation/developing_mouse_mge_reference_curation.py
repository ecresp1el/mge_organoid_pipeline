#!/usr/bin/env python3
"""Curate the first-checkpoint processed objects for three mouse MGE references.

This program deliberately stops before label transfer, integration, raw-read
download, reclustering, or de novo annotation. The supported commands are:

``source-audit``
    Validate the submitted registry, probe public resource URLs without
    downloading their bodies, and publish the starting source/raw manifests.
``inspect-study``
    Materialize the single P0 processed object for one candidate study and
    inspect its actual structure, annotations, and saved embeddings.
``checkpoint``
    Combine the three per-study inspection rows into the required early review
    table and mark the run package complete.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import importlib.metadata
import json
import os
import platform
import re
import shutil
import stat
import subprocess
import sys
import tempfile
import urllib.error
import urllib.parse
import urllib.request
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Sequence


STUDIES = ("LaManno2021", "Bandler2022", "Mayer2018")
REGISTRY_COLUMNS = (
    "paper", "dataset", "accession", "sample", "age", "region", "genotype",
    "technology", "resource_type", "filename", "format", "url",
    "size_if_known", "download_priority", "downloaded", "reason", "notes",
)
CHECKPOINT_COLUMNS = (
    "paper", "P0_file", "actual_object_type", "dimensions",
    "author_embedding_present", "annotation_columns_present", "MGE_selectable",
    "age_selectable", "cell_level_labels_immediately_usable", "next_minimal_action",
)
RAW_COLUMNS = (
    "paper", "accession", "sample", "raw_repository", "raw_accession",
    "raw_available", "assay", "reference_genome_if_known",
    "library_chemistry_if_known", "estimated_size_if_known",
    "can_reprocess_from_raw", "notes",
)
ANNOTATION_TERMS = (
    "cluster", "class", "subclass", "cell_type", "celltype", "type", "subtype",
    "identity", "annotation", "taxonomy", "development", "age", "embryonic",
    "tissue", "region", "forebrain", "telenceph", "ganglion", "mge", "progenitor",
    "radial", "neuroblast", "microglia", "lineage", "state",
)
MGE_TERMS = ("mge", "medial ganglionic", "ganglionic eminence", "ventral telenceph")
AGE_TERMS = ("age", "embryonic", "developmental", "stage", "timepoint")


class CurationError(RuntimeError):
    """Raised when a required curation invariant is not met."""


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def sha256_file(path: Path, chunk_size: int = 8 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while True:
            block = handle.read(chunk_size)
            if not block:
                break
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> List[Dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise CurationError(f"TSV lacks a header: {path}")
        return [dict(row) for row in reader]


def write_tsv(path: Path, columns: Sequence[str], rows: Iterable[Mapping[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(prefix=f".{path.name}.", dir=str(path.parent), text=True)
    temporary = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, delimiter="\t", fieldnames=list(columns), lineterminator="\n", extrasaction="ignore")
            writer.writeheader()
            for row in rows:
                writer.writerow({column: row.get(column, "") for column in columns})
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    finally:
        if temporary.exists():
            temporary.unlink()


def validate_registry(path: Path) -> List[Dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        observed = tuple(reader.fieldnames or ())
        if observed != REGISTRY_COLUMNS:
            raise CurationError(f"Registry schema mismatch: {observed!r}")
        rows = [dict(row) for row in reader]
    if not rows:
        raise CurationError("Source registry is empty")
    if set(row["paper"] for row in rows) != set(STUDIES):
        raise CurationError("Registry must contain exactly the three candidate papers")
    for study in STUDIES:
        p0 = [row for row in rows if row["paper"] == study and row["download_priority"] == "P0"]
        if len(p0) != 1:
            raise CurationError(f"{study} must have exactly one P0 object; observed {len(p0)}")
        if not p0[0]["filename"] or not p0[0]["url"]:
            raise CurationError(f"{study} P0 row requires filename and URL")
    allowed = {"P0", "P1", "P2", "RAW"}
    invalid = sorted(set(row["download_priority"] for row in rows) - allowed)
    if invalid:
        raise CurationError(f"Invalid download priorities: {invalid}")
    return rows


def probe_url(url: str, timeout: int = 30) -> Dict[str, str]:
    result = {"url": url, "access_status": "FAIL", "http_status": "", "content_length": "", "final_url": "", "error": ""}
    headers = {"User-Agent": "paper3-pcdh19-reference-curation/1.0"}
    request = urllib.request.Request(url, headers=headers, method="HEAD")
    try:
        with urllib.request.urlopen(request, timeout=timeout) as response:
            result.update({
                "access_status": "PASS",
                "http_status": str(getattr(response, "status", "")),
                "content_length": response.headers.get("Content-Length", ""),
                "final_url": response.geturl(),
            })
    except urllib.error.HTTPError as error:
        if error.code not in (403, 405):
            result["http_status"] = str(error.code)
            result["error"] = f"HTTPError: {error}"
            return result
        fallback = urllib.request.Request(url, headers={**headers, "Range": "bytes=0-0"}, method="GET")
        try:
            with urllib.request.urlopen(fallback, timeout=timeout) as response:
                result.update({
                    "access_status": "PASS",
                    "http_status": str(getattr(response, "status", "")),
                    "content_length": response.headers.get("Content-Range", response.headers.get("Content-Length", "")),
                    "final_url": response.geturl(),
                })
        except Exception as fallback_error:  # network diagnostics must be recorded, not hidden
            result["error"] = f"{type(fallback_error).__name__}: {fallback_error}"
    except Exception as error:
        result["error"] = f"{type(error).__name__}: {error}"
    return result


def raw_manifest_rows(registry: Sequence[Mapping[str, str]], probe_by_url: Mapping[str, Mapping[str, str]]) -> List[Dict[str, str]]:
    rows = []
    for row in registry:
        if row["download_priority"] != "RAW":
            continue
        probe = probe_by_url.get(row["url"], {})
        rows.append({
            "paper": row["paper"],
            "accession": row["accession"],
            "sample": row["sample"],
            "raw_repository": "NCBI SRA/BioProject or GEO",
            "raw_accession": row["accession"],
            "raw_available": "yes" if probe.get("access_status") == "PASS" else "public_record_registered_probe_failed",
            "assay": row["technology"],
            "reference_genome_if_known": "",
            "library_chemistry_if_known": "",
            "estimated_size_if_known": row["size_if_known"],
            "can_reprocess_from_raw": "yes_if_complete_run_accessions_validate",
            "notes": row["notes"],
        })
    return rows


def command_source_audit(args: argparse.Namespace) -> None:
    registry_path = Path(args.registry).resolve()
    run_dir = Path(args.run_dir).resolve()
    source_root = Path(args.source_root).resolve()
    rows = validate_registry(registry_path)
    source_root.mkdir(parents=True, exist_ok=True)
    probes = []
    for row in rows:
        probe = probe_url(row["url"], args.timeout)
        probe.update({"paper": row["paper"], "accession": row["accession"], "resource_type": row["resource_type"], "download_priority": row["download_priority"]})
        probes.append(probe)
    probe_by_url = {row["url"]: row for row in probes}
    manifest = []
    for row in rows:
        item = dict(row)
        cached = source_root / row["paper"] / "source" / row["filename"] if row["filename"] else None
        item["downloaded"] = "yes" if cached and cached.is_file() else "no"
        manifest.append(item)
    write_tsv(run_dir / "tables" / "manifest.tsv", REGISTRY_COLUMNS, manifest)
    probe_columns = ("paper", "accession", "resource_type", "download_priority", "url", "access_status", "http_status", "content_length", "final_url", "error")
    write_tsv(run_dir / "tables" / "source_url_access_audit.tsv", probe_columns, probes)
    write_tsv(run_dir / "tables" / "raw_access_manifest.tsv", RAW_COLUMNS, raw_manifest_rows(rows, probe_by_url))
    failures = sum(row["access_status"] != "PASS" for row in probes)
    write_tsv(run_dir / "tables" / "source_audit_summary.tsv", ("metric", "value"), (
        {"metric": "registry_rows", "value": len(rows)},
        {"metric": "P0_rows", "value": sum(row["download_priority"] == "P0" for row in rows)},
        {"metric": "RAW_rows", "value": sum(row["download_priority"] == "RAW" for row in rows)},
        {"metric": "url_probe_failures", "value": failures},
        {"metric": "audit_timestamp_utc", "value": utc_now()},
    ))


def p0_row(registry_path: Path, study: str) -> Dict[str, str]:
    rows = validate_registry(registry_path)
    matches = [row for row in rows if row["paper"] == study and row["download_priority"] == "P0"]
    return matches[0]


def materialize_p0(row: Mapping[str, str], source_root: Path) -> Path:
    directory = source_root / row["paper"] / "source"
    directory.mkdir(parents=True, exist_ok=True)
    path = directory / row["filename"]
    if path.is_file():
        return path
    if path.exists():
        raise CurationError(f"P0 target exists but is not a regular file: {path}")
    descriptor, temporary_name = tempfile.mkstemp(prefix=".download.", dir=str(directory))
    os.close(descriptor)
    temporary = Path(temporary_name)
    try:
        request = urllib.request.Request(row["url"], headers={"User-Agent": "paper3-pcdh19-reference-curation/1.0"})
        with urllib.request.urlopen(request, timeout=120) as response, temporary.open("wb") as output:
            shutil.copyfileobj(response, output, 8 * 1024 * 1024)
            output.flush()
            os.fsync(output.fileno())
        if temporary.stat().st_size == 0:
            raise CurationError(f"Downloaded P0 file is empty: {row['url']}")
        os.replace(temporary, path)
        path.chmod(stat.S_IRUSR | stat.S_IRGRP | stat.S_IROTH)
    finally:
        if temporary.exists():
            temporary.unlink()
    return path


def safe_label(value: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9._-]+", "_", value).strip("_")
    return cleaned[:140] or "unnamed"


def inspect_lamanno(path: Path, study_dir: Path) -> Dict[str, str]:
    try:
        import anndata as ad
        import pandas as pd
    except ImportError as error:
        raise CurationError(f"La Manno inspection requires anndata and pandas: {error}")
    audit_dir = study_dir / "audit"
    metadata_dir = study_dir / "metadata"
    counts_dir = audit_dir / "obs_value_counts"
    for directory in (audit_dir, metadata_dir, counts_dir):
        directory.mkdir(parents=True, exist_ok=True)
    data = ad.read_h5ad(path, backed="r")
    try:
        structure = [
            f"path={path}", f"shape={data.shape}", f"X_type={type(data.X).__module__}.{type(data.X).__name__}",
            f"raw_present={data.raw is not None}", f"layers={list(data.layers.keys())}",
            f"obsm={list(data.obsm.keys())}", f"uns={list(data.uns.keys())}",
            f"obs_columns={list(data.obs.columns)}", f"var_columns={list(data.var.columns)}",
        ]
        (audit_dir / "h5ad_structure.txt").write_text("\n".join(structure) + "\n", encoding="utf-8")
        inventory = []
        candidate_columns = []
        mge_columns = []
        age_columns = []
        for column in data.obs.columns:
            series = data.obs[column]
            lower = str(column).lower()
            n_unique = int(series.nunique(dropna=False))
            sample_values = [str(value) for value in series.dropna().astype(str).unique()[:12]]
            inventory.append({
                "column": column, "dtype": str(series.dtype), "n_unique_including_na": n_unique,
                "missing": int(series.isna().sum()), "sample_values": "|".join(sample_values),
            })
            name_candidate = any(term in lower for term in ANNOTATION_TERMS)
            values_lower = "|".join(sample_values).lower()
            value_candidate = any(term in values_lower for term in ANNOTATION_TERMS)
            if (name_candidate or value_candidate) and n_unique <= 5000:
                candidate_columns.append(str(column))
                value_counts = series.astype("string").fillna("NA").value_counts(dropna=False)
                write_tsv(counts_dir / f"{safe_label(str(column))}.tsv", ("published_label", "n_cells"), (
                    {"published_label": index, "n_cells": int(count)} for index, count in value_counts.items()
                ))
            if any(term in lower or term in values_lower for term in MGE_TERMS):
                mge_columns.append(str(column))
            if any(term in lower for term in AGE_TERMS):
                age_columns.append(str(column))
        write_tsv(audit_dir / "obs_columns.tsv", ("column", "dtype", "n_unique_including_na", "missing", "sample_values"), inventory)
        write_tsv(metadata_dir / "candidate_annotation_columns.tsv", ("annotation_column",), ({"annotation_column": column} for column in candidate_columns))
        embeddings = list(data.obsm.keys())
        author_embeddings = [key for key in embeddings if "umap" in key.lower() or "tsne" in key.lower()]
        write_tsv(audit_dir / "embedding_inventory.tsv", ("obsm_key", "coordinate_shape", "author_provided_not_recomputed"), (
            {"obsm_key": key, "coordinate_shape": str(getattr(data.obsm[key], "shape", "")), "author_provided_not_recomputed": "yes"}
            for key in embeddings
        ))
        return {
            "paper": "LaManno2021", "P0_file": path.name,
            "actual_object_type": "anndata.AnnData (backed H5AD)",
            "dimensions": f"{data.n_obs} cells x {data.n_vars} genes",
            "author_embedding_present": "yes:" + "|".join(author_embeddings) if author_embeddings else "no_saved_umap_or_tsne_key_detected",
            "annotation_columns_present": "yes:" + "|".join(candidate_columns) if candidate_columns else "no_plausible_columns_detected",
            "MGE_selectable": "yes_candidate_columns:" + "|".join(sorted(set(mge_columns))) if mge_columns else "not_yet_proven",
            "age_selectable": "yes_candidate_columns:" + "|".join(sorted(set(age_columns))) if age_columns else "not_yet_proven",
            "cell_level_labels_immediately_usable": "yes_pending_semantic_validation" if candidate_columns else "no",
            "next_minimal_action": "Review candidate obs values and author taxonomy; then verify MGE/age subset and plot only saved author embeddings.",
        }
    finally:
        try:
            data.file.close()
        except Exception:
            pass


def inspect_bandler(path: Path, study_dir: Path, r_helper: Path) -> Dict[str, str]:
    audit_dir = study_dir / "audit"
    audit_dir.mkdir(parents=True, exist_ok=True)
    summary_path = audit_dir / "checkpoint_row.tsv"
    command = ["Rscript", str(r_helper), "--input", str(path), "--study-dir", str(study_dir), "--summary", str(summary_path)]
    completed = subprocess.run(command, text=True, capture_output=True)
    (audit_dir / "r_inspection_stdout.txt").write_text(completed.stdout, encoding="utf-8")
    (audit_dir / "r_inspection_stderr.txt").write_text(completed.stderr, encoding="utf-8")
    if completed.returncode != 0:
        raise CurationError(f"Bandler RDS inspection failed with exit {completed.returncode}")
    rows = read_tsv(summary_path)
    if len(rows) != 1:
        raise CurationError("Bandler R helper did not publish exactly one checkpoint row")
    return rows[0]


def inspect_mayer(path: Path, study_dir: Path) -> Dict[str, str]:
    audit_dir = study_dir / "audit"
    metadata_dir = study_dir / "metadata"
    audit_dir.mkdir(parents=True, exist_ok=True)
    metadata_dir.mkdir(parents=True, exist_ok=True)
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", errors="strict", newline="") as handle:
        reader = csv.reader(handle)
        try:
            header = next(reader)
        except StopIteration:
            raise CurationError("Mayer processed matrix is empty")
        rows = 0
        first_fields = []
        for record in reader:
            rows += 1
            if len(record) != len(header):
                raise CurationError(f"Mayer CSV row {rows + 1} has {len(record)} columns; expected {len(header)}")
            if len(first_fields) < 12:
                first_fields.append(record[0] if record else "")
    cell_identifiers = header[1:]
    header_text = "|".join(cell_identifiers).lower()
    mge_evidence = any(term in header_text for term in ("gsm2790898", "mge_e13.5", "mge"))
    age_evidence = "e13.5" in header_text or "e13_5" in header_text
    write_tsv(metadata_dir / "matrix_column_inventory.tsv", ("column_index", "column_name"), (
        {"column_index": index, "column_name": value} for index, value in enumerate(header)
    ))
    (audit_dir / "csv_structure.txt").write_text(
        "\n".join((
            f"path={path}", "format=gzip-compressed CSV expression matrix", f"data_rows={rows}",
            f"columns={len(header)}", f"first_column={header[0] if header else ''}",
            f"first_row_identifiers={first_fields}",
        )) + "\n", encoding="utf-8",
    )
    return {
        "paper": "Mayer2018", "P0_file": path.name,
        "actual_object_type": "gzip-compressed CSV expression matrix",
        "dimensions": f"{rows} data rows x {max(0, len(header) - 1)} putative cell columns",
        "author_embedding_present": "no_embedding_in_expression_CSV",
        "annotation_columns_present": "no_cell_metadata_columns_in_expression_CSV_header" if len(header) > 1 else "no",
        "MGE_selectable": "yes_from_cell_identifier_evidence" if mge_evidence else "sample_known_MGE_but_cell_ID_linkage_not_yet_proven",
        "age_selectable": "yes_from_cell_identifier_evidence" if age_evidence else "sample_known_E13.5_but_cell_ID_linkage_not_yet_proven",
        "cell_level_labels_immediately_usable": "no",
        "next_minimal_action": "Locate the smallest author metadata/embedding artifact that maps these cell identifiers to published labels; preserve Lhx6-positive selection bias.",
    }


def command_inspect_study(args: argparse.Namespace) -> None:
    study = args.study
    if study not in STUDIES:
        raise CurationError(f"Unsupported study: {study}")
    run_dir = Path(args.run_dir).resolve()
    source_root = Path(args.source_root).resolve()
    registry_path = Path(args.registry).resolve()
    study_dir = run_dir / study
    for child in ("metadata", "figures", "audit"):
        (study_dir / child).mkdir(parents=True, exist_ok=True)
    row = p0_row(registry_path, study)
    path = materialize_p0(row, source_root)
    source_record = {
        "paper": study, "filename": path.name, "resolved_path": str(path),
        "url": row["url"], "downloaded_or_reused_at_utc": utc_now(),
        "bytes": path.stat().st_size, "sha256": sha256_file(path), "read_only": "yes" if not os.access(path, os.W_OK) else "filesystem_reports_writable",
    }
    write_tsv(study_dir / "audit" / "source_file.tsv", tuple(source_record), (source_record,))
    if study == "LaManno2021":
        summary = inspect_lamanno(path, study_dir)
    elif study == "Bandler2022":
        summary = inspect_bandler(path, study_dir, Path(args.r_helper).resolve())
    else:
        summary = inspect_mayer(path, study_dir)
    write_tsv(study_dir / "audit" / "checkpoint_row.tsv", CHECKPOINT_COLUMNS, (summary,))
    write_tsv(study_dir / "audit" / "inspection_status.tsv", ("paper", "status", "completed_utc"), (
        {"paper": study, "status": "PASS", "completed_utc": utc_now()},
    ))


def package_versions() -> List[Dict[str, str]]:
    rows = [
        {"component": "python", "version": platform.python_version()},
        {"component": "platform", "version": platform.platform()},
    ]
    for name in ("anndata", "h5py", "matplotlib", "pandas", "scipy"):
        try:
            version = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            version = "NOT_INSTALLED"
        rows.append({"component": name, "version": version})
    return rows


def command_checkpoint(args: argparse.Namespace) -> None:
    run_dir = Path(args.run_dir).resolve()
    source_root = Path(args.source_root).resolve()
    combined = []
    for study in STUDIES:
        status_path = run_dir / study / "audit" / "inspection_status.tsv"
        checkpoint_path = run_dir / study / "audit" / "checkpoint_row.tsv"
        if not status_path.is_file() or not checkpoint_path.is_file():
            raise CurationError(f"Missing completed P0 inspection for {study}")
        status = read_tsv(status_path)
        if len(status) != 1 or status[0].get("status") != "PASS":
            raise CurationError(f"P0 inspection did not pass for {study}")
        rows = read_tsv(checkpoint_path)
        if len(rows) != 1:
            raise CurationError(f"Expected one checkpoint row for {study}")
        combined.append(rows[0])
    write_tsv(run_dir / "tables" / "early_processed_object_checkpoint.tsv", CHECKPOINT_COLUMNS, combined)
    manifest_path = run_dir / "tables" / "manifest.tsv"
    manifest = read_tsv(manifest_path)
    for row in manifest:
        cached = source_root / row["paper"] / "source" / row["filename"] if row.get("filename") else None
        row["downloaded"] = "yes" if cached and cached.is_file() else "no"
    write_tsv(manifest_path, REGISTRY_COLUMNS, manifest)
    write_tsv(run_dir / "tables" / "software_versions.tsv", ("component", "version"), package_versions())
    success = (
        "PASS\n"
        f"completed={utc_now()}\n"
        "step=00_developing_mouse_mge_reference_curation\n"
        "scope=EARLY_PROCESSED_OBJECT_CHECKPOINT_ONLY\n"
        "stop_before=reconstruction,heavy_analysis,reference_selection,label_transfer\n"
    )
    (run_dir / "SUCCESS.txt").write_text(success, encoding="utf-8")
    failed = run_dir / "FAILED.txt"
    if failed.exists():
        failed.unlink()


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    source = subparsers.add_parser("source-audit", help="Validate and probe the source registry")
    source.add_argument("--registry", required=True)
    source.add_argument("--run-dir", required=True)
    source.add_argument("--source-root", required=True)
    source.add_argument("--timeout", type=int, default=30)
    source.set_defaults(func=command_source_audit)
    inspect = subparsers.add_parser("inspect-study", help="Download/reuse and inspect one P0 object")
    inspect.add_argument("--study", choices=STUDIES, required=True)
    inspect.add_argument("--registry", required=True)
    inspect.add_argument("--run-dir", required=True)
    inspect.add_argument("--source-root", required=True)
    inspect.add_argument("--r-helper", required=True)
    inspect.set_defaults(func=command_inspect_study)
    checkpoint = subparsers.add_parser("checkpoint", help="Combine the three P0 inspection rows")
    checkpoint.add_argument("--run-dir", required=True)
    checkpoint.add_argument("--source-root", required=True)
    checkpoint.set_defaults(func=command_checkpoint)
    return parser


def main(argv: Sequence[str] | None = None) -> None:
    args = build_parser().parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    try:
        main()
    except CurationError as error:
        print(f"ERROR: {error}", file=sys.stderr)
        raise SystemExit(2)
