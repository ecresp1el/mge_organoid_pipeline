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
from abc import ABC, abstractmethod
from collections import Counter, defaultdict
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
    "age_selectable", "cell_level_labels_immediately_usable",
    "published_sample_count", "published_sample_ids", "ages", "tissues_regions",
    "processed_object_sample_count", "processed_object_sample_ids",
    "published_samples_linked_to_P0_count", "unresolved_object_sample_count",
    "genotypes", "sc_data_type", "library_strategy", "library_source",
    "library_selection", "library_layout", "instrument_or_platform",
    "library_chemistry", "experimental_selection", "raw_accessions",
    "sample_inventory_status", "next_minimal_action",
)
RAW_COLUMNS = (
    "paper", "accession", "sample", "raw_repository", "raw_accession",
    "raw_available", "assay", "reference_genome_if_known",
    "library_chemistry_if_known", "estimated_size_if_known",
    "can_reprocess_from_raw", "notes",
)
SAMPLE_COLUMNS = (
    "paper", "dataset_accession", "sample_id", "sample_title", "sample_scope",
    "age", "tissue", "region", "genotype", "sex", "strain", "organism",
    "pool_id", "replicate", "quality_control_status", "project_or_cohort",
    "sc_data_type", "technology", "library_strategy", "library_source",
    "library_selection", "library_layout", "instrument_model", "platform_id",
    "library_chemistry", "sequencing_protocol", "reference_genome",
    "cell_isolation", "dissociation_method", "experimental_selection",
    "raw_accession", "biosample_accession", "processed_object_membership",
    "n_cells_in_processed_object", "metadata_source", "evidence_status", "notes",
)
SAMPLE_METADATA_LONG_COLUMNS = (
    "paper", "sample_id", "metadata_field", "metadata_value", "source_file",
)
STUDY_SUMMARY_COLUMNS = (
    "paper", "published_sample_count", "published_sample_ids", "ages",
    "processed_object_sample_count", "processed_object_sample_ids",
    "published_samples_linked_to_P0_count", "unresolved_object_sample_count",
    "tissues", "regions", "genotypes", "sexes", "strains", "organisms",
    "pool_ids", "replicates", "quality_control_statuses", "projects_or_cohorts",
    "sc_data_types", "technologies", "library_strategies", "library_sources",
    "library_selections", "library_layouts", "instrument_models", "platform_ids",
    "library_chemistries", "sequencing_protocols", "reference_genomes",
    "cell_isolation_methods", "dissociation_methods", "experimental_selections",
    "raw_accessions", "processed_object_membership", "sample_inventory_status",
    "remaining_sample_metadata_gaps",
)
CURATION_LEDGER_COLUMNS = (
    "paper", "requirement", "status", "evidence_file", "observed_value",
    "next_action",
)
SAMPLE_SCHEMA_COLUMNS = ("field", "meaning", "evidence_rule")
SAMPLE_SCHEMA = (
    ("sample_id", "Published or author-defined sample/library identifier", "Never synthesize an identifier from biology alone"),
    ("sample_scope", "Whether the row describes a published sample, author registry entry, or P0-object membership", "Keep paper metadata separate from object linkage"),
    ("age", "Developmental or postnatal age/stage", "Preserve exact submitted text; bounded ranges remain ranges"),
    ("tissue", "Submitted source tissue", "Preserve deposited or author value"),
    ("region", "Anatomical region or pool name", "Preserve deposited or author value"),
    ("genotype", "Genotype, condition, or reporter status", "Do not replace missing genotype with WT unless explicitly documented"),
    ("sex", "Submitted sex", "Blank when unavailable"),
    ("strain", "Mouse strain/background", "Blank when unavailable"),
    ("pool_id", "Pool or multiplex group", "Preserve author ID"),
    ("replicate", "Embryo, litter, donor, or replicate identifier", "Preserve its original granularity"),
    ("quality_control_status", "Author-reported sample QC state", "Retain failed and excluded samples visibly"),
    ("project_or_cohort", "Author project or cohort grouping", "Preserve author grouping"),
    ("sc_data_type", "Single-cell/single-nucleus and measured modality", "Do not infer modality from filename alone"),
    ("technology", "Capture or profiling technology", "Use deposited or author metadata"),
    ("library_strategy", "Repository library strategy", "Preserve repository vocabulary"),
    ("library_source", "Repository library source", "Preserve repository vocabulary"),
    ("library_selection", "Repository library selection", "Preserve repository vocabulary"),
    ("library_layout", "Single- or paired-end layout", "Blank when unavailable"),
    ("instrument_model", "Sequencing instrument", "Preserve repository value"),
    ("platform_id", "Repository platform accession", "Preserve repository value"),
    ("library_chemistry", "10x/Smart-seq version or kit", "Record only explicit evidence"),
    ("sequencing_protocol", "Protocol name or deposited protocol detail", "Record only explicit evidence"),
    ("reference_genome", "Genome build/reference used by authors", "Blank until verified"),
    ("cell_isolation", "Sorting, enrichment, capture, or reporter selection", "Selection bias must remain visible"),
    ("dissociation_method", "Tissue dissociation method", "Blank until verified"),
    ("experimental_selection", "Biological enrichment or exclusion relevant to reference suitability", "Never collapse into genotype"),
    ("raw_accession", "Raw run/experiment accession", "Verify without downloading reads"),
    ("biosample_accession", "BioSample accession", "Preserve repository identifier"),
    ("processed_object_membership", "Whether sample IDs are proven to occur in the P0 object", "Do not claim linkage without a stable cell/sample ID"),
    ("n_cells_in_processed_object", "Cells linked to the sample in the P0 object", "Blank if a sample-to-cell join is not proven"),
    ("metadata_source", "File/object supplying the value", "Every merged value retains source provenance"),
    ("evidence_status", "Strength or limitation of the sample mapping", "Use explicit unresolved states"),
)
ANNOTATION_TERMS = (
    "cluster", "class", "subclass", "cell_type", "celltype", "type", "subtype",
    "identity", "annotation", "taxonomy", "development", "age", "embryonic",
    "tissue", "region", "forebrain", "telenceph", "ganglion", "mge", "progenitor",
    "radial", "neuroblast", "microglia", "lineage", "state",
)
MGE_TERMS = ("mge", "medial ganglionic", "ganglionic eminence", "ventral telenceph")
AGE_TERMS = ("age", "embryonic", "developmental", "stage", "timepoint")
SAMPLE_ID_PRIORITY = (
    "sampleid", "sample_id", "sample", "sample_name", "samplename", "library_id",
    "library", "dataset", "pool", "orig.ident",
)


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


class SourceRegistry:
    """Validated registry with explicit study, P0, metadata, and raw views."""

    METADATA_RESOURCE_TYPES = {"geo_series_metadata", "author_sample_registry"}

    def __init__(self, path: Path):
        self.path = path.resolve()
        with self.path.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            observed = tuple(reader.fieldnames or ())
            if observed != REGISTRY_COLUMNS:
                raise CurationError(f"Registry schema mismatch: {observed!r}")
            self.rows = [dict(row) for row in reader]
        self._validate()

    def _validate(self) -> None:
        if not self.rows:
            raise CurationError("Source registry is empty")
        if set(row["paper"] for row in self.rows) != set(STUDIES):
            raise CurationError("Registry must contain exactly the three candidate papers")
        for study in STUDIES:
            p0 = [row for row in self.rows if row["paper"] == study and row["download_priority"] == "P0"]
            if len(p0) != 1:
                raise CurationError(f"{study} must have exactly one P0 object; observed {len(p0)}")
            if not p0[0]["filename"] or not p0[0]["url"]:
                raise CurationError(f"{study} P0 row requires filename and URL")
        allowed = {"P0", "P1", "P2", "RAW"}
        invalid = sorted(set(row["download_priority"] for row in self.rows) - allowed)
        if invalid:
            raise CurationError(f"Invalid download priorities: {invalid}")

    def for_study(self, study: str) -> List[Dict[str, str]]:
        return [dict(row) for row in self.rows if row["paper"] == study]

    def p0(self, study: str) -> Dict[str, str]:
        return next(dict(row) for row in self.rows if row["paper"] == study and row["download_priority"] == "P0")

    def metadata_sources(self, study: str) -> List[Dict[str, str]]:
        return [
            dict(row) for row in self.rows
            if row["paper"] == study and row["resource_type"] in self.METADATA_RESOURCE_TYPES
        ]

    def raw_sources(self, study: str | None = None) -> List[Dict[str, str]]:
        return [
            dict(row) for row in self.rows
            if row["download_priority"] == "RAW" and (study is None or row["paper"] == study)
        ]


def validate_registry(path: Path) -> List[Dict[str, str]]:
    """Compatibility helper used by tests and simple callers."""
    return SourceRegistry(path).rows


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


class ImmutableSourceCache:
    """Atomic, read-only cache for processed P0 objects and small metadata."""

    def __init__(self, root: Path):
        self.root = root.resolve()

    @staticmethod
    def category(row: Mapping[str, str]) -> str:
        return "metadata_sources" if row["resource_type"] in SourceRegistry.METADATA_RESOURCE_TYPES else "source"

    def path_for(self, row: Mapping[str, str]) -> Path:
        return self.root / row["paper"] / self.category(row) / row["filename"]

    def materialize(self, row: Mapping[str, str]) -> Path:
        path = self.path_for(row)
        path.parent.mkdir(parents=True, exist_ok=True)
        if path.is_file():
            return path
        if path.exists():
            raise CurationError(f"Cache target exists but is not a regular file: {path}")
        descriptor, temporary_name = tempfile.mkstemp(prefix=".download.", dir=str(path.parent))
        os.close(descriptor)
        temporary = Path(temporary_name)
        try:
            request = urllib.request.Request(row["url"], headers={"User-Agent": "paper3-pcdh19-reference-curation/1.0"})
            with urllib.request.urlopen(request, timeout=120) as response, temporary.open("wb") as output:
                shutil.copyfileobj(response, output, 8 * 1024 * 1024)
                output.flush()
                os.fsync(output.fileno())
            if temporary.stat().st_size == 0:
                raise CurationError(f"Downloaded file is empty: {row['url']}")
            os.replace(temporary, path)
            path.chmod(stat.S_IRUSR | stat.S_IRGRP | stat.S_IROTH)
        finally:
            if temporary.exists():
                temporary.unlink()
        return path


class SourceAuditWorkflow:
    """Publish the pre-download registry and public-access evidence."""

    def __init__(self, registry: SourceRegistry, run_dir: Path, cache: ImmutableSourceCache, timeout: int):
        self.registry = registry
        self.run_dir = run_dir.resolve()
        self.cache = cache
        self.timeout = timeout

    def run(self) -> None:
        self.cache.root.mkdir(parents=True, exist_ok=True)
        probes = []
        for row in self.registry.rows:
            probe = probe_url(row["url"], self.timeout)
            probe.update({"paper": row["paper"], "accession": row["accession"], "resource_type": row["resource_type"], "download_priority": row["download_priority"]})
            probes.append(probe)
        probe_by_url = {row["url"]: row for row in probes}
        manifest = []
        for row in self.registry.rows:
            item = dict(row)
            cached = self.cache.path_for(row) if row["filename"] else None
            item["downloaded"] = "yes" if cached and cached.is_file() else "no"
            manifest.append(item)
        write_tsv(self.run_dir / "tables" / "manifest.tsv", REGISTRY_COLUMNS, manifest)
        probe_columns = ("paper", "accession", "resource_type", "download_priority", "url", "access_status", "http_status", "content_length", "final_url", "error")
        write_tsv(self.run_dir / "tables" / "source_url_access_audit.tsv", probe_columns, probes)
        write_tsv(self.run_dir / "tables" / "raw_access_manifest.tsv", RAW_COLUMNS, raw_manifest_rows(self.registry.rows, probe_by_url))
        failures = sum(row["access_status"] != "PASS" for row in probes)
        write_tsv(self.run_dir / "tables" / "source_audit_summary.tsv", ("metric", "value"), (
            {"metric": "registry_rows", "value": len(self.registry.rows)},
            {"metric": "P0_rows", "value": sum(row["download_priority"] == "P0" for row in self.registry.rows)},
            {"metric": "metadata_source_rows", "value": sum(row["resource_type"] in SourceRegistry.METADATA_RESOURCE_TYPES for row in self.registry.rows)},
            {"metric": "RAW_rows", "value": len(self.registry.raw_sources())},
            {"metric": "url_probe_failures", "value": failures},
            {"metric": "audit_timestamp_utc", "value": utc_now()},
        ))


def command_source_audit(args: argparse.Namespace) -> None:
    registry = SourceRegistry(Path(args.registry))
    SourceAuditWorkflow(
        registry, Path(args.run_dir), ImmutableSourceCache(Path(args.source_root)), args.timeout,
    ).run()


def p0_row(registry_path: Path, study: str) -> Dict[str, str]:
    return SourceRegistry(registry_path).p0(study)


def materialize_p0(row: Mapping[str, str], source_root: Path) -> Path:
    return ImmutableSourceCache(source_root).materialize(row)


def unique_join(values: Iterable[object], separator: str = "|") -> str:
    cleaned = []
    for value in values:
        text = str(value).strip()
        if text and text.upper() not in {"NA", "NAN", "NONE"} and text not in cleaned:
            cleaned.append(text)
    return separator.join(cleaned)


def first_matching(mapping: Mapping[str, Sequence[str]], terms: Sequence[str]) -> str:
    values = []
    for key, entries in mapping.items():
        normalized = re.sub(r"[^a-z0-9]+", "_", key.lower()).strip("_")
        if any(term in normalized for term in terms):
            values.extend(entries)
    return unique_join(values)


class GeoSoftParser:
    """Preserve GEO sample fields and standardize technical sample metadata."""

    def __init__(self, paper: str, path: Path):
        self.paper = paper
        self.path = path
        self.long_rows: List[Dict[str, str]] = []

    @staticmethod
    def _clean(value: str) -> str:
        value = value.strip()
        if len(value) >= 2 and value[0] == value[-1] == '"':
            value = value[1:-1]
        return value

    def parse(self) -> List[Dict[str, str]]:
        opener = gzip.open if self.path.suffix == ".gz" else open
        samples = []
        current_id = ""
        fields: Dict[str, List[str]] = defaultdict(list)
        characteristics: Dict[str, List[str]] = defaultdict(list)

        def flush() -> None:
            nonlocal current_id, fields, characteristics
            if not current_id:
                return
            for key, values in fields.items():
                for value in values:
                    self.long_rows.append({
                        "paper": self.paper, "sample_id": current_id,
                        "metadata_field": key, "metadata_value": value,
                        "source_file": self.path.name,
                    })
            for key, values in characteristics.items():
                for value in values:
                    self.long_rows.append({
                        "paper": self.paper, "sample_id": current_id,
                        "metadata_field": f"characteristic:{key}", "metadata_value": value,
                        "source_file": self.path.name,
                    })
            relation_text = unique_join(fields.get("relation", []))
            raw_accessions = re.findall(r"\b(?:SRR|SRX|SRS|ERP|ERX|ERR)\d+\b", relation_text)
            biosamples = re.findall(r"\bSAM[NED][A-Z]?\d+\b", relation_text)
            title = unique_join(fields.get("title", []))
            source_name = unique_join(fields.get("source_name_ch1", []))
            selection_text = unique_join((
                first_matching(characteristics, ("selection", "enrichment", "sorting", "reporter")),
                title if any(term in title.lower() for term in ("lhx6", "gfp", "sorted", "positive", "negative")) else "",
            ))
            samples.append({
                "paper": self.paper,
                "dataset_accession": unique_join(fields.get("series_id", [])),
                "sample_id": current_id,
                "sample_title": title,
                "sample_scope": "published_GEO_sample",
                "age": first_matching(characteristics, ("age", "stage", "timepoint", "embryonic_day")),
                "tissue": unique_join((source_name, first_matching(characteristics, ("tissue", "source_name")))),
                "region": first_matching(characteristics, ("region", "anatomy", "brain_region")),
                "genotype": first_matching(characteristics, ("genotype", "treatment", "condition", "wild_type")),
                "sex": first_matching(characteristics, ("sex",)),
                "strain": first_matching(characteristics, ("strain", "background")),
                "organism": unique_join(fields.get("organism_ch1", [])),
                "pool_id": first_matching(characteristics, ("pool",)),
                "replicate": first_matching(characteristics, ("replicate", "embryo", "litter", "donor")),
                "quality_control_status": first_matching(characteristics, ("quality", "qc", "status")),
                "project_or_cohort": first_matching(characteristics, ("project", "cohort")),
                "sc_data_type": "single-cell RNA-seq",
                "technology": unique_join(fields.get("type", [])),
                "library_strategy": unique_join(fields.get("library_strategy", [])),
                "library_source": unique_join(fields.get("library_source", [])),
                "library_selection": unique_join(fields.get("library_selection", [])),
                "library_layout": unique_join(fields.get("library_layout", [])),
                "instrument_model": unique_join(fields.get("instrument_model", [])),
                "platform_id": unique_join(fields.get("platform_id", [])),
                "library_chemistry": first_matching(characteristics, ("chemistry", "10x_version", "kit")),
                "sequencing_protocol": first_matching(characteristics, ("sequencing_protocol", "protocol")),
                "reference_genome": first_matching(characteristics, ("genome", "reference")),
                "cell_isolation": first_matching(characteristics, ("isolation", "sorting", "capture")),
                "dissociation_method": first_matching(characteristics, ("dissociation",)),
                "experimental_selection": selection_text,
                "raw_accession": unique_join(raw_accessions),
                "biosample_accession": unique_join(biosamples),
                "processed_object_membership": "not_yet_joined_to_P0_object",
                "n_cells_in_processed_object": "",
                "metadata_source": self.path.name,
                "evidence_status": "published_GEO_sample_metadata",
                "notes": "All original GEO sample fields are preserved in published_sample_metadata_long.tsv.",
            })
            current_id = ""
            fields = defaultdict(list)
            characteristics = defaultdict(list)

        with opener(self.path, "rt", encoding="utf-8", errors="replace", newline="") as handle:
            for raw in handle:
                line = raw.rstrip("\r\n")
                if line.startswith("^SAMPLE = "):
                    flush()
                    current_id = self._clean(line.split("=", 1)[1])
                    continue
                if not current_id or not line.startswith("!Sample_") or " = " not in line:
                    continue
                key, value = line.split(" = ", 1)
                key = key[len("!Sample_"):]
                value = self._clean(value)
                if key == "characteristics_ch1" and ":" in value:
                    characteristic, payload = value.split(":", 1)
                    characteristics[characteristic.strip()].append(payload.strip())
                else:
                    fields[key].append(value)
        flush()
        return samples


class AuthorSampleRegistryParser:
    """Preserve author sample registries without inventing unavailable fields."""

    def __init__(self, paper: str, path: Path, definition: Mapping[str, str]):
        self.paper = paper
        self.path = path
        self.definition = definition
        self.long_rows: List[Dict[str, str]] = []

    @staticmethod
    def _value(record: Mapping[str, str], *names: str) -> str:
        normalized = {
            re.sub(r"[^a-z0-9]+", "_", key.lower()).strip("_"): value
            for key, value in record.items()
        }
        for name in names:
            key = re.sub(r"[^a-z0-9]+", "_", name.lower()).strip("_")
            if normalized.get(key):
                return normalized[key]
        return ""

    def parse(self) -> List[Dict[str, str]]:
        with self.path.open("r", encoding="utf-8", errors="replace", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if not reader.fieldnames:
                raise CurationError(f"Author sample registry lacks a header: {self.path}")
            rows = []
            for record in reader:
                sample_id = self._value(record, "SampleID", "sample_id", "sample")
                if not sample_id:
                    continue
                for field, value in record.items():
                    self.long_rows.append({
                        "paper": self.paper, "sample_id": sample_id,
                        "metadata_field": field, "metadata_value": value or "",
                        "source_file": self.path.name,
                    })
                pool = self._value(record, "Pool", "pool_id")
                age_match = re.search(
                    r"(?:^|[_-])((?:E|P)\d+(?:[.-]\d+)?(?:-\d+(?:[.-]\d+)?)?)",
                    pool, flags=re.IGNORECASE,
                )
                age = self._value(record, "Age", "stage", "timepoint") or (age_match.group(1) if age_match else "")
                tissue = re.sub(
                    r"[_-](?:E|P)\d+(?:[.-]\d+)?(?:-\d+(?:[.-]\d+)?)?.*$", "", pool,
                    flags=re.IGNORECASE,
                )
                qc = self._value(record, "QC", "quality_control", "status")
                project = self._value(record, "Project", "cohort")
                timepoint_pool = self._value(record, "TimepointPool")
                rows.append({
                    "paper": self.paper, "dataset_accession": "author_registry",
                    "sample_id": sample_id, "sample_title": "", "sample_scope": "author_sample_registry",
                    "age": age, "tissue": self._value(record, "Tissue") or tissue,
                    "region": pool, "genotype": self._value(record, "Genotype", "condition"),
                    "sex": self._value(record, "Sex"), "strain": self._value(record, "Strain"),
                    "organism": "Mus musculus", "sc_data_type": "single-cell RNA-seq",
                    "pool_id": pool, "replicate": self._value(record, "Replicate"),
                    "quality_control_status": qc, "project_or_cohort": project,
                    "technology": self._value(record, "Technology") or self.definition.get("technology", ""),
                    "library_strategy": "RNA-Seq",
                    "library_source": "TRANSCRIPTOMIC", "library_selection": "cDNA",
                    "library_layout": "", "instrument_model": "", "platform_id": "",
                    "library_chemistry": self._value(record, "Chemistry"), "sequencing_protocol": "",
                    "reference_genome": "", "cell_isolation": "", "dissociation_method": "",
                    "experimental_selection": "", "raw_accession": "", "biosample_accession": "",
                    "processed_object_membership": "not_yet_joined_to_P0_object",
                    "n_cells_in_processed_object": "", "metadata_source": self.path.name,
                    "evidence_status": "author_sample_registry",
                    "notes": unique_join((
                        "Blank fields were not present in the author registry and were not inferred.",
                        f"TimepointPool={timepoint_pool}" if timepoint_pool else "",
                    ), separator="; "),
                })
        return rows


class PublishedMetadataCollector:
    """Download small metadata resources and publish per-sample evidence."""

    def __init__(self, study: str, registry: SourceRegistry, cache: ImmutableSourceCache, study_dir: Path):
        self.study = study
        self.registry = registry
        self.cache = cache
        self.study_dir = study_dir

    def run(self) -> List[Dict[str, str]]:
        samples: List[Dict[str, str]] = []
        long_rows: List[Dict[str, str]] = []
        source_rows = []
        for definition in self.registry.metadata_sources(self.study):
            path = self.cache.materialize(definition)
            source_rows.append({
                "paper": self.study, "accession": definition["accession"], "filename": path.name,
                "resolved_path": str(path), "url": definition["url"], "bytes": path.stat().st_size,
                "sha256": sha256_file(path), "downloaded_or_reused_at_utc": utc_now(),
            })
            if "soft" in definition["format"].lower() or ".soft" in path.name.lower():
                parser = GeoSoftParser(self.study, path)
            else:
                parser = AuthorSampleRegistryParser(self.study, path, definition)
            samples.extend(parser.parse())
            long_rows.extend(parser.long_rows)
        metadata_dir = self.study_dir / "metadata"
        audit_dir = self.study_dir / "audit"
        write_tsv(metadata_dir / "published_sample_inventory.tsv", SAMPLE_COLUMNS, merge_sample_rows(samples))
        write_tsv(metadata_dir / "published_sample_metadata_long.tsv", SAMPLE_METADATA_LONG_COLUMNS, long_rows)
        write_tsv(audit_dir / "metadata_source_files.tsv", (
            "paper", "accession", "filename", "resolved_path", "url", "bytes", "sha256", "downloaded_or_reused_at_utc",
        ), source_rows)
        return merge_sample_rows(samples)


def merge_sample_rows(rows: Sequence[Mapping[str, str]]) -> List[Dict[str, str]]:
    """Deduplicate sample rows while retaining every nonempty observed value."""
    grouped: Dict[tuple[str, str], List[Mapping[str, str]]] = defaultdict(list)
    for row in rows:
        grouped[(row.get("paper", ""), row.get("sample_id", ""))].append(row)
    merged = []
    for key in sorted(grouped):
        records = grouped[key]
        item = {column: unique_join(record.get(column, "") for record in records) for column in SAMPLE_COLUMNS}
        item["paper"], item["sample_id"] = key
        merged.append(item)
    return merged


def safe_label(value: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9._-]+", "_", value).strip("_")
    return cleaned[:140] or "unnamed"


class AnnDataSampleInventory:
    """Derive a sample-level ledger from author-provided AnnData observations."""

    TECHNICAL_PATTERNS = (
        "sample", "library", "age", "stage", "timepoint", "tissue", "region",
        "genotype", "condition", "treatment", "sex", "strain", "batch", "replicate",
        "donor", "embryo", "litter", "chemistry", "platform", "technology", "accession",
    )

    def __init__(self, paper: str, obs, metadata_dir: Path):
        self.paper = paper
        self.obs = obs
        self.metadata_dir = metadata_dir

    def _sample_column(self) -> str | None:
        normalized = {re.sub(r"[^a-z0-9.]+", "_", str(column).lower()).strip("_"): str(column) for column in self.obs.columns}
        for candidate in SAMPLE_ID_PRIORITY:
            key = re.sub(r"[^a-z0-9.]+", "_", candidate.lower()).strip("_")
            if key in normalized:
                return normalized[key]
        for key, original in normalized.items():
            if "sample" in key and "sample" != "sample_count":
                return original
        return None

    def _matching_column(self, terms: Sequence[str]) -> str | None:
        for column in self.obs.columns:
            normalized = str(column).lower().replace(" ", "_")
            if any(term in normalized for term in terms):
                return str(column)
        return None

    def run(self) -> List[Dict[str, str]]:
        sample_column = self._sample_column()
        technical_columns = [
            str(column) for column in self.obs.columns
            if any(term in str(column).lower() for term in self.TECHNICAL_PATTERNS)
            and int(self.obs[column].nunique(dropna=False)) <= 5000
        ]
        if sample_column is None:
            groups = [("UNRESOLVED_SAMPLE_ID", self.obs)]
        else:
            groups = list(self.obs.groupby(sample_column, observed=True, dropna=False))
        field_map = {
            "age": self._matching_column(("age", "stage", "timepoint", "embryonic")),
            "tissue": self._matching_column(("tissue",)),
            "region": self._matching_column(("region", "anatom", "area", "tissue")),
            "genotype": self._matching_column(("genotype", "condition", "treatment")),
            "sex": self._matching_column(("sex",)),
            "strain": self._matching_column(("strain", "background")),
            "pool_id": self._matching_column(("pool",)),
            "replicate": self._matching_column(("replicate", "embryo", "litter", "donor")),
            "quality_control_status": self._matching_column(("quality", "qc", "status")),
            "project_or_cohort": self._matching_column(("project", "cohort")),
            "technology": self._matching_column(("technology", "assay")),
            "library_chemistry": self._matching_column(("chemistry", "10x_version", "kit")),
            "sequencing_protocol": self._matching_column(("sequencing_protocol", "protocol")),
            "platform_id": self._matching_column(("platform",)),
            "cell_isolation": self._matching_column(("isolation", "sorting", "capture")),
            "dissociation_method": self._matching_column(("dissociation",)),
            "raw_accession": self._matching_column(("sra", "run_accession", "raw_accession")),
            "biosample_accession": self._matching_column(("biosample",)),
        }
        rows = []
        long_rows = []
        for sample_value, frame in groups:
            sample_id = str(sample_value)
            standardized = {key: unique_join(frame[column].dropna().astype(str).unique()) if column else "" for key, column in field_map.items()}
            rows.append({
                "paper": self.paper, "dataset_accession": "processed_author_H5AD",
                "sample_id": sample_id, "sample_title": "", "sample_scope": "cells_in_P0_processed_object",
                "age": standardized["age"], "tissue": standardized["tissue"], "region": standardized["region"],
                "genotype": standardized["genotype"], "sex": standardized["sex"], "strain": standardized["strain"],
                "organism": "Mus musculus", "sc_data_type": "single-cell RNA-seq",
                "pool_id": standardized["pool_id"], "replicate": standardized["replicate"],
                "quality_control_status": standardized["quality_control_status"],
                "project_or_cohort": standardized["project_or_cohort"],
                "technology": standardized["technology"], "library_strategy": "RNA-Seq",
                "library_source": "TRANSCRIPTOMIC", "library_selection": "cDNA",
                "library_layout": "", "instrument_model": "", "platform_id": standardized["platform_id"],
                "library_chemistry": standardized["library_chemistry"],
                "sequencing_protocol": standardized["sequencing_protocol"], "reference_genome": "",
                "cell_isolation": standardized["cell_isolation"],
                "dissociation_method": standardized["dissociation_method"],
                "experimental_selection": "", "raw_accession": standardized["raw_accession"],
                "biosample_accession": standardized["biosample_accession"],
                "processed_object_membership": "yes",
                "n_cells_in_processed_object": len(frame), "metadata_source": "AnnData.obs",
                "evidence_status": "author_cell_metadata" if sample_column else "sample_id_not_identified_in_author_cell_metadata",
                "notes": f"sample_id_column={sample_column or 'NOT_IDENTIFIED'}",
            })
            for column in technical_columns:
                values = frame[column].astype("string").fillna("NA").value_counts(dropna=False)
                for value, count in values.items():
                    long_rows.append({
                        "paper": self.paper, "sample_id": sample_id,
                        "metadata_field": column, "metadata_value": f"{value} [n_cells={int(count)}]",
                        "source_file": "AnnData.obs",
                    })
        write_tsv(self.metadata_dir / "object_sample_inventory.tsv", SAMPLE_COLUMNS, rows)
        write_tsv(self.metadata_dir / "object_sample_metadata_long.tsv", SAMPLE_METADATA_LONG_COLUMNS, long_rows)
        write_tsv(self.metadata_dir / "object_sample_id_field.tsv", ("sample_id_column", "status"), ({
            "sample_id_column": sample_column or "", "status": "FOUND" if sample_column else "NOT_FOUND",
        },))
        return rows


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
        AnnDataSampleInventory("LaManno2021", data.obs, metadata_dir).run()
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
    sample_counts: Counter[str] = Counter()
    for identifier in cell_identifiers:
        match = re.search(r"\bGSM\d+\b", identifier)
        if match:
            sample_counts[match.group(0)] += 1
            continue
        match = re.search(r"(?:CGE|MGE|Cortex|Subcortex)[_.-][A-Za-z0-9_.-]+?(?=[_.-][ACGTN]{8,}|$)", identifier, flags=re.IGNORECASE)
        if match:
            sample_counts[match.group(0)] += 1
    object_samples = []
    if sample_counts:
        for sample_id, count in sorted(sample_counts.items()):
            object_samples.append({
                "paper": "Mayer2018", "dataset_accession": "GSE104156",
                "sample_id": sample_id, "sample_title": sample_id,
                "sample_scope": "cells_in_P0_expression_matrix", "age": "", "tissue": "", "region": "",
                "genotype": "", "sex": "", "strain": "", "organism": "Mus musculus",
                "sc_data_type": "single-cell RNA-seq", "technology": "10x Chromium",
                "library_strategy": "RNA-Seq", "library_source": "TRANSCRIPTOMIC",
                "library_selection": "cDNA", "library_layout": "", "instrument_model": "",
                "platform_id": "", "library_chemistry": "10x Chromium Single Cell 3-prime v2",
                "reference_genome": "", "experimental_selection": "",
                "raw_accession": "", "biosample_accession": "",
                "processed_object_membership": "yes", "n_cells_in_processed_object": count,
                "metadata_source": "processed_matrix_column_identifiers",
                "evidence_status": "sample_token_recovered_from_cell_identifier",
                "notes": "Must be joined to published GEO sample metadata before biological use.",
            })
    else:
        object_samples.append({
            "paper": "Mayer2018", "dataset_accession": "GSE104156",
            "sample_id": "UNRESOLVED_FROM_MATRIX_COLUMNS", "sample_title": "",
            "sample_scope": "P0_expression_matrix", "age": "", "tissue": "", "region": "",
            "genotype": "", "sex": "", "strain": "", "organism": "Mus musculus",
            "sc_data_type": "single-cell RNA-seq", "technology": "10x Chromium",
            "library_strategy": "RNA-Seq", "library_source": "TRANSCRIPTOMIC",
            "library_selection": "cDNA", "library_layout": "", "instrument_model": "",
            "platform_id": "", "library_chemistry": "10x Chromium Single Cell 3-prime v2",
            "reference_genome": "", "experimental_selection": "Lhx6-GFP selected study design",
            "raw_accession": "", "biosample_accession": "",
            "processed_object_membership": "yes_unresolved_by_sample", "n_cells_in_processed_object": len(cell_identifiers),
            "metadata_source": "processed_matrix_column_identifiers",
            "evidence_status": "sample_token_not_recovered",
            "notes": "Published sample metadata remains separate; no sample assignment was inferred.",
        })
    write_tsv(metadata_dir / "object_sample_inventory.tsv", SAMPLE_COLUMNS, object_samples)
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


class StudyObjectInspector(ABC):
    """Study-specific processed-object inspection interface."""

    def __init__(self, path: Path, study_dir: Path, r_helper: Path):
        self.path = path
        self.study_dir = study_dir
        self.r_helper = r_helper

    @abstractmethod
    def inspect(self) -> Dict[str, str]:
        raise NotImplementedError


class LaMannoObjectInspector(StudyObjectInspector):
    def inspect(self) -> Dict[str, str]:
        return inspect_lamanno(self.path, self.study_dir)


class BandlerObjectInspector(StudyObjectInspector):
    def inspect(self) -> Dict[str, str]:
        return inspect_bandler(self.path, self.study_dir, self.r_helper)


class MayerObjectInspector(StudyObjectInspector):
    def inspect(self) -> Dict[str, str]:
        return inspect_mayer(self.path, self.study_dir)


class StudyObjectInspectorFactory:
    INSPECTORS = {
        "LaManno2021": LaMannoObjectInspector,
        "Bandler2022": BandlerObjectInspector,
        "Mayer2018": MayerObjectInspector,
    }

    @classmethod
    def create(cls, study: str, path: Path, study_dir: Path, r_helper: Path) -> StudyObjectInspector:
        try:
            inspector_class = cls.INSPECTORS[study]
        except KeyError as error:
            raise CurationError(f"Unsupported study: {study}") from error
        return inspector_class(path, study_dir, r_helper)


def registry_declared_sample(study: str, p0: Mapping[str, str], summary: Mapping[str, str]) -> Dict[str, str]:
    sample_id_match = re.search(r"GSM\d+", p0.get("sample", ""))
    sample_id = sample_id_match.group(0) if sample_id_match else p0.get("sample", "") or "P0_SCOPE"
    return {
        "paper": study, "dataset_accession": p0.get("accession", ""), "sample_id": sample_id,
        "sample_title": p0.get("sample", ""), "sample_scope": "registry_declared_P0_scope",
        "age": p0.get("age", ""), "tissue": p0.get("region", ""), "region": p0.get("region", ""),
        "genotype": p0.get("genotype", ""), "sex": "", "strain": "", "organism": "Mus musculus",
        "sc_data_type": "single-cell RNA-seq", "technology": p0.get("technology", ""),
        "library_strategy": "RNA-Seq", "library_source": "TRANSCRIPTOMIC", "library_selection": "cDNA",
        "library_layout": "", "instrument_model": "", "platform_id": "", "library_chemistry": "",
        "reference_genome": "", "experimental_selection": "",
        "raw_accession": "", "biosample_accession": "", "processed_object_membership": "P0_scope",
        "n_cells_in_processed_object": "", "metadata_source": "submitted_source_registry",
        "evidence_status": "registry_scope_not_cell_barcode_join",
        "notes": unique_join((p0.get("notes", ""), "Object dimensions were not converted into a cell count without proving axis orientation."), separator="; "),
    }


def reconcile_object_samples(published: Sequence[Mapping[str, str]], object_rows: Sequence[Mapping[str, str]]) -> List[Dict[str, str]]:
    by_id = {row.get("sample_id", "").lower(): row.get("sample_id", "") for row in published if row.get("sample_id")}
    by_title = {row.get("sample_title", "").lower(): row.get("sample_id", "") for row in published if row.get("sample_title")}
    reconciled = []
    for source in object_rows:
        row = dict(source)
        token = row.get("sample_id", "").lower()
        if token in by_id:
            row["sample_id"] = by_id[token]
        elif token in by_title:
            row["sample_id"] = by_title[token]
        else:
            for title, sample_id in by_title.items():
                if title and (title in token or token in title):
                    row["sample_id"] = sample_id
                    break
        reconciled.append(row)
    return reconciled


def summarize_sample_inventory(study: str, samples: Sequence[Mapping[str, str]]) -> Dict[str, str]:
    published = [
        row for row in samples
        if any(scope.startswith(("published_", "author_")) for scope in row.get("sample_scope", "").split("|"))
        and not row.get("sample_id", "").startswith("UNRESOLVED")
    ]
    counted = published or [row for row in samples if not row.get("sample_id", "").startswith("UNRESOLVED")]
    published_ids = {row.get("sample_id", "") for row in counted if row.get("sample_id")}
    object_rows = [
        row for row in samples
        if any(
            token == "yes" or token.startswith("yes_") or token == "P0_scope"
            for token in row.get("processed_object_membership", "").split("|")
        )
    ]
    object_ids = {
        row.get("sample_id", "") for row in object_rows
        if row.get("sample_id") and not row.get("sample_id", "").startswith("UNRESOLVED")
    }
    unresolved_count = sum(
        row.get("sample_id", "").startswith("UNRESOLVED") for row in samples
    )
    field_map = {
        "published_sample_ids": "sample_id", "ages": "age", "tissues": "tissue", "regions": "region",
        "genotypes": "genotype", "sexes": "sex", "strains": "strain", "organisms": "organism",
        "pool_ids": "pool_id", "replicates": "replicate",
        "quality_control_statuses": "quality_control_status", "projects_or_cohorts": "project_or_cohort",
        "sc_data_types": "sc_data_type", "technologies": "technology", "library_strategies": "library_strategy",
        "library_sources": "library_source", "library_selections": "library_selection",
        "library_layouts": "library_layout", "instrument_models": "instrument_model",
        "platform_ids": "platform_id", "library_chemistries": "library_chemistry",
        "sequencing_protocols": "sequencing_protocol", "reference_genomes": "reference_genome",
        "cell_isolation_methods": "cell_isolation", "dissociation_methods": "dissociation_method",
        "experimental_selections": "experimental_selection",
        "raw_accessions": "raw_accession", "processed_object_membership": "processed_object_membership",
    }
    summary = {
        "paper": study,
        "published_sample_count": str(len(published_ids)),
        "processed_object_sample_count": str(len(object_ids)) if object_ids else "0",
        "processed_object_sample_ids": unique_join(sorted(object_ids)),
        "published_samples_linked_to_P0_count": str(len(published_ids & object_ids)),
        "unresolved_object_sample_count": str(unresolved_count),
    }
    for output, source in field_map.items():
        summary[output] = unique_join(row.get(source, "") for row in samples)
    summary["published_sample_ids"] = unique_join(sorted(published_ids))
    required = ("sample_id", "age", "tissue", "technology", "library_strategy", "library_source", "library_selection")
    gaps = [field for field in required if not unique_join(row.get(field, "") for row in samples)]
    summary["sample_inventory_status"] = "COMPLETE_FOR_CHECKPOINT" if not gaps else "PARTIAL_WITH_EXPLICIT_GAPS"
    summary["remaining_sample_metadata_gaps"] = unique_join(gaps)
    return summary


def enrich_checkpoint(summary: Mapping[str, str], sample_summary: Mapping[str, str]) -> Dict[str, str]:
    item = dict(summary)
    mapping = {
        "published_sample_count": "published_sample_count", "published_sample_ids": "published_sample_ids",
        "processed_object_sample_count": "processed_object_sample_count",
        "processed_object_sample_ids": "processed_object_sample_ids",
        "published_samples_linked_to_P0_count": "published_samples_linked_to_P0_count",
        "unresolved_object_sample_count": "unresolved_object_sample_count",
        "ages": "ages", "tissues_regions": "tissues", "genotypes": "genotypes",
        "sc_data_type": "sc_data_types", "library_strategy": "library_strategies",
        "library_source": "library_sources", "library_selection": "library_selections",
        "library_layout": "library_layouts", "instrument_or_platform": "instrument_models",
        "library_chemistry": "library_chemistries", "experimental_selection": "experimental_selections",
        "raw_accessions": "raw_accessions", "sample_inventory_status": "sample_inventory_status",
    }
    for destination, source in mapping.items():
        item[destination] = sample_summary.get(source, "")
    item["tissues_regions"] = unique_join((sample_summary.get("tissues", ""), sample_summary.get("regions", "")))
    item["instrument_or_platform"] = unique_join((sample_summary.get("instrument_models", ""), sample_summary.get("platform_ids", "")))
    return item


class StudyInspectionWorkflow:
    """Coordinate metadata collection, P0 inspection, and sample-ledger publication."""

    def __init__(self, study: str, registry: SourceRegistry, run_dir: Path, cache: ImmutableSourceCache, r_helper: Path):
        self.study = study
        self.registry = registry
        self.run_dir = run_dir.resolve()
        self.cache = cache
        self.r_helper = r_helper.resolve()
        self.study_dir = self.run_dir / study

    def run(self) -> None:
        for child in ("metadata", "figures", "audit"):
            (self.study_dir / child).mkdir(parents=True, exist_ok=True)
        p0 = self.registry.p0(self.study)
        path = self.cache.materialize(p0)
        source_record = {
            "paper": self.study, "filename": path.name, "resolved_path": str(path),
            "url": p0["url"], "downloaded_or_reused_at_utc": utc_now(),
            "bytes": path.stat().st_size, "sha256": sha256_file(path),
            "read_only": "yes" if not os.access(path, os.W_OK) else "filesystem_reports_writable",
        }
        write_tsv(self.study_dir / "audit" / "source_file.tsv", tuple(source_record), (source_record,))
        published_samples = PublishedMetadataCollector(self.study, self.registry, self.cache, self.study_dir).run()
        inspector = StudyObjectInspectorFactory.create(self.study, path, self.study_dir, self.r_helper)
        object_summary = inspector.inspect()
        object_inventory_path = self.study_dir / "metadata" / "object_sample_inventory.tsv"
        object_samples = read_tsv(object_inventory_path) if object_inventory_path.is_file() else []
        if not object_samples:
            object_samples = [registry_declared_sample(self.study, p0, object_summary)]
            write_tsv(object_inventory_path, SAMPLE_COLUMNS, object_samples)
        object_samples = reconcile_object_samples(published_samples, object_samples)
        combined_samples = merge_sample_rows([*published_samples, *object_samples])
        write_tsv(self.study_dir / "metadata" / "sample_inventory.tsv", SAMPLE_COLUMNS, combined_samples)
        sample_summary = summarize_sample_inventory(self.study, combined_samples)
        sample_summary["raw_accessions"] = unique_join((
            sample_summary.get("raw_accessions", ""),
            unique_join(row.get("accession", "") for row in self.registry.raw_sources(self.study)),
        ))
        write_tsv(self.study_dir / "metadata" / "study_sample_summary.tsv", STUDY_SUMMARY_COLUMNS, (sample_summary,))
        checkpoint = enrich_checkpoint(object_summary, sample_summary)
        write_tsv(self.study_dir / "audit" / "checkpoint_row.tsv", CHECKPOINT_COLUMNS, (checkpoint,))
        write_tsv(self.study_dir / "audit" / "inspection_status.tsv", ("paper", "status", "completed_utc"), (
            {"paper": self.study, "status": "PASS", "completed_utc": utc_now()},
        ))


def command_inspect_study(args: argparse.Namespace) -> None:
    registry = SourceRegistry(Path(args.registry))
    StudyInspectionWorkflow(
        args.study, registry, Path(args.run_dir), ImmutableSourceCache(Path(args.source_root)), Path(args.r_helper),
    ).run()


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


class CurationRequirementEvaluator:
    """Turn observed sample/object evidence into an explicit readiness ledger."""

    def __init__(self, run_dir: Path, study: str, checkpoint: Mapping[str, str], summary: Mapping[str, str]):
        self.run_dir = run_dir
        self.study = study
        self.checkpoint = checkpoint
        self.summary = summary
        self.rows: List[Dict[str, str]] = []

    def _add(self, requirement: str, status: str, evidence_file: str, observed_value: str, next_action: str) -> None:
        self.rows.append({
            "paper": self.study, "requirement": requirement, "status": status,
            "evidence_file": evidence_file, "observed_value": observed_value,
            "next_action": next_action,
        })

    @staticmethod
    def _presence_status(value: str) -> str:
        lower = value.lower()
        if not value:
            return "MISSING"
        if any(term in lower for term in ("unresolved", "not_yet", "not yet", "pending", "partial", "unknown")):
            return "PARTIAL"
        return "PASS"

    def evaluate(self) -> List[Dict[str, str]]:
        sample_evidence = f"{self.study}/metadata/sample_inventory.tsv"
        study_evidence = f"{self.study}/metadata/study_sample_summary.tsv"
        checkpoint_evidence = f"{self.study}/audit/checkpoint_row.tsv"
        sample_count = self.summary.get("published_sample_count", "")
        sample_ids = self.summary.get("published_sample_ids", "")
        count_status = "PASS" if sample_count and sample_ids else "MISSING"
        self._add(
            "published_sample_count_and_ids", count_status, sample_evidence,
            f"n={sample_count or 'UNKNOWN'}; ids={sample_ids or 'NOT_RECOVERED'}",
            "Recover an author/GEO sample registry" if count_status != "PASS" else "Retain exact IDs for object/barcode reconciliation",
        )
        object_sample_count = self.summary.get("processed_object_sample_count", "")
        object_sample_ids = self.summary.get("processed_object_sample_ids", "")
        unresolved = self.summary.get("unresolved_object_sample_count", "")
        object_sample_status = "PASS" if object_sample_count not in {"", "0"} and object_sample_ids else "MISSING"
        if unresolved not in {"", "0"}:
            object_sample_status = "PARTIAL"
        self._add(
            "processed_object_sample_count_and_ids", object_sample_status, sample_evidence,
            f"n={object_sample_count or 'UNKNOWN'}; ids={object_sample_ids or 'NOT_RECOVERED'}; unresolved_rows={unresolved or '0'}",
            "Identify and validate the sample field or barcode prefix in the P0 object" if object_sample_status != "PASS" else "Retain the proven object sample field",
        )
        sample_requirements = (
            ("developmental_age", "ages"), ("tissue", "tissues"), ("anatomical_region", "regions"),
            ("genotype_or_condition", "genotypes"), ("sex", "sexes"), ("strain", "strains"),
            ("pool_and_replicate_structure", "pool_ids"), ("sample_quality_control", "quality_control_statuses"),
            ("single_cell_data_type", "sc_data_types"), ("capture_technology", "technologies"),
            ("library_strategy", "library_strategies"), ("library_source", "library_sources"),
            ("library_selection", "library_selections"), ("library_layout", "library_layouts"),
            ("instrument_and_platform", "instrument_models"), ("library_chemistry", "library_chemistries"),
            ("sequencing_protocol", "sequencing_protocols"), ("reference_genome", "reference_genomes"),
            ("cell_isolation", "cell_isolation_methods"), ("dissociation_method", "dissociation_methods"),
            ("experimental_selection_bias", "experimental_selections"), ("raw_accessions", "raw_accessions"),
        )
        for requirement, field in sample_requirements:
            value = self.summary.get(field, "")
            if requirement == "instrument_and_platform":
                value = unique_join((value, self.summary.get("platform_ids", "")))
            status = self._presence_status(value)
            self._add(
                requirement, status, study_evidence, value or "NOT_RECORDED",
                "Verify from deposited/author methods metadata; do not infer" if status != "PASS" else "No early-checkpoint action",
            )
        membership = self.summary.get("processed_object_membership", "")
        membership_status = "PASS" if re.search(r"(?:^|\|)yes(?:\||$)", membership) else self._presence_status(membership)
        if "not_yet_joined" in membership or "unresolved" in membership.lower():
            membership_status = "PARTIAL"
        self._add(
            "sample_to_processed_object_linkage", membership_status, sample_evidence,
            membership or "NOT_RECORDED",
            "Prove the sample/cell-barcode join before using sample-level biology" if membership_status != "PASS" else "Retain the stable join key",
        )
        object_requirements = (
            ("processed_object_structure_inspected", "actual_object_type"),
            ("author_embedding_present", "author_embedding_present"),
            ("cell_level_annotation_columns", "annotation_columns_present"),
            ("MGE_directly_selectable", "MGE_selectable"),
            ("age_directly_selectable", "age_selectable"),
            ("cell_level_labels_immediately_usable", "cell_level_labels_immediately_usable"),
        )
        for requirement, field in object_requirements:
            value = self.checkpoint.get(field, "")
            lower = value.lower()
            if not value or lower == "no" or lower.startswith("no_"):
                status = "MISSING"
            elif any(term in lower for term in ("pending", "not_yet", "not yet", "not proven", "sample_known", "sample_level")):
                status = "PARTIAL"
            else:
                status = "PASS"
            self._add(requirement, status, checkpoint_evidence, value or "NOT_RECORDED", self.checkpoint.get("next_minimal_action", ""))
        for requirement, next_action in (
            ("all_published_annotation_levels_and_labels", "Create and verify annotation_dictionary.tsv from author labels"),
            ("required_broad_MGE_dissection_classes", "Audit RG, cycling progenitor, IPC, neuroblast, neuronal, glial, immune, and vascular author labels"),
            ("MGE_specific_developmental_states", "Quantify exact author states in MGE-relevant cells"),
            ("number_of_MGE_and_age_matched_cells", "Prove MGE and age filters, then count their intersection"),
            ("author_embedding_reproduction", "Render only saved author coordinates and record their exact key/source"),
            ("reference_readiness_classification", "Complete annotation/barcode/embedding audit before assigning readiness"),
        ):
            self._add(requirement, "NOT_ASSESSED", "", "EARLY_CHECKPOINT_STOP", next_action)
        return self.rows


class CheckpointPublisher:
    """Validate all study inspections and atomically publish combined checkpoint tables."""

    def __init__(self, run_dir: Path, registry: SourceRegistry, cache: ImmutableSourceCache):
        self.run_dir = run_dir.resolve()
        self.registry = registry
        self.cache = cache

    def _study_results(self, study: str) -> tuple[Dict[str, str], Dict[str, str], List[Dict[str, str]]]:
        status_path = self.run_dir / study / "audit" / "inspection_status.tsv"
        checkpoint_path = self.run_dir / study / "audit" / "checkpoint_row.tsv"
        summary_path = self.run_dir / study / "metadata" / "study_sample_summary.tsv"
        sample_path = self.run_dir / study / "metadata" / "sample_inventory.tsv"
        for path in (status_path, checkpoint_path, summary_path, sample_path):
            if not path.is_file():
                raise CurationError(f"Missing completed inspection output for {study}: {path}")
        status = read_tsv(status_path)
        if len(status) != 1 or status[0].get("status") != "PASS":
            raise CurationError(f"P0 inspection did not pass for {study}")
        checkpoint_rows = read_tsv(checkpoint_path)
        summary_rows = read_tsv(summary_path)
        if len(checkpoint_rows) != 1 or len(summary_rows) != 1:
            raise CurationError(f"Expected one checkpoint and one sample-summary row for {study}")
        return checkpoint_rows[0], summary_rows[0], read_tsv(sample_path)

    def run(self) -> None:
        checkpoints = []
        summaries = []
        samples = []
        ledger = []
        for study in STUDIES:
            checkpoint, summary, study_samples = self._study_results(study)
            checkpoints.append(checkpoint)
            summaries.append(summary)
            samples.extend(study_samples)
            ledger.extend(CurationRequirementEvaluator(self.run_dir, study, checkpoint, summary).evaluate())
        tables = self.run_dir / "tables"
        write_tsv(tables / "early_processed_object_checkpoint.tsv", CHECKPOINT_COLUMNS, checkpoints)
        write_tsv(tables / "all_candidate_reference_samples.tsv", SAMPLE_COLUMNS, samples)
        write_tsv(tables / "study_sample_summary.tsv", STUDY_SUMMARY_COLUMNS, summaries)
        write_tsv(tables / "reference_curation_requirements_ledger.tsv", CURATION_LEDGER_COLUMNS, ledger)
        write_tsv(tables / "sample_metadata_data_dictionary.tsv", SAMPLE_SCHEMA_COLUMNS, (
            {"field": field, "meaning": meaning, "evidence_rule": rule} for field, meaning, rule in SAMPLE_SCHEMA
        ))
        manifest_path = tables / "manifest.tsv"
        manifest = read_tsv(manifest_path)
        for row in manifest:
            cached = self.cache.path_for(row) if row.get("filename") else None
            row["downloaded"] = "yes" if cached and cached.is_file() else "no"
        write_tsv(manifest_path, REGISTRY_COLUMNS, manifest)
        write_tsv(tables / "software_versions.tsv", ("component", "version"), package_versions())
        success = (
            "PASS\n"
            f"completed={utc_now()}\n"
            "step=00_developing_mouse_mge_reference_curation\n"
            "scope=EARLY_PROCESSED_OBJECT_AND_SAMPLE_METADATA_CHECKPOINT_ONLY\n"
            "stop_before=reconstruction,heavy_analysis,reference_selection,label_transfer\n"
        )
        (self.run_dir / "SUCCESS.txt").write_text(success, encoding="utf-8")
        failed = self.run_dir / "FAILED.txt"
        if failed.exists():
            failed.unlink()


def command_checkpoint(args: argparse.Namespace) -> None:
    CheckpointPublisher(
        Path(args.run_dir), SourceRegistry(Path(args.registry)), ImmutableSourceCache(Path(args.source_root)),
    ).run()


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
    checkpoint.add_argument("--registry", required=True)
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
