#!/usr/bin/env python3
"""Download and validate the processed GSE94641 embryonic mouse MGE reference.

This independent Step 00 component stops after reference validation. It does
not load Paper 3 query cells, transfer labels, reconstruct clusters, manually
annotate cells, run differential expression, or modify frozen PCDH19 Steps
03--07. GEO's processed per-cell nRPKM/read-count files and published metadata
are cached under ``references/GSE94641``; deterministic validation tables are
atomically published under ``results/mge_reference_mapping_gse94641``.
"""

from __future__ import print_function

import argparse
import csv
import gzip
import hashlib
import json
import os
import re
import shutil
import sys
import tarfile
import tempfile
import urllib.request
from collections import Counter, OrderedDict


class ReferenceValidationError(RuntimeError):
    """Raised when a locked reference or scientific invariant fails."""


def sha256_file(path, chunk_size=8 * 1024 * 1024):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        while True:
            chunk = handle.read(chunk_size)
            if not chunk:
                break
            digest.update(chunk)
    return digest.hexdigest()


def write_tsv(path, header, rows):
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, delimiter="\t", fieldnames=header, lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in header})


class ValidationRecorder(object):
    """Collect explicit PASS checks and fail immediately on a mismatch."""

    def __init__(self):
        self.rows = []

    def require(self, name, observed, expected, details=""):
        passed = observed == expected
        self.rows.append({
            "check_name": name,
            "status": "PASS" if passed else "FAIL",
            "observed": observed,
            "expected": expected,
            "details": details,
        })
        if not passed:
            raise ReferenceValidationError(
                "{} failed: observed {!r}, expected {!r}".format(
                    name, observed, expected
                )
            )


class LockedReferenceDownloader(object):
    """Download checksum-locked GEO assets with atomic cache publication."""

    def __init__(self, reference_root, sources):
        self.reference_root = reference_root
        self.source_root = os.path.join(reference_root, "source_files")
        self.sources = sources

    def _verify(self, path, definition):
        if not os.path.isfile(path):
            raise ReferenceValidationError("Missing locked source: {}".format(path))
        size = os.path.getsize(path)
        if size != int(definition["bytes"]):
            raise ReferenceValidationError(
                "Locked size mismatch for {}: {} != {}".format(
                    path, size, definition["bytes"]
                )
            )
        observed = sha256_file(path)
        if observed != definition["sha256"]:
            raise ReferenceValidationError(
                "Locked checksum mismatch for {}".format(path)
            )

    def materialize(self):
        os.makedirs(self.source_root, exist_ok=True)
        for definition in self.sources:
            path = os.path.join(self.source_root, definition["name"])
            if os.path.exists(path):
                self._verify(path, definition)
                print("verified cached source: {}".format(definition["name"]))
                continue
            descriptor, temporary = tempfile.mkstemp(
                prefix=".download.", dir=self.source_root
            )
            os.close(descriptor)
            try:
                request = urllib.request.Request(
                    definition["url"],
                    headers={"User-Agent": "paper3-pcdh19-gse94641-validation/1.0"},
                )
                print("downloading: {}".format(definition["url"]))
                with urllib.request.urlopen(request) as response, open(temporary, "wb") as output:
                    shutil.copyfileobj(response, output, 1024 * 1024)
                    output.flush()
                    os.fsync(output.fileno())
                self._verify(temporary, definition)
                if os.path.exists(path):
                    raise ReferenceValidationError(
                        "Source appeared during download; refusing overwrite: {}".format(path)
                    )
                os.replace(temporary, path)
            finally:
                if os.path.exists(temporary):
                    os.unlink(temporary)
        return self.source_root


class GeoMetadataParser(object):
    """Recover every per-cell GEO field and named biological characteristic."""

    def __init__(self, series_matrix_path, expected_cells):
        self.path = series_matrix_path
        self.expected_cells = expected_cells
        self.field_sources = OrderedDict()
        self.characteristic_fields = []

    @staticmethod
    def _parse_row(raw):
        return next(csv.reader([raw.rstrip("\n")], delimiter="\t", quotechar='"'))

    def parse(self):
        columns = OrderedDict()
        repeated = Counter()
        with gzip.open(self.path, "rt", encoding="utf-8", errors="strict", newline="") as handle:
            for raw in handle:
                if not raw.startswith("!Sample_"):
                    continue
                parsed = self._parse_row(raw)
                values = parsed[1:]
                if len(values) != self.expected_cells:
                    raise ReferenceValidationError(
                        "Per-cell GEO field {} has {} values, expected {}".format(
                            parsed[0], len(values), self.expected_cells
                        )
                    )
                base = parsed[0][len("!Sample_"):]
                if base == "characteristics_ch1":
                    prefixes = []
                    stripped = []
                    for value in values:
                        if ":" not in value:
                            raise ReferenceValidationError(
                                "Malformed GEO characteristic: {}".format(value)
                            )
                        prefix, payload = value.split(":", 1)
                        prefixes.append(prefix.strip())
                        stripped.append(payload.strip())
                    if len(set(prefixes)) != 1:
                        raise ReferenceValidationError(
                            "A characteristics row mixes field names: {}".format(sorted(set(prefixes)))
                        )
                    field = prefixes[0]
                    if field in columns:
                        raise ReferenceValidationError("Duplicate characteristic field: {}".format(field))
                    columns[field] = stripped
                    self.field_sources[field] = "Sample_characteristics_ch1"
                    self.characteristic_fields.append(field)
                    continue
                repeated[base] += 1
                field = base if repeated[base] == 1 else "{}_{}".format(base, repeated[base])
                columns[field] = values
                self.field_sources[field] = "Sample_{}".format(base)

        required = ["geo_accession", "title", "age", "extendedphenotype", "subtype", "nonneuralcelltype"]
        missing = [field for field in required if field not in columns]
        if missing:
            raise ReferenceValidationError("Missing required GEO fields: {}".format(", ".join(missing)))
        rows = []
        for index in range(self.expected_cells):
            row = OrderedDict((field, values[index]) for field, values in columns.items())
            title_parts = row["title"].split()
            row["author_cell_id"] = title_parts[0] if title_parts else ""
            rows.append(row)
        return rows, list(columns.keys())


class ProcessedExpressionValidator(object):
    """Safely extract and validate the 225 GEO processed cell-expression files."""

    FILE_PATTERN = re.compile(r"^(GSM\d+)_sample(\d+)\.tsv\.gz$")

    def __init__(self, source_root, reference_root, expected_cells, expected_gene_rows):
        self.archive = os.path.join(source_root, "GSE94641_RAW.tar")
        self.filelist = os.path.join(source_root, "filelist.txt")
        self.expression_root = os.path.join(reference_root, "processed_expression")
        self.expected_cells = expected_cells
        self.expected_gene_rows = expected_gene_rows

    def _filelist_records(self):
        records = {}
        with open(self.filelist, "r", newline="") as handle:
            reader = csv.reader(handle, delimiter="\t")
            for row in reader:
                if len(row) >= 5 and row[0] == "File":
                    records[row[1]] = int(row[3])
        return records

    def extract(self):
        expected_sizes = self._filelist_records()
        if len(expected_sizes) != self.expected_cells:
            raise ReferenceValidationError("GEO filelist does not contain 225 processed files")
        with tarfile.open(self.archive, "r") as archive:
            members = archive.getmembers()
            names = [member.name for member in members]
            for member in members:
                if not member.isfile() or os.path.basename(member.name) != member.name:
                    raise ReferenceValidationError("Unsafe/non-file tar member: {}".format(member.name))
                if not self.FILE_PATTERN.match(member.name):
                    raise ReferenceValidationError("Unexpected tar member name: {}".format(member.name))
            if set(names) != set(expected_sizes):
                raise ReferenceValidationError("Tar members differ from GEO filelist")
            if os.path.isdir(self.expression_root):
                observed = sorted(
                    name for name in os.listdir(self.expression_root)
                    if os.path.isfile(os.path.join(self.expression_root, name))
                )
                if observed != sorted(names):
                    raise ReferenceValidationError(
                        "Existing processed-expression cache is incomplete or different"
                    )
                return
            parent = os.path.dirname(self.expression_root)
            os.makedirs(parent, exist_ok=True)
            staging = tempfile.mkdtemp(prefix=".processed-expression.", dir=parent)
            try:
                for member in members:
                    source = archive.extractfile(member)
                    destination = os.path.join(staging, member.name)
                    with source, open(destination, "wb") as output:
                        shutil.copyfileobj(source, output, 1024 * 1024)
                    if os.path.getsize(destination) != expected_sizes[member.name]:
                        raise ReferenceValidationError("Extracted size mismatch: {}".format(member.name))
                os.replace(staging, self.expression_root)
            finally:
                if os.path.isdir(staging):
                    shutil.rmtree(staging)

    def validate(self):
        filenames = sorted(os.listdir(self.expression_root))
        shared_ids = None
        manifest_rows = []
        for filename in filenames:
            match = self.FILE_PATTERN.match(filename)
            if not match:
                raise ReferenceValidationError("Unexpected extracted filename: {}".format(filename))
            path = os.path.join(self.expression_root, filename)
            identifiers = []
            nonzero_nrpkm = 0
            nonzero_counts = 0
            with gzip.open(path, "rt", encoding="utf-8", errors="strict", newline="") as handle:
                comments = [handle.readline().rstrip("\n") for _ in range(3)]
                expected_comments = [
                    "#ID_REF=ID in platform",
                    "#VALUE=Size-factor-adjusted RPKM (nRPKM)",
                    "#count=Number of reads uniquely aligning to gene model",
                ]
                if comments != expected_comments:
                    raise ReferenceValidationError("Unexpected expression metadata: {}".format(filename))
                reader = csv.DictReader(handle, delimiter="\t")
                if reader.fieldnames != ["ID_REF", "VALUE", "count"]:
                    raise ReferenceValidationError("Unexpected expression columns: {}".format(filename))
                for row in reader:
                    identifier = row["ID_REF"]
                    if not identifier.isdigit():
                        raise ReferenceValidationError("Non-numeric gene ID in {}".format(filename))
                    identifiers.append(identifier)
                    try:
                        value = float(row["VALUE"])
                        count = int(row["count"])
                    except ValueError:
                        raise ReferenceValidationError("Invalid expression value in {}".format(filename))
                    if value < 0 or count < 0:
                        raise ReferenceValidationError("Negative expression value in {}".format(filename))
                    nonzero_nrpkm += int(value > 0)
                    nonzero_counts += int(count > 0)
            if len(identifiers) != self.expected_gene_rows:
                raise ReferenceValidationError("Gene-row count mismatch in {}".format(filename))
            if len(set(identifiers)) != len(identifiers):
                raise ReferenceValidationError("Duplicate gene IDs in {}".format(filename))
            if shared_ids is None:
                shared_ids = identifiers
            elif identifiers != shared_ids:
                raise ReferenceValidationError("Gene ID/order mismatch in {}".format(filename))
            manifest_rows.append({
                "geo_accession": match.group(1),
                "sample_number": int(match.group(2)),
                "relative_path": os.path.join("processed_expression", filename),
                "gene_rows": len(identifiers),
                "nonzero_nrpkm_genes": nonzero_nrpkm,
                "nonzero_read_count_genes": nonzero_counts,
                "bytes": os.path.getsize(path),
                "sha256": sha256_file(path),
            })
        if len(manifest_rows) != self.expected_cells:
            raise ReferenceValidationError("Processed expression file count mismatch")
        manifest_rows.sort(key=lambda row: row["sample_number"])
        return shared_ids, manifest_rows


class PlatformAnnotationInspector(object):
    """Determine whether the GEO virtual platform directly supplies symbols."""

    def __init__(self, platform_soft_path):
        self.path = platform_soft_path

    def inspect(self):
        platform_fields = OrderedDict()
        table_header = None
        with gzip.open(self.path, "rt", encoding="utf-8", errors="replace") as handle:
            for raw in handle:
                if raw.startswith("!Platform_") and " = " in raw:
                    field, value = raw.rstrip("\n").split(" = ", 1)
                    if field not in platform_fields:
                        platform_fields[field] = value
                if raw.startswith("!platform_table_begin"):
                    table_header = next(handle).rstrip("\n").split("\t")
                    break
                if raw.startswith("^SAMPLE"):
                    break
        symbol_columns = []
        if table_header:
            symbol_columns = [field for field in table_header if "symbol" in field.lower()]
        return {
            "platform_table_present": table_header is not None,
            "platform_table_columns": "|".join(table_header or []),
            "gene_symbol_columns": "|".join(symbol_columns),
            "gene_symbols_directly_available": bool(symbol_columns),
            "platform_title": platform_fields.get("!Platform_title", ""),
            "platform_technology": platform_fields.get("!Platform_technology", ""),
            "platform_distribution": platform_fields.get("!Platform_distribution", ""),
        }


class ReferenceValidationPublisher(object):
    """Build and atomically publish the validation-only reference package."""

    def __init__(self, output_root, recorder):
        self.output_root = output_root
        self.recorder = recorder

    @staticmethod
    def _field_summary(metadata, fields, sources):
        rows = []
        for field in fields:
            values = [row[field] for row in metadata]
            distinct = sorted(set(values))
            display = " | ".join(distinct) if len(distinct) <= 20 else "[{} distinct; see cell metadata]".format(len(distinct))
            rows.append({
                "field": field,
                "geo_source": sources[field],
                "cells": len(values),
                "missing_or_blank": sum(value == "" for value in values),
                "distinct_values": len(distinct),
                "values_when_compact": display,
            })
        return rows

    @staticmethod
    def _annotation_tables(metadata, characteristic_fields):
        counts = []
        by_age = []
        ages = sorted(set(row["age"] for row in metadata))
        for field in characteristic_fields:
            labels = sorted(set(row[field] for row in metadata))
            for label in labels:
                cells = [row for row in metadata if row[field] == label]
                counts.append({
                    "annotation_field": field,
                    "annotation_label": label,
                    "all_cells": len(cells),
                    "e15_5_cells": sum(row["age"] == "E15.5" for row in cells),
                    "fraction_all_reference": "{:.6f}".format(len(cells) / float(len(metadata))),
                })
                for age in ages:
                    by_age.append({
                        "annotation_field": field,
                        "annotation_label": label,
                        "developmental_age": age,
                        "cells": sum(row["age"] == age for row in cells),
                    })
        return counts, by_age

    @staticmethod
    def _manifest(root):
        rows = []
        for directory, _, filenames in os.walk(root):
            for filename in sorted(filenames):
                path = os.path.join(directory, filename)
                relative = os.path.relpath(path, root)
                if relative == "output_manifest.tsv":
                    continue
                rows.append({
                    "relative_path": relative,
                    "bytes": os.path.getsize(path),
                    "sha256": sha256_file(path),
                })
        rows.sort(key=lambda row: row["relative_path"])
        return rows

    def publish(self, workflow_data):
        if os.path.isdir(self.output_root):
            manifest_path = os.path.join(self.output_root, "output_manifest.tsv")
            if not os.path.isfile(manifest_path):
                raise ReferenceValidationError("Existing validation output lacks a manifest")
            with open(manifest_path, "r", newline="") as handle:
                for row in csv.DictReader(handle, delimiter="\t"):
                    path = os.path.join(self.output_root, row["relative_path"])
                    if (not os.path.isfile(path) or os.path.getsize(path) != int(row["bytes"])
                            or sha256_file(path) != row["sha256"]):
                        raise ReferenceValidationError("Existing output fails manifest verification")
            print("verified existing validation package: {}".format(self.output_root))
            return

        parent = os.path.dirname(self.output_root)
        os.makedirs(parent, exist_ok=True)
        staging = tempfile.mkdtemp(prefix=".gse94641-validation.", dir=parent)
        try:
            metadata = workflow_data["metadata"]
            fields = workflow_data["metadata_fields"]
            characteristic_fields = workflow_data["characteristic_fields"]
            counts, by_age = self._annotation_tables(metadata, characteristic_fields)
            cell_header = ["author_cell_id"] + fields + ["expression_relative_path"]
            write_tsv(os.path.join(staging, "GSE94641_cell_metadata.tsv"), cell_header, metadata)
            write_tsv(
                os.path.join(staging, "GSE94641_published_cell_fields.tsv"),
                ["field", "geo_source", "cells", "missing_or_blank", "distinct_values", "values_when_compact"],
                self._field_summary(metadata, fields, workflow_data["field_sources"]),
            )
            write_tsv(
                os.path.join(staging, "GSE94641_annotation_counts.tsv"),
                ["annotation_field", "annotation_label", "all_cells", "e15_5_cells", "fraction_all_reference"],
                counts,
            )
            write_tsv(
                os.path.join(staging, "GSE94641_annotation_by_age.tsv"),
                ["annotation_field", "annotation_label", "developmental_age", "cells"],
                by_age,
            )
            write_tsv(
                os.path.join(staging, "GSE94641_expression_file_manifest.tsv"),
                ["geo_accession", "sample_number", "relative_path", "gene_rows", "nonzero_nrpkm_genes", "nonzero_read_count_genes", "bytes", "sha256"],
                workflow_data["expression_manifest"],
            )
            write_tsv(
                os.path.join(staging, "GSE94641_gene_identifiers.tsv"),
                ["ID_REF", "identifier_format", "gene_symbol", "gene_symbol_directly_available"],
                ({
                    "ID_REF": identifier,
                    "identifier_format": "numeric_gene_id_consistent_with_NCBI_Entrez",
                    "gene_symbol": "",
                    "gene_symbol_directly_available": "false",
                } for identifier in workflow_data["gene_ids"]),
            )
            source_rows = []
            for definition in workflow_data["lock"]["sources"]:
                source_rows.append({
                    "name": definition["name"],
                    "role": definition["role"],
                    "url": definition["url"],
                    "bytes": definition["bytes"],
                    "sha256": definition["sha256"],
                })
            write_tsv(
                os.path.join(staging, "source_manifest.tsv"),
                ["name", "role", "url", "bytes", "sha256"], source_rows,
            )
            summary = workflow_data["summary"]
            write_tsv(
                os.path.join(staging, "GSE94641_reference_summary.tsv"),
                ["property", "value", "interpretation"],
                ({"property": key, "value": value[0], "interpretation": value[1]} for key, value in summary.items()),
            )
            write_tsv(
                os.path.join(staging, "validation_checks.tsv"),
                ["check_name", "status", "observed", "expected", "details"],
                self.recorder.rows,
            )
            write_tsv(
                os.path.join(staging, "output_manifest.tsv"),
                ["relative_path", "bytes", "sha256"], self._manifest(staging),
            )
            os.replace(staging, self.output_root)
            print("published validation package: {}".format(self.output_root))
        finally:
            if os.path.isdir(staging):
                shutil.rmtree(staging)


class GSE94641ValidationWorkflow(object):
    """Coordinate independent download, metadata parsing, validation, and I/O."""

    def __init__(self, lock_path, reference_root, output_root):
        with open(lock_path, "r") as handle:
            self.lock = json.load(handle, object_pairs_hook=OrderedDict)
        self.lock_path = lock_path
        self.reference_root = reference_root
        self.output_root = output_root
        self.validation = ValidationRecorder()

    def run(self):
        downloader = LockedReferenceDownloader(self.reference_root, self.lock["sources"])
        source_root = downloader.materialize()
        expression = ProcessedExpressionValidator(
            source_root, self.reference_root,
            int(self.lock["expected_cells"]), int(self.lock["expected_gene_rows"]),
        )
        expression.extract()
        gene_ids, expression_manifest = expression.validate()

        parser = GeoMetadataParser(
            os.path.join(source_root, "GSE94641_series_matrix.txt.gz"),
            int(self.lock["expected_cells"]),
        )
        metadata, fields = parser.parse()
        expression_by_gsm = {row["geo_accession"]: row for row in expression_manifest}
        for row in metadata:
            if row["geo_accession"] not in expression_by_gsm:
                raise ReferenceValidationError("Metadata cell lacks processed expression")
            row["expression_relative_path"] = expression_by_gsm[row["geo_accession"]]["relative_path"]

        platform = PlatformAnnotationInspector(
            os.path.join(source_root, "GPL17021_family.soft.gz")
        ).inspect()
        ages = sorted(set(row["age"] for row in metadata))
        self.validation.require("accession", self.lock["accession"], "GSE94641")
        self.validation.require("cell_count", len(metadata), int(self.lock["expected_cells"]))
        self.validation.require("processed_file_count", len(expression_manifest), int(self.lock["expected_cells"]))
        self.validation.require("gene_rows", len(gene_ids), int(self.lock["expected_gene_rows"]))
        self.validation.require("developmental_ages", ages, sorted(self.lock["expected_ages"]))
        self.validation.require(
            "published_characteristic_fields", parser.characteristic_fields,
            self.lock["expected_characteristic_fields"],
        )
        self.validation.require("metadata_expression_gsm_identity", len(expression_by_gsm), len(metadata))
        self.validation.require("all_gene_ids_numeric", all(identifier.isdigit() for identifier in gene_ids), True)
        self.validation.require("platform_table_absent", platform["platform_table_present"], False)
        self.validation.require("gene_symbols_not_directly_available", platform["gene_symbols_directly_available"], False)

        summary = OrderedDict([
            ("reference_accession", (self.lock["accession"], "NCBI GEO series")),
            ("cells", (len(metadata), "one processed file and one GEO metadata column per cell")),
            ("developmental_ages", ("|".join(ages), "published age characteristic")),
            ("e15_5_cells", (sum(row["age"] == "E15.5" for row in metadata), "cells matching the query's approximate developmental stage")),
            ("expression_matrix_dimensions", ("{} genes x {} cells".format(len(gene_ids), len(metadata)), "distributed by GEO as one gene-by-measurement file per cell")),
            ("expression_measurements", ("VALUE=nRPKM|count=uniquely-aligning reads", "published processed quantities; not raw FASTQ")),
            ("gene_identifier_format", ("numeric_gene_id_consistent_with_NCBI_Entrez", "ID_REF values are numeric; GEO does not provide a platform annotation table here")),
            ("gene_symbols_directly_available", (str(platform["gene_symbols_directly_available"]).lower(), "GPL17021 is a virtual platform without a gene-symbol annotation table in this deposit")),
            ("published_characteristic_fields", ("|".join(parser.characteristic_fields), "all Sample_characteristics_ch1 fields")),
            ("label_transfer_status", ("NOT_RUN", "critical stopping point after download and validation")),
            ("query_data_status", ("NOT_LOADED", "PCDH19 query cells and classifications were not accessed")),
        ])
        publisher = ReferenceValidationPublisher(self.output_root, self.validation)
        publisher.publish({
            "lock": self.lock,
            "metadata": metadata,
            "metadata_fields": fields,
            "field_sources": parser.field_sources,
            "characteristic_fields": parser.characteristic_fields,
            "gene_ids": gene_ids,
            "expression_manifest": expression_manifest,
            "summary": summary,
        })
        print("cells={} genes={} E15.5={}".format(
            len(metadata), len(gene_ids), sum(row["age"] == "E15.5" for row in metadata)
        ))
        print("label_transfer=NOT_RUN")


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--lock", required=True)
    parser.add_argument("--reference-root", required=True)
    parser.add_argument("--output-root", required=True)
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    GSE94641ValidationWorkflow(
        os.path.abspath(args.lock),
        os.path.abspath(args.reference_root),
        os.path.abspath(args.output_root),
    ).run()


if __name__ == "__main__":
    try:
        main()
    except ReferenceValidationError as error:
        print("ERROR: {}".format(error), file=sys.stderr)
        sys.exit(2)
