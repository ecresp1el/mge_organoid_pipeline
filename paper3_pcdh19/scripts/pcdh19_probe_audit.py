#!/usr/bin/env python3
"""Locked, technical-only Pcdh19 Flex probe audit for AGC request 15662-JZ."""

import argparse
import csv
import datetime as dt
import gzip
import hashlib
import json
import os
import platform
import re
import shutil
import sys
import tempfile
import urllib.request
from collections import OrderedDict

import h5py
import numpy as np


PATTERN_ORDER = ["A+B+C", "A+B", "A+C", "A only", "B+C", "B only", "C only", "none"]
PATTERN_HEADER = ["barcode", "A_UMI", "B_UMI", "C_UMI", "Pcdh19_total_UMI", "detection_pattern"]
SUMMARY_HEADER = [
    "sample_id", "summary_type", "category", "probe_label", "probe_id",
    "total_umi", "cells", "detected_cells", "denominator_cells", "fraction",
    "percent", "percent_among_any_pcdh19", "mean_umi_among_detected",
    "max_umi_per_cell", "median_total_pcdh19_umi", "mean_total_pcdh19_umi",
    "max_total_pcdh19_umi",
]
VALIDATION_HEADER = ["sample_id", "check_name", "status", "observed", "expected", "details"]
COORDINATE_HEADER = [
    "probe_label", "probe_suffix", "probe_id", "probe_seq", "10x_region", "chrom",
    "chrom_start_0based", "chrom_end_0based", "strand_10x_bed", "bed_block_count",
    "bed_block_sizes", "bed_block_starts", "transcript_id_set",
    "target_transcript_exon_assignment", "other_pcdh19_exon_intersections",
    "pcdh19_gene_strand", "five_prime_order",
]


def fail(message):
    raise RuntimeError(message)


def sha256_file(path, chunk_size=8 * 1024 * 1024):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        while True:
            chunk = handle.read(chunk_size)
            if not chunk:
                break
            digest.update(chunk)
    return digest.hexdigest()


def require_sha256(path, expected, label):
    if not os.path.isfile(path):
        fail("{} is missing: {}".format(label, path))
    observed = sha256_file(path)
    if observed != expected:
        fail("{} checksum mismatch: observed {}, expected {}, path {}".format(label, observed, expected, path))
    return observed


def download_locked(url, path, expected_sha256):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    if os.path.exists(path):
        require_sha256(path, expected_sha256, "cached reference")
        return
    fd, temporary = tempfile.mkstemp(prefix=".download.", dir=os.path.dirname(path))
    os.close(fd)
    try:
        request = urllib.request.Request(
            url,
            headers={"User-Agent": "Mozilla/5.0 paper3-pcdh19-probe-audit/1.0"},
        )
        with urllib.request.urlopen(request) as response, open(temporary, "wb") as output:
            shutil.copyfileobj(response, output, 1024 * 1024)
            output.flush()
            os.fsync(output.fileno())
        require_sha256(temporary, expected_sha256, "downloaded reference")
        if os.path.exists(path):
            fail("Reference appeared during download and will not be overwritten: {}".format(path))
        os.replace(temporary, path)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


def write_tsv(path, header, rows):
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=header, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in header})


def publish_file_locked(temporary, destination):
    if os.path.exists(destination):
        if sha256_file(temporary) == sha256_file(destination):
            os.unlink(temporary)
            return "EXISTING_IDENTICAL"
        fail("Refusing to overwrite a different validated output: {}".format(destination))
    os.replace(temporary, destination)
    return "PUBLISHED"


def read_delivered_probe_set(path):
    metadata = OrderedDict()
    data_lines = []
    with open(path, "r", newline="") as handle:
        for line in handle:
            if line.startswith("#"):
                key, value = line[1:].rstrip("\r\n").split("=", 1)
                metadata[key] = value
            else:
                data_lines.append(line)
    return metadata, list(csv.DictReader(data_lines))


def parse_gtf_attributes(text):
    attributes = {}
    for item in text.split(";"):
        item = item.strip()
        if not item:
            continue
        fields = item.split(None, 1)
        if len(fields) == 2:
            attributes[fields[0]] = fields[1].strip().strip('"')
    return attributes


def read_pcdh19_gtf(path, gene_id):
    gene = None
    exons = []
    with gzip.open(path, "rt") as handle:
        for line in handle:
            if 'gene_id "{}"'.format(gene_id) not in line:
                continue
            fields = line.rstrip("\n").split("\t")
            attrs = parse_gtf_attributes(fields[8])
            if fields[2] == "gene":
                gene = {
                    "chrom": fields[0], "start": int(fields[3]) - 1, "end": int(fields[4]),
                    "strand": fields[6], "gene_name": attrs.get("gene_name", ""),
                }
            elif fields[2] == "exon":
                exon_number = attrs.get("exon_number")
                exons.append({
                    "chrom": fields[0], "start": int(fields[3]) - 1, "end": int(fields[4]),
                    "strand": fields[6], "transcript_id": attrs["transcript_id"],
                    "transcript_name": attrs.get("transcript_name", ""),
                    "exon_number": int(exon_number) if exon_number is not None else None,
                    "exon_id": attrs.get("exon_id", ""),
                })
    if gene is None:
        fail("Gene {} was not found in {}".format(gene_id, path))
    return gene, exons


def validate_reference_and_write(lock, cellranger_root, paper3_root):
    output_root = os.path.join(paper3_root, "results", "pcdh19_probe_audit")
    reference_output = os.path.join(output_root, "references")
    reference_cache = os.path.join(paper3_root, "inputs", "pcdh19_probe_audit", "references")
    os.makedirs(reference_output, exist_ok=True)
    os.makedirs(reference_cache, exist_ok=True)

    refs = lock["reference_inputs"]
    delivered = os.path.join(cellranger_root, refs["delivered_probe_set"]["basename"])
    gtf = os.environ.get("PAPER3_GRCM39_GTF", "")
    reference_json = os.environ.get("PAPER3_GRCM39_REFERENCE_JSON", "")
    if not gtf or not reference_json:
        fail("PAPER3_GRCM39_GTF and PAPER3_GRCM39_REFERENCE_JSON must be configured")
    bed = os.path.join(reference_cache, refs["tenx_bed"]["basename"])
    metadata_tsv = os.path.join(reference_cache, refs["tenx_metadata"]["basename"])

    download_locked(refs["tenx_bed"]["url"], bed, refs["tenx_bed"]["sha256"])
    download_locked(refs["tenx_metadata"]["url"], metadata_tsv, refs["tenx_metadata"]["sha256"])
    require_sha256(delivered, refs["delivered_probe_set"]["sha256"], "delivered probe_set.csv")
    require_sha256(gtf, refs["grcm39_gtf"]["sha256"], "GRCm39 GTF")
    require_sha256(reference_json, refs["grcm39_reference_json"]["sha256"], "GRCm39 reference.json")

    with open(reference_json, "r") as handle:
        reference_info = json.load(handle)
    reference_lock = refs["grcm39_reference_json"]
    if reference_info.get("version") != reference_lock["version"]:
        fail("Reference version mismatch")
    if reference_info.get("mkref_version") != reference_lock["mkref_version"]:
        fail("mkref version mismatch")
    if reference_info.get("input_gtf_files") != [reference_lock["input_gtf"]]:
        fail("Reference input GTF mismatch")

    panel_lock = lock["panel"]
    panel_header, delivered_rows = read_delivered_probe_set(delivered)
    expected_header = {
        "probe_set_file_format": panel_lock["probe_set_file_format"],
        "panel_name": panel_lock["panel_name"],
        "reference_genome": panel_lock["reference_genome"],
        "reference_version": panel_lock["reference_version"],
    }
    for key, expected in expected_header.items():
        if panel_header.get(key) != expected:
            fail("Delivered panel header {} mismatch: {} != {}".format(key, panel_header.get(key), expected))

    expected_ids = {probe["probe_id"] for probe in lock["probes"]}
    delivered_by_id = {row["probe_id"]: row for row in delivered_rows if row["gene_id"] == panel_lock["gene_id"]}
    if set(delivered_by_id) != expected_ids:
        fail("Delivered Pcdh19 probe IDs differ from the locked three-probe set")

    metadata_by_id = {}
    with open(metadata_tsv, "r", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["gene_id"] == panel_lock["gene_id"]:
                if row["probe_id"] in metadata_by_id:
                    fail("Duplicate metadata probe ID: {}".format(row["probe_id"]))
                metadata_by_id[row["probe_id"]] = row
    if set(metadata_by_id) != expected_ids:
        fail("10x metadata Pcdh19 probe IDs differ from the locked three-probe set")

    bed_by_id = {}
    with open(bed, "r") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 12 or fields[3] not in expected_ids:
                continue
            if fields[3] in bed_by_id:
                fail("Duplicate BED probe ID: {}".format(fields[3]))
            bed_by_id[fields[3]] = fields
    if set(bed_by_id) != expected_ids:
        fail("10x BED Pcdh19 probe IDs differ from the locked three-probe set")

    gene, exons = read_pcdh19_gtf(gtf, panel_lock["gene_id"])
    if gene["chrom"] != "chrX" or gene["strand"] != panel_lock["gene_strand"] or gene["gene_name"] != panel_lock["gene_name"]:
        fail("Locked Pcdh19 gene annotation does not match the GTF")

    coordinate_rows = []
    derived_order = sorted(lock["probes"], key=lambda p: p["start"], reverse=(gene["strand"] == "-"))
    if [p["label"] for p in derived_order] != [p["label"] for p in sorted(lock["probes"], key=lambda p: p["five_prime_order"])]:
        fail("Locked 5-prime probe order does not match the GTF gene strand")

    for probe in sorted(lock["probes"], key=lambda p: p["five_prime_order"]):
        probe_id = probe["probe_id"]
        delivered_row = delivered_by_id[probe_id]
        metadata_row = metadata_by_id[probe_id]
        bed_row = bed_by_id[probe_id]
        for source_name, value in [("delivered sequence", delivered_row["probe_seq"]), ("metadata sequence", metadata_row["probe_seq"])]:
            if value != probe["probe_seq"]:
                fail("{} mismatch for {}".format(source_name, probe_id))
        if delivered_row["region"] != probe["region"] or metadata_row["region"] != probe["region"]:
            fail("10x region mismatch for {}".format(probe_id))
        transcript_ids = metadata_row["transcript_id_set"].split(";")
        if transcript_ids != probe["transcript_id_set"]:
            fail("transcript_id_set mismatch for {}".format(probe_id))
        observed_bed = (bed_row[0], int(bed_row[1]), int(bed_row[2]), bed_row[5])
        expected_bed = (probe["chrom"], probe["start"], probe["end"], probe["bed_strand"])
        if observed_bed != expected_bed:
            fail("BED coordinate mismatch for {}: {} != {}".format(probe_id, observed_bed, expected_bed))
        if bed_row[9:] != ["1", "50", "0"]:
            fail("Expected one 50-base BED12 block for {}".format(probe_id))

        target_assignments = []
        for transcript_id in transcript_ids:
            full = [exon for exon in exons if exon["transcript_id"] == transcript_id and exon["start"] <= probe["start"] and exon["end"] >= probe["end"]]
            if len(full) != 1:
                fail("Probe {} is not entirely within exactly one exon of {}".format(probe_id, transcript_id))
            exon_number = full[0]["exon_number"]
            if exon_number != probe["expected_exons"][transcript_id]:
                fail("Exon assignment mismatch for {} in {}".format(probe_id, transcript_id))
            target_assignments.append("{}:exon{}".format(transcript_id, exon_number))

        other_intersections = []
        for exon in exons:
            if exon["transcript_id"] in transcript_ids:
                continue
            overlap = max(0, min(probe["end"], exon["end"]) - max(probe["start"], exon["start"]))
            if overlap:
                exon_label = exon["exon_number"] if exon["exon_number"] is not None else "NA"
                other_intersections.append("{}:exon{}:overlap_bp={}".format(exon["transcript_id"], exon_label, overlap))

        coordinate_rows.append({
            "probe_label": probe["label"], "probe_suffix": probe["suffix"], "probe_id": probe_id,
            "probe_seq": probe["probe_seq"], "10x_region": probe["region"], "chrom": probe["chrom"],
            "chrom_start_0based": probe["start"], "chrom_end_0based": probe["end"],
            "strand_10x_bed": probe["bed_strand"], "bed_block_count": bed_row[9],
            "bed_block_sizes": bed_row[10], "bed_block_starts": bed_row[11],
            "transcript_id_set": ";".join(transcript_ids),
            "target_transcript_exon_assignment": ";".join(target_assignments),
            "other_pcdh19_exon_intersections": ";".join(sorted(set(other_intersections))),
            "pcdh19_gene_strand": gene["strand"], "five_prime_order": probe["five_prime_order"],
        })

    coordinate_target = os.path.join(reference_output, "pcdh19_probe_coordinates.tsv")
    fd, coordinate_tmp = tempfile.mkstemp(prefix=".coordinates.", dir=reference_output, text=True)
    os.close(fd)
    write_tsv(coordinate_tmp, COORDINATE_HEADER, coordinate_rows)
    coordinate_sha = sha256_file(coordinate_tmp)
    publish_file_locked(coordinate_tmp, coordinate_target)

    manifest_rows = []
    inputs = [
        ("delivered_probe_set", "reference_input", panel_lock["panel_name"], "delivered", delivered, refs["delivered_probe_set"]["sha256"]),
        ("tenx_bed", "reference_input", panel_lock["panel_name"], refs["tenx_bed"]["url"], bed, refs["tenx_bed"]["sha256"]),
        ("tenx_metadata", "reference_input", panel_lock["panel_name"], refs["tenx_metadata"]["url"], metadata_tsv, refs["tenx_metadata"]["sha256"]),
        ("grcm39_gtf", "reference_input", "GRCm39-2024-A;GENCODE-vM33", "local_10x_reference", gtf, refs["grcm39_gtf"]["sha256"]),
        ("grcm39_reference_json", "reference_input", "GRCm39-2024-A;mkref-8.0.0", "local_10x_reference", reference_json, refs["grcm39_reference_json"]["sha256"]),
        ("pcdh19_probe_coordinates", "derived_reference", lock["pipeline_version"], "pipeline_derived", coordinate_target, coordinate_sha),
    ]
    for input_id, role, version, source, path, checksum in inputs:
        manifest_rows.append({
            "input_id": input_id, "role": role, "version": version, "source_uri": source,
            "resolved_path": path, "sha256": checksum, "bytes": os.path.getsize(path), "validation_status": "PASS",
        })
    manifest_header = ["input_id", "role", "version", "source_uri", "resolved_path", "sha256", "bytes", "validation_status"]
    manifest_target = os.path.join(reference_output, "reference_manifest.tsv")
    fd, manifest_tmp = tempfile.mkstemp(prefix=".manifest.", dir=reference_output, text=True)
    os.close(fd)
    write_tsv(manifest_tmp, manifest_header, manifest_rows)
    publish_file_locked(manifest_tmp, manifest_target)
    return {"bed": bed, "metadata": metadata_tsv, "gtf": gtf, "delivered": delivered}


def positions_for_rows(dataset, rows, chunk_size=5_000_000):
    positions = []
    matched_rows = []
    total = dataset.shape[0]
    for start in range(0, total, chunk_size):
        stop = min(start + chunk_size, total)
        indices = dataset[start:stop]
        mask = np.zeros(indices.shape, dtype=np.bool_)
        for row in rows:
            mask |= indices == row
        local = np.flatnonzero(mask)
        if local.size:
            positions.append(local.astype(np.int64) + start)
            matched_rows.append(indices[local])
    if not positions:
        return np.empty(0, dtype=np.int64), np.empty(0, dtype=np.int64)
    return np.concatenate(positions), np.concatenate(matched_rows)


def extract_sparse_rows(matrix, rows):
    indptr = matrix["indptr"][:]
    positions, matched = positions_for_rows(matrix["indices"], rows)
    values = matrix["data"][positions].astype(np.int64)
    columns = np.searchsorted(indptr, positions, side="right") - 1
    output = np.zeros((len(rows), indptr.size - 1), dtype=np.int64)
    for slot, row in enumerate(rows):
        selected = matched == row
        output[slot, columns[selected]] = values[selected]
    return output


def check_row(sample_id, name, observed, expected, details=""):
    status = "PASS" if observed == expected else "FAIL"
    return {
        "sample_id": sample_id, "check_name": name, "status": status,
        "observed": observed, "expected": expected, "details": details,
    }


def format_float(value):
    return "{:.12f}".format(float(value))


def sample_paths(cellranger_root, sample_id):
    directory = os.path.join(cellranger_root, "per_sample_outs", sample_id)
    return {
        "directory": directory,
        "probe_h5": os.path.join(directory, "sample_raw_probe_bc_matrix.h5"),
        "gene_h5": os.path.join(directory, "sample_filtered_feature_bc_matrix.h5"),
        "barcodes_csv": os.path.join(directory, "sample_filtered_barcodes.csv"),
    }


def verify_existing_sample(directory, prototype=None):
    required = ["pcdh19_probe_patterns.tsv", "pcdh19_probe_summary.tsv", "validation.tsv", "checksums.sha256"]
    if not os.path.isdir(directory):
        return False
    for name in required:
        if not os.path.isfile(os.path.join(directory, name)):
            fail("Existing sample output is incomplete and will not be overwritten: {}".format(directory))
    expected_checksums = {}
    with open(os.path.join(directory, "checksums.sha256"), "r") as handle:
        for line in handle:
            checksum, name = line.rstrip("\n").split("  ", 1)
            expected_checksums[name] = checksum
    for name in required[:-1]:
        observed = sha256_file(os.path.join(directory, name))
        if expected_checksums.get(name) != observed:
            fail("Existing output checksum failure: {}/{}".format(directory, name))
    with open(os.path.join(directory, "validation.tsv"), "r", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows or any(row["status"] != "PASS" for row in rows):
        fail("Existing validation is not fully PASS: {}".format(directory))
    if prototype:
        observed = sha256_file(os.path.join(directory, "pcdh19_probe_patterns.tsv"))
        if observed != prototype["pattern_table_sha256"]:
            fail("Existing JZ-1 prototype table is not byte-equivalent to the lock")
    return True


def process_sample(lock, cellranger_root, paper3_root, sample_id, prototype_gate=False):
    output_root = os.path.join(paper3_root, "results", "pcdh19_probe_audit")
    per_sample_root = os.path.join(output_root, "per_sample")
    os.makedirs(per_sample_root, exist_ok=True)
    destination = os.path.join(per_sample_root, sample_id)
    prototype = lock["prototype"] if prototype_gate else None
    if verify_existing_sample(destination, prototype):
        print("{}: existing validated output verified; no overwrite".format(sample_id), flush=True)
        return "EXISTING_VALIDATED"

    paths = sample_paths(cellranger_root, sample_id)
    for key in ["probe_h5", "gene_h5", "barcodes_csv"]:
        if not os.path.isfile(paths[key]):
            fail("{} input missing for {}: {}".format(key, sample_id, paths[key]))
    staging = tempfile.mkdtemp(prefix=".{}.tmp.".format(sample_id), dir=per_sample_root)
    try:
        probe_ids = [probe["probe_id"].encode() for probe in lock["probes"]]
        with h5py.File(paths["probe_h5"], "r") as handle:
            matrix = handle["matrix"]
            feature_ids = matrix["features/id"][:]
            matching = [int(i) for i, value in enumerate(feature_ids) if value in set(probe_ids)]
            if len(matching) != 3 or len(set(matching)) != 3:
                fail("{} does not contain exactly the three locked probe features".format(sample_id))
            lookup = {value: i for i, value in enumerate(feature_ids)}
            rows = np.array([lookup[probe_id] for probe_id in probe_ids], dtype=np.int64)
            target_sets = matrix["features/target_sets"]
            if list(target_sets.keys()) != [lock["panel"]["panel_name"]]:
                fail("{} probe matrix target-set version is not the locked panel".format(sample_id))
            for row, probe in zip(rows, lock["probes"]):
                if matrix["features/gene_id"][row].decode() != lock["panel"]["gene_id"]:
                    fail("{} probe feature gene ID mismatch".format(sample_id))
                if matrix["features/gene_name"][row].decode() != lock["panel"]["gene_name"]:
                    fail("{} probe feature gene name mismatch".format(sample_id))
                if matrix["features/probe_region"][row].decode() != probe["region"]:
                    fail("{} probe feature region mismatch for {}".format(sample_id, probe["probe_id"]))
            filtered_flags = matrix["filtered_barcodes"][:].astype(bool)
            filtered_columns = np.flatnonzero(filtered_flags)
            probe_barcodes = matrix["barcodes"][:][filtered_columns]
            probe_counts = extract_sparse_rows(matrix, rows)[:, filtered_columns]

        genomes = []
        csv_barcodes_list = []
        with open(paths["barcodes_csv"], "r", newline="") as handle:
            for row in csv.reader(handle):
                if len(row) != 2:
                    fail("Malformed filtered barcode row for {}".format(sample_id))
                genomes.append(row[0])
                csv_barcodes_list.append(row[1].encode())
        csv_barcodes = np.array(csv_barcodes_list, dtype=probe_barcodes.dtype)

        with h5py.File(paths["gene_h5"], "r") as handle:
            matrix = handle["matrix"]
            gene_barcodes = matrix["barcodes"][:]
            gene_ids = matrix["features/id"][:]
            hits = np.flatnonzero(gene_ids == lock["panel"]["gene_id"].encode())
            if hits.size != 1:
                fail("{} gene matrix does not contain exactly one Pcdh19 feature".format(sample_id))
            if matrix["features/genome"][int(hits[0])].decode() != lock["panel"]["reference_genome"]:
                fail("{} Pcdh19 gene feature is not GRCm39".format(sample_id))
            if handle.attrs.get("chemistry_description", b"").decode() != "GEM-X Flex v2":
                fail("{} gene matrix chemistry is not GEM-X Flex v2".format(sample_id))
            gene_counts = extract_sparse_rows(matrix, np.array([int(hits[0])], dtype=np.int64))[0]

        total_cells = len(csv_barcodes)
        A, B, C = probe_counts
        reconstructed = A + B + C
        a, b, c = A > 0, B > 0, C > 0
        patterns = np.empty(total_cells, dtype="U6")
        patterns[a & b & c] = "A+B+C"
        patterns[a & b & ~c] = "A+B"
        patterns[a & ~b & c] = "A+C"
        patterns[a & ~b & ~c] = "A only"
        patterns[~a & b & c] = "B+C"
        patterns[~a & b & ~c] = "B only"
        patterns[~a & ~b & c] = "C only"
        patterns[~a & ~b & ~c] = "none"
        pattern_counts = OrderedDict((pattern, int(np.count_nonzero(patterns == pattern))) for pattern in PATTERN_ORDER)
        any_detected = int(np.count_nonzero(reconstructed > 0))
        downstream = b | c
        downstream_cells = int(np.count_nonzero(downstream))
        downstream_A_positive = int(np.count_nonzero(downstream & a))
        downstream_A_negative = int(np.count_nonzero(downstream & ~a))
        bc_A_negative = int(np.count_nonzero(~a & b & c))

        validations = []
        validations.append(check_row(sample_id, "filtered_barcode_count_probe_vs_csv", len(probe_barcodes), total_cells))
        validations.append(check_row(sample_id, "filtered_barcode_count_gene_vs_csv", len(gene_barcodes), total_cells))
        validations.append(check_row(sample_id, "filtered_barcodes_probe_vs_csv", bool(np.array_equal(probe_barcodes, csv_barcodes)), True))
        validations.append(check_row(sample_id, "filtered_barcodes_gene_vs_csv", bool(np.array_equal(gene_barcodes, csv_barcodes)), True))
        validations.append(check_row(sample_id, "filtered_barcode_genome", sorted(set(genomes)), [lock["panel"]["reference_genome"]]))
        validations.append(check_row(sample_id, "duplicate_filtered_barcodes", total_cells - len(set(csv_barcodes_list)), 0))
        validations.append(check_row(sample_id, "per_barcode_probe_sum_matches_gene", int(np.count_nonzero(reconstructed != gene_counts)), 0))
        validations.append(check_row(sample_id, "reconstructed_total_matches_gene_total", int(reconstructed.sum()), int(gene_counts.sum())))
        validations.append(check_row(sample_id, "pattern_partition_matches_filtered_cells", sum(pattern_counts.values()), total_cells))
        for label, counts, detected in zip(["A", "B", "C"], [A, B, C], [a, b, c]):
            pattern_marginal = sum(count for pattern, count in pattern_counts.items() if label in pattern.replace("only", ""))
            validations.append(check_row(sample_id, "{}_pattern_marginal_detected".format(label), pattern_marginal, int(np.count_nonzero(detected))))
            validations.append(check_row(sample_id, "{}_summary_total_umi".format(label), int(counts.sum()), int(probe_counts[["A", "B", "C"].index(label)].sum())))

        pattern_path = os.path.join(staging, "pcdh19_probe_patterns.tsv")
        with open(pattern_path, "w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(PATTERN_HEADER)
            for barcode, av, bv, cv, total, pattern in zip(csv_barcodes, A, B, C, reconstructed, patterns):
                writer.writerow([barcode.decode(), int(av), int(bv), int(cv), int(total), pattern])

        serialized_rows = 0
        serialized_bad_sums = 0
        serialized_duplicates = 0
        serialized_seen = set()
        serialized_patterns = OrderedDict((pattern, 0) for pattern in PATTERN_ORDER)
        with open(pattern_path, "r", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if reader.fieldnames != PATTERN_HEADER:
                fail("Serialized pattern table header mismatch for {}".format(sample_id))
            for row in reader:
                serialized_rows += 1
                barcode = row["barcode"]
                if barcode in serialized_seen:
                    serialized_duplicates += 1
                serialized_seen.add(barcode)
                if int(row["A_UMI"]) + int(row["B_UMI"]) + int(row["C_UMI"]) != int(row["Pcdh19_total_UMI"]):
                    serialized_bad_sums += 1
                serialized_patterns[row["detection_pattern"]] += 1
        validations.append(check_row(sample_id, "serialized_pattern_rows", serialized_rows, total_cells))
        validations.append(check_row(sample_id, "serialized_per_barcode_sum_failures", serialized_bad_sums, 0))
        validations.append(check_row(sample_id, "serialized_duplicate_barcodes", serialized_duplicates, 0))
        validations.append(check_row(sample_id, "serialized_pattern_counts", dict(serialized_patterns), dict(pattern_counts)))

        summary_rows = []
        for probe, counts, detected in zip(lock["probes"], [A, B, C], [a, b, c]):
            detected_cells = int(np.count_nonzero(detected))
            summary_rows.append({
                "sample_id": sample_id, "summary_type": "probe", "category": probe["label"],
                "probe_label": probe["label"], "probe_id": probe["probe_id"], "total_umi": int(counts.sum()),
                "detected_cells": detected_cells, "denominator_cells": total_cells,
                "fraction": format_float(detected_cells / float(total_cells)),
                "percent": format_float(100.0 * detected_cells / float(total_cells)),
                "mean_umi_among_detected": format_float(counts.sum() / float(detected_cells)) if detected_cells else format_float(0),
                "max_umi_per_cell": int(counts.max()) if total_cells else 0,
            })
        summary_rows.append({
            "sample_id": sample_id, "summary_type": "gene_reconstruction", "category": "Pcdh19",
            "total_umi": int(reconstructed.sum()), "cells": any_detected, "detected_cells": any_detected,
            "denominator_cells": total_cells, "fraction": format_float(any_detected / float(total_cells)),
            "percent": format_float(100.0 * any_detected / float(total_cells)),
        })
        for pattern in PATTERN_ORDER:
            mask = patterns == pattern
            values = reconstructed[mask]
            cells = int(mask.sum())
            summary_rows.append({
                "sample_id": sample_id, "summary_type": "pattern", "category": pattern, "cells": cells,
                "denominator_cells": total_cells, "fraction": format_float(cells / float(total_cells)),
                "percent": format_float(100.0 * cells / float(total_cells)),
                "percent_among_any_pcdh19": "" if pattern == "none" or not any_detected else format_float(100.0 * cells / float(any_detected)),
                "median_total_pcdh19_umi": format_float(np.median(values)) if cells else "",
                "mean_total_pcdh19_umi": format_float(np.mean(values)) if cells else "",
                "max_total_pcdh19_umi": int(values.max()) if cells else "",
            })
        for category, cells in [("A_positive", downstream_A_positive), ("A_negative", downstream_A_negative)]:
            summary_rows.append({
                "sample_id": sample_id, "summary_type": "downstream_A_status", "category": category,
                "cells": cells, "denominator_cells": downstream_cells,
                "fraction": format_float(cells / float(downstream_cells)) if downstream_cells else format_float(0),
                "percent": format_float(100.0 * cells / float(downstream_cells)) if downstream_cells else format_float(0),
            })
        summary_rows.append({
            "sample_id": sample_id, "summary_type": "downstream_BC_A_negative", "category": "B+C_A-negative",
            "cells": bc_A_negative, "denominator_cells": downstream_cells,
            "fraction": format_float(bc_A_negative / float(downstream_cells)) if downstream_cells else format_float(0),
            "percent": format_float(100.0 * bc_A_negative / float(downstream_cells)) if downstream_cells else format_float(0),
        })
        summary_path = os.path.join(staging, "pcdh19_probe_summary.tsv")
        write_tsv(summary_path, SUMMARY_HEADER, summary_rows)

        pattern_sha = sha256_file(pattern_path)
        validations.append(check_row(sample_id, "pattern_table_sha256_recorded", pattern_sha, pattern_sha))
        if prototype_gate:
            expected = lock["prototype"]
            validations.append(check_row(sample_id, "prototype_filtered_cells", total_cells, expected["filtered_cells"]))
            validations.append(check_row(sample_id, "prototype_pcdh19_total_umi", int(reconstructed.sum()), expected["pcdh19_total_umi"]))
            validations.append(check_row(sample_id, "prototype_pattern_table_byte_equivalence", pattern_sha, expected["pattern_table_sha256"]))
            for label, counts, detected in zip(["A", "B", "C"], [A, B, C], [a, b, c]):
                validations.append(check_row(sample_id, "prototype_{}_total_umi".format(label), int(counts.sum()), expected["probe_total_umi"][label]))
                validations.append(check_row(sample_id, "prototype_{}_detected_cells".format(label), int(np.count_nonzero(detected)), expected["probe_detected_cells"][label]))
            for pattern in PATTERN_ORDER:
                validations.append(check_row(sample_id, "prototype_pattern_{}".format(pattern.replace(" ", "_").replace("+", "plus")), pattern_counts[pattern], expected["pattern_counts"][pattern]))

        failed = [row for row in validations if row["status"] != "PASS"]
        validation_path = os.path.join(staging, "validation.tsv")
        write_tsv(validation_path, VALIDATION_HEADER, validations)
        if failed:
            fail("{} validation failed: {}".format(sample_id, ", ".join(row["check_name"] for row in failed)))

        checksums = []
        for name in ["pcdh19_probe_patterns.tsv", "pcdh19_probe_summary.tsv", "validation.tsv"]:
            checksums.append((sha256_file(os.path.join(staging, name)), name))
        with open(os.path.join(staging, "checksums.sha256"), "w") as handle:
            for checksum, name in checksums:
                handle.write("{}  {}\n".format(checksum, name))
        if os.path.exists(destination):
            fail("Destination appeared during processing and will not be overwritten: {}".format(destination))
        os.replace(staging, destination)
        staging = None
        verify_existing_sample(destination, prototype)
        print("{}: published {} cells, {} Pcdh19 UMIs".format(sample_id, total_cells, int(reconstructed.sum())), flush=True)
        return "PUBLISHED"
    finally:
        if staging and os.path.isdir(staging):
            shutil.rmtree(staging)


def combine_outputs(lock, paper3_root, sample_ids):
    output_root = os.path.join(paper3_root, "results", "pcdh19_probe_audit")
    combined = os.path.join(output_root, "combined")
    os.makedirs(combined, exist_ok=True)
    probe_rows = []
    pattern_rows = []
    validation_rows = []
    for sample_id in sample_ids:
        directory = os.path.join(output_root, "per_sample", sample_id)
        verify_existing_sample(directory, lock["prototype"] if sample_id == lock["prototype"]["sample_id"] else None)
        with open(os.path.join(directory, "pcdh19_probe_summary.tsv"), "r", newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                if row["summary_type"] in ("probe", "gene_reconstruction"):
                    probe_rows.append(row)
                else:
                    pattern_rows.append(row)
        with open(os.path.join(directory, "validation.tsv"), "r", newline="") as handle:
            validation_rows.extend(csv.DictReader(handle, delimiter="\t"))
    outputs = [
        ("all_samples_pcdh19_probe_summary.tsv", SUMMARY_HEADER, probe_rows),
        ("all_samples_pcdh19_pattern_summary.tsv", SUMMARY_HEADER, pattern_rows),
        ("all_samples_validation.tsv", VALIDATION_HEADER, validation_rows),
    ]
    for name, header, rows in outputs:
        fd, temporary = tempfile.mkstemp(prefix=".{}.".format(name), dir=combined, text=True)
        os.close(fd)
        write_tsv(temporary, header, rows)
        publish_file_locked(temporary, os.path.join(combined, name))


def write_environment(lock_path, script_path, paper3_root):
    output_root = os.path.join(paper3_root, "results", "pcdh19_probe_audit")
    os.makedirs(output_root, exist_ok=True)
    rows = [
        {"key": "pipeline_version", "value": "1.0.0"},
        {"key": "python", "value": sys.version.replace("\n", " ")},
        {"key": "h5py", "value": h5py.__version__},
        {"key": "numpy", "value": np.__version__},
        {"key": "pipeline_script_sha256", "value": sha256_file(script_path)},
        {"key": "reference_lock_sha256", "value": sha256_file(lock_path)},
    ]
    target = os.path.join(output_root, "software_environment.tsv")
    fd, temporary = tempfile.mkstemp(prefix=".environment.", dir=output_root, text=True)
    os.close(fd)
    write_tsv(temporary, ["key", "value"], rows)
    publish_file_locked(temporary, target)


def write_output_manifest(paper3_root):
    output_root = os.path.join(paper3_root, "results", "pcdh19_probe_audit")
    target = os.path.join(output_root, "output_manifest.tsv")
    rows = []
    for root, directories, files in os.walk(output_root):
        directories[:] = sorted(name for name in directories if not name.startswith("."))
        for name in sorted(files):
            path = os.path.join(root, name)
            if path == target or name.startswith("."):
                continue
            rows.append({
                "relative_path": os.path.relpath(path, output_root), "bytes": os.path.getsize(path),
                "sha256": sha256_file(path),
            })
    rows.sort(key=lambda row: row["relative_path"])
    fd, temporary = tempfile.mkstemp(prefix=".output_manifest.", dir=output_root, text=True)
    os.close(fd)
    write_tsv(temporary, ["relative_path", "bytes", "sha256"], rows)
    publish_file_locked(temporary, target)


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=["run-all", "references", "sample", "combine"])
    parser.add_argument("--lock", required=True)
    parser.add_argument("--cellranger-root", required=True)
    parser.add_argument("--paper3-root", required=True)
    parser.add_argument("--sample-id")
    parser.add_argument("--prototype-only", action="store_true")
    return parser.parse_args()


def main():
    args = parse_args()
    with open(args.lock, "r") as handle:
        lock = json.load(handle, object_pairs_hook=OrderedDict)
    script_path = os.path.realpath(__file__)
    all_samples = ["15662-JZ-{}".format(i) for i in range(1, 13)]
    validate_reference_and_write(lock, args.cellranger_root, args.paper3_root)
    write_environment(args.lock, script_path, args.paper3_root)
    if args.command == "references":
        return
    if args.command == "sample":
        if not args.sample_id or args.sample_id not in all_samples:
            fail("--sample-id must be one of the 12 locked technical sample IDs")
        process_sample(lock, args.cellranger_root, args.paper3_root, args.sample_id, args.sample_id == lock["prototype"]["sample_id"])
        return
    if args.command == "combine":
        combine_outputs(lock, args.paper3_root, all_samples)
        write_output_manifest(args.paper3_root)
        return

    prototype_id = lock["prototype"]["sample_id"]
    process_sample(lock, args.cellranger_root, args.paper3_root, prototype_id, True)
    if args.prototype_only:
        return
    for sample_id in all_samples[1:]:
        process_sample(lock, args.cellranger_root, args.paper3_root, sample_id, False)
    combine_outputs(lock, args.paper3_root, all_samples)
    write_output_manifest(args.paper3_root)


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        sys.stderr.write("FATAL: {}\n".format(exc))
        raise
