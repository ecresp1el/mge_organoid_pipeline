#!/usr/bin/env python3
"""Locked raw three-probe EGFP audit for AGC request 15662-JZ.

Inputs are the Cell Ranger raw probe matrix, the vendor-filtered gene matrix,
and the vendor-filtered barcode CSV for each technical sample. Outputs remain
raw integer UMIs. This script does not normalize, call cell types, infer
genotype, or classify cells biologically.
"""

import argparse
import csv
import json
import os
import platform
import shutil
import sys
import tempfile
from collections import OrderedDict

import h5py
import numpy as np

import pcdh19_probe_audit as common


PATTERN_ORDER = [
    "E01+E02+E03", "E01+E02", "E01+E03", "E01 only",
    "E02+E03", "E02 only", "E03 only", "none",
]
PATTERN_HEADER = [
    "barcode", "probe01_UMI", "probe02_UMI", "probe03_UMI",
    "EGFP_total_UMI", "detection_pattern",
]
SUMMARY_HEADER = [
    "sample_id", "summary_type", "category", "probe_label", "probe_id",
    "total_umi", "cells", "detected_cells", "denominator_cells", "fraction",
    "percent", "percent_among_any_egfp", "mean_umi_among_detected",
    "max_umi_per_cell", "median_total_egfp_umi", "mean_total_egfp_umi",
    "max_total_egfp_umi",
]
VALIDATION_HEADER = ["sample_id", "check_name", "status", "observed", "expected", "details"]
REFERENCE_HEADER = [
    "probe_label", "probe_id", "probe_seq_5to3", "included", "10x_region",
    "alignment_reference", "reference_scope", "orientation_vs_egfp_mrna",
    "egfp_target_seq_5to3", "egfp_cds_start_1based", "egfp_cds_end_1based",
    "aligned_bases", "matches", "mismatches", "gaps", "percent_identity",
    "exact_match", "egfp_5prime_order",
]
DESIGN_HEADER = [
    "technical_sample_id", "submitted_sample_name", "genotype", "sex", "design_group",
    "filtered_cells", "egfp_total_umi", "any_egfp_cells", "fraction_any_egfp",
    "percent_any_egfp", "probe01_total_umi", "probe01_detected_cells",
    "probe02_total_umi", "probe02_detected_cells", "probe03_total_umi",
    "probe03_detected_cells",
]


def sample_paths(cellranger_root, sample_id):
    directory = os.path.join(cellranger_root, "per_sample_outs", sample_id)
    return {
        "probe_h5": os.path.join(directory, "sample_raw_probe_bc_matrix.h5"),
        "gene_h5": os.path.join(directory, "sample_filtered_feature_bc_matrix.h5"),
        "barcodes_csv": os.path.join(directory, "sample_filtered_barcodes.csv"),
    }


def locked_path(root, record):
    return os.path.join(root, record["relative_path"])


def validate_references(lock, cellranger_root, paper3_root, bundle_root):
    if lock.get("pipeline_name") != "paper3_egfp_probe_count_audit" or lock.get("pipeline_version") != "1.0.0":
        common.fail("Unexpected EGFP pipeline identity/version in lock")
    refs = lock["reference_inputs"]
    delivered = os.path.join(cellranger_root, refs["delivered_probe_set"]["basename"])
    common.require_sha256(delivered, refs["delivered_probe_set"]["sha256"], "delivered probe_set.csv")

    source_paths = {
        "xgfp_lock": locked_path(bundle_root, refs["xgfp_lock"]),
        "xgfp_alignment": locked_path(paper3_root, refs["xgfp_alignment"]),
        "xgfp_validation": locked_path(paper3_root, refs["xgfp_validation"]),
        "xgfp_output_manifest": locked_path(paper3_root, refs["xgfp_output_manifest"]),
        "sample_key": locked_path(bundle_root, refs["sample_key"]),
        "shared_helper": locked_path(bundle_root, refs["shared_helper"]),
    }
    for key, path in source_paths.items():
        common.require_sha256(path, refs[key]["sha256"], key)

    with open(source_paths["xgfp_validation"], "r", newline="") as handle:
        validation_rows = list(csv.DictReader(handle, delimiter="\t"))
    if not validation_rows or any(row["status"] != "PASS" for row in validation_rows):
        common.fail("Step 02b X-GFP compatibility validation is not fully PASS")

    panel_metadata, panel_rows = common.read_delivered_probe_set(delivered)
    panel = lock["panel"]
    expected_panel = {
        "probe_set_file_format": panel["probe_set_file_format"],
        "panel_name": panel["panel_name"],
        "reference_genome": panel["reference_genome"],
        "reference_version": panel["reference_version"],
    }
    for key, expected in expected_panel.items():
        if panel_metadata.get(key) != expected:
            common.fail("Delivered panel header {} mismatch: {} != {}".format(key, panel_metadata.get(key), expected))

    expected_ids = [probe["probe_id"] for probe in lock["probes"]]
    delivered_rows = [row for row in panel_rows if row["gene_id"] == panel["gene_id"]]
    delivered_by_id = {row["probe_id"]: row for row in delivered_rows}
    if len(delivered_rows) != 3 or set(delivered_by_id) != set(expected_ids):
        common.fail("Delivered panel does not contain exactly the three locked EGFP probe IDs")

    with open(source_paths["xgfp_alignment"], "r", newline="") as handle:
        alignment_rows = list(csv.DictReader(handle, delimiter="\t"))
    alignment_by_id = {row["probe_id"]: row for row in alignment_rows}
    if len(alignment_rows) != 3 or set(alignment_by_id) != set(expected_ids):
        common.fail("Step 02b alignment does not contain exactly the three locked EGFP probe IDs")

    reference_rows = []
    for probe in lock["probes"]:
        probe_id = probe["probe_id"]
        delivered_row = delivered_by_id[probe_id]
        alignment = alignment_by_id[probe_id]
        required = {
            "probe_seq": probe["probe_seq"], "included": "TRUE", "region": probe["region"],
        }
        for key, expected in required.items():
            if delivered_row[key] != expected:
                common.fail("Delivered {} mismatch for {}".format(key, probe_id))
        if alignment["probe_seq_5to3"] != probe["probe_seq"] or alignment["10x_region"] != probe["region"]:
            common.fail("Step 02b alignment identity mismatch for {}".format(probe_id))
        if alignment["exact_match"] != "TRUE" or alignment["matches"] != "50" or alignment["mismatches"] != "0" or alignment["gaps"] != "0":
            common.fail("Step 02b exact-alignment gate failed for {}".format(probe_id))
        if int(alignment["egfp_cds_start_1based"]) != probe["egfp_cds_start_1based"] or int(alignment["egfp_cds_end_1based"]) != probe["egfp_cds_end_1based"]:
            common.fail("Locked EGFP CDS interval mismatch for {}".format(probe_id))
        if int(alignment["egfp_5prime_order"]) != probe["egfp_5prime_order"]:
            common.fail("Locked EGFP 5-prime order mismatch for {}".format(probe_id))
        reference_rows.append({
            "probe_label": probe["label"], "probe_id": probe_id,
            "probe_seq_5to3": probe["probe_seq"], "included": "TRUE",
            "10x_region": probe["region"], "alignment_reference": alignment["alignment_reference"],
            "reference_scope": alignment["reference_scope"],
            "orientation_vs_egfp_mrna": alignment["orientation_vs_egfp_mrna"],
            "egfp_target_seq_5to3": alignment["egfp_target_seq_5to3"],
            "egfp_cds_start_1based": alignment["egfp_cds_start_1based"],
            "egfp_cds_end_1based": alignment["egfp_cds_end_1based"],
            "aligned_bases": alignment["aligned_bases"], "matches": alignment["matches"],
            "mismatches": alignment["mismatches"], "gaps": alignment["gaps"],
            "percent_identity": alignment["percent_identity"], "exact_match": alignment["exact_match"],
            "egfp_5prime_order": alignment["egfp_5prime_order"],
        })

    output_root = os.path.join(paper3_root, "results", "egfp_probe_audit")
    reference_root = os.path.join(output_root, "references")
    os.makedirs(reference_root, exist_ok=True)
    reference_target = os.path.join(reference_root, "egfp_probe_reference.tsv")
    fd, temporary = tempfile.mkstemp(prefix=".egfp_reference.", dir=reference_root, text=True)
    os.close(fd)
    common.write_tsv(temporary, REFERENCE_HEADER, reference_rows)
    reference_sha = common.sha256_file(temporary)
    common.publish_file_locked(temporary, reference_target)

    manifest_rows = [
        {"input_id": "delivered_probe_set", "role": "reference_input", "version": panel["panel_name"] + ";" + panel["reference_version"], "resolved_path": delivered, "sha256": refs["delivered_probe_set"]["sha256"], "validation_status": "PASS"},
        {"input_id": "xgfp_lock", "role": "upstream_lock", "version": "02b-1.0.0", "resolved_path": source_paths["xgfp_lock"], "sha256": refs["xgfp_lock"]["sha256"], "validation_status": "PASS"},
        {"input_id": "xgfp_alignment", "role": "upstream_validated_reference", "version": "U55762.1-EGFP-CDS", "resolved_path": source_paths["xgfp_alignment"], "sha256": refs["xgfp_alignment"]["sha256"], "validation_status": "PASS"},
        {"input_id": "xgfp_validation", "role": "upstream_validation", "version": "02b-1.0.0", "resolved_path": source_paths["xgfp_validation"], "sha256": refs["xgfp_validation"]["sha256"], "validation_status": "PASS"},
        {"input_id": "xgfp_output_manifest", "role": "upstream_manifest", "version": "02b-1.0.0", "resolved_path": source_paths["xgfp_output_manifest"], "sha256": refs["xgfp_output_manifest"]["sha256"], "validation_status": "PASS"},
        {"input_id": "sample_key", "role": "registered_sample_metadata_for_descriptive_join_only", "version": "2026-08-27", "resolved_path": source_paths["sample_key"], "sha256": refs["sample_key"]["sha256"], "validation_status": "PASS"},
        {"input_id": "shared_helper", "role": "locked_code_dependency", "version": "pcdh19-audit-1.0.0", "resolved_path": source_paths["shared_helper"], "sha256": refs["shared_helper"]["sha256"], "validation_status": "PASS"},
        {"input_id": "egfp_probe_reference", "role": "derived_reference", "version": lock["pipeline_version"], "resolved_path": reference_target, "sha256": reference_sha, "validation_status": "PASS"},
    ]
    for row in manifest_rows:
        row["bytes"] = os.path.getsize(row["resolved_path"])
    manifest_header = ["input_id", "role", "version", "resolved_path", "sha256", "bytes", "validation_status"]
    manifest_target = os.path.join(reference_root, "reference_manifest.tsv")
    fd, temporary = tempfile.mkstemp(prefix=".reference_manifest.", dir=reference_root, text=True)
    os.close(fd)
    common.write_tsv(temporary, manifest_header, manifest_rows)
    common.publish_file_locked(temporary, manifest_target)
    return source_paths


def verify_existing_sample(directory, prototype=None):
    required = ["egfp_probe_patterns.tsv", "egfp_probe_summary.tsv", "validation.tsv", "checksums.sha256"]
    if not os.path.isdir(directory):
        return False
    for name in required:
        if not os.path.isfile(os.path.join(directory, name)):
            common.fail("Existing EGFP sample output is incomplete and will not be overwritten: {}".format(directory))
    expected = {}
    with open(os.path.join(directory, "checksums.sha256"), "r") as handle:
        for line in handle:
            checksum, name = line.rstrip("\n").split("  ", 1)
            expected[name] = checksum
    for name in required[:-1]:
        observed = common.sha256_file(os.path.join(directory, name))
        if expected.get(name) != observed:
            common.fail("Existing EGFP output checksum failure: {}/{}".format(directory, name))
    with open(os.path.join(directory, "validation.tsv"), "r", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows or any(row["status"] != "PASS" for row in rows):
        common.fail("Existing EGFP validation is not fully PASS: {}".format(directory))
    if prototype and common.sha256_file(os.path.join(directory, "egfp_probe_patterns.tsv")) != prototype["pattern_table_sha256"]:
        common.fail("Existing JZ-1 EGFP barcode table is not byte-equivalent to the frozen prototype")
    return True


def process_sample(lock, cellranger_root, paper3_root, sample_id, prototype_gate=False):
    output_root = os.path.join(paper3_root, "results", "egfp_probe_audit")
    per_sample_root = os.path.join(output_root, "per_sample")
    os.makedirs(per_sample_root, exist_ok=True)
    destination = os.path.join(per_sample_root, sample_id)
    prototype = lock["prototype"] if prototype_gate else None
    if verify_existing_sample(destination, prototype):
        print("{}: existing validated EGFP output verified; no overwrite".format(sample_id), flush=True)
        return

    paths = sample_paths(cellranger_root, sample_id)
    for key, path in paths.items():
        if not os.path.isfile(path):
            common.fail("{} input missing for {}: {}".format(key, sample_id, path))
    staging = tempfile.mkdtemp(prefix=".{}.tmp.".format(sample_id), dir=per_sample_root)
    try:
        panel = lock["panel"]
        encoded_probe_ids = [probe["probe_id"].encode() for probe in lock["probes"]]
        with h5py.File(paths["probe_h5"], "r") as handle:
            matrix = handle["matrix"]
            feature_ids = matrix["features/id"][:]
            matching = np.flatnonzero(np.isin(feature_ids, np.array(encoded_probe_ids, dtype=feature_ids.dtype)))
            if matching.size != 3:
                common.fail("{} raw probe matrix does not contain exactly the three locked EGFP features".format(sample_id))
            lookup = {value: i for i, value in enumerate(feature_ids)}
            rows = np.array([lookup[probe_id] for probe_id in encoded_probe_ids], dtype=np.int64)
            if list(matrix["features/target_sets"].keys()) != [panel["panel_name"]]:
                common.fail("{} raw probe target set is not the locked panel".format(sample_id))
            for row, probe in zip(rows, lock["probes"]):
                checks = {
                    "gene_id": panel["gene_id"], "gene_name": panel["gene_name"],
                    "probe_region": probe["region"], "genome": panel["reference_genome"],
                }
                for field, expected in checks.items():
                    if matrix["features/" + field][row].decode() != expected:
                        common.fail("{} raw probe {} mismatch for {}".format(sample_id, field, probe["probe_id"]))
            filtered_columns = np.flatnonzero(matrix["filtered_barcodes"][:].astype(bool))
            probe_barcodes = matrix["barcodes"][:][filtered_columns]
            probe_counts = common.extract_sparse_rows(matrix, rows)[:, filtered_columns]

        genomes = []
        csv_barcode_list = []
        with open(paths["barcodes_csv"], "r", newline="") as handle:
            for row in csv.reader(handle):
                if len(row) != 2:
                    common.fail("Malformed filtered barcode CSV row for {}".format(sample_id))
                genomes.append(row[0])
                csv_barcode_list.append(row[1].encode())
        csv_barcodes = np.array(csv_barcode_list, dtype=probe_barcodes.dtype)

        with h5py.File(paths["gene_h5"], "r") as handle:
            matrix = handle["matrix"]
            gene_barcodes = matrix["barcodes"][:]
            hits = np.flatnonzero(matrix["features/id"][:] == panel["gene_id"].encode())
            if hits.size != 1:
                common.fail("{} filtered gene matrix does not contain exactly one EGFP feature".format(sample_id))
            hit = int(hits[0])
            if matrix["features/name"][hit].decode() != panel["gene_name"]:
                common.fail("{} filtered gene matrix EGFP name mismatch".format(sample_id))
            if matrix["features/genome"][hit].decode() != panel["reference_genome"]:
                common.fail("{} filtered gene matrix EGFP genome mismatch".format(sample_id))
            chemistry = handle.attrs.get("chemistry_description", b"")
            chemistry = chemistry.decode() if hasattr(chemistry, "decode") else chemistry
            if chemistry != panel["chemistry_description"]:
                common.fail("{} chemistry mismatch: {}".format(sample_id, chemistry))
            if list(matrix["features/target_sets"].keys()) != [panel["panel_name"]]:
                common.fail("{} filtered gene target set is not the locked panel".format(sample_id))
            gene_counts = common.extract_sparse_rows(matrix, np.array([hit], dtype=np.int64))[0]

        total_cells = len(csv_barcodes)
        E01, E02, E03 = probe_counts
        reconstructed = E01 + E02 + E03
        e01, e02, e03 = E01 > 0, E02 > 0, E03 > 0
        patterns = np.empty(total_cells, dtype="U13")
        masks = OrderedDict([
            ("E01+E02+E03", e01 & e02 & e03),
            ("E01+E02", e01 & e02 & ~e03),
            ("E01+E03", e01 & ~e02 & e03),
            ("E01 only", e01 & ~e02 & ~e03),
            ("E02+E03", ~e01 & e02 & e03),
            ("E02 only", ~e01 & e02 & ~e03),
            ("E03 only", ~e01 & ~e02 & e03),
            ("none", ~e01 & ~e02 & ~e03),
        ])
        for pattern, mask in masks.items():
            patterns[mask] = pattern
        pattern_counts = OrderedDict((pattern, int(mask.sum())) for pattern, mask in masks.items())
        any_detected = int(np.count_nonzero(reconstructed))

        validations = [
            common.check_row(sample_id, "filtered_barcode_count_probe_vs_csv", len(probe_barcodes), total_cells),
            common.check_row(sample_id, "filtered_barcode_count_gene_vs_csv", len(gene_barcodes), total_cells),
            common.check_row(sample_id, "filtered_barcodes_probe_vs_csv", bool(np.array_equal(probe_barcodes, csv_barcodes)), True),
            common.check_row(sample_id, "filtered_barcodes_gene_vs_csv", bool(np.array_equal(gene_barcodes, csv_barcodes)), True),
            common.check_row(sample_id, "filtered_barcode_genome", sorted(set(genomes)), [panel["reference_genome"]]),
            common.check_row(sample_id, "duplicate_filtered_barcodes", total_cells - len(set(csv_barcode_list)), 0),
            common.check_row(sample_id, "per_barcode_probe_sum_matches_gene", int(np.count_nonzero(reconstructed != gene_counts)), 0),
            common.check_row(sample_id, "reconstructed_total_matches_gene_total", int(reconstructed.sum()), int(gene_counts.sum())),
            common.check_row(sample_id, "pattern_partition_matches_filtered_cells", sum(pattern_counts.values()), total_cells),
        ]
        for label, detected in zip(["E01", "E02", "E03"], [e01, e02, e03]):
            marginal = sum(count for pattern, count in pattern_counts.items() if label in pattern)
            validations.append(common.check_row(sample_id, "{}_pattern_marginal_detected".format(label), marginal, int(np.count_nonzero(detected))))

        pattern_path = os.path.join(staging, "egfp_probe_patterns.tsv")
        with open(pattern_path, "w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(PATTERN_HEADER)
            for barcode, v1, v2, v3, total, pattern in zip(csv_barcodes, E01, E02, E03, reconstructed, patterns):
                writer.writerow([barcode.decode(), int(v1), int(v2), int(v3), int(total), pattern])

        serialized_rows = 0
        serialized_bad_sums = 0
        serialized_seen = set()
        serialized_patterns = OrderedDict((pattern, 0) for pattern in PATTERN_ORDER)
        with open(pattern_path, "r", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if reader.fieldnames != PATTERN_HEADER:
                common.fail("Serialized EGFP barcode header differs from the lock")
            for row in reader:
                serialized_rows += 1
                serialized_bad_sums += int(int(row["probe01_UMI"]) + int(row["probe02_UMI"]) + int(row["probe03_UMI"]) != int(row["EGFP_total_UMI"]))
                serialized_seen.add(row["barcode"])
                if row["detection_pattern"] not in serialized_patterns:
                    common.fail("Unexpected serialized EGFP pattern: {}".format(row["detection_pattern"]))
                serialized_patterns[row["detection_pattern"]] += 1
        validations.extend([
            common.check_row(sample_id, "serialized_barcode_rows", serialized_rows, total_cells),
            common.check_row(sample_id, "serialized_unique_barcodes", len(serialized_seen), total_cells),
            common.check_row(sample_id, "serialized_per_barcode_probe_sum", serialized_bad_sums, 0),
            common.check_row(sample_id, "serialized_pattern_counts", serialized_patterns, pattern_counts),
        ])

        summary_rows = []
        for probe, counts, detected in zip(lock["probes"], [E01, E02, E03], [e01, e02, e03]):
            detected_cells = int(np.count_nonzero(detected))
            summary_rows.append({
                "sample_id": sample_id, "summary_type": "probe", "category": probe["label"],
                "probe_label": probe["label"], "probe_id": probe["probe_id"],
                "total_umi": int(counts.sum()), "detected_cells": detected_cells,
                "denominator_cells": total_cells,
                "fraction": common.format_float(detected_cells / float(total_cells)),
                "percent": common.format_float(100.0 * detected_cells / float(total_cells)),
                "mean_umi_among_detected": common.format_float(np.mean(counts[detected])) if detected_cells else "",
                "max_umi_per_cell": int(counts.max()) if total_cells else 0,
            })
        summary_rows.append({
            "sample_id": sample_id, "summary_type": "gene_reconstruction", "category": "EGFP",
            "total_umi": int(reconstructed.sum()), "detected_cells": any_detected,
            "denominator_cells": total_cells, "fraction": common.format_float(any_detected / float(total_cells)),
            "percent": common.format_float(100.0 * any_detected / float(total_cells)),
            "mean_umi_among_detected": common.format_float(np.mean(reconstructed[reconstructed > 0])) if any_detected else "",
            "max_umi_per_cell": int(reconstructed.max()) if total_cells else 0,
        })
        for pattern, mask in masks.items():
            values = reconstructed[mask]
            cells = int(mask.sum())
            summary_rows.append({
                "sample_id": sample_id, "summary_type": "pattern", "category": pattern,
                "cells": cells, "denominator_cells": total_cells,
                "fraction": common.format_float(cells / float(total_cells)),
                "percent": common.format_float(100.0 * cells / float(total_cells)),
                "percent_among_any_egfp": "" if pattern == "none" or not any_detected else common.format_float(100.0 * cells / float(any_detected)),
                "median_total_egfp_umi": common.format_float(np.median(values)) if cells else "",
                "mean_total_egfp_umi": common.format_float(np.mean(values)) if cells else "",
                "max_total_egfp_umi": int(values.max()) if cells else "",
            })
        summary_path = os.path.join(staging, "egfp_probe_summary.tsv")
        common.write_tsv(summary_path, SUMMARY_HEADER, summary_rows)

        pattern_sha = common.sha256_file(pattern_path)
        validations.append(common.check_row(sample_id, "pattern_table_sha256_recorded", pattern_sha, pattern_sha))
        if prototype_gate:
            expected = lock["prototype"]
            validations.extend([
                common.check_row(sample_id, "prototype_filtered_cells", total_cells, expected["filtered_cells"]),
                common.check_row(sample_id, "prototype_egfp_total_umi", int(reconstructed.sum()), expected["egfp_total_umi"]),
                common.check_row(sample_id, "prototype_pattern_table_byte_equivalence", pattern_sha, expected["pattern_table_sha256"]),
            ])
            for label, counts, detected in zip(["E01", "E02", "E03"], [E01, E02, E03], [e01, e02, e03]):
                validations.append(common.check_row(sample_id, "prototype_{}_total_umi".format(label), int(counts.sum()), expected["probe_total_umi"][label]))
                validations.append(common.check_row(sample_id, "prototype_{}_detected_cells".format(label), int(np.count_nonzero(detected)), expected["probe_detected_cells"][label]))
            for pattern in PATTERN_ORDER:
                validations.append(common.check_row(sample_id, "prototype_pattern_{}".format(pattern.replace(" ", "_").replace("+", "plus")), pattern_counts[pattern], expected["pattern_counts"][pattern]))

        validation_path = os.path.join(staging, "validation.tsv")
        common.write_tsv(validation_path, VALIDATION_HEADER, validations)
        failed = [row for row in validations if row["status"] != "PASS"]
        if failed:
            common.fail("{} EGFP validation failed: {}".format(sample_id, ", ".join(row["check_name"] for row in failed)))

        with open(os.path.join(staging, "checksums.sha256"), "w") as handle:
            for name in ["egfp_probe_patterns.tsv", "egfp_probe_summary.tsv", "validation.tsv"]:
                handle.write("{}  {}\n".format(common.sha256_file(os.path.join(staging, name)), name))
        if os.path.exists(destination):
            common.fail("EGFP destination appeared during processing and will not be overwritten: {}".format(destination))
        os.replace(staging, destination)
        staging = None
        verify_existing_sample(destination, prototype)
        print("{}: published {} cells, {} EGFP UMIs".format(sample_id, total_cells, int(reconstructed.sum())), flush=True)
    finally:
        if staging and os.path.isdir(staging):
            shutil.rmtree(staging)


def read_sample_key(path, sample_ids):
    with open(path, "r", newline="") as handle:
        rows = list(csv.DictReader(handle))
    by_id = {row["technical_sample_id"]: row for row in rows}
    if len(rows) != len(sample_ids) or set(by_id) != set(sample_ids):
        common.fail("Registered sample key is not an exact one-to-one mapping of the 12 technical samples")
    return by_id


def combine_outputs(lock, paper3_root, bundle_root, sample_ids):
    output_root = os.path.join(paper3_root, "results", "egfp_probe_audit")
    combined = os.path.join(output_root, "combined")
    os.makedirs(combined, exist_ok=True)
    probe_rows, pattern_rows, validation_rows, design_rows = [], [], [], []
    sample_key_path = locked_path(bundle_root, lock["reference_inputs"]["sample_key"])
    sample_key = read_sample_key(sample_key_path, sample_ids)
    for sample_id in sample_ids:
        directory = os.path.join(output_root, "per_sample", sample_id)
        verify_existing_sample(directory, lock["prototype"] if sample_id == lock["prototype"]["sample_id"] else None)
        by_type = {}
        with open(os.path.join(directory, "egfp_probe_summary.tsv"), "r", newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                if row["summary_type"] in ("probe", "gene_reconstruction"):
                    probe_rows.append(row)
                else:
                    pattern_rows.append(row)
                by_type[(row["summary_type"], row["category"])] = row
        with open(os.path.join(directory, "validation.tsv"), "r", newline="") as handle:
            validation_rows.extend(csv.DictReader(handle, delimiter="\t"))
        key = sample_key[sample_id]
        gene = by_type[("gene_reconstruction", "EGFP")]
        p1, p2, p3 = [by_type[("probe", label)] for label in ["E01", "E02", "E03"]]
        design_rows.append({
            "technical_sample_id": sample_id, "submitted_sample_name": key["submitted_sample_name"],
            "genotype": key["genotype"], "sex": key["sex"], "design_group": key["design_group"],
            "filtered_cells": gene["denominator_cells"], "egfp_total_umi": gene["total_umi"],
            "any_egfp_cells": gene["detected_cells"], "fraction_any_egfp": gene["fraction"],
            "percent_any_egfp": gene["percent"], "probe01_total_umi": p1["total_umi"],
            "probe01_detected_cells": p1["detected_cells"], "probe02_total_umi": p2["total_umi"],
            "probe02_detected_cells": p2["detected_cells"], "probe03_total_umi": p3["total_umi"],
            "probe03_detected_cells": p3["detected_cells"],
        })
    outputs = [
        ("all_samples_egfp_probe_summary.tsv", SUMMARY_HEADER, probe_rows),
        ("all_samples_egfp_pattern_summary.tsv", SUMMARY_HEADER, pattern_rows),
        ("all_samples_validation.tsv", VALIDATION_HEADER, validation_rows),
        ("all_samples_egfp_design_summary.tsv", DESIGN_HEADER, design_rows),
    ]
    for name, header, rows in outputs:
        fd, temporary = tempfile.mkstemp(prefix=".{}.".format(name), dir=combined, text=True)
        os.close(fd)
        common.write_tsv(temporary, header, rows)
        common.publish_file_locked(temporary, os.path.join(combined, name))


def write_environment(lock_path, script_path, paper3_root):
    output_root = os.path.join(paper3_root, "results", "egfp_probe_audit")
    os.makedirs(output_root, exist_ok=True)
    rows = [
        {"key": "pipeline_version", "value": "1.0.0"},
        {"key": "python", "value": sys.version.replace("\n", " ")},
        {"key": "implementation", "value": platform.python_implementation()},
        {"key": "h5py", "value": h5py.__version__},
        {"key": "numpy", "value": np.__version__},
        {"key": "pipeline_script_sha256", "value": common.sha256_file(script_path)},
        {"key": "reference_lock_sha256", "value": common.sha256_file(lock_path)},
        {"key": "count_scale", "value": "raw_integer_UMI"},
        {"key": "barcode_denominator", "value": "Cell_Ranger_vendor_filtered_barcodes"},
        {"key": "biological_classification", "value": "none"},
    ]
    target = os.path.join(output_root, "software_environment.tsv")
    fd, temporary = tempfile.mkstemp(prefix=".environment.", dir=output_root, text=True)
    os.close(fd)
    common.write_tsv(temporary, ["key", "value"], rows)
    common.publish_file_locked(temporary, target)


def write_output_manifest(paper3_root):
    output_root = os.path.join(paper3_root, "results", "egfp_probe_audit")
    target = os.path.join(output_root, "output_manifest.tsv")
    rows = []
    for root, directories, files in os.walk(output_root):
        directories[:] = sorted(name for name in directories if not name.startswith("."))
        for name in sorted(files):
            path = os.path.join(root, name)
            if path == target or name.startswith("."):
                continue
            rows.append({
                "relative_path": os.path.relpath(path, output_root),
                "bytes": os.path.getsize(path), "sha256": common.sha256_file(path),
            })
    rows.sort(key=lambda row: row["relative_path"])
    fd, temporary = tempfile.mkstemp(prefix=".output_manifest.", dir=output_root, text=True)
    os.close(fd)
    common.write_tsv(temporary, ["relative_path", "bytes", "sha256"], rows)
    common.publish_file_locked(temporary, target)


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=["run-all", "references", "sample", "combine"])
    parser.add_argument("--lock", required=True)
    parser.add_argument("--cellranger-root", required=True)
    parser.add_argument("--paper3-root", required=True)
    parser.add_argument("--bundle-root", required=True)
    parser.add_argument("--sample-id")
    parser.add_argument("--prototype-only", action="store_true")
    return parser.parse_args()


def main():
    args = parse_args()
    with open(args.lock, "r") as handle:
        lock = json.load(handle, object_pairs_hook=OrderedDict)
    sample_ids = ["15662-JZ-{}".format(i) for i in range(1, 13)]
    validate_references(lock, args.cellranger_root, args.paper3_root, args.bundle_root)
    write_environment(args.lock, os.path.realpath(__file__), args.paper3_root)
    if args.command == "references":
        return
    if args.command == "sample":
        if not args.sample_id or args.sample_id not in sample_ids:
            common.fail("--sample-id must be one of the 12 locked technical sample IDs")
        process_sample(lock, args.cellranger_root, args.paper3_root, args.sample_id, args.sample_id == lock["prototype"]["sample_id"])
        return
    if args.command == "combine":
        combine_outputs(lock, args.paper3_root, args.bundle_root, sample_ids)
        write_output_manifest(args.paper3_root)
        return

    prototype_id = lock["prototype"]["sample_id"]
    process_sample(lock, args.cellranger_root, args.paper3_root, prototype_id, True)
    if args.prototype_only:
        return
    for sample_id in sample_ids[1:]:
        process_sample(lock, args.cellranger_root, args.paper3_root, sample_id, False)
    combine_outputs(lock, args.paper3_root, args.bundle_root, sample_ids)
    write_output_manifest(args.paper3_root)


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        sys.stderr.write("FATAL: {}\n".format(exc))
        raise
