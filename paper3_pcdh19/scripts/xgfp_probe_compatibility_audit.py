#!/usr/bin/env python3
"""Reproduce the locked construct-level X-GFP/Flex sequence audit for 15662-JZ."""

import argparse
import csv
import hashlib
import json
import os
import platform
import shutil
import sys
import tempfile
import urllib.request
from collections import OrderedDict


def fail(message):
    raise RuntimeError(message)


def sha256_bytes(value):
    return hashlib.sha256(value).hexdigest()


def sha256_file(path, chunk_size=8 * 1024 * 1024):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        while True:
            chunk = handle.read(chunk_size)
            if not chunk:
                return digest.hexdigest()
            digest.update(chunk)


def reverse_complement(sequence):
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def write_tsv(path, header, rows):
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=header, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in header})


def parse_panel(path):
    metadata = OrderedDict()
    rows = []
    data_lines = []
    with open(path, "r", newline="") as handle:
        for line in handle:
            if line.startswith("#"):
                key, value = line[1:].rstrip("\r\n").split("=", 1)
                metadata[key] = value
            else:
                data_lines.append(line)
    rows.extend(csv.DictReader(data_lines))
    return metadata, rows


def parse_fasta_bytes(raw):
    lines = raw.decode("ascii").splitlines()
    if not lines or not lines[0].startswith(">"):
        fail("Downloaded EGFP reference is not FASTA")
    sequence = "".join(line.strip() for line in lines[1:]).upper()
    if not sequence or set(sequence) - set("ACGTN"):
        fail("Downloaded EGFP reference contains invalid sequence")
    return lines[0][1:], sequence


def load_locked_reference(lock, cache_root):
    reference = lock["egfp_sequence_reference"]
    os.makedirs(cache_root, exist_ok=True)
    path = os.path.join(cache_root, reference["accession"] + ".fasta")
    if os.path.isfile(path):
        with open(path, "rb") as handle:
            raw = handle.read()
    else:
        request = urllib.request.Request(
            reference["efetch_url"],
            headers={"User-Agent": "paper3-xgfp-probe-audit/1.0"},
        )
        with urllib.request.urlopen(request) as response:
            raw = response.read()
        fd, temporary = tempfile.mkstemp(prefix=".U55762.1.", dir=cache_root)
        try:
            with os.fdopen(fd, "wb") as handle:
                handle.write(raw)
                handle.flush()
                os.fsync(handle.fileno())
            if os.path.exists(path):
                fail("Reference appeared while downloading: {}".format(path))
            os.replace(temporary, path)
        finally:
            if os.path.exists(temporary):
                os.unlink(temporary)
    header, sequence = parse_fasta_bytes(raw)
    if len(sequence) != reference["record_length"]:
        fail("U55762.1 length mismatch")
    if sha256_bytes(sequence.encode("ascii")) != reference["record_sequence_sha256"]:
        fail("U55762.1 sequence checksum mismatch")
    start = reference["egfp_cds_start_1based"] - 1
    end = reference["egfp_cds_end_1based"]
    cds = sequence[start:end]
    if len(cds) != reference["egfp_cds_length"]:
        fail("EGFP CDS length mismatch")
    if sha256_bytes(cds.encode("ascii")) != reference["egfp_cds_sha256"]:
        fail("EGFP CDS checksum mismatch")
    return path, header, sequence, cds


def verify_existing(output_root):
    manifest = os.path.join(output_root, "output_manifest.tsv")
    if not os.path.isdir(output_root):
        return False
    if not os.path.isfile(manifest):
        fail("Existing X-GFP audit output is incomplete: {}".format(output_root))
    with open(manifest, "r", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    for row in rows:
        path = os.path.join(output_root, row["relative_path"])
        if not os.path.isfile(path) or sha256_file(path) != row["sha256"]:
            fail("Existing X-GFP audit checksum failure: {}".format(path))
    validation = os.path.join(output_root, "validation.tsv")
    with open(validation, "r", newline="") as handle:
        checks = list(csv.DictReader(handle, delimiter="\t"))
    if not checks or any(row["status"] != "PASS" for row in checks):
        fail("Existing X-GFP audit is not fully PASS")
    return True


def add_check(rows, name, observed, expected, details=""):
    rows.append({
        "check_name": name,
        "status": "PASS" if observed == expected else "FAIL",
        "observed": observed,
        "expected": expected,
        "details": details,
    })


def run(args):
    with open(args.lock, "r") as handle:
        lock = json.load(handle, object_pairs_hook=OrderedDict)
    output_root = os.path.join(args.paper3_root, "results", "xgfp_probe_audit")
    if verify_existing(output_root):
        print("Existing validated X-GFP probe audit verified; no overwrite", flush=True)
        return

    panel_lock = lock["delivered_panel"]
    if os.path.basename(args.probe_set) != panel_lock["basename"]:
        fail("Delivered panel basename mismatch")
    if not os.path.isfile(args.probe_set):
        fail("Delivered probe_set.csv is missing: {}".format(args.probe_set))
    if sha256_file(args.probe_set) != panel_lock["sha256"]:
        fail("Delivered probe_set.csv checksum mismatch")

    panel_header, panel_rows = parse_panel(args.probe_set)
    for key in ["probe_set_file_format", "panel_name", "reference_genome", "reference_version"]:
        if panel_header.get(key) != panel_lock[key]:
            fail("Panel header mismatch for {}".format(key))
    egfp_rows = [row for row in panel_rows if row["gene_id"] == "EGFP" or row["gene_name"] == "EGFP"]
    expected = OrderedDict((row["probe_id"], row) for row in lock["expected_probes"])
    observed = OrderedDict((row["probe_id"], row) for row in egfp_rows)
    if len(observed) != len(egfp_rows):
        fail("Duplicate EGFP probe IDs in delivered panel")
    if set(observed) != set(expected):
        fail("Delivered EGFP probe set differs from locked three-probe set")
    for probe_id, expected_row in expected.items():
        row = observed[probe_id]
        for key in ["probe_seq", "included", "region"]:
            if row[key] != expected_row[key]:
                fail("{} mismatch for {}".format(key, probe_id))

    cache_root = os.path.join(args.paper3_root, "inputs", "xgfp_probe_audit", "references")
    fasta_path, fasta_header, record_sequence, cds = load_locked_reference(lock, cache_root)
    reference = lock["egfp_sequence_reference"]
    validations = []
    add_check(validations, "delivered_panel_sha256", sha256_file(args.probe_set), panel_lock["sha256"])
    add_check(validations, "exact_egfp_probe_id_set", sorted(observed), sorted(expected))
    add_check(validations, "egfp_probe_count", len(observed), 3)
    add_check(validations, "reference_record_length", len(record_sequence), reference["record_length"])
    add_check(validations, "reference_record_sequence_sha256", sha256_bytes(record_sequence.encode("ascii")), reference["record_sequence_sha256"])
    add_check(validations, "egfp_cds_length", len(cds), reference["egfp_cds_length"])
    add_check(validations, "egfp_cds_sha256", sha256_bytes(cds.encode("ascii")), reference["egfp_cds_sha256"])

    alignments = []
    for probe_id, expected_row in expected.items():
        probe = expected_row["probe_seq"]
        reverse_target = reverse_complement(probe)
        direct_positions = [i for i in range(len(cds) - len(probe) + 1) if cds[i:i + len(probe)] == probe]
        reverse_positions = [i for i in range(len(cds) - len(probe) + 1) if cds[i:i + len(probe)] == reverse_target]
        add_check(validations, probe_id + "_unique_reverse_complement_match", len(reverse_positions), 1)
        add_check(validations, probe_id + "_direct_match_count", len(direct_positions), 0)
        if len(reverse_positions) != 1:
            fail("{} does not have exactly one reverse-complement EGFP CDS match".format(probe_id))
        start0 = reverse_positions[0]
        cds_start = start0 + 1
        cds_end = start0 + len(probe)
        record_start = reference["egfp_cds_start_1based"] + start0
        record_end = record_start + len(probe) - 1
        alignments.append({
            "probe_id": probe_id,
            "probe_seq_5to3": probe,
            "included": expected_row["included"],
            "10x_region": expected_row["region"],
            "alignment_reference": reference["accession"],
            "reference_scope": "Clontech_EGFP_CDS",
            "orientation_vs_egfp_mrna": "reverse_complement",
            "egfp_target_seq_5to3": reverse_target,
            "egfp_cds_start_1based": cds_start,
            "egfp_cds_end_1based": cds_end,
            "u55762_start_1based": record_start,
            "u55762_end_1based": record_end,
            "aligned_bases": len(probe),
            "matches": len(probe),
            "mismatches": 0,
            "gaps": 0,
            "percent_identity": "100.000",
            "exact_match": "TRUE",
        })
    ordered = sorted(alignments, key=lambda row: row["egfp_cds_start_1based"])
    for rank, row in enumerate(ordered, 1):
        row["egfp_5prime_order"] = rank
    add_check(validations, "all_three_probes_exact_50bp_matches", sum(row["exact_match"] == "TRUE" and row["aligned_bases"] == 50 for row in alignments), 3)
    add_check(validations, "probe_order_on_egfp_mrna", [row["probe_id"] for row in ordered], ["EGFP|EGFP|probe03", "EGFP|EGFP|probe02", "EGFP|EGFP|probe01"])
    if any(row["status"] != "PASS" for row in validations):
        fail("One or more X-GFP audit validations failed")

    parent = os.path.dirname(output_root)
    os.makedirs(parent, exist_ok=True)
    staging = tempfile.mkdtemp(prefix=".xgfp_probe_audit.", dir=parent)
    try:
        references_dir = os.path.join(staging, "references")
        os.makedirs(references_dir)
        alignment_header = [
            "probe_id", "probe_seq_5to3", "included", "10x_region",
            "alignment_reference", "reference_scope", "orientation_vs_egfp_mrna",
            "egfp_target_seq_5to3", "egfp_cds_start_1based", "egfp_cds_end_1based",
            "u55762_start_1based", "u55762_end_1based", "aligned_bases", "matches",
            "mismatches", "gaps", "percent_identity", "exact_match", "egfp_5prime_order",
        ]
        write_tsv(os.path.join(staging, "gfp_probe_alignment.tsv"), alignment_header, ordered)
        write_tsv(
            os.path.join(staging, "validation.tsv"),
            ["check_name", "status", "observed", "expected", "details"],
            validations,
        )
        with open(os.path.join(references_dir, "U55762.1_EGFP_CDS.fasta"), "w") as handle:
            handle.write(">U55762.1:679-1398 Clontech EGFP CDS; construct-level X-GFP alignment reference\n")
            for i in range(0, len(cds), 60):
                handle.write(cds[i:i + 60] + "\n")

        source_rows = []
        source_rows.append({
            "source_id": "delivered_probe_set", "role": "exact Flex probe IDs and sequences",
            "citation": "Cell Ranger delivered probe_set.csv", "url": "",
            "resolved_path": args.probe_set, "sha256": sha256_file(args.probe_set),
            "status": "PASS",
        })
        source_rows.append({
            "source_id": "clontech_egfp_fasta", "role": "downloaded full sequence reference",
            "citation": "NCBI GenBank U55762.1", "url": reference["efetch_url"],
            "resolved_path": fasta_path, "sha256": sha256_file(fasta_path), "status": "PASS",
        })
        for row in lock["evidence_sources"]:
            source_rows.append({
                "source_id": row["source_id"], "role": row["role"],
                "citation": row["citation"], "url": row["url"],
                "resolved_path": "", "sha256": "", "status": "EVIDENCE_REGISTERED",
            })
        write_tsv(
            os.path.join(references_dir, "source_manifest.tsv"),
            ["source_id", "role", "citation", "url", "resolved_path", "sha256", "status"],
            source_rows,
        )
        conclusion_rows = [
            {"field": "reporter_line", "value": lock["reporter_line"]["symbol"]},
            {"field": "reporter_construct", "value": lock["reporter_line"]["construct"]},
            {"field": "alignment_sequence_scope", "value": "Clontech EGFP coding sequence specified for the original pCX-EGFP construct"},
            {"field": "exact_probe_matches", "value": "3_of_3"},
            {"field": "sequence_compatibility", "value": "PASS"},
            {"field": "theoretical_flex_capability", "value": "SUPPORTED"},
            {"field": "integrated_allele_full_sequence", "value": "NOT_DEPOSITED_IN_CITED_LINE_RECORDS"},
            {"field": "interpretation", "value": "The delivered Flex panel was sequence-compatible with the original X-GFP reporter coding sequence; this does not prove reporter transcription, RNA preservation, hybridization efficiency, or observed counts in any sample or cell."},
        ]
        write_tsv(os.path.join(staging, "conclusion.tsv"), ["field", "value"], conclusion_rows)
        environment = [
            {"key": "pipeline_name", "value": lock["pipeline_name"]},
            {"key": "pipeline_version", "value": lock["pipeline_version"]},
            {"key": "python", "value": sys.version.replace("\n", " ")},
            {"key": "platform", "value": platform.platform()},
            {"key": "pipeline_script_sha256", "value": sha256_file(os.path.realpath(__file__))},
            {"key": "reference_lock_sha256", "value": sha256_file(args.lock)},
            {"key": "fasta_header", "value": fasta_header},
        ]
        write_tsv(os.path.join(staging, "software_environment.tsv"), ["key", "value"], environment)

        manifest_rows = []
        for root, directories, files in os.walk(staging):
            directories.sort()
            for name in sorted(files):
                path = os.path.join(root, name)
                relative = os.path.relpath(path, staging)
                if relative == "output_manifest.tsv":
                    continue
                manifest_rows.append({
                    "relative_path": relative,
                    "bytes": os.path.getsize(path),
                    "sha256": sha256_file(path),
                })
        manifest_rows.sort(key=lambda row: row["relative_path"])
        write_tsv(os.path.join(staging, "output_manifest.tsv"), ["relative_path", "bytes", "sha256"], manifest_rows)
        if os.path.exists(output_root):
            fail("X-GFP output appeared during processing and will not be overwritten")
        os.replace(staging, output_root)
        staging = None
        verify_existing(output_root)
        print("Published validated X-GFP audit: 3/3 exact probe matches", flush=True)
    finally:
        if staging and os.path.isdir(staging):
            shutil.rmtree(staging)


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--lock", required=True)
    parser.add_argument("--probe-set", required=True)
    parser.add_argument("--paper3-root", required=True)
    return parser.parse_args()


if __name__ == "__main__":
    try:
        run(parse_args())
    except Exception as exc:
        sys.stderr.write("FATAL: {}\n".format(exc))
        raise
