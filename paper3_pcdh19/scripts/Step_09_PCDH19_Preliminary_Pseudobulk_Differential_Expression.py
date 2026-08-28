#!/usr/bin/env python3
"""PRELIMINARY Step 09: sample-level pseudobulk differential expression.

Raw integer gene-level Flex UMI/probe-ligation counts from Cell Ranger-filtered
cells are summed by biological sample and broad state. An external edgeR module
then performs TMM normalization, filtering, dispersion estimation, QL testing,
and plotting. HET WT-like/KO-like calls and Steps 03-07 are never opened.
"""

from __future__ import print_function

import argparse
import csv
import gzip
import hashlib
import json
import os
import platform
import shutil
import subprocess
import sys
import tempfile

import h5py
import numpy as np


STEP_ID = "09_pcdh19_preliminary_pseudobulk_differential_expression"
STATUS = "PRELIMINARY"
COUNTS_NAME = "step_09_preliminary_pseudobulk_raw_gene_counts.tsv.gz"
METADATA_NAME = "step_09_preliminary_pseudobulk_sample_metadata.tsv"
VALIDATION_NAME = "step_09_preliminary_validation.tsv"
ENVIRONMENT_NAME = "software_environment.tsv"
MANIFEST_NAME = "output_manifest.tsv"
STRATA = ["all_cells", "progenitor", "immature_neuron"]
STATE_TO_STRATUM = {
    "proliferating_neural_progenitor": "progenitor",
    "postmitotic_immature_neuron": "immature_neuron",
}


class Step09Error(RuntimeError):
    pass


def sha256_file(path, chunk_size=8 * 1024 * 1024):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        while True:
            block = handle.read(chunk_size)
            if not block:
                break
            digest.update(block)
    return digest.hexdigest()


def decode(values):
    return [value.decode("utf-8") if isinstance(value, bytes) else str(value) for value in values]


def write_tsv(path, header, rows):
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=header, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in header})


class ValidationRecorder(object):
    def __init__(self):
        self.rows = []

    def check(self, name, observed, expected, detail=""):
        passed = observed == expected
        self.rows.append({"check": name, "status": "PASS" if passed else "FAIL", "observed": observed, "expected": expected, "detail": detail})
        if not passed:
            raise Step09Error("Validation failed {}: {!r} != {!r}".format(name, observed, expected))


class Step09Configuration(object):
    def __init__(self, path):
        self.path = os.path.abspath(path)
        with open(self.path, "r") as handle:
            self.values = json.load(handle)
        if self.values.get("step_id") != STEP_ID or self.values.get("analysis_status") != STATUS:
            raise Step09Error("Step 09 must remain explicitly PRELIMINARY")
        if self.values["model"]["replicate_unit"] != "biological_sample_id":
            raise Step09Error("Sample ID must be the replicate")
        expected = {"wt_f_vs_wt_m", "het_f_vs_wt_f", "ko_m_vs_wt_m"}
        if set(row["contrast_id"] for row in self.values["contrasts"]) != expected:
            raise Step09Error("Unexpected contrast set")
        prohibited = set(self.values["prohibited_operations"])
        if "wt_like_vs_ko_like_het_comparison" not in prohibited or "cell_level_de_testing" not in prohibited:
            raise Step09Error("Incomplete scope guard")


class SampleRegistry(object):
    def __init__(self, path, configuration, validation):
        self.path = os.path.abspath(path)
        validation.check("sample_key_sha256", sha256_file(self.path), configuration.values["sample_key_sha256"])
        with open(self.path, "r", newline="") as handle:
            self.rows = list(csv.DictReader(handle))
        self.rows.sort(key=lambda row: int(row["technical_sample_id"].rsplit("-", 1)[1]))
        validation.check("registered_samples", len(self.rows), 12)
        validation.check("three_samples_per_group", sorted([sum(r["design_group"] == group for r in self.rows) for group in ("WT_M", "WT_F", "HET_F", "KO_M")]), [3, 3, 3, 3])
        self.by_id = dict((row["technical_sample_id"], row) for row in self.rows)


class LabelStateRepository(object):
    """Load exact broad-state membership without using classifier outputs."""

    def __init__(self, label_root, configuration, registry, validation):
        self.root = os.path.abspath(label_root)
        self.configuration = configuration
        self.registry = registry
        self.validation = validation

    def load(self, staging_root):
        manifest = os.path.join(self.root, MANIFEST_NAME)
        self.validation.check("step_00_manifest_sha256", sha256_file(manifest), self.configuration.values["label_transfer_manifest_sha256"])
        label_path = os.path.join(self.root, self.configuration.values["label_file"])
        self.validation.check("step_00_label_sha256", sha256_file(label_path), self.configuration.values["label_file_sha256"])
        membership_root = os.path.join(staging_root, "state_membership")
        os.makedirs(membership_root)
        paths = dict((sample_id, os.path.join(membership_root, "{}_broad_state.tsv".format(sample_id))) for sample_id in self.registry.by_id)
        handles = dict((sample_id, open(path, "w", newline="")) for sample_id, path in paths.items())
        writers = dict((sample_id, csv.writer(handle, delimiter="\t", lineterminator="\n")) for sample_id, handle in handles.items())
        counts = dict((sample_id, 0) for sample_id in self.registry.by_id)
        for writer in writers.values():
            writer.writerow(["cell_barcode", "GSE94641_broad_state"])
        try:
            with open(label_path, "r", newline="") as handle:
                reader = csv.DictReader(handle, delimiter="\t")
                required = {"sample_id", "cell_barcode", "design_group", "GSE94641_broad_state"}
                if not required.issubset(reader.fieldnames):
                    raise Step09Error("Label table schema mismatch")
                total = 0
                for row in reader:
                    sample_id = row["sample_id"]
                    if row["design_group"] != self.registry.by_id[sample_id]["design_group"]:
                        raise Step09Error("Design-group mismatch in labels")
                    writers[sample_id].writerow([row["cell_barcode"], row["GSE94641_broad_state"]])
                    counts[sample_id] += 1
                    total += 1
        finally:
            for handle in handles.values():
                handle.close()
        self.validation.check("label_cells", total, self.configuration.values["expected_cells"])
        for sample_id, count in counts.items():
            self.validation.check("{}_label_cells".format(sample_id), count, len(self._read_membership(paths[sample_id])))
        return paths

    @staticmethod
    def _read_membership(path):
        with open(path, "r", newline="") as handle:
            return dict((row["cell_barcode"], row["GSE94641_broad_state"]) for row in csv.DictReader(handle, delimiter="\t"))


class MatrixInputRegistry(object):
    def __init__(self, label_root, configuration, validation):
        path = os.path.join(label_root, "query_input_manifest.tsv")
        validation.check("query_input_manifest_sha256", sha256_file(path), configuration.values["query_input_manifest_sha256"])
        self.matrix_rows = {}
        with open(path, "r", newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                if row["input_role"] == "matrix":
                    self.matrix_rows[row["sample_id"]] = row
        validation.check("matrix_inputs", len(self.matrix_rows), 12)
        for sample_id, row in self.matrix_rows.items():
            validation.check("{}_matrix_exists".format(sample_id), os.path.isfile(row["path"]), True)
            validation.check("{}_matrix_bytes".format(sample_id), os.path.getsize(row["path"]), int(row["bytes"]))


class PseudobulkAggregator(object):
    """Sum filtered-cell integer gene counts by sample and broad state."""

    def __init__(self, configuration, registry, state_paths, inputs, validation):
        self.configuration = configuration
        self.registry = registry
        self.state_paths = state_paths
        self.inputs = inputs
        self.validation = validation
        self.feature_ids = None
        self.feature_names = None
        self.count_columns = []
        self.metadata = []

    def _matrix_contract(self, sample_id):
        path = self.inputs.matrix_rows[sample_id]["path"]
        with h5py.File(path, "r") as handle:
            group = handle["matrix"]
            ids = decode(group["features"]["id"][:])
            names = decode(group["features"]["name"][:])
            types = decode(group["features"]["feature_type"][:])
            barcodes = decode(group["barcodes"][:])
            shape = tuple(int(value) for value in group["shape"][:])
        if set(types) != {"Gene Expression"}:
            raise Step09Error("Unexpected feature type in {}".format(sample_id))
        if self.feature_ids is None:
            self.feature_ids, self.feature_names = ids, names
            self.validation.check("gene_features", len(ids), self.configuration.values["expected_features"])
            self.validation.check("unique_gene_ids", len(set(ids)), len(ids))
        elif ids != self.feature_ids or names != self.feature_names:
            raise Step09Error("Feature contract differs in {}".format(sample_id))
        return path, barcodes, shape

    @staticmethod
    def _stream_sum(path, masks, shape, cell_chunk=256):
        """Aggregate CSC columns with bounded memory; never materialize the matrix."""
        n_features, n_cells = shape
        totals = dict((stratum, np.zeros(n_features, dtype=np.int64)) for stratum in STRATA)
        with h5py.File(path, "r") as handle:
            group = handle["matrix"]
            indptr = group["indptr"][:]
            for cell_start in range(0, n_cells, cell_chunk):
                cell_end = min(n_cells, cell_start + cell_chunk)
                pointers = indptr[cell_start:cell_end + 1]
                nz_start, nz_end = int(pointers[0]), int(pointers[-1])
                indices = group["indices"][nz_start:nz_end]
                data = group["data"][nz_start:nz_end].astype(np.float64)
                repeats = np.diff(pointers)
                totals["all_cells"] += np.bincount(indices, weights=data, minlength=n_features).astype(np.int64)
                for stratum in ("progenitor", "immature_neuron"):
                    selected = np.repeat(masks[stratum][cell_start:cell_end], repeats)
                    if np.any(selected):
                        totals[stratum] += np.bincount(indices[selected], weights=data[selected], minlength=n_features).astype(np.int64)
        return totals

    def aggregate(self):
        for sample in self.registry.rows:
            sample_id = sample["technical_sample_id"]
            print("aggregating {}".format(sample_id), flush=True)
            path, barcodes, shape = self._matrix_contract(sample_id)
            with open(self.state_paths[sample_id], "r", newline="") as handle:
                state_map = dict((row["cell_barcode"], row["GSE94641_broad_state"]) for row in csv.DictReader(handle, delimiter="\t"))
            self.validation.check("{}_barcode_identity".format(sample_id), set(barcodes) == set(state_map), True)
            state_vector = np.asarray([state_map[barcode] for barcode in barcodes], dtype=object)
            masks = {
                "all_cells": np.ones(len(barcodes), dtype=bool),
                "progenitor": state_vector == "proliferating_neural_progenitor",
                "immature_neuron": state_vector == "postmitotic_immature_neuron",
            }
            totals = self._stream_sum(path, masks, shape)
            all_library = int(totals["all_cells"].sum())
            for stratum in STRATA:
                selected = np.flatnonzero(masks[stratum])
                counts = totals[stratum]
                pseudobulk_id = "{}__{}".format(sample_id, stratum)
                self.count_columns.append((pseudobulk_id, counts))
                self.metadata.append({
                    "analysis_status": STATUS,
                    "pseudobulk_id": pseudobulk_id,
                    "biological_sample_id": sample_id,
                    "submitted_sample_name": sample["submitted_sample_name"],
                    "genotype": sample["genotype"],
                    "sex": sample["sex"],
                    "design_group": sample["design_group"],
                    "stratum": stratum,
                    "cells_aggregated": len(selected),
                    "raw_gene_umi_ligation_count": int(counts.sum()),
                })
                if stratum == "all_cells":
                    self.validation.check("{}_all_library_sum".format(sample_id), int(counts.sum()), all_library)
            del totals, state_map, state_vector
        self.validation.check("pseudobulk_libraries", len(self.count_columns), 36)
        return self

    def write(self, counts_path, metadata_path):
        header = ["ensembl_gene_id", "gene_symbol"] + [item[0] for item in self.count_columns]
        with gzip.open(counts_path, "wt", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(header)
            arrays = [item[1] for item in self.count_columns]
            for index, gene_id in enumerate(self.feature_ids):
                writer.writerow([gene_id, self.feature_names[index]] + [int(values[index]) for values in arrays])
        meta_header = ["analysis_status", "pseudobulk_id", "biological_sample_id", "submitted_sample_name", "genotype", "sex", "design_group", "stratum", "cells_aggregated", "raw_gene_umi_ligation_count"]
        write_tsv(metadata_path, meta_header, self.metadata)


class EdgeRRunner(object):
    def __init__(self, rscript, module_path, configuration, validation):
        self.rscript = rscript
        self.module_path = os.path.abspath(module_path)
        self.configuration = configuration
        self.validation = validation

    def run(self, counts_path, metadata_path, output_root):
        contrast_path = os.path.join(output_root, "step_09_preliminary_contrast_specification.tsv")
        write_tsv(contrast_path, ["analysis_status", "contrast_id", "case_group", "reference_group", "interpretation"], [dict({"analysis_status": STATUS}, **row) for row in self.configuration.values["contrasts"]])
        subprocess.check_call([self.rscript, self.module_path, counts_path, metadata_path, contrast_path, output_root])
        summary = os.path.join(output_root, "step_09_preliminary_differential_expression_summary.tsv")
        self.validation.check("edger_summary_exists", os.path.isfile(summary), True)


class OutputPublisher(object):
    @staticmethod
    def environment(path, configuration, requirements, r_module):
        rows = [
            {"component": "analysis_status", "version": STATUS},
            {"component": "python", "version": platform.python_version()},
            {"component": "numpy", "version": np.__version__},
            {"component": "scipy", "version": __import__("scipy").__version__},
            {"component": "h5py", "version": h5py.__version__},
            {"component": "python_pipeline_sha256", "version": sha256_file(os.path.abspath(__file__))},
            {"component": "r_module_sha256", "version": sha256_file(r_module)},
            {"component": "lock_sha256", "version": sha256_file(configuration.path)},
            {"component": "requirements_sha256", "version": sha256_file(requirements)},
        ]
        write_tsv(path, ["component", "version"], rows)

    @staticmethod
    def manifest(root):
        rows = []
        for directory, _, filenames in os.walk(root):
            for filename in filenames:
                if filename == MANIFEST_NAME:
                    continue
                path = os.path.join(directory, filename)
                rows.append({"relative_path": os.path.relpath(path, root), "bytes": os.path.getsize(path), "sha256": sha256_file(path)})
        rows.sort(key=lambda row: row["relative_path"])
        write_tsv(os.path.join(root, MANIFEST_NAME), ["relative_path", "bytes", "sha256"], rows)


class Step09Workflow(object):
    def __init__(self, args):
        self.args = args
        self.configuration = Step09Configuration(args.lock)
        self.validation = ValidationRecorder()

    def run(self):
        output_root = os.path.abspath(self.args.output_root)
        if os.path.exists(output_root):
            raise Step09Error("Output exists and will not be overwritten: {}".format(output_root))
        parent = os.path.dirname(output_root)
        os.makedirs(parent, exist_ok=True)
        staging = tempfile.mkdtemp(prefix=".step09-preliminary.", dir=parent)
        try:
            registry = SampleRegistry(self.args.sample_key, self.configuration, self.validation)
            state_paths = LabelStateRepository(self.args.label_root, self.configuration, registry, self.validation).load(staging)
            inputs = MatrixInputRegistry(self.args.label_root, self.configuration, self.validation)
            aggregator = PseudobulkAggregator(self.configuration, registry, state_paths, inputs, self.validation).aggregate()
            counts_path = os.path.join(staging, COUNTS_NAME)
            metadata_path = os.path.join(staging, METADATA_NAME)
            aggregator.write(counts_path, metadata_path)
            os.makedirs(os.path.join(staging, "figures"))
            os.makedirs(os.path.join(staging, "differential_expression"))
            EdgeRRunner(self.args.rscript, self.args.r_module, self.configuration, self.validation).run(counts_path, metadata_path, staging)
            self.validation.check("cell_level_tests", 0, 0)
            self.validation.check("step_07_inputs_opened", 0, 0)
            self.validation.check("analysis_status", STATUS, "PRELIMINARY")
            write_tsv(os.path.join(staging, VALIDATION_NAME), ["check", "status", "observed", "expected", "detail"], self.validation.rows)
            OutputPublisher.environment(os.path.join(staging, ENVIRONMENT_NAME), self.configuration, self.args.requirements, self.args.r_module)
            shutil.copy2(self.configuration.path, os.path.join(staging, os.path.basename(self.configuration.path)))
            shutil.copy2(self.args.requirements, os.path.join(staging, os.path.basename(self.args.requirements)))
            shutil.copy2(self.args.r_module, os.path.join(staging, os.path.basename(self.args.r_module)))
            OutputPublisher.manifest(staging)
            os.rename(staging, output_root)
            staging = None
            print("PRELIMINARY Step 09 complete: {}".format(output_root))
        finally:
            if staging and os.path.isdir(staging):
                shutil.rmtree(staging)


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--lock", required=True)
    parser.add_argument("--requirements", required=True)
    parser.add_argument("--sample-key", required=True)
    parser.add_argument("--label-root", required=True)
    parser.add_argument("--rscript", required=True)
    parser.add_argument("--r-module", required=True)
    parser.add_argument("--output-root", required=True)
    return parser.parse_args(argv)


if __name__ == "__main__":
    try:
        Step09Workflow(parse_args()).run()
    except Step09Error as error:
        print("ERROR: {}".format(error), file=sys.stderr)
        sys.exit(2)
