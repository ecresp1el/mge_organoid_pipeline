#!/usr/bin/env python3
"""Build the Step 03 PCDH19 WT/KO-male classification-ready cell table.

This step is the first biological annotation layer after the locked probe
audits. It reads the registered sample key and the existing Step 02a
one-row-per-filtered-barcode Pcdh19 probe tables, identifies WT-male and
KO-male cells as the known ground-truth classes, and writes one deterministic
intermediate table for later model development.

The output retains technical sample identity, registered genotype and sex,
cell barcode, raw Probe A/B/C UMI counts, per-probe detection flags, the
three-probe total, and the existing detection pattern. It does not normalize
counts, split data, encode predictors beyond direct detection flags, fit a
classifier, score HET cells, or evaluate predictions.

The implementation deliberately separates configuration, metadata, probe
observation, input loading, and output publication responsibilities so future
steps can consume the stable table without extending this setup step into a
monolithic prediction program. It uses only the Python standard library and
is compatible with the existing Paper 3 Python 3.6 environment.
"""

import argparse
import csv
import datetime as dt
import hashlib
import json
import os
import platform
import shutil
import sys
import tempfile


STEP_OUTPUT_DIRECTORY = "step_03_pcdh19_genotype_classification_setup"
CLASSIFICATION_TABLE_NAME = (
    "pcdh19_wt_ko_male_genotype_classification_ready_cells.tsv"
)
SUMMARY_NAME = "pcdh19_genotype_classification_ground_truth_summary.tsv"
VALIDATION_NAME = "pcdh19_genotype_classification_setup_validation.tsv"
ENVIRONMENT_NAME = "software_environment.tsv"
MANIFEST_NAME = "output_manifest.tsv"

CLASSIFICATION_HEADER = [
    "technical_sample_id",
    "submitted_sample_name",
    "genotype",
    "sex",
    "design_group",
    "classification_role",
    "ground_truth_class",
    "ground_truth_label",
    "cell_barcode",
    "A_UMI",
    "A_detected",
    "B_UMI",
    "B_detected",
    "C_UMI",
    "C_detected",
    "Pcdh19_total_UMI",
    "detection_pattern",
]

SUMMARY_HEADER = [
    "technical_sample_id",
    "submitted_sample_name",
    "genotype",
    "sex",
    "design_group",
    "ground_truth_class",
    "ground_truth_label",
    "cells",
    "pcdh19_detected_cells",
    "A_detected_cells",
    "B_detected_cells",
    "C_detected_cells",
    "total_pcdh19_umi",
]

VALIDATION_HEADER = [
    "step_id", "check_name", "status", "observed", "expected", "details"
]


class ClassificationSetupError(RuntimeError):
    """Raised when a Step 03 input or output invariant is violated."""


def sha256_file(path, chunk_size=8 * 1024 * 1024):
    """Return a lowercase SHA-256 digest without loading a file into memory."""
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        while True:
            chunk = handle.read(chunk_size)
            if not chunk:
                break
            digest.update(chunk)
    return digest.hexdigest()


def write_tsv(path, header, rows):
    """Write dictionaries as a deterministic, LF-terminated TSV."""
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(
            handle, delimiter="\t", fieldnames=header, lineterminator="\n"
        )
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in header})


class StepConfiguration(object):
    """Load and validate the frozen scientific scope of Step 03."""

    REQUIRED_KEYS = {
        "step_id",
        "pipeline_version",
        "upstream_step_id",
        "upstream_result_directory",
        "sample_key_sha256",
        "expected_sample_ids",
        "ground_truth_classes",
        "required_probe_pattern_columns",
    }

    def __init__(self, lock_path):
        self.lock_path = os.path.abspath(lock_path)
        if not os.path.isfile(self.lock_path):
            raise ClassificationSetupError(
                "Step 03 lock is missing: {}".format(self.lock_path)
            )
        with open(self.lock_path, "r") as handle:
            self.values = json.load(handle)
        missing = sorted(self.REQUIRED_KEYS.difference(self.values))
        if missing:
            raise ClassificationSetupError(
                "Step 03 lock is missing keys: {}".format(", ".join(missing))
            )
        self._validate_ground_truth_classes()

    def _validate_ground_truth_classes(self):
        expected = {("WT", "M", "WT", 0), ("KO", "M", "KO", 1)}
        observed = set()
        for item in self.values["ground_truth_classes"]:
            observed.add(
                (
                    item.get("genotype"),
                    item.get("sex"),
                    item.get("ground_truth_class"),
                    item.get("ground_truth_label"),
                )
            )
        if observed != expected:
            raise ClassificationSetupError(
                "Ground-truth classes must be exactly WT-male=0 and KO-male=1"
            )

    @property
    def step_id(self):
        return self.values["step_id"]

    @property
    def pipeline_version(self):
        return self.values["pipeline_version"]

    @property
    def expected_sample_ids(self):
        return list(self.values["expected_sample_ids"])

    @property
    def required_probe_columns(self):
        return list(self.values["required_probe_pattern_columns"])

    def ground_truth_for(self, genotype, sex):
        """Return the configured target class for one genotype/sex pair."""
        for item in self.values["ground_truth_classes"]:
            if item["genotype"] == genotype and item["sex"] == sex:
                return item["ground_truth_class"], item["ground_truth_label"]
        return None


class SampleAnnotation(object):
    """Represent one registered technical sample and its biological labels."""

    REQUIRED_COLUMNS = {
        "technical_sample_id",
        "submitted_sample_name",
        "genotype",
        "sex",
        "design_group",
    }

    def __init__(self, row, ground_truth):
        self.technical_sample_id = row["technical_sample_id"]
        self.submitted_sample_name = row["submitted_sample_name"]
        self.genotype = row["genotype"]
        self.sex = row["sex"]
        self.design_group = row["design_group"]
        self.ground_truth = ground_truth

    @property
    def is_ground_truth(self):
        return self.ground_truth is not None

    @property
    def ground_truth_class(self):
        return self.ground_truth[0] if self.ground_truth else None

    @property
    def ground_truth_label(self):
        return self.ground_truth[1] if self.ground_truth else None


class ProbeObservation(object):
    """Validate and represent raw Probe A/B/C information for one cell."""

    PATTERNS = {
        (1, 1, 1): "A+B+C",
        (1, 1, 0): "A+B",
        (1, 0, 1): "A+C",
        (1, 0, 0): "A only",
        (0, 1, 1): "B+C",
        (0, 1, 0): "B only",
        (0, 0, 1): "C only",
        (0, 0, 0): "none",
    }

    def __init__(self, row):
        self.cell_barcode = row.get("barcode", "")
        if not self.cell_barcode:
            raise ClassificationSetupError("Encountered an empty cell barcode")
        try:
            self.a_umi = int(row["A_UMI"])
            self.b_umi = int(row["B_UMI"])
            self.c_umi = int(row["C_UMI"])
            self.total_umi = int(row["Pcdh19_total_UMI"])
        except (KeyError, TypeError, ValueError) as error:
            raise ClassificationSetupError(
                "Invalid probe count for barcode {}: {}".format(
                    self.cell_barcode, error
                )
            )
        if min(self.a_umi, self.b_umi, self.c_umi, self.total_umi) < 0:
            raise ClassificationSetupError(
                "Negative probe count for barcode {}".format(self.cell_barcode)
            )
        if self.a_umi + self.b_umi + self.c_umi != self.total_umi:
            raise ClassificationSetupError(
                "Probe counts do not sum for barcode {}".format(self.cell_barcode)
            )
        self.a_detected = int(self.a_umi > 0)
        self.b_detected = int(self.b_umi > 0)
        self.c_detected = int(self.c_umi > 0)
        expected_pattern = self.PATTERNS[
            (self.a_detected, self.b_detected, self.c_detected)
        ]
        self.detection_pattern = row.get("detection_pattern", "")
        if self.detection_pattern != expected_pattern:
            raise ClassificationSetupError(
                "Detection pattern mismatch for barcode {}: observed {!r}, "
                "expected {!r}".format(
                    self.cell_barcode, self.detection_pattern, expected_pattern
                )
            )

    @property
    def pcdh19_detected(self):
        return int(self.total_umi > 0)


class PCDH19ProbeCellLoader(object):
    """Load registered sample metadata and validated Step 02a cell tables."""

    def __init__(self, configuration, sample_key_path, probe_audit_root):
        self.configuration = configuration
        self.sample_key_path = os.path.abspath(sample_key_path)
        self.probe_audit_root = os.path.abspath(probe_audit_root)
        self.manifest_path = os.path.join(
            self.probe_audit_root, "output_manifest.tsv"
        )
        self.manifest = self._read_upstream_manifest()

    def _read_upstream_manifest(self):
        if not os.path.isfile(self.manifest_path):
            raise ClassificationSetupError(
                "Step 02a output manifest is missing: {}".format(
                    self.manifest_path
                )
            )
        manifest = {}
        with open(self.manifest_path, "r", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if reader.fieldnames != ["relative_path", "bytes", "sha256"]:
                raise ClassificationSetupError(
                    "Unexpected Step 02a output manifest schema"
                )
            for row in reader:
                manifest[row["relative_path"]] = row
        return manifest

    def load_sample_annotations(self):
        expected_sha = self.configuration.values["sample_key_sha256"]
        if not os.path.isfile(self.sample_key_path):
            raise ClassificationSetupError(
                "Registered sample key is missing: {}".format(
                    self.sample_key_path
                )
            )
        observed_sha = sha256_file(self.sample_key_path)
        if observed_sha != expected_sha:
            raise ClassificationSetupError(
                "Registered sample-key checksum mismatch: observed {}, "
                "expected {}".format(observed_sha, expected_sha)
            )
        annotations = []
        with open(self.sample_key_path, "r", newline="") as handle:
            reader = csv.DictReader(handle)
            missing = sorted(
                SampleAnnotation.REQUIRED_COLUMNS.difference(
                    reader.fieldnames or []
                )
            )
            if missing:
                raise ClassificationSetupError(
                    "Sample key is missing columns: {}".format(", ".join(missing))
                )
            seen = set()
            for row in reader:
                sample_id = row["technical_sample_id"]
                if sample_id in seen:
                    raise ClassificationSetupError(
                        "Duplicate sample-key ID: {}".format(sample_id)
                    )
                seen.add(sample_id)
                ground_truth = self.configuration.ground_truth_for(
                    row["genotype"], row["sex"]
                )
                annotations.append(SampleAnnotation(row, ground_truth))
        observed_ids = [item.technical_sample_id for item in annotations]
        if observed_ids != self.configuration.expected_sample_ids:
            raise ClassificationSetupError(
                "Sample-key IDs/order differ from the Step 03 lock"
            )
        ground_truth = [item for item in annotations if item.is_ground_truth]
        observed_groups = {
            (item.genotype, item.sex, item.ground_truth_class)
            for item in ground_truth
        }
        expected_groups = {("WT", "M", "WT"), ("KO", "M", "KO")}
        if observed_groups != expected_groups:
            raise ClassificationSetupError(
                "Ground-truth samples are not exactly WT-male and KO-male"
            )
        return annotations

    def _validated_pattern_path(self, sample_id):
        relative_path = os.path.join(
            "per_sample", sample_id, "pcdh19_probe_patterns.tsv"
        )
        manifest_row = self.manifest.get(relative_path)
        if manifest_row is None:
            raise ClassificationSetupError(
                "Step 02a manifest does not list {}".format(relative_path)
            )
        path = os.path.join(self.probe_audit_root, relative_path)
        if not os.path.isfile(path):
            raise ClassificationSetupError(
                "Step 02a per-cell input is missing: {}".format(path)
            )
        observed_size = os.path.getsize(path)
        if observed_size != int(manifest_row["bytes"]):
            raise ClassificationSetupError(
                "Step 02a per-cell input size mismatch: {}".format(path)
            )
        observed_sha = sha256_file(path)
        if observed_sha != manifest_row["sha256"]:
            raise ClassificationSetupError(
                "Step 02a per-cell input checksum mismatch: {}".format(path)
            )
        return path

    def iter_probe_observations(self, sample):
        """Yield validated observations for one ground-truth sample."""
        if not sample.is_ground_truth:
            raise ClassificationSetupError(
                "Non-ground-truth sample requested: {}".format(
                    sample.technical_sample_id
                )
            )
        path = self._validated_pattern_path(sample.technical_sample_id)
        with open(path, "r", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if reader.fieldnames != self.configuration.required_probe_columns:
                raise ClassificationSetupError(
                    "Unexpected Step 02a pattern-table schema for {}".format(
                        sample.technical_sample_id
                    )
                )
            seen_barcodes = set()
            for row in reader:
                observation = ProbeObservation(row)
                if observation.cell_barcode in seen_barcodes:
                    raise ClassificationSetupError(
                        "Duplicate cell barcode within sample {}: {}".format(
                            sample.technical_sample_id,
                            observation.cell_barcode,
                        )
                    )
                seen_barcodes.add(observation.cell_barcode)
                yield observation


class ClassificationReadyTableWriter(object):
    """Serialize the stable WT/KO ground-truth table and compact provenance."""

    def __init__(self, configuration, loader, output_root, script_path):
        self.configuration = configuration
        self.loader = loader
        self.output_root = os.path.abspath(output_root)
        self.script_path = os.path.abspath(script_path)
        self.validation_rows = []

    def _validation(self, name, observed, expected, details=""):
        status = "PASS" if str(observed) == str(expected) else "FAIL"
        row = {
            "step_id": self.configuration.step_id,
            "check_name": name,
            "status": status,
            "observed": observed,
            "expected": expected,
            "details": details,
        }
        self.validation_rows.append(row)
        if status != "PASS":
            raise ClassificationSetupError(
                "Validation failed for {}: observed {}, expected {}".format(
                    name, observed, expected
                )
            )

    @staticmethod
    def _classification_row(sample, observation):
        return {
            "technical_sample_id": sample.technical_sample_id,
            "submitted_sample_name": sample.submitted_sample_name,
            "genotype": sample.genotype,
            "sex": sample.sex,
            "design_group": sample.design_group,
            "classification_role": "ground_truth_training_candidate",
            "ground_truth_class": sample.ground_truth_class,
            "ground_truth_label": sample.ground_truth_label,
            "cell_barcode": observation.cell_barcode,
            "A_UMI": observation.a_umi,
            "A_detected": observation.a_detected,
            "B_UMI": observation.b_umi,
            "B_detected": observation.b_detected,
            "C_UMI": observation.c_umi,
            "C_detected": observation.c_detected,
            "Pcdh19_total_UMI": observation.total_umi,
            "detection_pattern": observation.detection_pattern,
        }

    def _write_classification_table(self, path, annotations):
        summaries = []
        total_cells = 0
        combined_keys = set()
        with open(path, "w", newline="") as handle:
            writer = csv.DictWriter(
                handle,
                delimiter="\t",
                fieldnames=CLASSIFICATION_HEADER,
                lineterminator="\n",
            )
            writer.writeheader()
            for sample in annotations:
                if not sample.is_ground_truth:
                    continue
                counts = {
                    "cells": 0,
                    "pcdh19_detected_cells": 0,
                    "A_detected_cells": 0,
                    "B_detected_cells": 0,
                    "C_detected_cells": 0,
                    "total_pcdh19_umi": 0,
                }
                for observation in self.loader.iter_probe_observations(sample):
                    key = (sample.technical_sample_id, observation.cell_barcode)
                    if key in combined_keys:
                        raise ClassificationSetupError(
                            "Duplicate sample/barcode key: {} {}".format(*key)
                        )
                    combined_keys.add(key)
                    writer.writerow(self._classification_row(sample, observation))
                    counts["cells"] += 1
                    counts["pcdh19_detected_cells"] += observation.pcdh19_detected
                    counts["A_detected_cells"] += observation.a_detected
                    counts["B_detected_cells"] += observation.b_detected
                    counts["C_detected_cells"] += observation.c_detected
                    counts["total_pcdh19_umi"] += observation.total_umi
                if counts["cells"] == 0:
                    raise ClassificationSetupError(
                        "Ground-truth sample has no cells: {}".format(
                            sample.technical_sample_id
                        )
                    )
                total_cells += counts["cells"]
                summary = {
                    "technical_sample_id": sample.technical_sample_id,
                    "submitted_sample_name": sample.submitted_sample_name,
                    "genotype": sample.genotype,
                    "sex": sample.sex,
                    "design_group": sample.design_group,
                    "ground_truth_class": sample.ground_truth_class,
                    "ground_truth_label": sample.ground_truth_label,
                }
                summary.update(counts)
                summaries.append(summary)
        self._validation(
            "classification_table_rows_equal_unique_sample_barcode_keys",
            total_cells,
            len(combined_keys),
        )
        self._validation(
            "ground_truth_samples_are_three_wt_male_and_three_ko_male",
            ",".join(
                "{}:{}".format(key, sum(1 for item in summaries if item["ground_truth_class"] == key))
                for key in ("WT", "KO")
            ),
            "WT:3,KO:3",
        )
        return summaries, total_cells

    def _write_environment(self, path):
        rows = [
            {"key": "step_id", "value": self.configuration.step_id},
            {"key": "pipeline_version", "value": self.configuration.pipeline_version},
            {"key": "python_version", "value": platform.python_version()},
            {"key": "platform", "value": platform.platform()},
            {"key": "python_script", "value": self.script_path},
            {"key": "python_script_sha256", "value": sha256_file(self.script_path)},
            {"key": "lock_path", "value": self.configuration.lock_path},
            {"key": "lock_sha256", "value": sha256_file(self.configuration.lock_path)},
            {"key": "sample_key_path", "value": self.loader.sample_key_path},
            {"key": "sample_key_sha256", "value": sha256_file(self.loader.sample_key_path)},
            {"key": "upstream_step", "value": self.configuration.values["upstream_step_id"]},
            {"key": "upstream_manifest_path", "value": self.loader.manifest_path},
            {"key": "upstream_manifest_sha256", "value": sha256_file(self.loader.manifest_path)},
            {"key": "classification_fitted", "value": "no"},
        ]
        write_tsv(path, ["key", "value"], rows)

    @staticmethod
    def _write_manifest(directory):
        rows = []
        for filename in sorted(os.listdir(directory)):
            path = os.path.join(directory, filename)
            if filename == MANIFEST_NAME or not os.path.isfile(path):
                continue
            rows.append(
                {
                    "relative_path": filename,
                    "bytes": os.path.getsize(path),
                    "sha256": sha256_file(path),
                }
            )
        write_tsv(
            os.path.join(directory, MANIFEST_NAME),
            ["relative_path", "bytes", "sha256"],
            rows,
        )

    @staticmethod
    def _verify_manifest(directory):
        manifest_path = os.path.join(directory, MANIFEST_NAME)
        if not os.path.isfile(manifest_path):
            raise ClassificationSetupError(
                "Existing Step 03 output is incomplete: {}".format(directory)
            )
        listed = set()
        with open(manifest_path, "r", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if reader.fieldnames != ["relative_path", "bytes", "sha256"]:
                raise ClassificationSetupError("Unexpected Step 03 manifest schema")
            for row in reader:
                listed.add(row["relative_path"])
                path = os.path.join(directory, row["relative_path"])
                if not os.path.isfile(path):
                    raise ClassificationSetupError(
                        "Manifest-listed Step 03 file is missing: {}".format(path)
                    )
                if os.path.getsize(path) != int(row["bytes"]):
                    raise ClassificationSetupError(
                        "Step 03 output size mismatch: {}".format(path)
                    )
                if sha256_file(path) != row["sha256"]:
                    raise ClassificationSetupError(
                        "Step 03 output checksum mismatch: {}".format(path)
                    )
        expected = {
            CLASSIFICATION_TABLE_NAME,
            SUMMARY_NAME,
            VALIDATION_NAME,
            ENVIRONMENT_NAME,
        }
        if listed != expected:
            raise ClassificationSetupError(
                "Existing Step 03 manifest does not list the expected files"
            )

    def _verify_existing_provenance(self):
        self._verify_manifest(self.output_root)
        environment_path = os.path.join(self.output_root, ENVIRONMENT_NAME)
        with open(environment_path, "r", newline="") as handle:
            environment = {
                row["key"]: row["value"]
                for row in csv.DictReader(handle, delimiter="\t")
            }
        expected = {
            "python_script_sha256": sha256_file(self.script_path),
            "lock_sha256": sha256_file(self.configuration.lock_path),
            "sample_key_sha256": sha256_file(self.loader.sample_key_path),
            "upstream_manifest_sha256": sha256_file(self.loader.manifest_path),
        }
        for key, value in expected.items():
            if environment.get(key) != value:
                raise ClassificationSetupError(
                    "Existing Step 03 provenance differs for {}".format(key)
                )

    def build(self):
        """Build or checksum-verify the complete Step 03 result package."""
        if os.path.exists(self.output_root):
            self._verify_existing_provenance()
            print("EXISTING_VALIDATED\t{}".format(self.output_root))
            return

        annotations = self.loader.load_sample_annotations()
        parent = os.path.dirname(self.output_root)
        os.makedirs(parent, exist_ok=True)
        staging = tempfile.mkdtemp(
            prefix=".step_03_pcdh19_genotype_classification_setup.", dir=parent
        )
        try:
            table_path = os.path.join(staging, CLASSIFICATION_TABLE_NAME)
            summaries, total_cells = self._write_classification_table(
                table_path, annotations
            )
            self._validation(
                "classification_table_contains_ground_truth_cells",
                int(total_cells > 0),
                1,
            )
            self._validation(
                "classifier_fitted_in_step_03",
                "no",
                "no",
                "Step 03 prepares data only",
            )
            write_tsv(os.path.join(staging, SUMMARY_NAME), SUMMARY_HEADER, summaries)
            write_tsv(
                os.path.join(staging, VALIDATION_NAME),
                VALIDATION_HEADER,
                self.validation_rows,
            )
            self._write_environment(os.path.join(staging, ENVIRONMENT_NAME))
            self._write_manifest(staging)
            self._verify_manifest(staging)
            if os.path.exists(self.output_root):
                raise ClassificationSetupError(
                    "Step 03 output appeared during publication: {}".format(
                        self.output_root
                    )
                )
            os.replace(staging, self.output_root)
            staging = None
            print(
                "PUBLISHED\t{}\tcells={}".format(self.output_root, total_cells)
            )
        finally:
            if staging and os.path.exists(staging):
                shutil.rmtree(staging)


class PCDH19GenotypeClassificationSetup(object):
    """Coordinate the narrowly scoped Step 03 loading and publication flow."""

    def __init__(
        self,
        lock_path,
        sample_key_path,
        probe_audit_root,
        paper3_root,
        script_path,
    ):
        configuration = StepConfiguration(lock_path)
        loader = PCDH19ProbeCellLoader(
            configuration, sample_key_path, probe_audit_root
        )
        output_root = os.path.join(
            os.path.abspath(paper3_root), "results", STEP_OUTPUT_DIRECTORY
        )
        self.writer = ClassificationReadyTableWriter(
            configuration, loader, output_root, script_path
        )

    def run(self):
        self.writer.build()


def build_parser():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--lock", required=True, help="Step 03 JSON lock")
    parser.add_argument(
        "--sample-key", required=True, help="Registered Paper 3 sample_key.csv"
    )
    parser.add_argument(
        "--probe-audit-root",
        required=True,
        help="Validated results/pcdh19_probe_audit directory",
    )
    parser.add_argument(
        "--paper3-root",
        required=True,
        help="Paper 3 generated-output root containing results/ and logs/",
    )
    return parser


def main(argv=None):
    args = build_parser().parse_args(argv)
    setup = PCDH19GenotypeClassificationSetup(
        lock_path=args.lock,
        sample_key_path=args.sample_key,
        probe_audit_root=args.probe_audit_root,
        paper3_root=args.paper3_root,
        script_path=__file__,
    )
    setup.run()


if __name__ == "__main__":
    try:
        main()
    except ClassificationSetupError as error:
        sys.stderr.write("ERROR: {}\n".format(error))
        sys.exit(1)
