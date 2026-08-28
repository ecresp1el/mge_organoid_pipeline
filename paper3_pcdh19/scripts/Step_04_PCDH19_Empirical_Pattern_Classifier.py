#!/usr/bin/env python3
"""Fit the Step 04 empirical PCDH19 A/B/C probe-pattern classifier.

Scientific question
-------------------
Given one of the eight binary Probe A/B/C detection patterns in a known
WT-male or KO-male cell, how much empirical evidence does that pattern provide
for WT versus KO genotype?

This step reads only the manifested Step 03 classification-ready table. The
model uses ``A_detected``, ``B_detected``, and ``C_detected`` as predictors and
the registered WT/KO ground-truth class as the outcome. It computes pattern
counts, posterior class probabilities, class-conditional pattern
probabilities, and likelihood ratios without smoothing.

No train/test split, hard genotype call, confusion matrix, logistic
regression, HET-cell input, or HET-cell prediction is produced. In particular,
``000`` remains an observed no-probe-detection state and is not forced into a
WT or KO call.

Architecture
------------
Input reading, pattern encoding, empirical estimation, probability lookup,
descriptive evaluation, plotting, and atomic output publication are separate
classes. ``ProbabilisticClassifier`` supplies the small prediction interface
that later validated classifiers can share.
"""

import argparse
import csv
import hashlib
import json
import math
import os
import platform
import shutil
import sys
import tempfile

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


STEP_OUTPUT_DIRECTORY = "step_04_pcdh19_empirical_pattern_classifier"
MODEL_TABLE_NAME = "step_04_pcdh19_empirical_pattern_classifier.tsv"
DISTRIBUTION_TABLE_NAME = "step_04_pcdh19_pattern_distribution.tsv"
VALIDATION_NAME = "step_04_pcdh19_empirical_pattern_classifier_validation.tsv"
ENVIRONMENT_NAME = "software_environment.tsv"
MANIFEST_NAME = "output_manifest.tsv"

FREQUENCY_PLOT_NAME = "step_04_wt_vs_ko_pattern_frequency.png"
POSTERIOR_PLOT_NAME = "step_04_wt_ko_conditional_probability_by_pattern.png"
ABUNDANCE_PLOT_NAME = "step_04_pattern_cell_count_and_proportion.png"

MODEL_HEADER = [
    "pattern_code",
    "pattern_label",
    "A_detected",
    "B_detected",
    "C_detected",
    "wt_cells",
    "ko_cells",
    "total_cells",
    "p_wt_given_pattern",
    "p_ko_given_pattern",
    "wt_to_ko_cell_ratio",
    "ko_to_wt_cell_ratio",
    "p_pattern_given_wt",
    "p_pattern_given_ko",
    "wt_to_ko_likelihood_ratio",
    "ko_to_wt_likelihood_ratio",
    "hard_call",
]

DISTRIBUTION_HEADER = [
    "pattern_code",
    "pattern_label",
    "total_cells",
    "fraction_of_all_cells",
    "percent_of_all_cells",
    "wt_cells",
    "fraction_of_wt_cells",
    "percent_of_wt_cells",
    "ko_cells",
    "fraction_of_ko_cells",
    "percent_of_ko_cells",
]

VALIDATION_HEADER = [
    "step_id", "check_name", "status", "observed", "expected", "details"
]


class Step04Error(RuntimeError):
    """Raised when a Step 04 scientific or publication invariant fails."""


def sha256_file(path, chunk_size=8 * 1024 * 1024):
    """Return the SHA-256 digest of a file using bounded-memory reads."""
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        while True:
            chunk = handle.read(chunk_size)
            if not chunk:
                break
            digest.update(chunk)
    return digest.hexdigest()


def write_tsv(path, header, rows):
    """Write dictionaries as a deterministic LF-terminated TSV."""
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(
            handle, delimiter="\t", fieldnames=header, lineterminator="\n"
        )
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in header})


class ValidationRecorder(object):
    """Collect machine-readable checks and fail immediately on violations."""

    def __init__(self, step_id):
        self.step_id = step_id
        self.rows = []

    def require_equal(self, name, observed, expected, details=""):
        status = "PASS" if str(observed) == str(expected) else "FAIL"
        self.rows.append(
            {
                "step_id": self.step_id,
                "check_name": name,
                "status": status,
                "observed": observed,
                "expected": expected,
                "details": details,
            }
        )
        if status != "PASS":
            raise Step04Error(
                "Validation failed for {}: observed {}, expected {}".format(
                    name, observed, expected
                )
            )

    def require_close(self, name, observed, expected, tolerance=1e-12, details=""):
        status = "PASS" if abs(observed - expected) <= tolerance else "FAIL"
        self.rows.append(
            {
                "step_id": self.step_id,
                "check_name": name,
                "status": status,
                "observed": "{:.15g}".format(observed),
                "expected": "{:.15g}".format(expected),
                "details": details,
            }
        )
        if status != "PASS":
            raise Step04Error(
                "Validation failed for {}: observed {}, expected {}".format(
                    name, observed, expected
                )
            )


class Step04Configuration(object):
    """Load the locked Step 04 input, model, and output scope."""

    REQUIRED_KEYS = {
        "step_id",
        "pipeline_version",
        "upstream_step_id",
        "upstream_output_manifest_sha256",
        "classification_table",
        "required_input_columns",
        "pattern_order",
        "ground_truth_classes",
        "expected_total_cells",
        "probability_decimal_places",
        "plot_dpi",
        "hard_calls_permitted",
        "het_cells_permitted",
    }

    def __init__(self, lock_path):
        self.lock_path = os.path.abspath(lock_path)
        if not os.path.isfile(self.lock_path):
            raise Step04Error("Step 04 lock is missing: {}".format(self.lock_path))
        with open(self.lock_path, "r") as handle:
            self.values = json.load(handle)
        missing = sorted(self.REQUIRED_KEYS.difference(self.values))
        if missing:
            raise Step04Error(
                "Step 04 lock is missing keys: {}".format(", ".join(missing))
            )
        expected_patterns = [
            "000", "001", "010", "011", "100", "101", "110", "111"
        ]
        if self.values["pattern_order"] != expected_patterns:
            raise Step04Error("Step 04 must contain the eight A/B/C patterns")
        if self.values["hard_calls_permitted"] is not False:
            raise Step04Error("Step 04 lock must prohibit hard genotype calls")
        if self.values["het_cells_permitted"] is not False:
            raise Step04Error("Step 04 lock must prohibit HET-cell input")

    @property
    def step_id(self):
        return self.values["step_id"]

    @property
    def pattern_order(self):
        return list(self.values["pattern_order"])

    @property
    def decimal_places(self):
        return int(self.values["probability_decimal_places"])

    @property
    def plot_dpi(self):
        return int(self.values["plot_dpi"])


class ProbePatternEncoder(object):
    """Encode and label the only eight permitted binary A/B/C states."""

    LABELS = {
        "000": "none",
        "001": "C only",
        "010": "B only",
        "011": "B+C",
        "100": "A only",
        "101": "A+C",
        "110": "A+B",
        "111": "A+B+C",
    }

    def __init__(self, pattern_order):
        self.pattern_order = list(pattern_order)
        if set(self.pattern_order) != set(self.LABELS):
            raise Step04Error("Pattern encoder requires exactly eight states")

    @staticmethod
    def _binary(value, probe_label):
        try:
            result = int(value)
        except (TypeError, ValueError):
            raise Step04Error(
                "{} detection is not binary: {!r}".format(probe_label, value)
            )
        if result not in (0, 1):
            raise Step04Error(
                "{} detection is not binary: {!r}".format(probe_label, value)
            )
        return result

    def encode(self, a_detected, b_detected, c_detected):
        """Return a three-character state ordered strictly as A, B, C."""
        values = (
            self._binary(a_detected, "A"),
            self._binary(b_detected, "B"),
            self._binary(c_detected, "C"),
        )
        code = "{}{}{}".format(*values)
        if code not in self.LABELS:
            raise Step04Error("Unexpected probe pattern: {}".format(code))
        return code

    def decode(self, pattern_code):
        if pattern_code not in self.LABELS:
            raise Step04Error("Unknown probe pattern: {}".format(pattern_code))
        return tuple(int(value) for value in pattern_code)

    def label(self, pattern_code):
        if pattern_code not in self.LABELS:
            raise Step04Error("Unknown probe pattern: {}".format(pattern_code))
        return self.LABELS[pattern_code]


class LabeledPatternObservation(object):
    """Minimal fit record containing only binary pattern and known class."""

    def __init__(self, pattern_code, ground_truth_class):
        self.pattern_code = pattern_code
        self.ground_truth_class = ground_truth_class


class Step03ClassificationTableReader(object):
    """Validate and stream the manifested Step 03 ground-truth table."""

    def __init__(
        self, configuration, encoder, validation, step03_root
    ):
        self.configuration = configuration
        self.encoder = encoder
        self.validation = validation
        self.step03_root = os.path.abspath(step03_root)
        self.manifest_path = os.path.join(self.step03_root, MANIFEST_NAME)
        table_relative = self.configuration.values["classification_table"][
            "relative_path"
        ]
        self.table_path = os.path.join(self.step03_root, table_relative)
        self._validate_upstream_package()

    def _validate_upstream_package(self):
        if not os.path.isfile(self.manifest_path):
            raise Step04Error(
                "Step 03 output manifest is missing: {}".format(self.manifest_path)
            )
        manifest_sha = sha256_file(self.manifest_path)
        expected_manifest_sha = self.configuration.values[
            "upstream_output_manifest_sha256"
        ]
        self.validation.require_equal(
            "step_03_output_manifest_sha256",
            manifest_sha,
            expected_manifest_sha,
        )
        entries = {}
        with open(self.manifest_path, "r", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if reader.fieldnames != ["relative_path", "bytes", "sha256"]:
                raise Step04Error("Unexpected Step 03 output manifest schema")
            for row in reader:
                if row["relative_path"] in entries:
                    raise Step04Error(
                        "Duplicate Step 03 manifest path: {}".format(
                            row["relative_path"]
                        )
                    )
                entries[row["relative_path"]] = row
        expected = self.configuration.values["classification_table"]
        manifest_row = entries.get(expected["relative_path"])
        if manifest_row is None:
            raise Step04Error("Step 03 manifest does not list the input table")
        self.validation.require_equal(
            "step_03_table_manifest_bytes",
            manifest_row["bytes"],
            expected["bytes"],
        )
        self.validation.require_equal(
            "step_03_table_manifest_sha256",
            manifest_row["sha256"],
            expected["sha256"],
        )
        if not os.path.isfile(self.table_path):
            raise Step04Error(
                "Step 03 classification table is missing: {}".format(
                    self.table_path
                )
            )
        self.validation.require_equal(
            "step_03_table_observed_bytes",
            os.path.getsize(self.table_path),
            expected["bytes"],
        )
        self.validation.require_equal(
            "step_03_table_observed_sha256",
            sha256_file(self.table_path),
            expected["sha256"],
        )

    def iter_observations(self):
        """Yield only encoded A/B/C state and registered WT/KO target."""
        class_counts = {"WT": 0, "KO": 0}
        seen_keys = set()
        with open(self.table_path, "r", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            required_header = self.configuration.values["required_input_columns"]
            if reader.fieldnames != required_header:
                raise Step04Error("Unexpected Step 03 classification-table schema")
            for row in reader:
                ground_truth = row["ground_truth_class"]
                if ground_truth not in class_counts:
                    raise Step04Error(
                        "Step 04 input contains a non-WT/KO class: {}".format(
                            ground_truth
                        )
                    )
                expected_class = self.configuration.values["ground_truth_classes"][
                    ground_truth
                ]
                if row["genotype"] != expected_class["genotype"]:
                    raise Step04Error("Genotype/ground-truth mismatch in Step 03")
                if row["sex"] != expected_class["sex"] or row["sex"] != "M":
                    raise Step04Error("Step 04 input must contain male cells only")
                if row["ground_truth_label"] != str(expected_class["label"]):
                    raise Step04Error("Label/ground-truth mismatch in Step 03")
                if row["classification_role"] != "ground_truth_training_candidate":
                    raise Step04Error("Unexpected Step 03 classification role")
                key = (row["technical_sample_id"], row["cell_barcode"])
                if key in seen_keys:
                    raise Step04Error(
                        "Duplicate Step 03 sample/barcode key: {} {}".format(*key)
                    )
                seen_keys.add(key)
                pattern_code = self.encoder.encode(
                    row["A_detected"], row["B_detected"], row["C_detected"]
                )
                if row["detection_pattern"] != self.encoder.label(pattern_code):
                    raise Step04Error(
                        "Step 03 binary state and pattern label disagree for {} {}".format(
                            *key
                        )
                    )
                class_counts[ground_truth] += 1
                yield LabeledPatternObservation(pattern_code, ground_truth)

        expected_classes = self.configuration.values["ground_truth_classes"]
        for class_name in ("WT", "KO"):
            self.validation.require_equal(
                "{}_input_cells".format(class_name.lower()),
                class_counts[class_name],
                expected_classes[class_name]["cells"],
            )
        self.validation.require_equal(
            "step_04_input_cells",
            sum(class_counts.values()),
            self.configuration.values["expected_total_cells"],
        )
        self.validation.require_equal(
            "unique_step_03_sample_barcode_keys",
            len(seen_keys),
            sum(class_counts.values()),
        )


class ProbabilisticClassifier(object):
    """Minimal probability interface shared by present and future models."""

    @property
    def model_name(self):
        raise NotImplementedError

    def predict_proba(self, features):
        """Return a mapping containing WT and KO probabilities."""
        raise NotImplementedError


class EmpiricalPatternEstimator(object):
    """Compute eight pattern-by-ground-truth counts from fit observations."""

    def __init__(self, encoder):
        self.encoder = encoder

    def fit(self, observations):
        counts = {
            pattern: {"WT": 0, "KO": 0}
            for pattern in self.encoder.pattern_order
        }
        class_totals = {"WT": 0, "KO": 0}
        for observation in observations:
            if observation.pattern_code not in counts:
                raise Step04Error(
                    "Estimator received an unknown pattern: {}".format(
                        observation.pattern_code
                    )
                )
            if observation.ground_truth_class not in class_totals:
                raise Step04Error(
                    "Estimator received an unknown class: {}".format(
                        observation.ground_truth_class
                    )
                )
            counts[observation.pattern_code][observation.ground_truth_class] += 1
            class_totals[observation.ground_truth_class] += 1
        return EmpiricalPatternClassifier(self.encoder, counts, class_totals)


class EmpiricalPatternClassifier(ProbabilisticClassifier):
    """Reusable probability lookup model estimated from pattern counts."""

    def __init__(self, encoder, counts, class_totals):
        self.encoder = encoder
        self.counts = counts
        self.class_totals = class_totals
        if min(class_totals.values()) <= 0:
            raise Step04Error("Both WT and KO ground-truth classes are required")

    @property
    def model_name(self):
        return "empirical_pcdh19_probe_pattern_classifier"

    @staticmethod
    def _ratio(numerator, denominator):
        if denominator == 0:
            if numerator == 0:
                return None
            return float("inf")
        return numerator / denominator

    def evidence(self, pattern_code):
        """Return empirical probability evidence without making a hard call."""
        if pattern_code not in self.counts:
            raise Step04Error("Unknown model pattern: {}".format(pattern_code))
        wt_cells = self.counts[pattern_code]["WT"]
        ko_cells = self.counts[pattern_code]["KO"]
        total_cells = wt_cells + ko_cells
        if total_cells == 0:
            p_wt = None
            p_ko = None
        else:
            p_wt = wt_cells / total_cells
            p_ko = ko_cells / total_cells
        p_pattern_wt = wt_cells / self.class_totals["WT"]
        p_pattern_ko = ko_cells / self.class_totals["KO"]
        return {
            "pattern_code": pattern_code,
            "pattern_label": self.encoder.label(pattern_code),
            "wt_cells": wt_cells,
            "ko_cells": ko_cells,
            "total_cells": total_cells,
            "p_wt_given_pattern": p_wt,
            "p_ko_given_pattern": p_ko,
            "wt_to_ko_cell_ratio": self._ratio(wt_cells, ko_cells),
            "ko_to_wt_cell_ratio": self._ratio(ko_cells, wt_cells),
            "p_pattern_given_wt": p_pattern_wt,
            "p_pattern_given_ko": p_pattern_ko,
            "wt_to_ko_likelihood_ratio": self._ratio(
                p_pattern_wt, p_pattern_ko
            ),
            "ko_to_wt_likelihood_ratio": self._ratio(
                p_pattern_ko, p_pattern_wt
            ),
            "hard_call": None,
        }

    def predict_proba(self, features):
        pattern_code = features.get("pattern_code")
        evidence = self.evidence(pattern_code)
        return {
            "WT": evidence["p_wt_given_pattern"],
            "KO": evidence["p_ko_given_pattern"],
        }

    def evidence_rows(self):
        return [self.evidence(pattern) for pattern in self.encoder.pattern_order]


class ClassificationEvaluator(object):
    """Perform Step 04 descriptive checks without hard-call evaluation."""

    def __init__(self, validation):
        self.validation = validation

    def evaluate(self, classifier):
        evidence_rows = classifier.evidence_rows()
        total_cells = sum(row["total_cells"] for row in evidence_rows)
        self.validation.require_equal(
            "empirical_model_patterns", len(evidence_rows), 8
        )
        self.validation.require_equal(
            "empirical_model_cells",
            total_cells,
            sum(classifier.class_totals.values()),
        )
        self.validation.require_equal(
            "hard_genotype_calls_created",
            sum(row["hard_call"] is not None for row in evidence_rows),
            0,
            "Step 04 estimates evidence only",
        )
        self.validation.require_equal(
            "pattern_000_preserved",
            sum(row["pattern_code"] == "000" for row in evidence_rows),
            1,
            "No-probe-detection state remains in the model",
        )
        for row in evidence_rows:
            pattern = row["pattern_code"]
            if row["total_cells"] > 0:
                self.validation.require_close(
                    "{}_posterior_probability_sum".format(pattern),
                    row["p_wt_given_pattern"] + row["p_ko_given_pattern"],
                    1.0,
                )
            predicted = classifier.predict_proba({"pattern_code": pattern})
            self.validation.require_equal(
                "{}_probability_interface_matches_model".format(pattern),
                predicted,
                {
                    "WT": row["p_wt_given_pattern"],
                    "KO": row["p_ko_given_pattern"],
                },
            )

        total_fraction = sum(
            row["total_cells"] / total_cells for row in evidence_rows
        )
        self.validation.require_close(
            "pattern_distribution_fraction_sum", total_fraction, 1.0
        )
        return evidence_rows

    @staticmethod
    def distribution_rows(evidence_rows, class_totals):
        total_cells = sum(row["total_cells"] for row in evidence_rows)
        rows = []
        for evidence in evidence_rows:
            total_fraction = evidence["total_cells"] / total_cells
            wt_fraction = evidence["wt_cells"] / class_totals["WT"]
            ko_fraction = evidence["ko_cells"] / class_totals["KO"]
            rows.append(
                {
                    "pattern_code": evidence["pattern_code"],
                    "pattern_label": evidence["pattern_label"],
                    "total_cells": evidence["total_cells"],
                    "fraction_of_all_cells": total_fraction,
                    "percent_of_all_cells": 100.0 * total_fraction,
                    "wt_cells": evidence["wt_cells"],
                    "fraction_of_wt_cells": wt_fraction,
                    "percent_of_wt_cells": 100.0 * wt_fraction,
                    "ko_cells": evidence["ko_cells"],
                    "fraction_of_ko_cells": ko_fraction,
                    "percent_of_ko_cells": 100.0 * ko_fraction,
                }
            )
        return rows


class ClassificationPlotter(object):
    """Render Step 04 diagnostics from computed model rows only."""

    WT_COLOR = "#3B6FB6"
    KO_COLOR = "#D95F43"
    TOTAL_COLOR = "#667788"

    def __init__(self, dpi):
        self.dpi = dpi
        plt.rcParams.update(
            {
                "font.size": 9,
                "axes.spines.top": False,
                "axes.spines.right": False,
                "figure.facecolor": "white",
                "axes.facecolor": "white",
            }
        )

    @staticmethod
    def _axis_labels(rows):
        return [
            "{}\n{}".format(row["pattern_code"], row["pattern_label"])
            for row in rows
        ]

    def plot_wt_vs_ko_frequency(self, rows, path):
        x_values = list(range(len(rows)))
        width = 0.38
        figure, axis = plt.subplots(figsize=(10, 5.5))
        axis.bar(
            [value - width / 2 for value in x_values],
            [100.0 * row["p_pattern_given_wt"] for row in rows],
            width,
            color=self.WT_COLOR,
            label="WT male",
        )
        axis.bar(
            [value + width / 2 for value in x_values],
            [100.0 * row["p_pattern_given_ko"] for row in rows],
            width,
            color=self.KO_COLOR,
            label="KO male",
        )
        axis.set_xticks(x_values)
        axis.set_xticklabels(self._axis_labels(rows))
        axis.set_ylabel("Cells within genotype (%)")
        axis.set_xlabel("Binary Probe A/B/C pattern")
        axis.set_title("Step 04: WT versus KO frequency of each PCDH19 probe pattern")
        axis.legend(frameon=False)
        axis.grid(axis="y", color="#DDDDDD", linewidth=0.7)
        figure.tight_layout()
        figure.savefig(
            path,
            dpi=self.dpi,
            metadata={"Title": "Step 04 WT versus KO pattern frequency"},
        )
        plt.close(figure)

    def plot_conditional_probabilities(self, rows, path):
        x_values = list(range(len(rows)))
        width = 0.38
        figure, axis = plt.subplots(figsize=(10, 5.5))
        axis.bar(
            [value - width / 2 for value in x_values],
            [100.0 * row["p_wt_given_pattern"] for row in rows],
            width,
            color=self.WT_COLOR,
            label="P(WT | pattern)",
        )
        axis.bar(
            [value + width / 2 for value in x_values],
            [100.0 * row["p_ko_given_pattern"] for row in rows],
            width,
            color=self.KO_COLOR,
            label="P(KO | pattern)",
        )
        axis.set_xticks(x_values)
        axis.set_xticklabels(self._axis_labels(rows))
        axis.set_ylim(0, 105)
        axis.set_ylabel("Empirical conditional probability (%)")
        axis.set_xlabel("Binary Probe A/B/C pattern")
        axis.set_title("Step 04: Empirical genotype probability given probe pattern")
        axis.legend(frameon=False)
        axis.grid(axis="y", color="#DDDDDD", linewidth=0.7)
        figure.tight_layout()
        figure.savefig(
            path,
            dpi=self.dpi,
            metadata={"Title": "Step 04 conditional genotype probabilities"},
        )
        plt.close(figure)

    def plot_pattern_abundance(self, rows, path):
        x_values = list(range(len(rows)))
        counts = [row["total_cells"] for row in rows]
        figure, axis = plt.subplots(figsize=(10, 5.5))
        bars = axis.bar(x_values, counts, color=self.TOTAL_COLOR, width=0.66)
        axis.set_yscale("log")
        axis.set_xticks(x_values)
        axis.set_xticklabels(self._axis_labels(rows))
        axis.set_ylabel("Ground-truth cells (log scale)")
        axis.set_xlabel("Binary Probe A/B/C pattern")
        axis.set_title("Step 04: Cell count and proportion in each probe pattern")
        axis.grid(axis="y", color="#DDDDDD", linewidth=0.7, which="both")
        for bar, row in zip(bars, rows):
            axis.annotate(
                "{:,}\n{:.2f}%".format(
                    row["total_cells"], 100.0 * row["fraction_of_all_cells"]
                ),
                xy=(bar.get_x() + bar.get_width() / 2, bar.get_height()),
                xytext=(0, 4),
                textcoords="offset points",
                ha="center",
                va="bottom",
                fontsize=7.5,
            )
        figure.tight_layout()
        figure.savefig(
            path,
            dpi=self.dpi,
            metadata={"Title": "Step 04 probe-pattern cell abundance"},
        )
        plt.close(figure)

    def render_all(self, evidence_rows, distribution_rows, output_directory):
        paths = [
            os.path.join(output_directory, FREQUENCY_PLOT_NAME),
            os.path.join(output_directory, POSTERIOR_PLOT_NAME),
            os.path.join(output_directory, ABUNDANCE_PLOT_NAME),
        ]
        self.plot_wt_vs_ko_frequency(evidence_rows, paths[0])
        self.plot_conditional_probabilities(evidence_rows, paths[1])
        self.plot_pattern_abundance(distribution_rows, paths[2])
        return paths


class Step04OutputPublisher(object):
    """Serialize, manifest, verify, and atomically publish Step 04 outputs."""

    EXPECTED_MANIFEST_FILES = {
        MODEL_TABLE_NAME,
        DISTRIBUTION_TABLE_NAME,
        VALIDATION_NAME,
        ENVIRONMENT_NAME,
        FREQUENCY_PLOT_NAME,
        POSTERIOR_PLOT_NAME,
        ABUNDANCE_PLOT_NAME,
    }

    def __init__(
        self,
        configuration,
        validation,
        plotter,
        output_root,
        script_path,
        requirements_path,
        step03_manifest_path,
        step03_table_path,
    ):
        self.configuration = configuration
        self.validation = validation
        self.plotter = plotter
        self.output_root = os.path.abspath(output_root)
        self.script_path = os.path.abspath(script_path)
        self.requirements_path = os.path.abspath(requirements_path)
        self.step03_manifest_path = os.path.abspath(step03_manifest_path)
        self.step03_table_path = os.path.abspath(step03_table_path)

    def _format_float(self, value):
        if value is None:
            return ""
        if math.isinf(value):
            return "Inf"
        return ("{:.%df}" % self.configuration.decimal_places).format(value)

    def _serialized_model_rows(self, evidence_rows):
        rows = []
        float_columns = {
            "p_wt_given_pattern",
            "p_ko_given_pattern",
            "wt_to_ko_cell_ratio",
            "ko_to_wt_cell_ratio",
            "p_pattern_given_wt",
            "p_pattern_given_ko",
            "wt_to_ko_likelihood_ratio",
            "ko_to_wt_likelihood_ratio",
        }
        for evidence in evidence_rows:
            row = dict(evidence)
            a_detected, b_detected, c_detected = tuple(
                int(value) for value in evidence["pattern_code"]
            )
            row.update(
                {
                    "A_detected": a_detected,
                    "B_detected": b_detected,
                    "C_detected": c_detected,
                    "hard_call": "",
                }
            )
            for key in float_columns:
                row[key] = self._format_float(row[key])
            rows.append(row)
        return rows

    def _serialized_distribution_rows(self, distribution_rows):
        float_columns = {
            "fraction_of_all_cells",
            "percent_of_all_cells",
            "fraction_of_wt_cells",
            "percent_of_wt_cells",
            "fraction_of_ko_cells",
            "percent_of_ko_cells",
        }
        rows = []
        for distribution in distribution_rows:
            row = dict(distribution)
            for key in float_columns:
                row[key] = self._format_float(row[key])
            rows.append(row)
        return rows

    def _write_environment(self, path):
        rows = [
            {"key": "step_id", "value": self.configuration.step_id},
            {
                "key": "pipeline_version",
                "value": self.configuration.values["pipeline_version"],
            },
            {"key": "python_version", "value": platform.python_version()},
            {"key": "platform", "value": platform.platform()},
            {"key": "numpy_version", "value": np.__version__},
            {"key": "matplotlib_version", "value": matplotlib.__version__},
            {"key": "python_script", "value": self.script_path},
            {"key": "python_script_sha256", "value": sha256_file(self.script_path)},
            {"key": "lock_path", "value": self.configuration.lock_path},
            {
                "key": "lock_sha256",
                "value": sha256_file(self.configuration.lock_path),
            },
            {"key": "requirements_path", "value": self.requirements_path},
            {
                "key": "requirements_sha256",
                "value": sha256_file(self.requirements_path),
            },
            {
                "key": "upstream_step",
                "value": self.configuration.values["upstream_step_id"],
            },
            {"key": "step_03_manifest_path", "value": self.step03_manifest_path},
            {
                "key": "step_03_manifest_sha256",
                "value": sha256_file(self.step03_manifest_path),
            },
            {"key": "step_03_table_path", "value": self.step03_table_path},
            {
                "key": "step_03_table_sha256",
                "value": sha256_file(self.step03_table_path),
            },
            {"key": "predictors", "value": "A_detected,B_detected,C_detected"},
            {"key": "smoothing", "value": "none"},
            {"key": "hard_calls", "value": "none"},
            {"key": "het_cells", "value": "not_loaded"},
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

    @classmethod
    def _verify_manifest(cls, directory):
        manifest_path = os.path.join(directory, MANIFEST_NAME)
        if not os.path.isfile(manifest_path):
            raise Step04Error("Step 04 output manifest is missing")
        listed = set()
        with open(manifest_path, "r", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if reader.fieldnames != ["relative_path", "bytes", "sha256"]:
                raise Step04Error("Unexpected Step 04 output manifest schema")
            for row in reader:
                if row["relative_path"] in listed:
                    raise Step04Error(
                        "Duplicate Step 04 manifest path: {}".format(
                            row["relative_path"]
                        )
                    )
                listed.add(row["relative_path"])
                path = os.path.join(directory, row["relative_path"])
                if not os.path.isfile(path):
                    raise Step04Error(
                        "Manifest-listed Step 04 file is missing: {}".format(path)
                    )
                if os.path.getsize(path) != int(row["bytes"]):
                    raise Step04Error(
                        "Step 04 output size mismatch: {}".format(path)
                    )
                if sha256_file(path) != row["sha256"]:
                    raise Step04Error(
                        "Step 04 output checksum mismatch: {}".format(path)
                    )
        if listed != cls.EXPECTED_MANIFEST_FILES:
            raise Step04Error("Step 04 manifest file set is incomplete or unexpected")

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
            "requirements_sha256": sha256_file(self.requirements_path),
            "step_03_manifest_sha256": sha256_file(self.step03_manifest_path),
            "step_03_table_sha256": sha256_file(self.step03_table_path),
        }
        for key, value in expected.items():
            if environment.get(key) != value:
                raise Step04Error(
                    "Existing Step 04 provenance differs for {}".format(key)
                )

    def publish(self, evidence_rows, distribution_rows):
        if os.path.exists(self.output_root):
            self._verify_existing_provenance()
            print("EXISTING_VALIDATED\t{}".format(self.output_root))
            return
        parent = os.path.dirname(self.output_root)
        os.makedirs(parent, exist_ok=True)
        staging = tempfile.mkdtemp(
            prefix=".step_04_pcdh19_empirical_pattern_classifier.", dir=parent
        )
        try:
            write_tsv(
                os.path.join(staging, MODEL_TABLE_NAME),
                MODEL_HEADER,
                self._serialized_model_rows(evidence_rows),
            )
            write_tsv(
                os.path.join(staging, DISTRIBUTION_TABLE_NAME),
                DISTRIBUTION_HEADER,
                self._serialized_distribution_rows(distribution_rows),
            )
            plot_paths = self.plotter.render_all(
                evidence_rows, distribution_rows, staging
            )
            self.validation.require_equal(
                "diagnostic_plots_created",
                sum(os.path.getsize(path) > 0 for path in plot_paths),
                3,
            )
            write_tsv(
                os.path.join(staging, VALIDATION_NAME),
                VALIDATION_HEADER,
                self.validation.rows,
            )
            self._write_environment(os.path.join(staging, ENVIRONMENT_NAME))
            self._write_manifest(staging)
            self._verify_manifest(staging)
            if os.path.exists(self.output_root):
                raise Step04Error(
                    "Step 04 output appeared during publication: {}".format(
                        self.output_root
                    )
                )
            os.replace(staging, self.output_root)
            staging = None
            print(
                "PUBLISHED\t{}\tpatterns=8\tcells={}".format(
                    self.output_root,
                    sum(row["total_cells"] for row in evidence_rows),
                )
            )
        finally:
            if staging and os.path.exists(staging):
                shutil.rmtree(staging)


class PCDH19EmpiricalPatternClassifierStep(object):
    """Orchestrate the narrow Step 04 fit, diagnostics, and publication."""

    def __init__(
        self,
        lock_path,
        requirements_path,
        step03_root,
        paper3_root,
        script_path,
    ):
        self.configuration = Step04Configuration(lock_path)
        self.validation = ValidationRecorder(self.configuration.step_id)
        self.encoder = ProbePatternEncoder(self.configuration.pattern_order)
        self.reader = Step03ClassificationTableReader(
            self.configuration,
            self.encoder,
            self.validation,
            step03_root,
        )
        self.estimator = EmpiricalPatternEstimator(self.encoder)
        self.evaluator = ClassificationEvaluator(self.validation)
        self.plotter = ClassificationPlotter(self.configuration.plot_dpi)
        output_root = os.path.join(
            os.path.abspath(paper3_root), "results", STEP_OUTPUT_DIRECTORY
        )
        self.publisher = Step04OutputPublisher(
            self.configuration,
            self.validation,
            self.plotter,
            output_root,
            script_path,
            requirements_path,
            self.reader.manifest_path,
            self.reader.table_path,
        )

    def run(self):
        if os.path.exists(self.publisher.output_root):
            self.publisher.publish([], [])
            return
        classifier = self.estimator.fit(self.reader.iter_observations())
        evidence_rows = self.evaluator.evaluate(classifier)
        distribution_rows = self.evaluator.distribution_rows(
            evidence_rows, classifier.class_totals
        )
        self.publisher.publish(evidence_rows, distribution_rows)


def build_parser():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--lock", required=True, help="Step 04 JSON lock")
    parser.add_argument(
        "--requirements", required=True, help="Pinned Step 04 requirements"
    )
    parser.add_argument(
        "--step03-root",
        required=True,
        help="Validated Step 03 result directory",
    )
    parser.add_argument(
        "--paper3-root",
        required=True,
        help="Paper 3 generated-output root containing results/ and logs/",
    )
    return parser


def main(argv=None):
    args = build_parser().parse_args(argv)
    step = PCDH19EmpiricalPatternClassifierStep(
        lock_path=args.lock,
        requirements_path=args.requirements,
        step03_root=args.step03_root,
        paper3_root=args.paper3_root,
        script_path=__file__,
    )
    step.run()


if __name__ == "__main__":
    try:
        main()
    except Step04Error as error:
        sys.stderr.write("ERROR: {}\n".format(error))
        sys.exit(1)
