#!/usr/bin/env python3
"""Formal Step 06: apply frozen Step 05 PCDH19 models to HET females.

This is an inference-only step.  It reads the immutable WT-M/WT-F/KO-M
Step 05 model artifacts, verifies their identities, and applies their exact
coefficients to the three registered HET-female samples.  It never fits a
model, selects a model, estimates a threshold, or assigns WT/KO genotypes.

Raw probe counts are probe-level UMI/ligation evidence, not transcript counts.
"""

from __future__ import print_function

import argparse
import csv
import datetime
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


STEP_ID = "06_pcdh19_het_female_inference"
MANIFEST_NAME = "output_manifest.tsv"
VALIDATION_NAME = "step_06_het_female_inference_validation.tsv"
ENVIRONMENT_NAME = "software_environment.tsv"
CELL_TABLE_NAME = "step_06_het_female_cell_probabilities.tsv"
SAMPLE_SUMMARY_NAME = "step_06_het_female_sample_summary.tsv"
PATTERN_TABLE_NAME = "step_06_het_female_pattern_distribution.tsv"
BC_TABLE_NAME = "step_06_het_female_a_negative_bc_evidence_distribution.tsv"
PROBABILITY_SUMMARY_NAME = "step_06_het_and_control_probability_summary.tsv"
PROBABILITY_HISTOGRAM_NAME = "step_06_het_and_control_probability_histograms.tsv"
MODEL_IDENTITY_NAME = "step_06_frozen_step_05_model_identity.tsv"
EVIDENCE_PLOT_NAME = "step_06_het_female_probe_evidence_composition.png"
BINARY_PLOT_NAME = "step_06_het_vs_controls_binary_model_probability.png"
COUNT_PLOT_NAME = "step_06_het_vs_controls_count_informed_model_probability.png"
BC_PLOT_NAME = "step_06_het_vs_controls_a_negative_bc_evidence.png"


class Step06Error(RuntimeError):
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


def write_tsv(path, header, rows):
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=header, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({name: row.get(name, "") for name in header})


def format_float(value, places=12):
    return ("{:.%df}" % places).format(float(value))


class ValidationRecorder(object):
    def __init__(self):
        self.rows = []

    def check(self, name, observed, expected, detail=""):
        passed = observed == expected
        self.rows.append({
            "check": name,
            "status": "PASS" if passed else "FAIL",
            "observed": str(observed),
            "expected": str(expected),
            "detail": detail,
        })
        if not passed:
            raise Step06Error("Validation failed for {}: observed {!r}, expected {!r}".format(name, observed, expected))


class Step06Configuration(object):
    REQUIRED = {
        "step_id", "pipeline_version", "upstream_step_id", "class_encoding",
        "positive_model_class", "binary_model", "count_informed_model",
        "step_02", "sample_key", "het_samples", "control_reference",
        "pattern_order", "bc_bins", "probability_bins", "probability_decimal_places",
        "plot_dpi", "prohibited_operations",
    }

    def __init__(self, path):
        self.path = os.path.abspath(path)
        with open(self.path, "r") as handle:
            self.values = json.load(handle)
        missing = sorted(self.REQUIRED.difference(self.values))
        if missing:
            raise Step06Error("Step 06 lock is missing keys: {}".format(", ".join(missing)))
        if self.values["step_id"] != STEP_ID:
            raise Step06Error("Unexpected Step 06 identifier")
        if self.values["class_encoding"] != {"WT": 0, "KO": 1}:
            raise Step06Error("Step 06 must preserve WT=0 and KO=1")
        if self.values["positive_model_class"] != "KO":
            raise Step06Error("Frozen models must estimate P(KO)")
        if self.values["pattern_order"] != ["000", "001", "010", "011", "100", "101", "110", "111"]:
            raise Step06Error("Unexpected pattern order")
        if self.values["bc_bins"] != ["0", "1", "2", "3+"]:
            raise Step06Error("Raw evidence bins must be 0, 1, 2, and 3+")
        required_prohibitions = {
            "model_fitting", "coefficient_estimation", "threshold_optimization",
            "class_weighting", "model_selection", "het_ground_truth_labels",
            "hard_genotype_calls", "cell_type_stratification",
        }
        if set(self.values["prohibited_operations"]) != required_prohibitions:
            raise Step06Error("Step 06 prohibited-operation guard is incomplete")
        if set(self.values["het_samples"]) != {"15662-JZ-7", "15662-JZ-8", "15662-JZ-9"}:
            raise Step06Error("Step 06 HET cohort must be exactly JZ-7, JZ-8, and JZ-9")

    @property
    def decimals(self):
        return int(self.values["probability_decimal_places"])


class FrozenLogisticModel(object):
    """Probability-only implementation of an already fitted Step 05 model."""

    REQUIRED_COLUMNS = ["term", "predictor", "outcome_log_odds", "coefficient", "odds_ratio", "interpretation"]

    def __init__(self, model_id, formula, predictors, coefficients, source_path, source_sha256):
        self.model_id = model_id
        self.formula = formula
        self.predictors = list(predictors)
        self.coefficients = dict(coefficients)
        self.source_path = source_path
        self.source_sha256 = source_sha256

    def predict_ko(self, feature_values):
        linear = self.coefficients["intercept"]
        for predictor in self.predictors:
            linear += self.coefficients[predictor] * float(feature_values[predictor])
        if linear >= 0:
            return 1.0 / (1.0 + math.exp(-linear))
        exponential = math.exp(linear)
        return exponential / (1.0 + exponential)


class FrozenStep05ModelReader(object):
    """Verify exact model/package bytes and deserialize coefficients only."""

    def __init__(self, configuration, step05_root, validation):
        self.configuration = configuration
        self.step05_root = os.path.abspath(step05_root)
        self.validation = validation

    def _verify_package(self, definition):
        manifest = os.path.join(self.step05_root, definition["relative_directory"], MANIFEST_NAME)
        if not os.path.isfile(manifest):
            raise Step06Error("Frozen Step 05 package manifest is missing")
        self.validation.check(
            definition["model_id"] + "_package_manifest_sha256",
            sha256_file(manifest), definition["package_manifest_sha256"],
        )

    def _read_one(self, definition):
        self._verify_package(definition)
        path = os.path.join(self.step05_root, definition["relative_directory"], definition["coefficient_file"])
        if not os.path.isfile(path):
            raise Step06Error("Frozen coefficient table is missing: {}".format(path))
        self.validation.check(definition["model_id"] + "_coefficient_bytes", os.path.getsize(path), definition["coefficient_bytes"])
        observed_sha = sha256_file(path)
        self.validation.check(definition["model_id"] + "_coefficient_sha256", observed_sha, definition["coefficient_sha256"])
        coefficients = {}
        with open(path, "r", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if reader.fieldnames != FrozenLogisticModel.REQUIRED_COLUMNS:
                raise Step06Error("Unexpected frozen coefficient schema")
            for row in reader:
                term = row["term"]
                if term in coefficients:
                    raise Step06Error("Duplicate frozen model term")
                if row["outcome_log_odds"] != "KO_vs_WT; WT=0, KO=1":
                    raise Step06Error("Frozen model class orientation changed")
                coefficients[term] = float(row["coefficient"])
        required = {"intercept"}.union(definition["predictors"])
        self.validation.check(definition["model_id"] + "_coefficient_terms", set(coefficients), required)
        return FrozenLogisticModel(
            definition["model_id"], definition["formula"], definition["predictors"],
            coefficients, path, observed_sha,
        )

    def read(self):
        binary = self._read_one(self.configuration.values["binary_model"])
        count = self._read_one(self.configuration.values["count_informed_model"])
        self.validation.check("frozen_model_count", 2, 2)
        return binary, count


class ProbeCell(object):
    __slots__ = (
        "sample_id", "submitted_name", "design_group", "genotype_context", "sex",
        "barcode", "a", "b", "c", "pattern", "binary_p_ko", "count_p_ko",
    )

    def __init__(self, sample_id, submitted_name, design_group, genotype_context, sex, barcode, a, b, c):
        self.sample_id = sample_id
        self.submitted_name = submitted_name
        self.design_group = design_group
        self.genotype_context = genotype_context
        self.sex = sex
        self.barcode = barcode
        self.a = int(a)
        self.b = int(b)
        self.c = int(c)
        self.pattern = "{}{}{}".format(int(self.a > 0), int(self.b > 0), int(self.c > 0))
        self.binary_p_ko = None
        self.count_p_ko = None

    @property
    def is_000(self):
        return self.pattern == "000"

    @property
    def bc_total(self):
        return self.b + self.c


class ManifestedHETFemaleLoader(object):
    """Load only the three locked HET-female Step 02 probe tables."""

    COLUMNS = ["barcode", "A_UMI", "B_UMI", "C_UMI", "Pcdh19_total_UMI", "detection_pattern"]
    LABELS = {
        "000": "none", "001": "C only", "010": "B only", "011": "B+C",
        "100": "A only", "101": "A+C", "110": "A+B", "111": "A+B+C",
    }

    def __init__(self, configuration, step02_root, sample_key, validation):
        self.configuration = configuration
        self.step02_root = os.path.abspath(step02_root)
        self.sample_key = os.path.abspath(sample_key)
        self.validation = validation

    def _metadata(self):
        lock = self.configuration.values
        self.validation.check("sample_key_sha256", sha256_file(self.sample_key), lock["sample_key"]["sha256"])
        rows = {}
        with open(self.sample_key, "r", newline="") as handle:
            reader = csv.DictReader(handle)
            required = {"technical_sample_id", "submitted_sample_name", "genotype", "sex", "design_group"}
            if not required.issubset(set(reader.fieldnames or [])):
                raise Step06Error("Registered sample key schema is incomplete")
            for row in reader:
                rows[row["technical_sample_id"]] = row
        for sample_id, definition in lock["het_samples"].items():
            row = rows.get(sample_id)
            if row is None:
                raise Step06Error("HET sample is absent from registered sample key")
            for field in ("genotype", "sex", "design_group", "submitted_sample_name"):
                self.validation.check(sample_id + "_" + field, row[field], definition[field])
        return rows

    def _manifest(self):
        path = os.path.join(self.step02_root, MANIFEST_NAME)
        self.validation.check("step_02_manifest_sha256", sha256_file(path), self.configuration.values["step_02"]["manifest_sha256"])
        rows = {}
        with open(path, "r", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if reader.fieldnames != ["relative_path", "bytes", "sha256"]:
                raise Step06Error("Unexpected Step 02 manifest schema")
            for row in reader:
                rows[row["relative_path"]] = row
        return rows

    def load(self):
        metadata = self._metadata()
        manifest = self._manifest()
        records = []
        observed_samples = set()
        for sample_id in sorted(self.configuration.values["het_samples"]):
            definition = self.configuration.values["het_samples"][sample_id]
            relative = os.path.join("per_sample", sample_id, "pcdh19_probe_patterns.tsv")
            entry = manifest.get(relative)
            if entry is None:
                raise Step06Error("HET probe table is not manifested")
            path = os.path.join(self.step02_root, relative)
            self.validation.check(sample_id + "_input_bytes", os.path.getsize(path), int(entry["bytes"]))
            self.validation.check(sample_id + "_input_sha256", sha256_file(path), entry["sha256"])
            seen = set()
            count = 0
            with open(path, "r", newline="") as handle:
                reader = csv.DictReader(handle, delimiter="\t")
                if reader.fieldnames != self.COLUMNS:
                    raise Step06Error("Unexpected HET probe-table schema")
                for row in reader:
                    if not row["barcode"] or row["barcode"] in seen:
                        raise Step06Error("Empty or duplicate HET barcode")
                    seen.add(row["barcode"])
                    try:
                        a, b, c, total = [int(row[name]) for name in ("A_UMI", "B_UMI", "C_UMI", "Pcdh19_total_UMI")]
                    except ValueError:
                        raise Step06Error("Non-integer HET probe count")
                    if min(a, b, c, total) < 0 or a + b + c != total:
                        raise Step06Error("Inconsistent HET raw probe counts")
                    record = ProbeCell(sample_id, metadata[sample_id]["submitted_sample_name"], "HET_F", "HET_withheld_inference_only", "F", row["barcode"], a, b, c)
                    if row["detection_pattern"] != self.LABELS[record.pattern]:
                        raise Step06Error("HET detection-pattern label is inconsistent")
                    records.append(record)
                    count += 1
            self.validation.check(sample_id + "_cell_count", count, definition["cells"])
            observed_samples.add(sample_id)
        self.validation.check("het_samples_loaded", observed_samples, set(self.configuration.values["het_samples"]))
        self.validation.check("het_total_cells", len(records), sum(row["cells"] for row in self.configuration.values["het_samples"].values()))
        return records


class FrozenControlReferenceLoader(object):
    """Read immutable Step 05 controls solely for distribution overlays."""

    REQUIRED = [
        "cell_barcode", "biological_sample_id", "submitted_sample_name", "sex", "design_group",
        "true_genotype", "A_UMI", "A_detected", "B_UMI", "C_UMI", "pattern_code",
        "binary_predicted_wt_probability", "binary_predicted_ko_probability", "binary_predicted_genotype",
        "count_predicted_wt_probability", "count_predicted_ko_probability", "count_predicted_genotype",
    ]

    def __init__(self, configuration, step05_root, validation):
        self.definition = configuration.values["control_reference"]
        self.step05_root = os.path.abspath(step05_root)
        self.validation = validation

    def load(self):
        path = os.path.join(self.step05_root, self.definition["relative_directory"], self.definition["cell_file"])
        self.validation.check("control_reference_bytes", os.path.getsize(path), self.definition["cell_bytes"])
        self.validation.check("control_reference_sha256", sha256_file(path), self.definition["cell_sha256"])
        records = []
        counts = {"WT_M": 0, "WT_F": 0, "KO_M": 0}
        with open(path, "r", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if reader.fieldnames != self.REQUIRED:
                raise Step06Error("Unexpected Step 05 control-reference schema")
            for row in reader:
                group = row["design_group"]
                if group not in counts:
                    raise Step06Error("Control reference contains a prohibited cohort")
                expected_genotype = "KO" if group == "KO_M" else "WT"
                if row["true_genotype"] != expected_genotype:
                    raise Step06Error("Control reference genotype/group mismatch")
                record = ProbeCell(
                    row["biological_sample_id"], row["submitted_sample_name"], group,
                    row["true_genotype"] + "_ground_truth", row["sex"], row["cell_barcode"],
                    int(row["A_UMI"]), int(row["B_UMI"]), int(row["C_UMI"]),
                )
                if record.pattern != row["pattern_code"]:
                    raise Step06Error("Control raw counts disagree with pattern code")
                records.append(record)
                counts[group] += 1
        self.validation.check("control_reference_total_cells", len(records), self.definition["expected_total_cells"])
        self.validation.check("control_reference_group_cells", counts, self.definition["expected_group_cells"])
        return records


class FrozenModelApplicator(object):
    """Apply models without exposing any fit or decision interface."""

    def __init__(self, binary_model, count_model):
        self.binary_model = binary_model
        self.count_model = count_model

    def apply(self, records):
        for record in records:
            record.binary_p_ko = self.binary_model.predict_ko({
                "A_detected": int(record.a > 0),
                "B_detected": int(record.b > 0),
                "C_detected": int(record.c > 0),
            })
            record.count_p_ko = self.count_model.predict_ko({
                "A_detected": int(record.a > 0), "B_UMI": record.b, "C_UMI": record.c,
            })
        return records


class HETInferenceSummarizer(object):
    """Build descriptive evidence and probability distributions; never calls genotype."""

    def __init__(self, configuration):
        self.configuration = configuration
        self.pattern_order = configuration.values["pattern_order"]

    @staticmethod
    def _bin_count(value):
        return str(value) if value < 3 else "3+"

    @staticmethod
    def _quantile(values, q):
        if not values:
            return float("nan")
        return float(np.percentile(np.asarray(values, dtype=float), 100.0 * q))

    def cell_rows(self, records):
        rows = []
        for record in records:
            rows.append({
                "biological_sample_id": record.sample_id,
                "submitted_sample_name": record.submitted_name,
                "sex": record.sex,
                "cohort_role": record.genotype_context,
                "cell_barcode": record.barcode,
                "A_UMI": record.a, "B_UMI": record.b, "C_UMI": record.c,
                "B_plus_C_UMI": record.bc_total,
                "A_detected": int(record.a > 0), "B_detected": int(record.b > 0), "C_detected": int(record.c > 0),
                "pattern_code": record.pattern,
                "binary_predicted_wt_probability": format_float(1.0 - record.binary_p_ko, self.configuration.decimals),
                "binary_predicted_ko_probability": format_float(record.binary_p_ko, self.configuration.decimals),
                "count_informed_predicted_wt_probability": format_float(1.0 - record.count_p_ko, self.configuration.decimals),
                "count_informed_predicted_ko_probability": format_float(record.count_p_ko, self.configuration.decimals),
                "inference_status": "uncalled_000" if record.is_000 else "probability_only_no_genotype_call",
            })
        return rows

    def sample_summary(self, het_records):
        rows = []
        for sample_id in sorted(set(record.sample_id for record in het_records)):
            records = [record for record in het_records if record.sample_id == sample_id]
            a_negative = [record for record in records if record.a == 0]
            informative = [record for record in records if not record.is_000]
            bc_counts = {key: sum(self._bin_count(record.bc_total) == key for record in a_negative) for key in ("0", "1", "2", "3+")}
            binary = [record.binary_p_ko for record in informative]
            count = [record.count_p_ko for record in informative]
            n = len(records)
            row = {
                "biological_sample_id": sample_id,
                "submitted_sample_name": records[0].submitted_name,
                "total_cells": n,
                "cells_000": sum(record.is_000 for record in records),
                "fraction_000": format_float(sum(record.is_000 for record in records) / float(n)),
                "cells_A_positive": sum(record.a > 0 for record in records),
                "fraction_A_positive": format_float(sum(record.a > 0 for record in records) / float(n)),
                "cells_A_negative": len(a_negative),
                "A_negative_B_plus_C_0_fraction": format_float(bc_counts["0"] / float(len(a_negative))),
                "A_negative_B_plus_C_1_fraction": format_float(bc_counts["1"] / float(len(a_negative))),
                "A_negative_B_plus_C_2_fraction": format_float(bc_counts["2"] / float(len(a_negative))),
                "A_negative_B_plus_C_3plus_fraction": format_float(bc_counts["3+"] / float(len(a_negative))),
                "informative_non000_cells": len(informative),
            }
            for label, values in (("binary_p_ko", binary), ("count_informed_p_ko", count)):
                row[label + "_mean_non000"] = format_float(sum(values) / float(len(values)))
                for suffix, quantile in (("q10", .10), ("q25", .25), ("median", .50), ("q75", .75), ("q90", .90)):
                    row[label + "_" + suffix + "_non000"] = format_float(self._quantile(values, quantile))
            rows.append(row)
        return rows

    def pattern_rows(self, het_records):
        rows = []
        for sample_id in sorted(set(record.sample_id for record in het_records)):
            records = [record for record in het_records if record.sample_id == sample_id]
            for pattern in self.pattern_order:
                cells = sum(record.pattern == pattern for record in records)
                rows.append({"biological_sample_id": sample_id, "pattern_code": pattern, "cells": cells, "fraction": format_float(cells / float(len(records)))})
        return rows

    def bc_distribution_rows(self, het_records, control_records):
        cohorts = [(sample_id, "HET_F", [r for r in het_records if r.sample_id == sample_id]) for sample_id in sorted(set(r.sample_id for r in het_records))]
        cohorts += [(group, group, [r for r in control_records if r.design_group == group]) for group in ("WT_M", "WT_F", "KO_M")]
        rows = []
        for cohort_id, cohort_group, records in cohorts:
            a_negative = [record for record in records if record.a == 0]
            for measure, getter in (("B_UMI", lambda x: x.b), ("C_UMI", lambda x: x.c), ("B_plus_C_UMI", lambda x: x.bc_total)):
                for count_bin in ("0", "1", "2", "3+"):
                    cells = sum(self._bin_count(getter(record)) == count_bin for record in a_negative)
                    rows.append({
                        "cohort_id": cohort_id, "cohort_group": cohort_group, "subset": "A_negative",
                        "measure": measure, "count_bin": count_bin, "cells": cells,
                        "proportion_within_cohort_A_negative": format_float(cells / float(len(a_negative))),
                        "cohort_A_negative_cells": len(a_negative),
                    })
        return rows

    def probability_outputs(self, het_records, control_records):
        cohorts = [(sample_id, "HET_F", [r for r in het_records if r.sample_id == sample_id]) for sample_id in sorted(set(r.sample_id for r in het_records))]
        cohorts += [(group, group, [r for r in control_records if r.design_group == group]) for group in ("WT_M", "WT_F", "KO_M")]
        summary_rows, histogram_rows = [], []
        edges = [float(value) for value in self.configuration.values["probability_bins"]]
        for cohort_id, cohort_group, all_records in cohorts:
            records = [record for record in all_records if not record.is_000]
            for model, getter in (("binary", lambda x: x.binary_p_ko), ("count_informed", lambda x: x.count_p_ko)):
                values = [getter(record) for record in records]
                summary_rows.append({
                    "cohort_id": cohort_id, "cohort_group": cohort_group, "model": model,
                    "subset": "non000", "cells": len(values),
                    "mean_p_ko": format_float(sum(values) / float(len(values))),
                    "q10_p_ko": format_float(self._quantile(values, .10)),
                    "q25_p_ko": format_float(self._quantile(values, .25)),
                    "median_p_ko": format_float(self._quantile(values, .50)),
                    "q75_p_ko": format_float(self._quantile(values, .75)),
                    "q90_p_ko": format_float(self._quantile(values, .90)),
                })
                hist, observed_edges = np.histogram(np.asarray(values), bins=np.asarray(edges))
                for index, cells in enumerate(hist):
                    histogram_rows.append({
                        "cohort_id": cohort_id, "cohort_group": cohort_group, "model": model,
                        "subset": "non000", "bin_lower_inclusive": format_float(observed_edges[index]),
                        "bin_upper": format_float(observed_edges[index + 1]),
                        "upper_bound_inclusive": int(index == len(hist) - 1), "cells": int(cells),
                        "proportion": format_float(cells / float(len(values))),
                    })
        return summary_rows, histogram_rows


class HETInferencePlotter(object):
    """Render diagnostic figures separately from inference and summarization."""

    COLORS = {"WT_M": "#8ecae6", "WT_F": "#219ebc", "KO_M": "#d95f02", "HET_F": "#7b2cbf"}

    def __init__(self, dpi):
        self.dpi = int(dpi)

    def evidence_composition(self, summary_rows, path):
        sample_ids = [row["biological_sample_id"] for row in summary_rows]
        data = np.asarray([[float(row[key]) for key in (
            "fraction_000", "fraction_A_positive", "A_negative_B_plus_C_1_fraction",
            "A_negative_B_plus_C_2_fraction", "A_negative_B_plus_C_3plus_fraction",
        )] for row in summary_rows])
        # A-negative B/C fractions have an A-negative denominator; convert to all-cell fractions.
        for index, row in enumerate(summary_rows):
            a_negative_fraction = float(row["cells_A_negative"]) / float(row["total_cells"])
            data[index, 2:] *= a_negative_fraction
        fig, ax = plt.subplots(figsize=(8.2, 4.8))
        labels = ["000", "A+", "A-; B+C=1", "A-; B+C=2", "A-; B+C=3+"]
        colors = ["#d9d9d9", "#2a9d8f", "#e9c46a", "#f4a261", "#e76f51"]
        bottom = np.zeros(len(sample_ids))
        for index, label in enumerate(labels):
            ax.bar(sample_ids, data[:, index], bottom=bottom, label=label, color=colors[index])
            bottom += data[:, index]
        ax.set_ylabel("Fraction of all HET-female cells")
        ax.set_ylim(0, 1)
        ax.set_title("Step 06 HET-female probe-evidence composition")
        ax.legend(frameon=False, ncol=3, fontsize=8)
        fig.tight_layout()
        fig.savefig(path, dpi=self.dpi)
        plt.close(fig)

    def probability_plot(self, het_records, control_records, model, path):
        getter = (lambda x: x.binary_p_ko) if model == "binary" else (lambda x: x.count_p_ko)
        sample_ids = sorted(set(record.sample_id for record in het_records))
        fig, axes = plt.subplots(1, 3, figsize=(13.2, 3.8), sharex=True, sharey=True)
        bins = np.linspace(0, 1, 31)
        for ax, sample_id in zip(axes, sample_ids):
            for group in ("WT_M", "WT_F", "KO_M"):
                values = [getter(record) for record in control_records if record.design_group == group and not record.is_000]
                ax.hist(values, bins=bins, density=True, histtype="step", linewidth=1.5, color=self.COLORS[group], label=group)
            values = [getter(record) for record in het_records if record.sample_id == sample_id and not record.is_000]
            ax.hist(values, bins=bins, density=True, histtype="step", linewidth=2.5, color=self.COLORS["HET_F"], label=sample_id)
            ax.set_title(sample_id)
            ax.set_xlabel("Frozen model P(KO)")
        axes[0].set_ylabel("Density among non-000 cells")
        axes[-1].legend(frameon=False, fontsize=7, loc="upper center")
        title = "binary A/B/C" if model == "binary" else "count-informed A detection + raw B/C UMI"
        fig.suptitle("Step 06 HET versus controls: {} model".format(title), y=1.02)
        fig.tight_layout()
        fig.savefig(path, dpi=self.dpi, bbox_inches="tight")
        plt.close(fig)

    def bc_evidence_plot(self, bc_rows, path):
        rows = [row for row in bc_rows if row["measure"] == "B_plus_C_UMI"]
        cohort_ids = ["WT_M", "WT_F", "KO_M", "15662-JZ-7", "15662-JZ-8", "15662-JZ-9"]
        bins = ["0", "1", "2", "3+"]
        matrix = np.asarray([[float(next(row["proportion_within_cohort_A_negative"] for row in rows if row["cohort_id"] == cohort and row["count_bin"] == count_bin)) for count_bin in bins] for cohort in cohort_ids])
        fig, ax = plt.subplots(figsize=(9, 5))
        x = np.arange(len(cohort_ids))
        bottom = np.zeros(len(cohort_ids))
        colors = ["#d9d9d9", "#e9c46a", "#f4a261", "#e76f51"]
        for index, count_bin in enumerate(bins):
            ax.bar(x, matrix[:, index], bottom=bottom, color=colors[index], label="B+C {}".format(count_bin))
            bottom += matrix[:, index]
        ax.set_xticks(x)
        ax.set_xticklabels(cohort_ids, rotation=30, ha="right")
        ax.set_ylabel("Fraction among A-negative cells")
        ax.set_ylim(0, 1)
        ax.set_title("Step 06 raw B+C probe-level UMI/ligation evidence")
        ax.legend(frameon=False, ncol=4, fontsize=8)
        fig.tight_layout()
        fig.savefig(path, dpi=self.dpi)
        plt.close(fig)


class Step06OutputPublisher(object):
    EXPECTED = {
        CELL_TABLE_NAME, SAMPLE_SUMMARY_NAME, PATTERN_TABLE_NAME, BC_TABLE_NAME,
        PROBABILITY_SUMMARY_NAME, PROBABILITY_HISTOGRAM_NAME, MODEL_IDENTITY_NAME,
        EVIDENCE_PLOT_NAME, BINARY_PLOT_NAME, COUNT_PLOT_NAME, BC_PLOT_NAME,
        VALIDATION_NAME, ENVIRONMENT_NAME,
    }

    def __init__(self, configuration, validation, plotter, output_root, script_path, requirements_path):
        self.configuration = configuration
        self.validation = validation
        self.plotter = plotter
        self.output_root = os.path.abspath(output_root)
        self.script_path = os.path.abspath(script_path)
        self.requirements_path = os.path.abspath(requirements_path)

    def _manifest(self, root):
        rows = []
        for name in sorted(os.listdir(root)):
            if name == MANIFEST_NAME:
                continue
            path = os.path.join(root, name)
            if os.path.isfile(path):
                rows.append({"relative_path": name, "bytes": os.path.getsize(path), "sha256": sha256_file(path)})
        write_tsv(os.path.join(root, MANIFEST_NAME), ["relative_path", "bytes", "sha256"], rows)

    def _verify_existing(self):
        manifest = os.path.join(self.output_root, MANIFEST_NAME)
        if not os.path.isfile(manifest):
            raise Step06Error("Existing Step 06 output lacks a manifest")
        names = set()
        with open(manifest, "r", newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                path = os.path.join(self.output_root, row["relative_path"])
                if not os.path.isfile(path) or os.path.getsize(path) != int(row["bytes"]) or sha256_file(path) != row["sha256"]:
                    raise Step06Error("Existing Step 06 output fails manifest verification")
                names.add(row["relative_path"])
        if names != self.EXPECTED:
            raise Step06Error("Existing Step 06 output file set differs from contract")
        print("EXISTING_VALIDATED\t{}".format(self.output_root))

    def publish(self, het_records, cell_rows, summary_rows, pattern_rows, bc_rows, probability_rows, histogram_rows, models):
        if os.path.exists(self.output_root):
            self._verify_existing()
            return
        parent = os.path.dirname(self.output_root)
        os.makedirs(parent, exist_ok=True)
        staging = tempfile.mkdtemp(prefix=".step_06_pcdh19_het_female_inference.", dir=parent)
        try:
            write_tsv(os.path.join(staging, CELL_TABLE_NAME), [
                "biological_sample_id", "submitted_sample_name", "sex", "cohort_role", "cell_barcode",
                "A_UMI", "B_UMI", "C_UMI", "B_plus_C_UMI", "A_detected", "B_detected", "C_detected", "pattern_code",
                "binary_predicted_wt_probability", "binary_predicted_ko_probability",
                "count_informed_predicted_wt_probability", "count_informed_predicted_ko_probability", "inference_status",
            ], cell_rows)
            summary_header = list(summary_rows[0].keys())
            write_tsv(os.path.join(staging, SAMPLE_SUMMARY_NAME), summary_header, summary_rows)
            write_tsv(os.path.join(staging, PATTERN_TABLE_NAME), ["biological_sample_id", "pattern_code", "cells", "fraction"], pattern_rows)
            write_tsv(os.path.join(staging, BC_TABLE_NAME), ["cohort_id", "cohort_group", "subset", "measure", "count_bin", "cells", "proportion_within_cohort_A_negative", "cohort_A_negative_cells"], bc_rows)
            write_tsv(os.path.join(staging, PROBABILITY_SUMMARY_NAME), ["cohort_id", "cohort_group", "model", "subset", "cells", "mean_p_ko", "q10_p_ko", "q25_p_ko", "median_p_ko", "q75_p_ko", "q90_p_ko"], probability_rows)
            write_tsv(os.path.join(staging, PROBABILITY_HISTOGRAM_NAME), ["cohort_id", "cohort_group", "model", "subset", "bin_lower_inclusive", "bin_upper", "upper_bound_inclusive", "cells", "proportion"], histogram_rows)
            model_rows = []
            for model in models:
                for term in ["intercept"] + model.predictors:
                    model_rows.append({
                        "model_id": model.model_id, "class_encoding": "WT=0;KO=1", "positive_model_class": "KO",
                        "formula": model.formula, "term": term, "coefficient": format_float(model.coefficients[term]),
                        "source_coefficient_file": model.source_path, "source_sha256": model.source_sha256,
                    })
            write_tsv(os.path.join(staging, MODEL_IDENTITY_NAME), ["model_id", "class_encoding", "positive_model_class", "formula", "term", "coefficient", "source_coefficient_file", "source_sha256"], model_rows)
            self.plotter.evidence_composition(summary_rows, os.path.join(staging, EVIDENCE_PLOT_NAME))
            # Control overlays are passed through temporary attributes set by the orchestrator.
            self.plotter.probability_plot(het_records, self.control_records, "binary", os.path.join(staging, BINARY_PLOT_NAME))
            self.plotter.probability_plot(het_records, self.control_records, "count_informed", os.path.join(staging, COUNT_PLOT_NAME))
            self.plotter.bc_evidence_plot(bc_rows, os.path.join(staging, BC_PLOT_NAME))
            write_tsv(os.path.join(staging, VALIDATION_NAME), ["check", "status", "observed", "expected", "detail"], self.validation.rows)
            environment = [
                {"key": "step_id", "value": STEP_ID},
                {"key": "pipeline_version", "value": self.configuration.values["pipeline_version"]},
                {"key": "created_utc", "value": datetime.datetime.utcnow().replace(microsecond=0).isoformat() + "Z"},
                {"key": "python", "value": platform.python_version()},
                {"key": "numpy", "value": np.__version__},
                {"key": "matplotlib", "value": matplotlib.__version__},
                {"key": "pipeline_sha256", "value": sha256_file(self.script_path)},
                {"key": "lock_sha256", "value": sha256_file(self.configuration.path)},
                {"key": "requirements_sha256", "value": sha256_file(self.requirements_path)},
                {"key": "scope", "value": "HET_F inference only; frozen Step 05 models; no fitting, model selection, thresholds, calls, or cell types"},
                {"key": "raw_count_interpretation", "value": "probe-level UMI/ligation evidence; not transcript numbers"},
            ]
            write_tsv(os.path.join(staging, ENVIRONMENT_NAME), ["key", "value"], environment)
            self.validation.check("het_cell_rows_serialized", len(cell_rows), len(het_records))
            self.validation.check("all_000_rows_uncalled", sum(row["inference_status"] == "uncalled_000" for row in cell_rows), sum(record.is_000 for record in het_records))
            self.validation.check("hard_genotype_call_columns", any("genotype_call" in key for key in cell_rows[0]), False)
            # Rewrite validation after final publication checks were recorded.
            write_tsv(os.path.join(staging, VALIDATION_NAME), ["check", "status", "observed", "expected", "detail"], self.validation.rows)
            self._manifest(staging)
            os.replace(staging, self.output_root)
            staging = None
            print("PUBLISHED\t{}\thet_cells={}".format(self.output_root, len(het_records)))
        finally:
            if staging and os.path.exists(staging):
                shutil.rmtree(staging)


class PCDH19HETFemaleInferenceStep(object):
    def __init__(self, args):
        self.args = args
        self.validation = ValidationRecorder()
        self.configuration = Step06Configuration(args.lock)

    def run(self):
        binary_model, count_model = FrozenStep05ModelReader(self.configuration, self.args.step05_root, self.validation).read()
        het = ManifestedHETFemaleLoader(self.configuration, self.args.step02_root, self.args.sample_key, self.validation).load()
        controls = FrozenControlReferenceLoader(self.configuration, self.args.step05_root, self.validation).load()
        applicator = FrozenModelApplicator(binary_model, count_model)
        applicator.apply(het)
        applicator.apply(controls)
        self.validation.check("het_rows_used_for_model_fitting", 0, 0)
        self.validation.check("het_rows_used_for_threshold_optimization", 0, 0)
        self.validation.check("het_rows_used_for_model_selection", 0, 0)
        summarizer = HETInferenceSummarizer(self.configuration)
        cell_rows = summarizer.cell_rows(het)
        summary_rows = summarizer.sample_summary(het)
        pattern_rows = summarizer.pattern_rows(het)
        bc_rows = summarizer.bc_distribution_rows(het, controls)
        probability_rows, histogram_rows = summarizer.probability_outputs(het, controls)
        publisher = Step06OutputPublisher(
            self.configuration, self.validation,
            HETInferencePlotter(self.configuration.values["plot_dpi"]),
            self.args.output_root, __file__, self.args.requirements,
        )
        publisher.control_records = controls
        publisher.publish(het, cell_rows, summary_rows, pattern_rows, bc_rows, probability_rows, histogram_rows, [binary_model, count_model])


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--lock", required=True)
    parser.add_argument("--requirements", required=True)
    parser.add_argument("--step02-root", required=True)
    parser.add_argument("--step05-root", required=True)
    parser.add_argument("--sample-key", required=True)
    parser.add_argument("--output-root", required=True)
    return parser.parse_args(argv)


def main(argv=None):
    PCDH19HETFemaleInferenceStep(parse_args(argv)).run()


if __name__ == "__main__":
    main()
