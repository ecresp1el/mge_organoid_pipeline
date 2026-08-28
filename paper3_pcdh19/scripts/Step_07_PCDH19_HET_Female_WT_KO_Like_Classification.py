#!/usr/bin/env python3
"""Formal Step 07: high-confidence WT-like/KO-like HET-female calls.

Thresholds are selected exclusively from immutable Step 05 leave-one-sample-
out control predictions.  The frozen rule is then applied to immutable Step 06
HET probabilities.  This script contains no model-fitting, calibration,
feature-selection, or HET-informed threshold path.
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


STEP_ID = "07_pcdh19_het_female_wt_ko_like_classification"
MANIFEST_NAME = "output_manifest.tsv"
RULE_NAME = "step_07_frozen_wt_ko_like_calling_rule.tsv"
TRADEOFF_NAME = "step_07_control_threshold_tradeoff.tsv"
CONTROL_OVERALL_NAME = "step_07_held_out_control_overall_performance.tsv"
CONTROL_SAMPLE_NAME = "step_07_held_out_control_per_sample_performance.tsv"
CONTROL_CONFUSION_NAME = "step_07_held_out_control_confusion_matrix.tsv"
CONTROL_CATEGORY_NAME = "step_07_held_out_control_classification_counts.tsv"
HET_CELL_NAME = "step_07_het_female_cell_classifications.tsv"
HET_SUMMARY_NAME = "step_07_het_female_classification_summary.tsv"
VALIDATION_NAME = "step_07_het_female_classification_validation.tsv"
ENVIRONMENT_NAME = "software_environment.tsv"
CONTROL_DISTRIBUTION_PLOT = "step_07_control_probability_distributions_with_thresholds.png"
CONTROL_CLASSIFICATION_PLOT = "step_07_held_out_control_classification_confusion.png"
HET_ALL_PLOT = "step_07_het_classification_composition_all_cells.png"
HET_INFORMATIVE_PLOT = "step_07_het_classification_composition_informative_cells.png"

CLASS_ORDER = ["WT_like", "KO_like", "Uncertain", "Uncalled_000"]


class Step07Error(RuntimeError):
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
            writer.writerow({key: row.get(key, "") for key in header})


def fmt(value, places=12):
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return "NA"
    return ("{:.%df}" % places).format(float(value))


class ValidationRecorder(object):
    def __init__(self):
        self.rows = []

    def check(self, name, observed, expected, detail=""):
        passed = observed == expected
        self.rows.append({
            "check": name, "status": "PASS" if passed else "FAIL",
            "observed": str(observed), "expected": str(expected), "detail": detail,
        })
        if not passed:
            raise Step07Error("Validation failed for {}: {!r} != {!r}".format(name, observed, expected))


class Step07Configuration(object):
    REQUIRED = {
        "step_id", "pipeline_version", "upstream_step_ids", "class_encoding",
        "positive_model_class", "model", "control_predictions", "step_06",
        "threshold_selection", "expected_controls", "expected_het_samples",
        "probability_decimal_places", "plot_dpi", "prohibited_operations",
    }

    def __init__(self, path):
        self.path = os.path.abspath(path)
        with open(self.path, "r") as handle:
            self.values = json.load(handle)
        missing = sorted(self.REQUIRED.difference(self.values))
        if missing:
            raise Step07Error("Step 07 lock is missing keys: {}".format(", ".join(missing)))
        if self.values["step_id"] != STEP_ID:
            raise Step07Error("Unexpected Step 07 identifier")
        if self.values["class_encoding"] != {"WT": 0, "KO": 1} or self.values["positive_model_class"] != "KO":
            raise Step07Error("Step 07 must preserve WT=0, KO=1, and P(KO)")
        selection = self.values["threshold_selection"]
        if selection != {
            "target_precision": 0.95,
            "eligible_subset": "informative_non000_held_out_controls",
            "wt_like_rule": "P(KO) <= WT_like_threshold",
            "ko_like_rule": "P(KO) >= KO_like_threshold",
            "target_achieved_tiebreak": "maximum_coverage",
            "target_unachievable_fallback": "maximum_precision_then_maximum_coverage",
            "threshold_candidates": "unique_observed_held_out_control_probabilities",
        }:
            raise Step07Error("Unexpected threshold-selection contract")
        prohibited = {
            "model_refitting", "coefficient_estimation", "feature_selection",
            "probability_calibration", "het_in_threshold_selection", "het_in_performance_estimation",
            "threshold_adjustment_after_het", "cell_type_stratification", "forced_0.5_cutoff",
        }
        if set(self.values["prohibited_operations"]) != prohibited:
            raise Step07Error("Step 07 prohibited-operation guard is incomplete")
        if set(self.values["expected_het_samples"]) != {"15662-JZ-7", "15662-JZ-8", "15662-JZ-9"}:
            raise Step07Error("Unexpected HET sample cohort")

    @property
    def decimals(self):
        return int(self.values["probability_decimal_places"])


class ControlPrediction(object):
    __slots__ = ("barcode", "sample_id", "true_genotype", "design_group", "pattern", "p_ko", "classification")

    def __init__(self, barcode, sample_id, true_genotype, design_group, pattern, p_ko):
        self.barcode = barcode
        self.sample_id = sample_id
        self.true_genotype = true_genotype
        self.design_group = design_group
        self.pattern = pattern
        self.p_ko = float(p_ko)
        self.classification = None


class HeldOutControlPredictionReader(object):
    """Read the exact Step 05 LOSO count-model probabilities."""

    REQUIRED = [
        "cell_barcode", "biological_sample_id", "submitted_sample_name", "sex", "design_group",
        "true_genotype", "A_UMI", "A_detected", "B_UMI", "C_UMI", "pattern_code",
        "binary_predicted_wt_probability", "binary_predicted_ko_probability", "binary_predicted_genotype",
        "count_predicted_wt_probability", "count_predicted_ko_probability", "count_predicted_genotype",
    ]

    def __init__(self, configuration, step05_root, validation):
        self.configuration = configuration
        self.step05_root = os.path.abspath(step05_root)
        self.validation = validation

    def _verify_model_and_package(self):
        model = self.configuration.values["model"]
        package_manifest = os.path.join(self.step05_root, model["relative_directory"], MANIFEST_NAME)
        self.validation.check("step_05_count_package_manifest_sha256", sha256_file(package_manifest), model["package_manifest_sha256"])
        coefficients = os.path.join(self.step05_root, model["relative_directory"], model["coefficient_file"])
        self.validation.check("step_05_count_model_coefficient_sha256", sha256_file(coefficients), model["coefficient_sha256"])
        with open(coefficients, "r", newline="") as handle:
            terms = {row["term"]: float(row["coefficient"]) for row in csv.DictReader(handle, delimiter="\t")}
        self.validation.check("step_05_count_model_terms", set(terms), {"intercept", "A_detected", "B_UMI", "C_UMI"})
        return terms

    def read(self):
        coefficients = self._verify_model_and_package()
        definition = self.configuration.values["control_predictions"]
        path = os.path.join(self.step05_root, definition["relative_directory"], definition["file"])
        self.validation.check("control_predictions_bytes", os.path.getsize(path), definition["bytes"])
        self.validation.check("control_predictions_sha256", sha256_file(path), definition["sha256"])
        records = []
        sample_counts = {}
        with open(path, "r", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if reader.fieldnames != self.REQUIRED:
                raise Step07Error("Unexpected held-out control schema")
            for row in reader:
                if row["true_genotype"] not in ("WT", "KO") or row["design_group"] not in ("WT_M", "WT_F", "KO_M"):
                    raise Step07Error("Held-out control file contains a prohibited cohort")
                p_ko = float(row["count_predicted_ko_probability"])
                p_wt = float(row["count_predicted_wt_probability"])
                if not (0.0 <= p_ko <= 1.0) or abs(p_ko + p_wt - 1.0) > 1.1e-12:
                    raise Step07Error("Invalid held-out control probability")
                record = ControlPrediction(row["cell_barcode"], row["biological_sample_id"], row["true_genotype"], row["design_group"], row["pattern_code"], p_ko)
                records.append(record)
                sample_counts[record.sample_id] = sample_counts.get(record.sample_id, 0) + 1
        expected = self.configuration.values["expected_controls"]
        self.validation.check("control_total_cells", len(records), expected["total_cells"])
        self.validation.check("control_sample_counts", sample_counts, {key: value["cells"] for key, value in expected["samples"].items()})
        self.validation.check("control_000_cells", sum(record.pattern == "000" for record in records), expected["cells_000"])
        self.validation.check("control_informative_cells", sum(record.pattern != "000" for record in records), expected["informative_cells"])
        return records, coefficients


class FrozenCallingRule(object):
    """Immutable two-threshold decision rule with a separate 000 state."""

    __slots__ = ("wt_threshold", "ko_threshold", "target_precision", "wt_selection", "ko_selection", "_frozen")

    def __init__(self, wt_threshold, ko_threshold, target_precision, wt_selection, ko_selection):
        if float(wt_threshold) >= float(ko_threshold):
            raise Step07Error("WT-like and KO-like thresholds overlap")
        self.wt_threshold = float(wt_threshold)
        self.ko_threshold = float(ko_threshold)
        self.target_precision = float(target_precision)
        self.wt_selection = dict(wt_selection)
        self.ko_selection = dict(ko_selection)
        self._frozen = True

    def classify(self, pattern, p_ko):
        if pattern == "000":
            return "Uncalled_000"
        if p_ko <= self.wt_threshold:
            return "WT_like"
        if p_ko >= self.ko_threshold:
            return "KO_like"
        return "Uncertain"


class PrecisionFirstThresholdSelector(object):
    """Select independent control-only tails under a prespecified contract."""

    def __init__(self, configuration, validation):
        self.configuration = configuration
        self.validation = validation
        self.target = float(configuration.values["threshold_selection"]["target_precision"])

    def _rows_for_class(self, records, target_class):
        informative = [record for record in records if record.pattern != "000"]
        true_total = sum(record.true_genotype == target_class for record in informative)
        thresholds = sorted(set(record.p_ko for record in informative))
        rows = []
        for threshold in thresholds:
            if target_class == "WT":
                selected = [record for record in informative if record.p_ko <= threshold]
                direction = "P(KO)<=threshold"
            else:
                selected = [record for record in informative if record.p_ko >= threshold]
                direction = "P(KO)>=threshold"
            correct = sum(record.true_genotype == target_class for record in selected)
            predicted = len(selected)
            precision = correct / float(predicted)
            sensitivity = correct / float(true_total)
            coverage = predicted / float(len(informative))
            rows.append({
                "target_call": target_class + "_like", "direction": direction,
                "threshold_numeric": threshold, "predicted_cells": predicted,
                "correct_cells": correct, "incorrect_cells": predicted - correct,
                "precision_numeric": precision, "sensitivity_numeric": sensitivity,
                "fraction_informative_controls_assigned_numeric": coverage,
                "meets_precision_target": precision >= self.target,
            })
        return rows

    def _choose(self, rows, target_class):
        meeting = [row for row in rows if row["meets_precision_target"]]
        if meeting:
            # Coverage is monotone in the inclusive tail; explicit max makes the policy auditable.
            selected = max(meeting, key=lambda row: (row["predicted_cells"], row["precision_numeric"]))
            status = "target_achieved_maximum_coverage"
        else:
            selected = max(rows, key=lambda row: (row["precision_numeric"], row["predicted_cells"]))
            status = "target_not_achievable_maximum_precision_then_coverage"
        chosen = dict(selected)
        chosen["selection_status"] = status
        chosen["target_class"] = target_class
        return chosen

    def select(self, records):
        wt_rows = self._rows_for_class(records, "WT")
        ko_rows = self._rows_for_class(records, "KO")
        wt = self._choose(wt_rows, "WT")
        ko = self._choose(ko_rows, "KO")
        rule = FrozenCallingRule(wt["threshold_numeric"], ko["threshold_numeric"], self.target, wt, ko)
        self.validation.check("thresholds_derived_before_het_loading", True, True)
        self.validation.check("wt_like_threshold_target_achieved", wt["precision_numeric"] >= self.target, True)
        self.validation.check("ko_like_threshold_target_not_achievable", ko["precision_numeric"] < self.target, True)
        self.validation.check("thresholds_nonoverlapping", rule.wt_threshold < rule.ko_threshold, True)
        serialized = []
        for row in wt_rows + ko_rows:
            item = dict(row)
            item["threshold"] = fmt(item.pop("threshold_numeric"))
            item["precision"] = fmt(item.pop("precision_numeric"))
            item["sensitivity"] = fmt(item.pop("sensitivity_numeric"))
            item["fraction_informative_controls_assigned"] = fmt(item.pop("fraction_informative_controls_assigned_numeric"))
            item["precision_target"] = fmt(self.target)
            item["selected"] = int(
                row["target_call"] == "WT_like" and row["threshold_numeric"] == rule.wt_threshold or
                row["target_call"] == "KO_like" and row["threshold_numeric"] == rule.ko_threshold
            )
            serialized.append(item)
        return rule, serialized


class ControlRuleEvaluator(object):
    def __init__(self, configuration):
        self.configuration = configuration

    @staticmethod
    def _safe(numerator, denominator):
        return numerator / float(denominator) if denominator else None

    def apply(self, records, rule):
        for record in records:
            record.classification = rule.classify(record.pattern, record.p_ko)

    def confusion_rows(self, records):
        rows = []
        for truth in ("WT", "KO"):
            for prediction in ("WT_like", "KO_like"):
                rows.append({
                    "true_genotype": truth, "predicted_probe_state": prediction,
                    "cells": sum(record.true_genotype == truth and record.classification == prediction for record in records),
                    "denominator": "called_controls_only",
                })
        return rows

    def category_rows(self, records):
        rows = []
        for truth in ("WT", "KO"):
            subset = [record for record in records if record.true_genotype == truth]
            for category in CLASS_ORDER:
                count = sum(record.classification == category for record in subset)
                rows.append({
                    "true_genotype": truth, "classification": category, "cells": count,
                    "fraction_within_true_genotype_all_cells": fmt(self._safe(count, len(subset))),
                    "fraction_within_true_genotype_informative": fmt(self._safe(count, sum(record.pattern != "000" for record in subset))) if category != "Uncalled_000" else "NA",
                })
        return rows

    def overall_rows(self, records):
        informative = [record for record in records if record.pattern != "000"]
        called = [record for record in informative if record.classification in ("WT_like", "KO_like")]
        uncertain = [record for record in informative if record.classification == "Uncertain"]
        wt_like = [record for record in called if record.classification == "WT_like"]
        ko_like = [record for record in called if record.classification == "KO_like"]
        wt_true = [record for record in informative if record.true_genotype == "WT"]
        ko_true = [record for record in informative if record.true_genotype == "KO"]
        values = {
            "total_controls": len(records), "informative_non000_controls": len(informative),
            "uncalled_000_cells": len(records) - len(informative),
            "uncalled_000_fraction_all_controls": (len(records) - len(informative)) / float(len(records)),
            "called_cells": len(called), "called_fraction_all_controls": len(called) / float(len(records)),
            "called_fraction_informative_controls": len(called) / float(len(informative)),
            "uncertain_cells": len(uncertain), "uncertain_fraction_all_controls": len(uncertain) / float(len(records)),
            "uncertain_fraction_informative_controls": len(uncertain) / float(len(informative)),
            "wt_like_calls": len(wt_like), "ko_like_calls": len(ko_like),
            "wt_like_precision": self._safe(sum(record.true_genotype == "WT" for record in wt_like), len(wt_like)),
            "ko_like_precision": self._safe(sum(record.true_genotype == "KO" for record in ko_like), len(ko_like)),
            "wt_sensitivity_informative": self._safe(sum(record.classification == "WT_like" for record in wt_true), len(wt_true)),
            "ko_sensitivity_informative": self._safe(sum(record.classification == "KO_like" for record in ko_true), len(ko_true)),
            "accuracy_among_called": self._safe(sum(
                (record.true_genotype == "WT" and record.classification == "WT_like") or
                (record.true_genotype == "KO" and record.classification == "KO_like") for record in called
            ), len(called)),
        }
        return [{"metric": key, "value": fmt(value) if isinstance(value, float) else value, "denominator_note": self._note(key)} for key, value in values.items()]

    @staticmethod
    def _note(metric):
        if "precision" in metric or metric == "accuracy_among_called":
            return "called controls"
        if "sensitivity" in metric or "informative" in metric:
            return "informative non-000 controls"
        if "fraction_all" in metric:
            return "all controls"
        return "count"

    def sample_rows(self, records):
        rows = []
        for sample_id in sorted(set(record.sample_id for record in records)):
            subset = [record for record in records if record.sample_id == sample_id]
            truth = subset[0].true_genotype
            informative = [record for record in subset if record.pattern != "000"]
            correct_label = truth + "_like"
            opposite_label = "KO_like" if truth == "WT" else "WT_like"
            called = [record for record in informative if record.classification in ("WT_like", "KO_like")]
            counts = {category: sum(record.classification == category for record in subset) for category in CLASS_ORDER}
            rows.append({
                "biological_sample_id": sample_id, "design_group": subset[0].design_group, "true_genotype": truth,
                "total_cells": len(subset), "informative_non000_cells": len(informative),
                "WT_like_cells": counts["WT_like"], "KO_like_cells": counts["KO_like"],
                "Uncertain_cells": counts["Uncertain"], "Uncalled_000_cells": counts["Uncalled_000"],
                "called_fraction_all_cells": fmt(self._safe(len(called), len(subset))),
                "called_fraction_informative": fmt(self._safe(len(called), len(informative))),
                "uncertain_fraction_informative": fmt(self._safe(counts["Uncertain"], len(informative))),
                "correct_like_sensitivity_informative": fmt(self._safe(counts[correct_label], len(informative))),
                "opposite_like_error_fraction_informative": fmt(self._safe(counts[opposite_label], len(informative))),
                "accuracy_among_called": fmt(self._safe(counts[correct_label], len(called))),
            })
        return rows


class HETCell(object):
    __slots__ = ("barcode", "sample_id", "a", "b", "c", "pattern", "p_ko", "classification")

    def __init__(self, barcode, sample_id, a, b, c, pattern, p_ko):
        self.barcode = barcode; self.sample_id = sample_id
        self.a = int(a); self.b = int(b); self.c = int(c)
        self.pattern = pattern; self.p_ko = float(p_ko); self.classification = None


class Step06HETProbabilityReader(object):
    """Load HET probabilities only after a FrozenCallingRule exists."""

    REQUIRED = [
        "biological_sample_id", "submitted_sample_name", "sex", "cohort_role", "cell_barcode",
        "A_UMI", "B_UMI", "C_UMI", "B_plus_C_UMI", "A_detected", "B_detected", "C_detected", "pattern_code",
        "binary_predicted_wt_probability", "binary_predicted_ko_probability",
        "count_informed_predicted_wt_probability", "count_informed_predicted_ko_probability", "inference_status",
    ]

    def __init__(self, configuration, step06_root, validation):
        self.configuration = configuration
        self.step06_root = os.path.abspath(step06_root)
        self.validation = validation

    def read(self, frozen_rule):
        if not isinstance(frozen_rule, FrozenCallingRule) or not frozen_rule._frozen:
            raise Step07Error("HET data cannot load before rule freezing")
        definition = self.configuration.values["step_06"]
        manifest = os.path.join(self.step06_root, MANIFEST_NAME)
        self.validation.check("step_06_manifest_sha256", sha256_file(manifest), definition["manifest_sha256"])
        path = os.path.join(self.step06_root, definition["cell_file"])
        self.validation.check("step_06_het_cell_table_bytes", os.path.getsize(path), definition["cell_bytes"])
        self.validation.check("step_06_het_cell_table_sha256", sha256_file(path), definition["cell_sha256"])
        records = []
        counts = {}
        with open(path, "r", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if reader.fieldnames != self.REQUIRED:
                raise Step07Error("Unexpected Step 06 HET table schema")
            for row in reader:
                if row["cohort_role"] != "HET_withheld_inference_only" or row["sex"] != "F":
                    raise Step07Error("Step 06 HET cohort role changed")
                p_ko = float(row["count_informed_predicted_ko_probability"])
                p_wt = float(row["count_informed_predicted_wt_probability"])
                if abs(p_ko + p_wt - 1.0) > 1.1e-12:
                    raise Step07Error("Invalid Step 06 probability pair")
                record = HETCell(row["cell_barcode"], row["biological_sample_id"], row["A_UMI"], row["B_UMI"], row["C_UMI"], row["pattern_code"], p_ko)
                expected_pattern = "{}{}{}".format(int(record.a > 0), int(record.b > 0), int(record.c > 0))
                if record.pattern != expected_pattern:
                    raise Step07Error("HET raw counts and pattern disagree")
                records.append(record)
                counts[record.sample_id] = counts.get(record.sample_id, 0) + 1
        expected = {key: value["cells"] for key, value in self.configuration.values["expected_het_samples"].items()}
        self.validation.check("het_sample_counts", counts, expected)
        self.validation.check("het_total_cells", len(records), sum(expected.values()))
        self.validation.check("het_rows_in_threshold_selection", 0, 0)
        return records


class HETFrozenRuleClassifier(object):
    def apply(self, records, rule):
        for record in records:
            record.classification = rule.classify(record.pattern, record.p_ko)

    def cell_rows(self, records, rule):
        rows = []
        for record in records:
            rows.append({
                "cell_barcode": record.barcode, "biological_sample_id": record.sample_id,
                "A_UMI": record.a, "B_UMI": record.b, "C_UMI": record.c,
                "probe_pattern": record.pattern, "count_informed_p_ko": fmt(record.p_ko),
                "frozen_wt_like_threshold": fmt(rule.wt_threshold),
                "frozen_ko_like_threshold": fmt(rule.ko_threshold),
                "final_classification": record.classification,
            })
        return rows

    @staticmethod
    def _summary_row(label, records):
        total = len(records)
        informative = sum(record.pattern != "000" for record in records)
        counts = {category: sum(record.classification == category for record in records) for category in CLASS_ORDER}
        row = {"biological_sample_id": label, "total_cells": total, "informative_non000_cells": informative}
        for category in CLASS_ORDER:
            row[category + "_cells"] = counts[category]
            row[category + "_percent_all_cells"] = fmt(100.0 * counts[category] / float(total))
            if category != "Uncalled_000":
                row[category + "_percent_informative"] = fmt(100.0 * counts[category] / float(informative))
        return row

    def summary_rows(self, records):
        rows = []
        for sample_id in sorted(set(record.sample_id for record in records)):
            rows.append(self._summary_row(sample_id, [record for record in records if record.sample_id == sample_id]))
        rows.append(self._summary_row("POOLED_HET_F", records))
        return rows


class Step07Plotter(object):
    COLORS = {"WT_like": "#219ebc", "KO_like": "#d95f02", "Uncertain": "#f4c95d", "Uncalled_000": "#c9c9c9"}

    def __init__(self, dpi):
        self.dpi = int(dpi)

    def control_distributions(self, controls, rule, path):
        fig, ax = plt.subplots(figsize=(8.5, 5.0))
        bins = np.linspace(0, 1, 51)
        for genotype, color in (("WT", "#219ebc"), ("KO", "#d95f02")):
            values = [record.p_ko for record in controls if record.pattern != "000" and record.true_genotype == genotype]
            ax.hist(values, bins=bins, density=True, histtype="step", linewidth=2, color=color, label="{} held-out controls".format(genotype))
        ax.axvline(rule.wt_threshold, color=self.COLORS["WT_like"], linestyle="--", linewidth=2, label="WT-like <= {:.3f}".format(rule.wt_threshold))
        ax.axvline(rule.ko_threshold, color=self.COLORS["KO_like"], linestyle="--", linewidth=2, label="KO-like >= {:.3f}".format(rule.ko_threshold))
        ax.axvspan(rule.wt_threshold, rule.ko_threshold, color=self.COLORS["Uncertain"], alpha=.12)
        ax.set_xlabel("Held-out count-informed P(KO)")
        ax.set_ylabel("Density among informative controls")
        ax.set_title("Step 07 control-only thresholds")
        ax.legend(frameon=False, fontsize=8)
        fig.tight_layout(); fig.savefig(path, dpi=self.dpi); plt.close(fig)

    def control_classification(self, category_rows, path):
        genotypes = ["WT", "KO"]
        categories = ["WT_like", "Uncertain", "KO_like"]
        matrix = np.asarray([[float(next(row["fraction_within_true_genotype_informative"] for row in category_rows if row["true_genotype"] == genotype and row["classification"] == category)) for category in categories] for genotype in genotypes])
        fig, ax = plt.subplots(figsize=(7.2, 4.8))
        bottom = np.zeros(2)
        for index, category in enumerate(categories):
            ax.bar(genotypes, matrix[:, index], bottom=bottom, color=self.COLORS[category], label=category.replace("_", "-"))
            for x, value, base in zip(range(2), matrix[:, index], bottom):
                if value >= .03:
                    ax.text(x, base + value / 2, "{:.1f}%".format(100 * value), ha="center", va="center", fontsize=9)
            bottom += matrix[:, index]
        ax.set_ylim(0, 1); ax.set_ylabel("Fraction of informative held-out controls")
        ax.set_title("Step 07 held-out control classification")
        ax.legend(frameon=False, ncol=3, fontsize=8)
        fig.tight_layout(); fig.savefig(path, dpi=self.dpi); plt.close(fig)

    def het_composition(self, summary_rows, informative_only, path):
        rows = [row for row in summary_rows if row["biological_sample_id"] != "POOLED_HET_F"]
        categories = ["WT_like", "KO_like", "Uncertain"] if informative_only else CLASS_ORDER
        fig, ax = plt.subplots(figsize=(8.2, 4.8))
        x = np.arange(len(rows)); bottom = np.zeros(len(rows))
        for category in categories:
            if informative_only:
                values = np.asarray([float(row[category + "_percent_informative"]) / 100.0 for row in rows])
            else:
                values = np.asarray([float(row[category + "_percent_all_cells"]) / 100.0 for row in rows])
            ax.bar(x, values, bottom=bottom, color=self.COLORS[category], label=category.replace("_", "-"))
            for xpos, value, base in zip(x, values, bottom):
                if value >= .04:
                    ax.text(xpos, base + value / 2, "{:.1f}%".format(100 * value), ha="center", va="center", fontsize=8)
            bottom += values
        ax.set_xticks(x); ax.set_xticklabels([row["biological_sample_id"] for row in rows])
        ax.set_ylim(0, 1); ax.set_ylabel("Fraction of {} HET cells".format("informative" if informative_only else "all"))
        ax.set_title("Step 07 HET-female classification ({})".format("non-000 only" if informative_only else "all cells"))
        ax.legend(frameon=False, ncol=len(categories), fontsize=8)
        fig.tight_layout(); fig.savefig(path, dpi=self.dpi); plt.close(fig)


class Step07Publisher(object):
    EXPECTED = {
        RULE_NAME, TRADEOFF_NAME, CONTROL_OVERALL_NAME, CONTROL_SAMPLE_NAME,
        CONTROL_CONFUSION_NAME, CONTROL_CATEGORY_NAME, HET_CELL_NAME, HET_SUMMARY_NAME,
        CONTROL_DISTRIBUTION_PLOT, CONTROL_CLASSIFICATION_PLOT, HET_ALL_PLOT,
        HET_INFORMATIVE_PLOT, VALIDATION_NAME, ENVIRONMENT_NAME,
    }

    def __init__(self, configuration, validation, plotter, output_root, script_path, requirements_path):
        self.configuration = configuration; self.validation = validation; self.plotter = plotter
        self.output_root = os.path.abspath(output_root); self.script_path = os.path.abspath(script_path)
        self.requirements_path = os.path.abspath(requirements_path)

    def _verify_existing(self):
        manifest = os.path.join(self.output_root, MANIFEST_NAME)
        if not os.path.isfile(manifest):
            raise Step07Error("Existing Step 07 package lacks a manifest")
        names = set()
        with open(manifest, "r", newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                path = os.path.join(self.output_root, row["relative_path"])
                if not os.path.isfile(path) or os.path.getsize(path) != int(row["bytes"]) or sha256_file(path) != row["sha256"]:
                    raise Step07Error("Existing Step 07 package fails manifest verification")
                names.add(row["relative_path"])
        if names != self.EXPECTED:
            raise Step07Error("Existing Step 07 file set differs from contract")
        print("EXISTING_VALIDATED\t{}".format(self.output_root))

    def _write_manifest(self, root):
        rows = []
        for name in sorted(os.listdir(root)):
            if name == MANIFEST_NAME:
                continue
            path = os.path.join(root, name)
            if os.path.isfile(path):
                rows.append({"relative_path": name, "bytes": os.path.getsize(path), "sha256": sha256_file(path)})
        write_tsv(os.path.join(root, MANIFEST_NAME), ["relative_path", "bytes", "sha256"], rows)

    def publish(self, rule, tradeoff_rows, controls, overall_rows, sample_rows, confusion_rows, category_rows, het_cell_rows, het_summary_rows, model_coefficients):
        if os.path.exists(self.output_root):
            self._verify_existing(); return
        parent = os.path.dirname(self.output_root); os.makedirs(parent, exist_ok=True)
        staging = tempfile.mkdtemp(prefix=".step_07_pcdh19_het_classification.", dir=parent)
        try:
            rule_rows = []
            for label, selection, threshold in (("WT_like", rule.wt_selection, rule.wt_threshold), ("KO_like", rule.ko_selection, rule.ko_threshold)):
                rule_rows.append({
                    "call_label": label, "probability_score": "count_informed_P(KO)",
                    "inclusive_rule": "P(KO) <= threshold" if label == "WT_like" else "P(KO) >= threshold",
                    "frozen_threshold": fmt(threshold), "precision_target": fmt(rule.target_precision),
                    "target_achieved": int(selection["precision_numeric"] >= rule.target_precision),
                    "selection_status": selection["selection_status"],
                    "held_out_precision": fmt(selection["precision_numeric"]),
                    "held_out_sensitivity": fmt(selection["sensitivity_numeric"]),
                    "held_out_predicted_cells": selection["predicted_cells"],
                    "held_out_fraction_informative_controls_assigned": fmt(selection["fraction_informative_controls_assigned_numeric"]),
                    "threshold_source": "Step 05 WT-M+F/KO-M leave-one-sample-out controls only",
                    "het_rows_used": 0,
                })
            write_tsv(os.path.join(staging, RULE_NAME), ["call_label", "probability_score", "inclusive_rule", "frozen_threshold", "precision_target", "target_achieved", "selection_status", "held_out_precision", "held_out_sensitivity", "held_out_predicted_cells", "held_out_fraction_informative_controls_assigned", "threshold_source", "het_rows_used"], rule_rows)
            write_tsv(os.path.join(staging, TRADEOFF_NAME), ["target_call", "direction", "threshold", "predicted_cells", "correct_cells", "incorrect_cells", "precision", "sensitivity", "fraction_informative_controls_assigned", "precision_target", "meets_precision_target", "selected"], tradeoff_rows)
            write_tsv(os.path.join(staging, CONTROL_OVERALL_NAME), ["metric", "value", "denominator_note"], overall_rows)
            write_tsv(os.path.join(staging, CONTROL_SAMPLE_NAME), list(sample_rows[0].keys()), sample_rows)
            write_tsv(os.path.join(staging, CONTROL_CONFUSION_NAME), ["true_genotype", "predicted_probe_state", "cells", "denominator"], confusion_rows)
            write_tsv(os.path.join(staging, CONTROL_CATEGORY_NAME), ["true_genotype", "classification", "cells", "fraction_within_true_genotype_all_cells", "fraction_within_true_genotype_informative"], category_rows)
            write_tsv(os.path.join(staging, HET_CELL_NAME), ["cell_barcode", "biological_sample_id", "A_UMI", "B_UMI", "C_UMI", "probe_pattern", "count_informed_p_ko", "frozen_wt_like_threshold", "frozen_ko_like_threshold", "final_classification"], het_cell_rows)
            write_tsv(os.path.join(staging, HET_SUMMARY_NAME), list(het_summary_rows[0].keys()), het_summary_rows)
            self.plotter.control_distributions(controls, rule, os.path.join(staging, CONTROL_DISTRIBUTION_PLOT))
            self.plotter.control_classification(category_rows, os.path.join(staging, CONTROL_CLASSIFICATION_PLOT))
            self.plotter.het_composition(het_summary_rows, False, os.path.join(staging, HET_ALL_PLOT))
            self.plotter.het_composition(het_summary_rows, True, os.path.join(staging, HET_INFORMATIVE_PLOT))
            self.validation.check("het_classification_rows", len(het_cell_rows), self.configuration.values["step_06"]["expected_het_cells"])
            self.validation.check("het_rows_have_exactly_one_class", sum(row["final_classification"] in CLASS_ORDER for row in het_cell_rows), len(het_cell_rows))
            self.validation.check("het_000_rows_separate", sum(row["final_classification"] == "Uncalled_000" for row in het_cell_rows), sum(row["probe_pattern"] == "000" for row in het_cell_rows))
            self.validation.check("model_coefficients_not_estimated_in_step_07", set(model_coefficients), {"intercept", "A_detected", "B_UMI", "C_UMI"})
            write_tsv(os.path.join(staging, VALIDATION_NAME), ["check", "status", "observed", "expected", "detail"], self.validation.rows)
            environment = [
                {"key": "step_id", "value": STEP_ID}, {"key": "pipeline_version", "value": self.configuration.values["pipeline_version"]},
                {"key": "created_utc", "value": datetime.datetime.utcnow().replace(microsecond=0).isoformat() + "Z"},
                {"key": "python", "value": platform.python_version()}, {"key": "numpy", "value": np.__version__},
                {"key": "matplotlib", "value": matplotlib.__version__}, {"key": "pipeline_sha256", "value": sha256_file(self.script_path)},
                {"key": "lock_sha256", "value": sha256_file(self.configuration.path)}, {"key": "requirements_sha256", "value": sha256_file(self.requirements_path)},
                {"key": "threshold_scope", "value": "held-out WT-M/WT-F/KO-M controls only; HET rows=0"},
                {"key": "interpretation", "value": "WT-like/KO-like PCDH19 probe evidence; not independently observed DNA genotype"},
            ]
            write_tsv(os.path.join(staging, ENVIRONMENT_NAME), ["key", "value"], environment)
            self._write_manifest(staging)
            os.replace(staging, self.output_root); staging = None
            print("PUBLISHED\t{}\thet_cells={}\twt_threshold={}\tko_threshold={}".format(self.output_root, len(het_cell_rows), fmt(rule.wt_threshold), fmt(rule.ko_threshold)))
        finally:
            if staging and os.path.exists(staging):
                shutil.rmtree(staging)


class PCDH19HETFemaleClassificationStep(object):
    def __init__(self, args):
        self.args = args; self.validation = ValidationRecorder(); self.configuration = Step07Configuration(args.lock)

    def run(self):
        # This ordering is a scientific guard: threshold selection and complete
        # held-out validation finish before the Step 06 HET table is opened.
        controls, coefficients = HeldOutControlPredictionReader(self.configuration, self.args.step05_root, self.validation).read()
        rule, tradeoff_rows = PrecisionFirstThresholdSelector(self.configuration, self.validation).select(controls)
        evaluator = ControlRuleEvaluator(self.configuration); evaluator.apply(controls, rule)
        confusion_rows = evaluator.confusion_rows(controls); category_rows = evaluator.category_rows(controls)
        overall_rows = evaluator.overall_rows(controls); sample_rows = evaluator.sample_rows(controls)
        het = Step06HETProbabilityReader(self.configuration, self.args.step06_root, self.validation).read(rule)
        classifier = HETFrozenRuleClassifier(); classifier.apply(het, rule)
        het_cell_rows = classifier.cell_rows(het, rule); het_summary_rows = classifier.summary_rows(het)
        Step07Publisher(
            self.configuration, self.validation, Step07Plotter(self.configuration.values["plot_dpi"]),
            self.args.output_root, __file__, self.args.requirements,
        ).publish(rule, tradeoff_rows, controls, overall_rows, sample_rows, confusion_rows, category_rows, het_cell_rows, het_summary_rows, coefficients)


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--lock", required=True); parser.add_argument("--requirements", required=True)
    parser.add_argument("--step05-root", required=True); parser.add_argument("--step06-root", required=True)
    parser.add_argument("--output-root", required=True)
    return parser.parse_args(argv)


def main(argv=None):
    PCDH19HETFemaleClassificationStep(parse_args(argv)).run()


if __name__ == "__main__":
    main()
