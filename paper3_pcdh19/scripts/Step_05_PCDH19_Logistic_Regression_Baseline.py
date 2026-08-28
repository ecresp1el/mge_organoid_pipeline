#!/usr/bin/env python3
"""Fit the Step 05 three-feature PCDH19 logistic-regression baseline.

The model is an unpenalized binary logistic regression with the exact formula

``logit(P(KO)) = intercept + beta_A*A + beta_B*B + beta_C*C``.

Step 03 encoding is retained: WT=0 and KO=1, so positive coefficients increase
KO log odds and negative coefficients decrease KO log odds. The only
predictors are binary ``A_detected``, ``B_detected``, and ``C_detected``. The
fit retains all ``000`` cells; the intercept is their fitted KO log odds.

No UMI counts, transcriptome features, cell types, interactions, nonlinear
terms, penalty, or HET cells are introduced. Step 05 compares fitted
eight-pattern probabilities with the Step 04 empirical probabilities and
performs leave-one-registered-sample-out validation of the same model. The
held-out decision policy uses a fixed 0.5 probability threshold for
informative patterns, always leaves ``000`` uncalled, and never optimizes a
threshold on held-out data. The original WT-male/KO-male benchmark is retained
unchanged. A second Step 05 cohort includes WT males plus WT females as WT
ground truth and KO males as KO ground truth; all HET females remain excluded.
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

from Step_04_PCDH19_Empirical_Pattern_Classifier import (
    EmpiricalPatternEstimator,
    LabeledPatternObservation,
    ProbePatternEncoder,
    ProbabilisticClassifier,
    Step03ClassificationTableReader,
    Step04Configuration,
    Step04Error,
    ValidationRecorder,
)


STEP_OUTPUT_DIRECTORY = "step_05_pcdh19_logistic_regression_baseline"
COEFFICIENT_TABLE_NAME = "step_05_pcdh19_logistic_regression_coefficients.tsv"
PROBABILITY_TABLE_NAME = "step_05_pcdh19_logistic_pattern_probabilities.tsv"
COMPARISON_TABLE_NAME = "step_05_pcdh19_empirical_vs_logistic_comparison.tsv"
DIAGNOSTICS_NAME = "step_05_pcdh19_logistic_regression_diagnostics.tsv"
VALIDATION_NAME = "step_05_pcdh19_logistic_regression_validation.tsv"
ENVIRONMENT_NAME = "software_environment.tsv"
MANIFEST_NAME = "output_manifest.tsv"
HELD_OUT_VALIDATION_DIRECTORY = "sample_level_held_out_validation"
EXPANDED_VALIDATION_DIRECTORY = "wt_male_female_vs_ko_male_validation"

HELD_OUT_PREDICTIONS_NAME = "step_05_loso_held_out_cell_predictions.tsv"
HELD_OUT_CONFUSION_NAME = "step_05_loso_confusion_matrix.tsv"
HELD_OUT_PER_SAMPLE_NAME = "step_05_loso_per_sample_metrics.tsv"
HELD_OUT_OVERALL_NAME = "step_05_loso_overall_metrics.tsv"
HELD_OUT_FOLD_COEFFICIENTS_NAME = "step_05_loso_fold_model_coefficients.tsv"
HELD_OUT_CHECKS_NAME = "step_05_loso_validation_checks.tsv"
HELD_OUT_ENVIRONMENT_NAME = "software_environment.tsv"
HELD_OUT_MANIFEST_NAME = "output_manifest.tsv"

HELD_OUT_CONFUSION_PLOT_NAME = "step_05_loso_held_out_confusion_matrix.png"
HELD_OUT_ACCURACY_PLOT_NAME = "step_05_loso_per_sample_accuracy_called.png"
HELD_OUT_CALLED_PLOT_NAME = "step_05_loso_per_sample_percent_called.png"

EXPANDED_PREDICTIONS_NAME = "step_05_wt_mf_ko_m_held_out_cell_predictions.tsv"
EXPANDED_CONFUSION_NAME = "step_05_wt_mf_ko_m_confusion_matrix.tsv"
EXPANDED_PER_SAMPLE_NAME = "step_05_wt_mf_ko_m_per_sample_metrics.tsv"
EXPANDED_GROUP_METRICS_NAME = "step_05_wt_mf_ko_m_metrics_by_group.tsv"
EXPANDED_OVERALL_NAME = "step_05_wt_mf_ko_m_overall_metrics.tsv"
EXPANDED_FOLD_COEFFICIENTS_NAME = "step_05_wt_mf_ko_m_fold_model_coefficients.tsv"
EXPANDED_FULL_COEFFICIENTS_NAME = "step_05_wt_mf_ko_m_full_fit_coefficients.tsv"
EXPANDED_FULL_PROBABILITIES_NAME = "step_05_wt_mf_ko_m_full_fit_pattern_probabilities.tsv"
MALE_PATTERN_ERRORS_NAME = "step_05_male_only_wt_false_ko_errors_by_pattern.tsv"
EXPANDED_PATTERN_ERRORS_NAME = "step_05_wt_mf_ko_m_wt_false_ko_errors_by_pattern_and_sex.tsv"
COHORT_COMPARISON_NAME = "step_05_ground_truth_cohort_validation_comparison.tsv"
EXPANDED_CHECKS_NAME = "step_05_wt_mf_ko_m_validation_checks.tsv"
EXPANDED_ENVIRONMENT_NAME = "software_environment.tsv"
EXPANDED_MANIFEST_NAME = "output_manifest.tsv"
EXPANDED_CONFUSION_PLOT_NAME = "step_05_wt_mf_ko_m_held_out_confusion_matrix.png"
EXPANDED_ACCURACY_PLOT_NAME = "step_05_wt_mf_ko_m_per_sample_accuracy_called.png"
EXPANDED_CALLED_PLOT_NAME = "step_05_wt_mf_ko_m_per_sample_percent_called.png"
PATTERN_ERROR_PLOT_NAME = "step_05_wt_false_ko_errors_by_pattern_cohort.png"
COHORT_COMPARISON_PLOT_NAME = "step_05_ground_truth_cohort_validation_comparison.png"

PROBABILITY_PLOT_NAME = (
    "step_05_logistic_predicted_genotype_probability_by_pattern.png"
)
COMPARISON_PLOT_NAME = (
    "step_05_empirical_vs_logistic_probability_by_pattern.png"
)
COEFFICIENT_PLOT_NAME = "step_05_logistic_coefficient_odds_ratio.png"

COEFFICIENT_HEADER = [
    "term",
    "predictor",
    "outcome_log_odds",
    "coefficient",
    "odds_ratio",
    "interpretation",
]

PROBABILITY_HEADER = [
    "pattern_code",
    "pattern_label",
    "A_detected",
    "B_detected",
    "C_detected",
    "total_ground_truth_cells",
    "linear_predictor_log_odds_ko",
    "p_wt_logistic",
    "p_ko_logistic",
    "hard_call",
]

COMPARISON_HEADER = [
    "pattern_code",
    "pattern_label",
    "total_ground_truth_cells",
    "p_wt_empirical_step_04",
    "p_wt_logistic_step_05",
    "logistic_minus_empirical_p_wt",
    "absolute_difference_p_wt",
    "p_ko_empirical_step_04",
    "p_ko_logistic_step_05",
    "logistic_minus_empirical_p_ko",
    "absolute_difference_p_ko",
    "hard_call",
]

VALIDATION_HEADER = [
    "step_id", "check_name", "status", "observed", "expected", "details"
]

HELD_OUT_PREDICTION_HEADER = [
    "cell_barcode",
    "biological_sample_id",
    "submitted_sample_name",
    "true_genotype",
    "A_detected",
    "B_detected",
    "C_detected",
    "pattern_code",
    "predicted_wt_probability",
    "predicted_ko_probability",
    "predicted_genotype",
]

HELD_OUT_CONFUSION_HEADER = [
    "true_genotype", "predicted_wt_cells", "predicted_ko_cells", "total_called_cells"
]

HELD_OUT_PER_SAMPLE_HEADER = [
    "biological_sample_id",
    "submitted_sample_name",
    "true_genotype",
    "total_cells",
    "called_cells",
    "uncalled_cells",
    "percent_called",
    "accuracy_among_called",
    "ko_sensitivity_where_defined",
    "ko_specificity_where_defined",
]

HELD_OUT_OVERALL_HEADER = ["metric", "value", "definition"]

HELD_OUT_FOLD_COEFFICIENT_HEADER = [
    "held_out_biological_sample_id",
    "held_out_true_genotype",
    "training_samples",
    "training_cells",
    "term",
    "coefficient",
    "odds_ratio",
    "converged",
    "iterations",
]

EXPANDED_PREDICTION_HEADER = [
    "cell_barcode",
    "biological_sample_id",
    "submitted_sample_name",
    "sex",
    "design_group",
    "true_genotype",
    "A_detected",
    "B_detected",
    "C_detected",
    "pattern_code",
    "predicted_wt_probability",
    "predicted_ko_probability",
    "predicted_genotype",
]

EXPANDED_PER_SAMPLE_HEADER = [
    "biological_sample_id",
    "submitted_sample_name",
    "sex",
    "design_group",
    "true_genotype",
    "total_cells",
    "called_cells",
    "uncalled_cells",
    "percent_called",
    "accuracy_among_called",
    "ko_sensitivity_where_defined",
    "ko_specificity_where_defined",
]

GROUP_METRICS_HEADER = [
    "evaluation_group",
    "true_genotype",
    "sex",
    "total_cells",
    "called_cells",
    "uncalled_cells",
    "percent_called",
    "accuracy_among_called",
    "ko_sensitivity_where_defined",
    "ko_specificity_where_defined",
    "false_ko_cells",
]

EXPANDED_FOLD_COEFFICIENT_HEADER = [
    "held_out_biological_sample_id",
    "held_out_true_genotype",
    "held_out_sex",
    "held_out_design_group",
    "training_samples",
    "training_cells",
    "term",
    "coefficient",
    "odds_ratio",
    "converged",
    "iterations",
]

PATTERN_ERROR_HEADER = [
    "cohort",
    "wt_sex",
    "pattern_code",
    "true_wt_cells",
    "predicted_wt_cells",
    "predicted_ko_cells_false_ko",
    "uncalled_cells",
    "false_ko_rate_among_called",
    "percent_of_cohort_false_ko_errors",
]

COHORT_COMPARISON_HEADER = [
    "metric",
    "male_only_value",
    "wt_male_female_ko_male_value",
    "expanded_minus_male_only",
    "interpretation_denominator",
]


class Step05Error(RuntimeError):
    """Raised when a Step 05 model or output invariant fails."""


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
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(
            handle, delimiter="\t", fieldnames=header, lineterminator="\n"
        )
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in header})


class Step05Configuration(object):
    """Load and enforce the deliberately narrow Step 05 model contract."""

    REQUIRED_KEYS = {
        "step_id",
        "pipeline_version",
        "upstream_step_ids",
        "shared_framework_script",
        "step_04_output_manifest_sha256",
        "step_04_empirical_model",
        "pattern_order",
        "predictors",
        "outcome",
        "class_encoding",
        "positive_model_class",
        "model_formula",
        "fit_method",
        "penalty",
        "interaction_terms",
        "nonlinear_transformations",
        "umi_count_predictors",
        "whole_transcriptome_predictors",
        "cell_type_predictors",
        "het_cells_permitted",
        "hard_calls_permitted",
        "expected_class_cells",
        "expected_total_cells",
        "optimizer",
        "probability_decimal_places",
        "plot_dpi",
        "existing_base_package",
        "sample_level_validation",
        "existing_male_only_validation_package",
        "expanded_ground_truth_validation",
    }

    def __init__(self, lock_path, bundle_root):
        self.lock_path = os.path.abspath(lock_path)
        self.bundle_root = os.path.abspath(bundle_root)
        if not os.path.isfile(self.lock_path):
            raise Step05Error("Step 05 lock is missing: {}".format(self.lock_path))
        with open(self.lock_path, "r") as handle:
            self.values = json.load(handle)
        missing = sorted(self.REQUIRED_KEYS.difference(self.values))
        if missing:
            raise Step05Error(
                "Step 05 lock is missing keys: {}".format(", ".join(missing))
            )
        expected_predictors = ["A_detected", "B_detected", "C_detected"]
        if self.values["predictors"] != expected_predictors:
            raise Step05Error("Step 05 predictors must be exactly A, B, and C detection")
        if self.values["class_encoding"] != {"WT": 0, "KO": 1}:
            raise Step05Error("Step 05 class encoding must be WT=0 and KO=1")
        if self.values["positive_model_class"] != "KO":
            raise Step05Error("Step 05 logistic positive class must be KO")
        prohibited_flags = [
            "interaction_terms",
            "nonlinear_transformations",
            "umi_count_predictors",
            "whole_transcriptome_predictors",
            "cell_type_predictors",
            "het_cells_permitted",
            "hard_calls_permitted",
        ]
        if any(self.values[key] is not False for key in prohibited_flags):
            raise Step05Error("Step 05 lock enables a prohibited model feature")
        if self.values["penalty"] != "none":
            raise Step05Error("Step 05 logistic regression must be unpenalized")
        validation = self.values["sample_level_validation"]
        if validation["method"] != "leave_one_registered_sample_out":
            raise Step05Error("Step 05 validation must be leave-one-sample-out")
        if validation["holdout_unit_field"] != "technical_sample_id":
            raise Step05Error("Step 05 holdout unit must be technical_sample_id")
        calling = validation["calling_rule"]
        if calling["uncalled_pattern_codes"] != ["000"]:
            raise Step05Error("Step 05 validation must leave pattern 000 uncalled")
        if float(calling["probability_threshold"]) != 0.5:
            raise Step05Error("Step 05 validation threshold must be fixed at 0.5")
        if calling["threshold_optimized"] is not False:
            raise Step05Error("Held-out data must not optimize the call threshold")
        expanded = self.values["expanded_ground_truth_validation"]
        if expanded["method"] != "leave_one_registered_sample_out":
            raise Step05Error("Expanded Step 05 validation must be leave-one-sample-out")
        if expanded["holdout_unit_field"] != "technical_sample_id":
            raise Step05Error("Expanded Step 05 holdout unit must be technical_sample_id")
        expected_groups = {
            (row["genotype"], row["sex"], row["design_group"])
            for row in expanded["expected_samples"].values()
        }
        if expected_groups != {("WT", "M", "WT_M"), ("WT", "F", "WT_F"), ("KO", "M", "KO_M")}:
            raise Step05Error("Expanded ground truth must be WT-M, WT-F, and KO-M only")
        if sum(row["cells"] for row in expanded["expected_samples"].values()) != expanded["expected_total_cells"]:
            raise Step05Error("Expanded ground-truth cell total is inconsistent")
        if any(row != {"genotype": "HET", "sex": "F"} for row in expanded["excluded_samples"].values()):
            raise Step05Error("Expanded validation exclusions must be HET females")
        expanded_calling = expanded["calling_rule"]
        if expanded_calling != calling:
            raise Step05Error("Male-only and expanded validation calling rules must match")
        self.framework_script_path = os.path.join(
            self.bundle_root,
            self.values["shared_framework_script"]["relative_path"],
        )
        expected_sha = self.values["shared_framework_script"]["sha256"]
        if not os.path.isfile(self.framework_script_path):
            raise Step05Error("Shared Step 04 framework script is missing")
        if sha256_file(self.framework_script_path) != expected_sha:
            raise Step05Error("Shared Step 04 framework script checksum mismatch")

    @property
    def step_id(self):
        return self.values["step_id"]

    @property
    def pattern_order(self):
        return list(self.values["pattern_order"])

    @property
    def decimal_places(self):
        return int(self.values["probability_decimal_places"])


class Step04EmpiricalModelReader(object):
    """Validate and load the serialized Step 04 model used for comparison."""

    def __init__(self, configuration, validation, step04_root):
        self.configuration = configuration
        self.validation = validation
        self.step04_root = os.path.abspath(step04_root)
        self.manifest_path = os.path.join(self.step04_root, MANIFEST_NAME)
        model_relative = self.configuration.values["step_04_empirical_model"][
            "relative_path"
        ]
        self.model_path = os.path.join(self.step04_root, model_relative)
        self._validate_files()

    def _validate_files(self):
        if not os.path.isfile(self.manifest_path):
            raise Step05Error("Step 04 output manifest is missing")
        self.validation.require_equal(
            "step_04_output_manifest_sha256",
            sha256_file(self.manifest_path),
            self.configuration.values["step_04_output_manifest_sha256"],
        )
        entries = {}
        with open(self.manifest_path, "r", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if reader.fieldnames != ["relative_path", "bytes", "sha256"]:
                raise Step05Error("Unexpected Step 04 manifest schema")
            for row in reader:
                if row["relative_path"] in entries:
                    raise Step05Error("Duplicate Step 04 manifest path")
                entries[row["relative_path"]] = row
        expected = self.configuration.values["step_04_empirical_model"]
        manifest_row = entries.get(expected["relative_path"])
        if manifest_row is None:
            raise Step05Error("Step 04 manifest does not list empirical model")
        self.validation.require_equal(
            "step_04_model_manifest_bytes",
            manifest_row["bytes"],
            expected["bytes"],
        )
        self.validation.require_equal(
            "step_04_model_manifest_sha256",
            manifest_row["sha256"],
            expected["sha256"],
        )
        if not os.path.isfile(self.model_path):
            raise Step05Error("Step 04 empirical model is missing")
        self.validation.require_equal(
            "step_04_model_observed_bytes",
            os.path.getsize(self.model_path),
            expected["bytes"],
        )
        self.validation.require_equal(
            "step_04_model_observed_sha256",
            sha256_file(self.model_path),
            expected["sha256"],
        )

    def read_rows(self):
        rows = []
        with open(self.model_path, "r", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            required = {
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
                "hard_call",
            }
            if not required.issubset(set(reader.fieldnames or [])):
                raise Step05Error("Step 04 empirical model schema is incomplete")
            for row in reader:
                if row["hard_call"] != "":
                    raise Step05Error("Step 04 empirical model contains a hard call")
                rows.append(
                    {
                        "pattern_code": row["pattern_code"],
                        "pattern_label": row["pattern_label"],
                        "A_detected": int(row["A_detected"]),
                        "B_detected": int(row["B_detected"]),
                        "C_detected": int(row["C_detected"]),
                        "wt_cells": int(row["wt_cells"]),
                        "ko_cells": int(row["ko_cells"]),
                        "total_cells": int(row["total_cells"]),
                        "p_wt_given_pattern": float(row["p_wt_given_pattern"]),
                        "p_ko_given_pattern": float(row["p_ko_given_pattern"]),
                    }
                )
        self.validation.require_equal(
            "step_04_empirical_patterns",
            [row["pattern_code"] for row in rows],
            self.configuration.pattern_order,
        )
        return rows


class LogisticRegressionClassifier(ProbabilisticClassifier):
    """Unpenalized main-effects probability model with WT=0 and KO=1."""

    TERM_NAMES = ["intercept", "A_detected", "B_detected", "C_detected"]

    def __init__(self, encoder, coefficients, fit_diagnostics):
        self.encoder = encoder
        self.coefficients = np.asarray(coefficients, dtype=float)
        self.fit_diagnostics = dict(fit_diagnostics)
        if self.coefficients.shape != (4,):
            raise Step05Error("Logistic model requires four coefficients")

    @property
    def model_name(self):
        return "pcdh19_three_feature_logistic_regression_baseline"

    @staticmethod
    def _expit(value):
        if value >= 0:
            return 1.0 / (1.0 + math.exp(-value))
        exp_value = math.exp(value)
        return exp_value / (1.0 + exp_value)

    def _feature_vector(self, features):
        if "pattern_code" in features:
            a_value, b_value, c_value = self.encoder.decode(
                features["pattern_code"]
            )
        else:
            a_value = int(features["A_detected"])
            b_value = int(features["B_detected"])
            c_value = int(features["C_detected"])
            self.encoder.encode(a_value, b_value, c_value)
        return np.asarray([1.0, a_value, b_value, c_value], dtype=float)

    def linear_predictor(self, features):
        return float(np.dot(self._feature_vector(features), self.coefficients))

    def predict_proba(self, features):
        p_ko = self._expit(self.linear_predictor(features))
        return {"WT": 1.0 - p_ko, "KO": p_ko}

    def coefficient_rows(self):
        rows = []
        for index, (term, coefficient) in enumerate(
            zip(self.TERM_NAMES, self.coefficients)
        ):
            if index == 0:
                interpretation = "KO:WT odds when A=B=C=0"
                predictor = "baseline_000"
            else:
                interpretation = (
                    "Multiplicative change in KO:WT odds for 0-to-1 detection, "
                    "holding other probes fixed"
                )
                predictor = term
            rows.append(
                {
                    "term": term,
                    "predictor": predictor,
                    "outcome_log_odds": "KO_vs_WT; WT=0, KO=1",
                    "coefficient": float(coefficient),
                    "odds_ratio": float(math.exp(coefficient)),
                    "interpretation": interpretation,
                }
            )
        return rows


class LogisticRegressionEstimator(object):
    """Fit grouped-binomial MLE by Newton/IRLS with no penalty or interactions."""

    def __init__(self, configuration, encoder):
        self.configuration = configuration
        self.encoder = encoder

    @staticmethod
    def _log_likelihood(design, trials, ko_cells, coefficients):
        linear = np.dot(design, coefficients)
        return float(
            np.sum(ko_cells * linear - trials * np.logaddexp(0.0, linear))
        )

    def fit(self, empirical_classifier):
        evidence = empirical_classifier.evidence_rows()
        design = np.asarray(
            [[1.0] + list(self.encoder.decode(row["pattern_code"])) for row in evidence],
            dtype=float,
        )
        wt_cells = np.asarray([row["wt_cells"] for row in evidence], dtype=float)
        ko_cells = np.asarray([row["ko_cells"] for row in evidence], dtype=float)
        trials = wt_cells + ko_cells
        coefficients = np.zeros(4, dtype=float)
        overall_ko = float(np.sum(ko_cells))
        overall_wt = float(np.sum(wt_cells))
        coefficients[0] = math.log(overall_ko / overall_wt)
        optimizer = self.configuration.values["optimizer"]
        maximum_iterations = int(optimizer["maximum_iterations"])
        tolerance = float(optimizer["coefficient_tolerance"])
        maximum_halvings = int(optimizer["maximum_step_halvings"])
        converged = False
        step_halvings = 0

        for iteration in range(1, maximum_iterations + 1):
            linear = np.dot(design, coefficients)
            probabilities = 1.0 / (1.0 + np.exp(-linear))
            weights = trials * probabilities * (1.0 - probabilities)
            score = np.dot(design.T, ko_cells - trials * probabilities)
            information = np.dot(design.T, weights[:, None] * design)
            try:
                delta = np.linalg.solve(information, score)
            except np.linalg.LinAlgError as error:
                raise Step05Error(
                    "Logistic information matrix is singular: {}".format(error)
                )
            old_log_likelihood = self._log_likelihood(
                design, trials, ko_cells, coefficients
            )
            scale = 1.0
            candidate = coefficients + delta
            candidate_log_likelihood = self._log_likelihood(
                design, trials, ko_cells, candidate
            )
            halvings_this_iteration = 0
            while (
                candidate_log_likelihood < old_log_likelihood
                and halvings_this_iteration < maximum_halvings
            ):
                scale *= 0.5
                candidate = coefficients + scale * delta
                candidate_log_likelihood = self._log_likelihood(
                    design, trials, ko_cells, candidate
                )
                halvings_this_iteration += 1
            step_halvings += halvings_this_iteration
            if candidate_log_likelihood < old_log_likelihood:
                raise Step05Error("IRLS failed to find a likelihood-increasing step")
            maximum_change = float(np.max(np.abs(candidate - coefficients)))
            coefficients = candidate
            if maximum_change <= tolerance:
                converged = True
                break

        linear = np.dot(design, coefficients)
        probabilities = 1.0 / (1.0 + np.exp(-linear))
        weights = trials * probabilities * (1.0 - probabilities)
        final_score = np.dot(design.T, ko_cells - trials * probabilities)
        information = np.dot(design.T, weights[:, None] * design)
        diagnostics = {
            "converged": converged,
            "iterations": iteration,
            "step_halvings": step_halvings,
            "log_likelihood": self._log_likelihood(
                design, trials, ko_cells, coefficients
            ),
            "maximum_absolute_score": float(np.max(np.abs(final_score))),
            "information_matrix_rank": int(np.linalg.matrix_rank(information)),
            "information_matrix_condition_number": float(
                np.linalg.cond(information)
            ),
            "fit_cells": int(np.sum(trials)),
            "fit_patterns": len(evidence),
        }
        if not converged:
            raise Step05Error("Logistic regression did not converge")
        return LogisticRegressionClassifier(
            self.encoder, coefficients, diagnostics
        )


class LogisticModelEvaluator(object):
    """Compare probability behavior without thresholds or hard-call metrics."""

    def __init__(self, configuration, encoder, validation):
        self.configuration = configuration
        self.encoder = encoder
        self.validation = validation

    def evaluate(self, classifier, empirical_rows, recomputed_empirical):
        recomputed_by_pattern = {
            row["pattern_code"]: row
            for row in recomputed_empirical.evidence_rows()
        }
        probability_rows = []
        comparison_rows = []
        squared_differences = []
        weighted_squared_differences = []
        total_cells = sum(row["total_cells"] for row in empirical_rows)
        for empirical in empirical_rows:
            pattern = empirical["pattern_code"]
            recomputed = recomputed_by_pattern[pattern]
            self.validation.require_equal(
                "{}_step_04_counts_match_step_03".format(pattern),
                (
                    empirical["wt_cells"],
                    empirical["ko_cells"],
                    empirical["total_cells"],
                ),
                (
                    recomputed["wt_cells"],
                    recomputed["ko_cells"],
                    recomputed["total_cells"],
                ),
            )
            self.validation.require_close(
                "{}_step_04_p_ko_matches_step_03".format(pattern),
                empirical["p_ko_given_pattern"],
                recomputed["p_ko_given_pattern"],
                tolerance=5e-12,
            )
            features = {"pattern_code": pattern}
            probabilities = classifier.predict_proba(features)
            self.validation.require_close(
                "{}_logistic_probability_sum".format(pattern),
                probabilities["WT"] + probabilities["KO"],
                1.0,
            )
            a_value, b_value, c_value = self.encoder.decode(pattern)
            probability_rows.append(
                {
                    "pattern_code": pattern,
                    "pattern_label": empirical["pattern_label"],
                    "A_detected": a_value,
                    "B_detected": b_value,
                    "C_detected": c_value,
                    "total_ground_truth_cells": empirical["total_cells"],
                    "linear_predictor_log_odds_ko": classifier.linear_predictor(
                        features
                    ),
                    "p_wt_logistic": probabilities["WT"],
                    "p_ko_logistic": probabilities["KO"],
                    "hard_call": None,
                }
            )
            difference_ko = probabilities["KO"] - empirical["p_ko_given_pattern"]
            difference_wt = probabilities["WT"] - empirical["p_wt_given_pattern"]
            squared_differences.append(difference_ko ** 2)
            weighted_squared_differences.append(
                empirical["total_cells"] * difference_ko ** 2
            )
            comparison_rows.append(
                {
                    "pattern_code": pattern,
                    "pattern_label": empirical["pattern_label"],
                    "total_ground_truth_cells": empirical["total_cells"],
                    "p_wt_empirical_step_04": empirical["p_wt_given_pattern"],
                    "p_wt_logistic_step_05": probabilities["WT"],
                    "logistic_minus_empirical_p_wt": difference_wt,
                    "absolute_difference_p_wt": abs(difference_wt),
                    "p_ko_empirical_step_04": empirical["p_ko_given_pattern"],
                    "p_ko_logistic_step_05": probabilities["KO"],
                    "logistic_minus_empirical_p_ko": difference_ko,
                    "absolute_difference_p_ko": abs(difference_ko),
                    "hard_call": None,
                }
            )
        self.validation.require_equal(
            "logistic_probability_patterns", len(probability_rows), 8
        )
        self.validation.require_equal(
            "pattern_000_retained_in_logistic_fit",
            sum(row["pattern_code"] == "000" for row in probability_rows),
            1,
        )
        self.validation.require_equal(
            "hard_genotype_calls_created",
            sum(row["hard_call"] is not None for row in probability_rows),
            0,
        )
        metrics = {
            "unweighted_pattern_rmse_p_ko": math.sqrt(
                sum(squared_differences) / len(squared_differences)
            ),
            "cell_weighted_rmse_p_ko": math.sqrt(
                sum(weighted_squared_differences) / total_cells
            ),
            "maximum_absolute_pattern_difference_p_ko": max(
                abs(row["logistic_minus_empirical_p_ko"])
                for row in comparison_rows
            ),
        }
        return probability_rows, comparison_rows, metrics


class SampleAwarePatternRecord(object):
    """One WT/KO cell retaining the registered sample used for holdout."""

    __slots__ = (
        "cell_barcode",
        "biological_sample_id",
        "submitted_sample_name",
        "sex",
        "design_group",
        "true_genotype",
        "pattern_code",
    )

    def __init__(
        self,
        cell_barcode,
        biological_sample_id,
        submitted_sample_name,
        sex,
        design_group,
        true_genotype,
        pattern_code,
    ):
        self.cell_barcode = cell_barcode
        self.biological_sample_id = biological_sample_id
        self.submitted_sample_name = submitted_sample_name
        self.sex = sex
        self.design_group = design_group
        self.true_genotype = true_genotype
        self.pattern_code = pattern_code


class Step03SampleAwarePatternReader(object):
    """Load Step 03 cells while retaining the registered sample holdout key."""

    def __init__(self, configuration, shared_configuration, encoder, validation, table_path):
        self.configuration = configuration
        self.shared_configuration = shared_configuration
        self.encoder = encoder
        self.validation = validation
        self.table_path = table_path

    def read_records(self):
        expected_samples = self.configuration.values["sample_level_validation"][
            "expected_samples"
        ]
        observed_counts = {sample_id: 0 for sample_id in expected_samples}
        records = []
        with open(self.table_path, "r", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if reader.fieldnames != self.shared_configuration.values["required_input_columns"]:
                raise Step05Error("Unexpected Step 03 schema for sample-aware validation")
            for row in reader:
                sample_id = row["technical_sample_id"]
                if sample_id not in expected_samples:
                    raise Step05Error(
                        "Unexpected sample in Step 05 validation: {}".format(sample_id)
                    )
                expected_genotype = expected_samples[sample_id]["genotype"]
                if row["genotype"] != expected_genotype or row["ground_truth_class"] != expected_genotype:
                    raise Step05Error(
                        "Sample/genotype mismatch in held-out validation: {}".format(sample_id)
                    )
                if row["sex"] != "M":
                    raise Step05Error("Held-out validation input contains a non-male cell")
                pattern = self.encoder.encode(
                    row["A_detected"], row["B_detected"], row["C_detected"]
                )
                records.append(
                    SampleAwarePatternRecord(
                        row["cell_barcode"],
                        sample_id,
                        row["submitted_sample_name"],
                        row["sex"],
                        row["design_group"],
                        row["ground_truth_class"],
                        pattern,
                    )
                )
                observed_counts[sample_id] += 1
        for sample_id in expected_samples:
            self.validation.require_equal(
                "{}_held_out_unit_cells".format(sample_id),
                observed_counts[sample_id],
                expected_samples[sample_id]["cells"],
            )
        self.validation.require_equal(
            "sample_aware_validation_cells",
            len(records),
            self.configuration.values["expected_total_cells"],
        )
        self.validation.require_equal(
            "leave_one_sample_out_units",
            len(observed_counts),
            6,
        )
        return records


class ManifestedProbePatternCohortReader(object):
    """Extend the Step 03 male contract with manifested WT-female tables."""

    REQUIRED_COLUMNS = [
        "barcode", "A_UMI", "B_UMI", "C_UMI", "Pcdh19_total_UMI", "detection_pattern"
    ]
    PATTERN_LABELS = {
        "000": "none", "001": "C only", "010": "B only", "011": "B+C",
        "100": "A only", "101": "A+C", "110": "A+B", "111": "A+B+C",
    }

    def __init__(self, configuration, encoder, validation, sample_key_path, step02_root):
        self.configuration = configuration
        self.encoder = encoder
        self.validation = validation
        self.cohort = configuration.values["expanded_ground_truth_validation"]
        self.sample_key_path = os.path.abspath(sample_key_path)
        self.step02_root = os.path.abspath(step02_root)
        self.manifest_path = os.path.join(self.step02_root, "output_manifest.tsv")

    def _manifest(self):
        if not os.path.isfile(self.manifest_path):
            raise Step05Error("Expanded-cohort Step 02a manifest is missing")
        self.validation.require_equal(
            "expanded_step_02a_manifest_sha256",
            sha256_file(self.manifest_path),
            self.cohort["step_02a_output_manifest_sha256"],
        )
        rows = {}
        with open(self.manifest_path, "r", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if reader.fieldnames != ["relative_path", "bytes", "sha256"]:
                raise Step05Error("Unexpected Step 02a manifest schema")
            for row in reader:
                rows[row["relative_path"]] = row
        return rows

    def _sample_metadata(self):
        if not os.path.isfile(self.sample_key_path):
            raise Step05Error("Expanded-cohort sample key is missing")
        self.validation.require_equal(
            "expanded_sample_key_sha256",
            sha256_file(self.sample_key_path),
            self.cohort["sample_key_sha256"],
        )
        metadata = {}
        with open(self.sample_key_path, "r", newline="") as handle:
            reader = csv.DictReader(handle)
            required = {"technical_sample_id", "submitted_sample_name", "genotype", "sex", "design_group"}
            if not required.issubset(set(reader.fieldnames or [])):
                raise Step05Error("Expanded-cohort sample key schema is incomplete")
            for row in reader:
                sample_id = row["technical_sample_id"]
                if sample_id in metadata:
                    raise Step05Error("Duplicate sample in registered sample key")
                metadata[sample_id] = row
        expected = self.cohort["expected_samples"]
        excluded = self.cohort["excluded_samples"]
        for sample_id, definition in expected.items():
            row = metadata.get(sample_id)
            if row is None:
                raise Step05Error("Expanded ground-truth sample is not registered")
            for field in ("genotype", "sex", "design_group"):
                if row[field] != definition[field]:
                    raise Step05Error("Expanded sample metadata differs from lock: {}".format(sample_id))
        for sample_id, definition in excluded.items():
            row = metadata.get(sample_id)
            if row is None or row["genotype"] != definition["genotype"] or row["sex"] != definition["sex"]:
                raise Step05Error("HET-female exclusion metadata differs from lock")
        self.validation.require_equal("expanded_het_female_samples_excluded", len(excluded), 3)
        return metadata

    def read_records(self, male_records):
        manifest = self._manifest()
        metadata = self._sample_metadata()
        expected = self.cohort["expected_samples"]
        records = list(male_records)
        class_counts = {"WT": 0, "KO": 0}
        group_counts = {"WT_M": 0, "WT_F": 0, "KO_M": 0}
        male_sample_counts = {sample_id: 0 for sample_id, row in expected.items() if row["sex"] == "M"}
        for record in records:
            definition = expected.get(record.biological_sample_id)
            if definition is None or definition["sex"] != "M":
                raise Step05Error("Expanded cohort received a non-male Step 03 base record")
            if record.true_genotype != definition["genotype"] or record.sex != "M":
                raise Step05Error("Expanded cohort Step 03 male metadata mismatch")
            male_sample_counts[record.biological_sample_id] += 1
            class_counts[record.true_genotype] += 1
            group_counts[record.design_group] += 1
        for sample_id, observed in male_sample_counts.items():
            self.validation.require_equal(
                "expanded_{}_step_03_male_cells".format(sample_id),
                observed,
                expected[sample_id]["cells"],
            )
        for sample_id, definition in expected.items():
            if definition["sex"] == "M":
                continue
            relative_path = os.path.join("per_sample", sample_id, "pcdh19_probe_patterns.tsv")
            manifest_row = manifest.get(relative_path)
            if manifest_row is None:
                raise Step05Error("Step 02a manifest omits expanded sample {}".format(sample_id))
            path = os.path.join(self.step02_root, relative_path)
            if not os.path.isfile(path):
                raise Step05Error("Expanded sample probe-pattern table is missing")
            if os.path.getsize(path) != int(manifest_row["bytes"]) or sha256_file(path) != manifest_row["sha256"]:
                raise Step05Error("Expanded sample probe-pattern identity mismatch")
            seen = set()
            row_count = 0
            with open(path, "r", newline="") as handle:
                reader = csv.DictReader(handle, delimiter="\t")
                if reader.fieldnames != self.REQUIRED_COLUMNS:
                    raise Step05Error("Unexpected expanded sample pattern schema")
                for row in reader:
                    barcode = row["barcode"]
                    if not barcode or barcode in seen:
                        raise Step05Error("Empty or duplicate expanded sample barcode")
                    seen.add(barcode)
                    try:
                        counts = [int(row[name]) for name in ("A_UMI", "B_UMI", "C_UMI")]
                        total = int(row["Pcdh19_total_UMI"])
                    except ValueError:
                        raise Step05Error("Invalid expanded sample probe count")
                    if min(counts + [total]) < 0 or sum(counts) != total:
                        raise Step05Error("Expanded sample probe counts are inconsistent")
                    pattern = self.encoder.encode(*(int(value > 0) for value in counts))
                    if row["detection_pattern"] != self.PATTERN_LABELS[pattern]:
                        raise Step05Error("Expanded sample detection pattern is inconsistent")
                    registered = metadata[sample_id]
                    records.append(
                        SampleAwarePatternRecord(
                            barcode,
                            sample_id,
                            registered["submitted_sample_name"],
                            definition["sex"],
                            definition["design_group"],
                            definition["genotype"],
                            pattern,
                        )
                    )
                    row_count += 1
            self.validation.require_equal(
                "expanded_{}_cells".format(sample_id), row_count, definition["cells"]
            )
            class_counts[definition["genotype"]] += row_count
            group_counts[definition["design_group"]] += row_count
        self.validation.require_equal(
            "expanded_ground_truth_cells", len(records), self.cohort["expected_total_cells"]
        )
        for genotype, expected_count in self.cohort["expected_class_cells"].items():
            self.validation.require_equal(
                "expanded_{}_class_cells".format(genotype.lower()), class_counts[genotype], expected_count
            )
        for group, expected_count in self.cohort["expected_sex_class_cells"].items():
            self.validation.require_equal(
                "expanded_{}_cells".format(group.lower()), group_counts[group], expected_count
            )
        self.validation.require_equal(
            "expanded_het_cells_loaded",
            sum(record.biological_sample_id in self.cohort["excluded_samples"] for record in records),
            0,
        )
        return records


class HeldOutCallingPolicy(object):
    """Apply the predeclared 0.5 rule while always leaving 000 uncalled."""

    def __init__(self, calling_rule):
        self.uncalled_patterns = set(calling_rule["uncalled_pattern_codes"])
        self.threshold = float(calling_rule["probability_threshold"])

    def call(self, pattern_code, probabilities):
        if pattern_code in self.uncalled_patterns:
            return "uncalled"
        if probabilities["KO"] > self.threshold:
            return "KO"
        if probabilities["WT"] > self.threshold:
            return "WT"
        return "uncalled"


class LeaveOneSampleOutLogisticValidator(object):
    """Fit one logistic model per registered-sample holdout."""

    def __init__(
        self,
        configuration,
        encoder,
        validation,
        calling_policy,
        validation_config=None,
        check_prefix="",
    ):
        self.configuration = configuration
        self.encoder = encoder
        self.validation = validation
        self.calling_policy = calling_policy
        self.validation_config = validation_config or configuration.values["sample_level_validation"]
        self.check_prefix = check_prefix
        self.estimator = LogisticRegressionEstimator(configuration, encoder)

    def _check(self, name):
        return "{}{}".format(self.check_prefix, name)

    def run(self, records):
        sample_config = self.validation_config["expected_samples"]
        sample_order = list(sample_config.keys())
        predictions = []
        fold_coefficients = []
        for held_out_sample in sample_order:
            training_records = [
                record
                for record in records
                if record.biological_sample_id != held_out_sample
            ]
            test_records = [
                record
                for record in records
                if record.biological_sample_id == held_out_sample
            ]
            training_samples = sorted(
                {record.biological_sample_id for record in training_records}
            )
            self.validation.require_equal(
                self._check("{}_training_samples_exclude_holdout".format(held_out_sample)),
                held_out_sample in training_samples,
                False,
            )
            self.validation.require_equal(
                self._check("{}_training_sample_count".format(held_out_sample)),
                len(training_samples),
                len(sample_config) - 1,
            )
            self.validation.require_equal(
                self._check("{}_held_out_cells".format(held_out_sample)),
                len(test_records),
                sample_config[held_out_sample]["cells"],
            )
            observations = (
                LabeledPatternObservation(
                    record.pattern_code, record.true_genotype
                )
                for record in training_records
            )
            empirical_training = EmpiricalPatternEstimator(self.encoder).fit(
                observations
            )
            classifier = self.estimator.fit(empirical_training)
            self.validation.require_equal(
                self._check("{}_fold_model_converged".format(held_out_sample)),
                classifier.fit_diagnostics["converged"],
                True,
            )
            for coefficient_row in classifier.coefficient_rows():
                fold_coefficients.append(
                    {
                        "held_out_biological_sample_id": held_out_sample,
                        "held_out_true_genotype": sample_config[held_out_sample]["genotype"],
                        "held_out_sex": sample_config[held_out_sample].get("sex", "M"),
                        "held_out_design_group": sample_config[held_out_sample].get(
                            "design_group", "{}_M".format(sample_config[held_out_sample]["genotype"])
                        ),
                        "training_samples": ",".join(training_samples),
                        "training_cells": len(training_records),
                        "term": coefficient_row["term"],
                        "coefficient": coefficient_row["coefficient"],
                        "odds_ratio": coefficient_row["odds_ratio"],
                        "converged": classifier.fit_diagnostics["converged"],
                        "iterations": classifier.fit_diagnostics["iterations"],
                    }
                )
            for record in test_records:
                probabilities = classifier.predict_proba(
                    {"pattern_code": record.pattern_code}
                )
                a_value, b_value, c_value = self.encoder.decode(
                    record.pattern_code
                )
                predictions.append(
                    {
                        "cell_barcode": record.cell_barcode,
                        "biological_sample_id": record.biological_sample_id,
                        "submitted_sample_name": record.submitted_sample_name,
                        "sex": record.sex,
                        "design_group": record.design_group,
                        "true_genotype": record.true_genotype,
                        "A_detected": a_value,
                        "B_detected": b_value,
                        "C_detected": c_value,
                        "pattern_code": record.pattern_code,
                        "predicted_wt_probability": probabilities["WT"],
                        "predicted_ko_probability": probabilities["KO"],
                        "predicted_genotype": self.calling_policy.call(
                            record.pattern_code, probabilities
                        ),
                    }
                )
        self.validation.require_equal(
            self._check("held_out_prediction_rows"), len(predictions), len(records)
        )
        self.validation.require_equal(
            self._check("each_cell_predicted_once"),
            len(
                {
                    (row["biological_sample_id"], row["cell_barcode"])
                    for row in predictions
                }
            ),
            len(records),
        )
        self.validation.require_equal(
            self._check("fold_coefficient_rows"), len(fold_coefficients), 4 * len(sample_config)
        )
        return predictions, fold_coefficients


class HeldOutValidationEvaluator(object):
    """Compute called-cell confusion and sample/overall validation metrics."""

    def __init__(self, configuration, validation, validation_config=None, check_prefix=""):
        self.configuration = configuration
        self.validation = validation
        self.validation_config = validation_config or configuration.values["sample_level_validation"]
        self.check_prefix = check_prefix

    def _check(self, name):
        return "{}{}".format(self.check_prefix, name)

    @staticmethod
    def _safe_fraction(numerator, denominator):
        if denominator == 0:
            return None
        return numerator / denominator

    @staticmethod
    def _confusion(predictions):
        counts = {("WT", "WT"): 0, ("WT", "KO"): 0, ("KO", "WT"): 0, ("KO", "KO"): 0}
        for row in predictions:
            predicted = row["predicted_genotype"]
            if predicted in ("WT", "KO"):
                counts[(row["true_genotype"], predicted)] += 1
        return counts

    def evaluate(self, predictions):
        sample_config = self.validation_config["expected_samples"]
        per_sample = []
        for sample_id in sample_config:
            sample_rows = [row for row in predictions if row["biological_sample_id"] == sample_id]
            called = [row for row in sample_rows if row["predicted_genotype"] in ("WT", "KO")]
            correct = sum(row["predicted_genotype"] == row["true_genotype"] for row in called)
            true_genotype = sample_config[sample_id]["genotype"]
            sensitivity = None
            specificity = None
            if true_genotype == "KO":
                sensitivity = self._safe_fraction(correct, len(called))
            else:
                specificity = self._safe_fraction(correct, len(called))
            per_sample.append(
                {
                    "biological_sample_id": sample_id,
                    "submitted_sample_name": sample_rows[0]["submitted_sample_name"],
                    "sex": sample_rows[0].get("sex", sample_config[sample_id].get("sex", "M")),
                    "design_group": sample_rows[0].get(
                        "design_group", sample_config[sample_id].get("design_group", "")
                    ),
                    "true_genotype": true_genotype,
                    "total_cells": len(sample_rows),
                    "called_cells": len(called),
                    "uncalled_cells": len(sample_rows) - len(called),
                    "percent_called": 100.0 * self._safe_fraction(len(called), len(sample_rows)),
                    "accuracy_among_called": self._safe_fraction(correct, len(called)),
                    "ko_sensitivity_where_defined": sensitivity,
                    "ko_specificity_where_defined": specificity,
                }
            )
        confusion = self._confusion(predictions)
        tn = confusion[("WT", "WT")]
        fp = confusion[("WT", "KO")]
        fn = confusion[("KO", "WT")]
        tp = confusion[("KO", "KO")]
        total = len(predictions)
        called = tn + fp + fn + tp
        uncalled = total - called
        cohort_id = self.validation_config.get("cohort_id", "wt_male_vs_ko_male")
        overall_values = [
            ("total_cells", total, "All held-out cells in {}".format(cohort_id)),
            ("called_cells", called, "Non-000 cells receiving WT or KO at fixed threshold"),
            ("uncalled_cells", uncalled, "Pattern 000 plus any exact probability tie"),
            ("percent_called", 100.0 * called / total, "100 * called / total"),
            ("percent_uncalled", 100.0 * uncalled / total, "100 * uncalled / total"),
            ("accuracy_among_called", self._safe_fraction(tp + tn, called), "(TP+TN) / called"),
            ("ko_sensitivity_among_called", self._safe_fraction(tp, tp + fn), "TP / (TP+FN); KO is positive"),
            ("ko_specificity_among_called", self._safe_fraction(tn, tn + fp), "TN / (TN+FP); KO is positive"),
            ("true_ko_predicted_ko", tp, "TP"),
            ("true_ko_predicted_wt", fn, "FN"),
            ("true_wt_predicted_ko", fp, "FP"),
            ("true_wt_predicted_wt", tn, "TN"),
        ]
        overall = [
            {"metric": key, "value": value, "definition": definition}
            for key, value, definition in overall_values
        ]
        confusion_rows = [
            {"true_genotype": "WT", "predicted_wt_cells": tn, "predicted_ko_cells": fp, "total_called_cells": tn + fp},
            {"true_genotype": "KO", "predicted_wt_cells": fn, "predicted_ko_cells": tp, "total_called_cells": fn + tp},
        ]
        self.validation.require_equal(
            self._check("pattern_000_predictions_uncalled"),
            sum(row["pattern_code"] == "000" and row["predicted_genotype"] != "uncalled" for row in predictions),
            0,
        )
        self.validation.require_equal(
            self._check("called_confusion_total"), called, sum(row["called_cells"] for row in per_sample)
        )
        self.validation.require_equal(
            self._check("fixed_call_threshold"), self.validation_config["calling_rule"]["probability_threshold"], 0.5
        )
        self.validation.require_equal(
            self._check("held_out_threshold_optimized"), self.validation_config["calling_rule"]["threshold_optimized"], False
        )
        return confusion_rows, per_sample, overall


class ExistingMaleOnlyValidationReader(object):
    """Verify and read the immutable male-only validation benchmark."""

    def __init__(self, configuration, validation, base_root):
        self.configuration = configuration
        self.validation = validation
        definition = configuration.values["existing_male_only_validation_package"]
        self.root = os.path.join(os.path.abspath(base_root), definition["relative_directory"])
        self.manifest_path = os.path.join(self.root, HELD_OUT_MANIFEST_NAME)
        self.predictions_path = os.path.join(self.root, HELD_OUT_PREDICTIONS_NAME)
        self.overall_path = os.path.join(self.root, HELD_OUT_OVERALL_NAME)

    def verify(self):
        definition = self.configuration.values["existing_male_only_validation_package"]
        if not os.path.isfile(self.manifest_path):
            raise Step05Error("Male-only validation benchmark is missing")
        self.validation.require_equal(
            "male_only_validation_manifest_unchanged",
            sha256_file(self.manifest_path),
            definition["output_manifest_sha256"],
        )
        with open(self.manifest_path, "r", newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        for row in rows:
            path = os.path.join(self.root, row["relative_path"])
            if not os.path.isfile(path) or os.path.getsize(path) != int(row["bytes"]):
                raise Step05Error("Male-only validation benchmark is incomplete")
            if sha256_file(path) != row["sha256"]:
                raise Step05Error("Male-only validation artifact identity mismatch")
        self.validation.require_equal("male_only_validation_manifest_files", len(rows), 10)
        self.validation.require_equal(
            "male_only_prediction_table_unchanged",
            sha256_file(self.predictions_path),
            definition["held_out_predictions_sha256"],
        )

    def read_predictions(self):
        self.verify()
        rows = []
        with open(self.predictions_path, "r", newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                row["sex"] = "M"
                row["design_group"] = "{}_M".format(row["true_genotype"])
                rows.append(row)
        self.validation.require_equal("male_only_benchmark_prediction_rows", len(rows), 230269)
        return rows

    def read_overall(self):
        values = {}
        with open(self.overall_path, "r", newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                values[row["metric"]] = float(row["value"])
        return values


class PatternErrorAnalyzer(object):
    """Attribute held-out WT-to-KO errors to binary probe patterns."""

    PATTERNS = ["001", "010", "011", "100", "101", "110", "111"]

    @staticmethod
    def _rows_for_group(predictions, cohort, wt_sex, sex_filter=None):
        subset = [
            row for row in predictions
            if row["true_genotype"] == "WT"
            and row["pattern_code"] != "000"
            and (sex_filter is None or row.get("sex") == sex_filter)
        ]
        total_errors = sum(row["predicted_genotype"] == "KO" for row in subset)
        rows = []
        for pattern in PatternErrorAnalyzer.PATTERNS:
            pattern_rows = [row for row in subset if row["pattern_code"] == pattern]
            predicted_wt = sum(row["predicted_genotype"] == "WT" for row in pattern_rows)
            predicted_ko = sum(row["predicted_genotype"] == "KO" for row in pattern_rows)
            uncalled = sum(row["predicted_genotype"] == "uncalled" for row in pattern_rows)
            called = predicted_wt + predicted_ko
            rows.append(
                {
                    "cohort": cohort,
                    "wt_sex": wt_sex,
                    "pattern_code": pattern,
                    "true_wt_cells": len(pattern_rows),
                    "predicted_wt_cells": predicted_wt,
                    "predicted_ko_cells_false_ko": predicted_ko,
                    "uncalled_cells": uncalled,
                    "false_ko_rate_among_called": predicted_ko / called if called else None,
                    "percent_of_cohort_false_ko_errors": 100.0 * predicted_ko / total_errors if total_errors else None,
                }
            )
        return rows

    def analyze_male_only(self, predictions):
        return self._rows_for_group(predictions, "male_only", "M", "M")

    def analyze_expanded(self, predictions):
        rows = self._rows_for_group(predictions, "wt_male_female_vs_ko_male", "ALL")
        rows.extend(self._rows_for_group(predictions, "wt_male_female_vs_ko_male", "M", "M"))
        rows.extend(self._rows_for_group(predictions, "wt_male_female_vs_ko_male", "F", "F"))
        return rows


class ExpandedCohortEvaluator(object):
    """Create sex/genotype group metrics and cohort-benchmark comparisons."""

    GROUPS = [
        ("WT_M", "WT", "M"),
        ("WT_F", "WT", "F"),
        ("KO_M", "KO", "M"),
    ]

    @staticmethod
    def _group_row(predictions, label, genotype, sex):
        rows = [row for row in predictions if row["true_genotype"] == genotype and row.get("sex") == sex]
        called = [row for row in rows if row["predicted_genotype"] in ("WT", "KO")]
        correct = sum(row["predicted_genotype"] == genotype for row in called)
        false_ko = sum(row["predicted_genotype"] == "KO" for row in called) if genotype == "WT" else 0
        accuracy = correct / len(called) if called else None
        return {
            "evaluation_group": label,
            "true_genotype": genotype,
            "sex": sex,
            "total_cells": len(rows),
            "called_cells": len(called),
            "uncalled_cells": len(rows) - len(called),
            "percent_called": 100.0 * len(called) / len(rows) if rows else None,
            "accuracy_among_called": accuracy,
            "ko_sensitivity_where_defined": accuracy if genotype == "KO" else None,
            "ko_specificity_where_defined": accuracy if genotype == "WT" else None,
            "false_ko_cells": false_ko,
        }

    def group_metrics(self, predictions):
        return [self._group_row(predictions, *definition) for definition in self.GROUPS]

    @staticmethod
    def _overall_dict(overall_rows):
        return {row["metric"]: float(row["value"]) for row in overall_rows}

    @staticmethod
    def _specificity(predictions, sex=None):
        rows = [
            row for row in predictions
            if row["true_genotype"] == "WT"
            and row["predicted_genotype"] in ("WT", "KO")
            and (sex is None or row.get("sex") == sex)
        ]
        return sum(row["predicted_genotype"] == "WT" for row in rows) / len(rows) if rows else None

    @staticmethod
    def _sensitivity(predictions):
        rows = [row for row in predictions if row["true_genotype"] == "KO" and row["predicted_genotype"] in ("WT", "KO")]
        return sum(row["predicted_genotype"] == "KO" for row in rows) / len(rows) if rows else None

    def comparison(self, male_predictions, male_overall, expanded_predictions, expanded_overall_rows):
        expanded_overall = self._overall_dict(expanded_overall_rows)
        definitions = [
            ("percent_called_overall", male_overall["percent_called"], expanded_overall["percent_called"], "all cohort cells"),
            ("accuracy_among_called_overall", male_overall["accuracy_among_called"], expanded_overall["accuracy_among_called"], "called cells in each cohort"),
            ("ko_male_sensitivity_among_called", self._sensitivity(male_predictions), self._sensitivity(expanded_predictions), "called KO-male cells"),
            ("wt_male_specificity_among_called", self._specificity(male_predictions, "M"), self._specificity(expanded_predictions, "M"), "called WT-male cells; directly comparable"),
            ("all_wt_specificity_among_called", self._specificity(male_predictions), self._specificity(expanded_predictions), "called WT cells; expanded includes WT females"),
            ("wt_female_specificity_among_called", None, self._specificity(expanded_predictions, "F"), "called WT-female cells; no male-only counterpart"),
            ("wt_male_false_ko_cells", float(sum(row["true_genotype"] == "WT" and row["predicted_genotype"] == "KO" for row in male_predictions)), float(sum(row["true_genotype"] == "WT" and row.get("sex") == "M" and row["predicted_genotype"] == "KO" for row in expanded_predictions)), "held-out WT-male false-KO cells"),
        ]
        rows = []
        for metric, male, expanded, denominator in definitions:
            rows.append(
                {
                    "metric": metric,
                    "male_only_value": male,
                    "wt_male_female_ko_male_value": expanded,
                    "expanded_minus_male_only": expanded - male if male is not None and expanded is not None else None,
                    "interpretation_denominator": denominator,
                }
            )
        return rows


class ExpandedFullFitBuilder(object):
    """Fit the expanded full-data reference model without making cell calls."""

    def __init__(self, configuration, encoder):
        self.encoder = encoder
        self.estimator = LogisticRegressionEstimator(configuration, encoder)

    def fit(self, records):
        empirical = EmpiricalPatternEstimator(self.encoder).fit(
            LabeledPatternObservation(record.pattern_code, record.true_genotype) for record in records
        )
        classifier = self.estimator.fit(empirical)
        probability_rows = []
        for pattern in self.encoder.pattern_order:
            probabilities = classifier.predict_proba({"pattern_code": pattern})
            a_value, b_value, c_value = self.encoder.decode(pattern)
            evidence = empirical.evidence(pattern)
            probability_rows.append(
                {
                    "pattern_code": pattern,
                    "pattern_label": self.encoder.label(pattern),
                    "A_detected": a_value,
                    "B_detected": b_value,
                    "C_detected": c_value,
                    "total_ground_truth_cells": evidence["total_cells"],
                    "linear_predictor_log_odds_ko": classifier.linear_predictor({"pattern_code": pattern}),
                    "p_wt_logistic": probabilities["WT"],
                    "p_ko_logistic": probabilities["KO"],
                    "hard_call": None,
                }
            )
        return classifier, probability_rows


class LogisticRegressionPlotter(object):
    """Render Step 05 model diagnostics without fitting or reading cells."""

    WT_COLOR = "#3B6FB6"
    KO_COLOR = "#D95F43"
    EMPIRICAL_COLOR = "#555555"
    LOGISTIC_COLOR = "#8E5BA6"

    def __init__(self, dpi):
        self.dpi = int(dpi)
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
    def _labels(rows):
        return [
            "{}\n{}".format(row["pattern_code"], row["pattern_label"])
            for row in rows
        ]

    def plot_probabilities(self, rows, path):
        x_values = list(range(len(rows)))
        width = 0.38
        figure, axis = plt.subplots(figsize=(10, 5.5))
        axis.bar(
            [value - width / 2 for value in x_values],
            [100.0 * row["p_wt_logistic"] for row in rows],
            width,
            color=self.WT_COLOR,
            label="Logistic P(WT)",
        )
        axis.bar(
            [value + width / 2 for value in x_values],
            [100.0 * row["p_ko_logistic"] for row in rows],
            width,
            color=self.KO_COLOR,
            label="Logistic P(KO)",
        )
        axis.set_xticks(x_values)
        axis.set_xticklabels(self._labels(rows))
        axis.set_ylim(0, 105)
        axis.set_ylabel("Predicted probability (%)")
        axis.set_xlabel("Binary Probe A/B/C pattern")
        axis.set_title("Step 05: Logistic-regression genotype probability by pattern")
        axis.legend(frameon=False)
        axis.grid(axis="y", color="#DDDDDD", linewidth=0.7)
        figure.tight_layout()
        figure.savefig(path, dpi=self.dpi, metadata={"Title": "Step 05 logistic probabilities"})
        plt.close(figure)

    def plot_comparison(self, rows, path):
        x_values = list(range(len(rows)))
        figure, axes = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
        for axis, class_name in zip(axes, ("WT", "KO")):
            empirical_key = "p_{}_empirical_step_04".format(class_name.lower())
            logistic_key = "p_{}_logistic_step_05".format(class_name.lower())
            axis.plot(
                x_values,
                [100.0 * row[empirical_key] for row in rows],
                marker="o",
                color=self.EMPIRICAL_COLOR,
                label="Step 04 empirical",
            )
            axis.plot(
                x_values,
                [100.0 * row[logistic_key] for row in rows],
                marker="s",
                linestyle="--",
                color=self.LOGISTIC_COLOR,
                label="Step 05 logistic",
            )
            axis.set_ylim(0, 105)
            axis.set_ylabel("P({}) (%)".format(class_name))
            axis.grid(axis="y", color="#DDDDDD", linewidth=0.7)
            axis.legend(frameon=False, loc="best")
        axes[0].set_title("Step 05: Empirical versus logistic probability by pattern")
        axes[1].set_xticks(x_values)
        axes[1].set_xticklabels(self._labels(rows))
        axes[1].set_xlabel("Binary Probe A/B/C pattern")
        figure.tight_layout()
        figure.savefig(path, dpi=self.dpi, metadata={"Title": "Step 05 empirical logistic comparison"})
        plt.close(figure)

    def plot_coefficients(self, rows, path):
        terms = [row["term"] for row in rows]
        coefficients = [row["coefficient"] for row in rows]
        odds_ratios = [row["odds_ratio"] for row in rows]
        colors = [self.KO_COLOR if value > 0 else self.WT_COLOR for value in coefficients]
        figure, axes = plt.subplots(1, 2, figsize=(11, 5.5))
        axes[0].barh(terms, coefficients, color=colors)
        axes[0].axvline(0.0, color="black", linewidth=0.8)
        axes[0].set_xlabel("Coefficient: change in KO log odds")
        axes[0].set_title("Logistic coefficients")
        axes[0].grid(axis="x", color="#DDDDDD", linewidth=0.7)
        axes[1].barh(terms, odds_ratios, color=colors)
        axes[1].axvline(1.0, color="black", linewidth=0.8)
        axes[1].set_xscale("log")
        axes[1].set_xlabel("KO:WT odds ratio (log scale)")
        axes[1].set_title("Exponentiated coefficients")
        axes[1].grid(axis="x", color="#DDDDDD", linewidth=0.7, which="both")
        figure.suptitle("Step 05: Logistic-regression coefficients and odds ratios")
        figure.tight_layout()
        figure.savefig(path, dpi=self.dpi, metadata={"Title": "Step 05 logistic coefficients"})
        plt.close(figure)

    def render_all(self, probability_rows, comparison_rows, coefficient_rows, output_directory):
        paths = [
            os.path.join(output_directory, PROBABILITY_PLOT_NAME),
            os.path.join(output_directory, COMPARISON_PLOT_NAME),
            os.path.join(output_directory, COEFFICIENT_PLOT_NAME),
        ]
        self.plot_probabilities(probability_rows, paths[0])
        self.plot_comparison(comparison_rows, paths[1])
        self.plot_coefficients(coefficient_rows, paths[2])
        return paths


class HeldOutValidationPlotter(object):
    """Render only the three plots needed to assess sample-held-out behavior."""

    WT_COLOR = "#3B6FB6"
    KO_COLOR = "#D95F43"
    BAR_COLOR = "#667788"

    def __init__(self, dpi):
        self.dpi = int(dpi)

    def plot_confusion(self, confusion_rows, path):
        matrix = np.asarray(
            [
                [confusion_rows[0]["predicted_wt_cells"], confusion_rows[0]["predicted_ko_cells"]],
                [confusion_rows[1]["predicted_wt_cells"], confusion_rows[1]["predicted_ko_cells"]],
            ],
            dtype=int,
        )
        figure, axis = plt.subplots(figsize=(6.5, 5.5))
        image = axis.imshow(matrix, cmap="Blues")
        for row_index in range(2):
            for column_index in range(2):
                value = int(matrix[row_index, column_index])
                axis.text(
                    column_index,
                    row_index,
                    "{:,}".format(value),
                    ha="center",
                    va="center",
                    color="white" if value > matrix.max() / 2 else "black",
                    fontsize=12,
                )
        axis.set_xticks([0, 1])
        axis.set_xticklabels(["Predicted WT", "Predicted KO"])
        axis.set_yticks([0, 1])
        axis.set_yticklabels(["True WT", "True KO"])
        axis.set_title("Step 05 LOSO: Held-out confusion matrix (called cells)")
        figure.colorbar(image, ax=axis, label="Held-out cells")
        figure.tight_layout()
        figure.savefig(path, dpi=self.dpi, metadata={"Title": "Step 05 LOSO confusion matrix"})
        plt.close(figure)

    @staticmethod
    def _sample_colors(rows):
        return ["#3B6FB6" if row["true_genotype"] == "WT" else "#D95F43" for row in rows]

    def plot_per_sample_accuracy(self, rows, path):
        x_values = list(range(len(rows)))
        values = [100.0 * row["accuracy_among_called"] for row in rows]
        figure, axis = plt.subplots(figsize=(9, 5.5))
        bars = axis.bar(x_values, values, color=self._sample_colors(rows))
        axis.set_xticks(x_values)
        axis.set_xticklabels([row["biological_sample_id"] for row in rows], rotation=25, ha="right")
        axis.set_ylim(0, 105)
        axis.set_ylabel("Accuracy among called cells (%)")
        axis.set_title("Step 05 LOSO: Per-sample held-out accuracy")
        axis.grid(axis="y", color="#DDDDDD", linewidth=0.7)
        for bar, value in zip(bars, values):
            axis.text(bar.get_x() + bar.get_width() / 2, value + 1.0, "{:.1f}%".format(value), ha="center", va="bottom", fontsize=8)
        figure.tight_layout()
        figure.savefig(path, dpi=self.dpi, metadata={"Title": "Step 05 LOSO per-sample accuracy"})
        plt.close(figure)

    def plot_per_sample_called(self, rows, path):
        x_values = list(range(len(rows)))
        values = [row["percent_called"] for row in rows]
        figure, axis = plt.subplots(figsize=(9, 5.5))
        bars = axis.bar(x_values, values, color=self._sample_colors(rows))
        axis.set_xticks(x_values)
        axis.set_xticklabels([row["biological_sample_id"] for row in rows], rotation=25, ha="right")
        axis.set_ylim(0, max(values) * 1.2)
        axis.set_ylabel("Held-out cells called (%)")
        axis.set_title("Step 05 LOSO: Per-sample call rate (000 remains uncalled)")
        axis.grid(axis="y", color="#DDDDDD", linewidth=0.7)
        for bar, value in zip(bars, values):
            axis.text(bar.get_x() + bar.get_width() / 2, value + 0.3, "{:.1f}%".format(value), ha="center", va="bottom", fontsize=8)
        figure.tight_layout()
        figure.savefig(path, dpi=self.dpi, metadata={"Title": "Step 05 LOSO per-sample percent called"})
        plt.close(figure)

    def render_all(self, confusion_rows, per_sample_rows, output_directory):
        paths = [
            os.path.join(output_directory, HELD_OUT_CONFUSION_PLOT_NAME),
            os.path.join(output_directory, HELD_OUT_ACCURACY_PLOT_NAME),
            os.path.join(output_directory, HELD_OUT_CALLED_PLOT_NAME),
        ]
        self.plot_confusion(confusion_rows, paths[0])
        self.plot_per_sample_accuracy(per_sample_rows, paths[1])
        self.plot_per_sample_called(per_sample_rows, paths[2])
        return paths


class ExpandedValidationPlotter(HeldOutValidationPlotter):
    """Render expanded-cohort generalization and cohort-comparison diagnostics."""

    def plot_pattern_errors(self, male_rows, expanded_rows, path):
        male = {row["pattern_code"]: row for row in male_rows if row["wt_sex"] == "M"}
        expanded = {
            row["pattern_code"]: row
            for row in expanded_rows
            if row["wt_sex"] == "M"
        }
        patterns = PatternErrorAnalyzer.PATTERNS
        x_values = np.arange(len(patterns))
        width = 0.38
        figure, axis = plt.subplots(figsize=(10, 5.5))
        axis.bar(
            x_values - width / 2,
            [male[pattern]["predicted_ko_cells_false_ko"] for pattern in patterns],
            width,
            label="Male-only training",
            color="#667788",
        )
        axis.bar(
            x_values + width / 2,
            [expanded[pattern]["predicted_ko_cells_false_ko"] for pattern in patterns],
            width,
            label="WT-M+F / KO-M training",
            color="#8E5BA6",
        )
        axis.set_xticks(x_values)
        axis.set_xticklabels(patterns)
        axis.set_xlabel("Binary Probe A/B/C pattern")
        axis.set_ylabel("Held-out WT-male cells called KO")
        axis.set_title("Step 05: WT-male false-KO errors by probe pattern")
        axis.grid(axis="y", color="#DDDDDD", linewidth=0.7)
        axis.legend(frameon=False)
        figure.tight_layout()
        figure.savefig(path, dpi=self.dpi, metadata={"Title": "Step 05 WT false-KO pattern errors"})
        plt.close(figure)

    def plot_cohort_comparison(self, rows, path):
        wanted = [
            "wt_male_specificity_among_called",
            "ko_male_sensitivity_among_called",
            "accuracy_among_called_overall",
        ]
        indexed = {row["metric"]: row for row in rows}
        x_values = np.arange(len(wanted))
        width = 0.38
        figure, axis = plt.subplots(figsize=(9, 5.5))
        axis.bar(
            x_values - width / 2,
            [100.0 * indexed[key]["male_only_value"] for key in wanted],
            width,
            label="Male-only",
            color="#667788",
        )
        axis.bar(
            x_values + width / 2,
            [100.0 * indexed[key]["wt_male_female_ko_male_value"] for key in wanted],
            width,
            label="WT-M+F / KO-M",
            color="#8E5BA6",
        )
        axis.set_xticks(x_values)
        axis.set_xticklabels(["WT-M specificity", "KO-M sensitivity", "Overall accuracy"])
        axis.set_ylim(0, 105)
        axis.set_ylabel("Called-cell metric (%)")
        axis.set_title("Step 05: Ground-truth cohort validation comparison")
        axis.grid(axis="y", color="#DDDDDD", linewidth=0.7)
        axis.legend(frameon=False)
        figure.tight_layout()
        figure.savefig(path, dpi=self.dpi, metadata={"Title": "Step 05 cohort validation comparison"})
        plt.close(figure)

    def render_expanded(
        self,
        confusion_rows,
        per_sample_rows,
        male_pattern_rows,
        expanded_pattern_rows,
        comparison_rows,
        output_directory,
    ):
        paths = [
            os.path.join(output_directory, EXPANDED_CONFUSION_PLOT_NAME),
            os.path.join(output_directory, EXPANDED_ACCURACY_PLOT_NAME),
            os.path.join(output_directory, EXPANDED_CALLED_PLOT_NAME),
            os.path.join(output_directory, PATTERN_ERROR_PLOT_NAME),
            os.path.join(output_directory, COHORT_COMPARISON_PLOT_NAME),
        ]
        self.plot_confusion(confusion_rows, paths[0])
        self.plot_per_sample_accuracy(per_sample_rows, paths[1])
        self.plot_per_sample_called(per_sample_rows, paths[2])
        self.plot_pattern_errors(male_pattern_rows, expanded_pattern_rows, paths[3])
        self.plot_cohort_comparison(comparison_rows, paths[4])
        return paths


class Step05OutputPublisher(object):
    """Own deterministic Step 05 serialization, plots, and atomic publication."""

    EXPECTED_FILES = {
        COEFFICIENT_TABLE_NAME,
        PROBABILITY_TABLE_NAME,
        COMPARISON_TABLE_NAME,
        DIAGNOSTICS_NAME,
        VALIDATION_NAME,
        ENVIRONMENT_NAME,
        PROBABILITY_PLOT_NAME,
        COMPARISON_PLOT_NAME,
        COEFFICIENT_PLOT_NAME,
    }

    def __init__(self, configuration, validation, plotter, output_root, paths):
        self.configuration = configuration
        self.validation = validation
        self.plotter = plotter
        self.output_root = os.path.abspath(output_root)
        self.paths = {key: os.path.abspath(value) for key, value in paths.items()}

    def _format_float(self, value):
        return ("{:.%df}" % self.configuration.decimal_places).format(value)

    def _serialize_rows(self, rows, float_keys, blank_hard_call=False):
        serialized = []
        for source in rows:
            row = dict(source)
            for key in float_keys:
                row[key] = self._format_float(row[key])
            if blank_hard_call:
                row["hard_call"] = ""
            serialized.append(row)
        return serialized

    def _diagnostic_rows(self, classifier, metrics):
        baseline = classifier.predict_proba({"pattern_code": "000"})
        values = {
            "model_name": classifier.model_name,
            "class_encoding": "WT=0;KO=1",
            "positive_model_class": "KO=1",
            "model_formula": self.configuration.values["model_formula"],
            "predictors": ",".join(self.configuration.values["predictors"]),
            "penalty": "none",
            "interaction_terms": "none",
            "nonlinear_transformations": "none",
            "fit_representation": "8 grouped binomial patterns; equivalent to 230269 cell rows",
            "converged": str(classifier.fit_diagnostics["converged"]).lower(),
            "iterations": classifier.fit_diagnostics["iterations"],
            "step_halvings": classifier.fit_diagnostics["step_halvings"],
            "log_likelihood": self._format_float(classifier.fit_diagnostics["log_likelihood"]),
            "maximum_absolute_score": self._format_float(classifier.fit_diagnostics["maximum_absolute_score"]),
            "information_matrix_rank": classifier.fit_diagnostics["information_matrix_rank"],
            "information_matrix_condition_number": self._format_float(classifier.fit_diagnostics["information_matrix_condition_number"]),
            "fit_cells": classifier.fit_diagnostics["fit_cells"],
            "fit_patterns": classifier.fit_diagnostics["fit_patterns"],
            "pattern_000_included": "yes",
            "p_wt_logistic_000": self._format_float(baseline["WT"]),
            "p_ko_logistic_000": self._format_float(baseline["KO"]),
            "pattern_000_hard_call": "none",
            "unweighted_pattern_rmse_p_ko": self._format_float(metrics["unweighted_pattern_rmse_p_ko"]),
            "cell_weighted_rmse_p_ko": self._format_float(metrics["cell_weighted_rmse_p_ko"]),
            "maximum_absolute_pattern_difference_p_ko": self._format_float(metrics["maximum_absolute_pattern_difference_p_ko"]),
            "held_out_validation": "not_performed",
            "het_cells": "not_loaded",
        }
        return [{"key": key, "value": value} for key, value in values.items()]

    def _write_environment(self, path):
        rows = [
            {"key": "step_id", "value": self.configuration.step_id},
            {"key": "pipeline_version", "value": self.configuration.values["pipeline_version"]},
            {"key": "python_version", "value": platform.python_version()},
            {"key": "platform", "value": platform.platform()},
            {"key": "numpy_version", "value": np.__version__},
            {"key": "matplotlib_version", "value": matplotlib.__version__},
        ]
        for key in sorted(self.paths):
            rows.append({"key": "{}_path".format(key), "value": self.paths[key]})
            rows.append({"key": "{}_sha256".format(key), "value": sha256_file(self.paths[key])})
        write_tsv(path, ["key", "value"], rows)

    @staticmethod
    def _write_manifest(directory):
        rows = []
        for filename in sorted(os.listdir(directory)):
            path = os.path.join(directory, filename)
            if filename == MANIFEST_NAME or not os.path.isfile(path):
                continue
            rows.append({"relative_path": filename, "bytes": os.path.getsize(path), "sha256": sha256_file(path)})
        write_tsv(os.path.join(directory, MANIFEST_NAME), ["relative_path", "bytes", "sha256"], rows)

    @classmethod
    def _verify_manifest(cls, directory):
        manifest_path = os.path.join(directory, MANIFEST_NAME)
        if not os.path.isfile(manifest_path):
            raise Step05Error("Step 05 output manifest is missing")
        listed = set()
        with open(manifest_path, "r", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if reader.fieldnames != ["relative_path", "bytes", "sha256"]:
                raise Step05Error("Unexpected Step 05 manifest schema")
            for row in reader:
                if row["relative_path"] in listed:
                    raise Step05Error("Duplicate Step 05 manifest path")
                listed.add(row["relative_path"])
                path = os.path.join(directory, row["relative_path"])
                if not os.path.isfile(path):
                    raise Step05Error("Manifest-listed Step 05 file is missing")
                if os.path.getsize(path) != int(row["bytes"]) or sha256_file(path) != row["sha256"]:
                    raise Step05Error("Step 05 output identity mismatch: {}".format(path))
        if listed != cls.EXPECTED_FILES:
            raise Step05Error("Step 05 manifest file set is incomplete or unexpected")

    def _verify_existing_provenance(self):
        self._verify_manifest(self.output_root)
        with open(os.path.join(self.output_root, ENVIRONMENT_NAME), "r", newline="") as handle:
            environment = {row["key"]: row["value"] for row in csv.DictReader(handle, delimiter="\t")}
        for key, path in self.paths.items():
            if environment.get("{}_sha256".format(key)) != sha256_file(path):
                raise Step05Error("Existing Step 05 provenance differs for {}".format(key))

    def publish(self, classifier, probability_rows, comparison_rows, metrics):
        if os.path.exists(self.output_root):
            self._verify_existing_provenance()
            print("EXISTING_VALIDATED\t{}".format(self.output_root))
            return
        parent = os.path.dirname(self.output_root)
        os.makedirs(parent, exist_ok=True)
        staging = tempfile.mkdtemp(prefix=".step_05_pcdh19_logistic_regression_baseline.", dir=parent)
        try:
            coefficient_rows = classifier.coefficient_rows()
            write_tsv(
                os.path.join(staging, COEFFICIENT_TABLE_NAME),
                COEFFICIENT_HEADER,
                self._serialize_rows(coefficient_rows, {"coefficient", "odds_ratio"}),
            )
            write_tsv(
                os.path.join(staging, PROBABILITY_TABLE_NAME),
                PROBABILITY_HEADER,
                self._serialize_rows(
                    probability_rows,
                    {"linear_predictor_log_odds_ko", "p_wt_logistic", "p_ko_logistic"},
                    blank_hard_call=True,
                ),
            )
            write_tsv(
                os.path.join(staging, COMPARISON_TABLE_NAME),
                COMPARISON_HEADER,
                self._serialize_rows(
                    comparison_rows,
                    {
                        "p_wt_empirical_step_04", "p_wt_logistic_step_05",
                        "logistic_minus_empirical_p_wt", "absolute_difference_p_wt",
                        "p_ko_empirical_step_04", "p_ko_logistic_step_05",
                        "logistic_minus_empirical_p_ko", "absolute_difference_p_ko",
                    },
                    blank_hard_call=True,
                ),
            )
            write_tsv(os.path.join(staging, DIAGNOSTICS_NAME), ["key", "value"], self._diagnostic_rows(classifier, metrics))
            plot_paths = self.plotter.render_all(probability_rows, comparison_rows, coefficient_rows, staging)
            self.validation.require_equal("diagnostic_plots_created", sum(os.path.getsize(path) > 0 for path in plot_paths), 3)
            self.validation.require_equal("logistic_regression_converged", classifier.fit_diagnostics["converged"], True)
            self.validation.require_equal("information_matrix_full_rank", classifier.fit_diagnostics["information_matrix_rank"], 4)
            write_tsv(os.path.join(staging, VALIDATION_NAME), VALIDATION_HEADER, self.validation.rows)
            self._write_environment(os.path.join(staging, ENVIRONMENT_NAME))
            self._write_manifest(staging)
            self._verify_manifest(staging)
            if os.path.exists(self.output_root):
                raise Step05Error("Step 05 output appeared during publication")
            os.replace(staging, self.output_root)
            staging = None
            print("PUBLISHED\t{}\tpatterns=8\tcells={}".format(self.output_root, classifier.fit_diagnostics["fit_cells"]))
        finally:
            if staging and os.path.exists(staging):
                shutil.rmtree(staging)


class ExistingStep05BasePackageVerifier(object):
    """Protect the already published Step 05 fit while adding validation."""

    def __init__(self, configuration, base_root, validation):
        self.configuration = configuration
        self.base_root = os.path.abspath(base_root)
        self.validation = validation
        self.manifest_path = os.path.join(self.base_root, MANIFEST_NAME)

    def verify(self):
        if not os.path.isfile(self.manifest_path):
            raise Step05Error("Existing Step 05 base manifest is missing")
        base_config = self.configuration.values["existing_base_package"]
        observed_manifest_sha = sha256_file(self.manifest_path)
        base_environment_path = os.path.join(self.base_root, ENVIRONMENT_NAME)
        if not os.path.isfile(base_environment_path):
            raise Step05Error("Existing Step 05 base environment is missing")
        with open(base_environment_path, "r", newline="") as handle:
            base_environment = {
                row["key"]: row["value"]
                for row in csv.DictReader(handle, delimiter="\t")
            }
        if base_environment.get("pipeline_version") == "0.1.0":
            self.validation.require_equal(
                "historical_step_05_base_manifest_unchanged",
                observed_manifest_sha,
                base_config["output_manifest_sha256"],
            )
        entries = {}
        with open(self.manifest_path, "r", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for row in reader:
                path = os.path.join(self.base_root, row["relative_path"])
                if not os.path.isfile(path):
                    raise Step05Error("Existing Step 05 base artifact is missing")
                if os.path.getsize(path) != int(row["bytes"]) or sha256_file(path) != row["sha256"]:
                    raise Step05Error("Existing Step 05 base artifact identity mismatch")
                entries[row["relative_path"]] = row["sha256"]
        for filename, expected_sha in base_config["core_artifacts"].items():
            self.validation.require_equal(
                "existing_base_{}_sha256".format(filename),
                entries.get(filename),
                expected_sha,
            )
        self.validation.require_equal(
            "existing_step_05_base_package_accepted",
            True,
            True,
            "historical_manifest={} observed_manifest={}; core model artifacts are exact".format(
                base_config["output_manifest_sha256"], observed_manifest_sha
            ),
        )


class HeldOutValidationOutputPublisher(object):
    """Atomically publish the LOSO extension without altering base outputs."""

    EXPECTED_FILES = {
        HELD_OUT_PREDICTIONS_NAME,
        HELD_OUT_CONFUSION_NAME,
        HELD_OUT_PER_SAMPLE_NAME,
        HELD_OUT_OVERALL_NAME,
        HELD_OUT_FOLD_COEFFICIENTS_NAME,
        HELD_OUT_CHECKS_NAME,
        HELD_OUT_ENVIRONMENT_NAME,
        HELD_OUT_CONFUSION_PLOT_NAME,
        HELD_OUT_ACCURACY_PLOT_NAME,
        HELD_OUT_CALLED_PLOT_NAME,
    }

    def __init__(self, configuration, validation, plotter, output_root, paths):
        self.configuration = configuration
        self.validation = validation
        self.plotter = plotter
        self.output_root = os.path.abspath(output_root)
        self.paths = {key: os.path.abspath(value) for key, value in paths.items()}

    def _format_float(self, value):
        if value is None:
            return ""
        return ("{:.%df}" % self.configuration.decimal_places).format(value)

    def _prediction_rows(self, predictions):
        rows = []
        for source in predictions:
            row = dict(source)
            row["predicted_wt_probability"] = self._format_float(row["predicted_wt_probability"])
            row["predicted_ko_probability"] = self._format_float(row["predicted_ko_probability"])
            rows.append(row)
        return rows

    def _per_sample_rows(self, metrics):
        rows = []
        for source in metrics:
            row = dict(source)
            for key in [
                "percent_called",
                "accuracy_among_called",
                "ko_sensitivity_where_defined",
                "ko_specificity_where_defined",
            ]:
                row[key] = self._format_float(row[key])
            rows.append(row)
        return rows

    def _overall_rows(self, metrics):
        rows = []
        for source in metrics:
            row = dict(source)
            if isinstance(row["value"], float):
                row["value"] = self._format_float(row["value"])
            rows.append(row)
        return rows

    def _fold_rows(self, rows):
        serialized = []
        for source in rows:
            row = dict(source)
            row["coefficient"] = self._format_float(row["coefficient"])
            row["odds_ratio"] = self._format_float(row["odds_ratio"])
            row["converged"] = str(row["converged"]).lower()
            serialized.append(row)
        return serialized

    def _write_environment(self, path):
        calling = self.configuration.values["sample_level_validation"]["calling_rule"]
        rows = [
            {"key": "step_id", "value": self.configuration.step_id},
            {"key": "pipeline_version", "value": self.configuration.values["pipeline_version"]},
            {"key": "validation_method", "value": "leave_one_registered_sample_out"},
            {"key": "holdout_unit_field", "value": "technical_sample_id"},
            {"key": "holdout_unit_interpretation", "value": self.configuration.values["sample_level_validation"]["holdout_unit_interpretation"]},
            {"key": "class_encoding", "value": "WT=0;KO=1"},
            {"key": "probability_threshold", "value": calling["probability_threshold"]},
            {"key": "uncalled_patterns", "value": ",".join(calling["uncalled_pattern_codes"])},
            {"key": "tie_behavior", "value": calling["tie_behavior"]},
            {"key": "threshold_optimized", "value": str(calling["threshold_optimized"]).lower()},
            {"key": "python_version", "value": platform.python_version()},
            {"key": "numpy_version", "value": np.__version__},
            {"key": "matplotlib_version", "value": matplotlib.__version__},
        ]
        for key in sorted(self.paths):
            rows.append({"key": "{}_path".format(key), "value": self.paths[key]})
            rows.append({"key": "{}_sha256".format(key), "value": sha256_file(self.paths[key])})
        write_tsv(path, ["key", "value"], rows)

    @staticmethod
    def _write_manifest(directory):
        rows = []
        for filename in sorted(os.listdir(directory)):
            path = os.path.join(directory, filename)
            if filename == HELD_OUT_MANIFEST_NAME or not os.path.isfile(path):
                continue
            rows.append({"relative_path": filename, "bytes": os.path.getsize(path), "sha256": sha256_file(path)})
        write_tsv(os.path.join(directory, HELD_OUT_MANIFEST_NAME), ["relative_path", "bytes", "sha256"], rows)

    @classmethod
    def _verify_manifest(cls, directory):
        manifest_path = os.path.join(directory, HELD_OUT_MANIFEST_NAME)
        if not os.path.isfile(manifest_path):
            raise Step05Error("Held-out validation manifest is missing")
        listed = set()
        with open(manifest_path, "r", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for row in reader:
                path = os.path.join(directory, row["relative_path"])
                if row["relative_path"] in listed or not os.path.isfile(path):
                    raise Step05Error("Held-out validation manifest is invalid")
                listed.add(row["relative_path"])
                if os.path.getsize(path) != int(row["bytes"]) or sha256_file(path) != row["sha256"]:
                    raise Step05Error("Held-out validation artifact identity mismatch")
        if listed != cls.EXPECTED_FILES:
            raise Step05Error("Held-out validation manifest file set is unexpected")

    def _verify_existing_provenance(self):
        self._verify_manifest(self.output_root)
        with open(os.path.join(self.output_root, HELD_OUT_ENVIRONMENT_NAME), "r", newline="") as handle:
            environment = {row["key"]: row["value"] for row in csv.DictReader(handle, delimiter="\t")}
        for key, path in self.paths.items():
            if environment.get("{}_sha256".format(key)) != sha256_file(path):
                raise Step05Error("Held-out validation provenance differs for {}".format(key))

    def publish(self, predictions, confusion_rows, per_sample_rows, overall_rows, fold_rows):
        if os.path.exists(self.output_root):
            self._verify_existing_provenance()
            print("EXISTING_VALIDATED\t{}".format(self.output_root))
            return
        parent = os.path.dirname(self.output_root)
        staging = tempfile.mkdtemp(prefix=".sample_level_held_out_validation.", dir=parent)
        try:
            write_tsv(os.path.join(staging, HELD_OUT_PREDICTIONS_NAME), HELD_OUT_PREDICTION_HEADER, self._prediction_rows(predictions))
            write_tsv(os.path.join(staging, HELD_OUT_CONFUSION_NAME), HELD_OUT_CONFUSION_HEADER, confusion_rows)
            write_tsv(os.path.join(staging, HELD_OUT_PER_SAMPLE_NAME), HELD_OUT_PER_SAMPLE_HEADER, self._per_sample_rows(per_sample_rows))
            write_tsv(os.path.join(staging, HELD_OUT_OVERALL_NAME), HELD_OUT_OVERALL_HEADER, self._overall_rows(overall_rows))
            write_tsv(os.path.join(staging, HELD_OUT_FOLD_COEFFICIENTS_NAME), HELD_OUT_FOLD_COEFFICIENT_HEADER, self._fold_rows(fold_rows))
            plot_paths = self.plotter.render_all(confusion_rows, per_sample_rows, staging)
            self.validation.require_equal("held_out_validation_plots_created", sum(os.path.getsize(path) > 0 for path in plot_paths), 3)
            write_tsv(os.path.join(staging, HELD_OUT_CHECKS_NAME), VALIDATION_HEADER, self.validation.rows)
            self._write_environment(os.path.join(staging, HELD_OUT_ENVIRONMENT_NAME))
            self._write_manifest(staging)
            self._verify_manifest(staging)
            if os.path.exists(self.output_root):
                raise Step05Error("Held-out validation output appeared during publication")
            os.replace(staging, self.output_root)
            staging = None
            print("PUBLISHED\t{}\tfolds=6\tcells={}".format(self.output_root, len(predictions)))
        finally:
            if staging and os.path.exists(staging):
                shutil.rmtree(staging)


class ExpandedCohortOutputPublisher(object):
    """Publish the expanded Step 05 cohort comparison as a separate package."""

    EXPECTED_FILES = {
        EXPANDED_PREDICTIONS_NAME,
        EXPANDED_CONFUSION_NAME,
        EXPANDED_PER_SAMPLE_NAME,
        EXPANDED_GROUP_METRICS_NAME,
        EXPANDED_OVERALL_NAME,
        EXPANDED_FOLD_COEFFICIENTS_NAME,
        EXPANDED_FULL_COEFFICIENTS_NAME,
        EXPANDED_FULL_PROBABILITIES_NAME,
        MALE_PATTERN_ERRORS_NAME,
        EXPANDED_PATTERN_ERRORS_NAME,
        COHORT_COMPARISON_NAME,
        EXPANDED_CHECKS_NAME,
        EXPANDED_ENVIRONMENT_NAME,
        EXPANDED_CONFUSION_PLOT_NAME,
        EXPANDED_ACCURACY_PLOT_NAME,
        EXPANDED_CALLED_PLOT_NAME,
        PATTERN_ERROR_PLOT_NAME,
        COHORT_COMPARISON_PLOT_NAME,
    }

    def __init__(self, configuration, validation, plotter, output_root, paths):
        self.configuration = configuration
        self.validation = validation
        self.plotter = plotter
        self.output_root = os.path.abspath(output_root)
        self.paths = {key: os.path.abspath(value) for key, value in paths.items()}

    def _format(self, value):
        if value is None or value == "":
            return ""
        return ("{:.%df}" % self.configuration.decimal_places).format(float(value))

    def _serialize(self, rows, float_keys, blank_hard_call=False):
        serialized = []
        for source in rows:
            row = dict(source)
            for key in float_keys:
                row[key] = self._format(row.get(key))
            if blank_hard_call:
                row["hard_call"] = ""
            serialized.append(row)
        return serialized

    def _overall_rows(self, rows):
        serialized = []
        for source in rows:
            row = dict(source)
            if isinstance(row["value"], float):
                row["value"] = self._format(row["value"])
            serialized.append(row)
        return serialized

    def _fold_rows(self, rows):
        serialized = self._serialize(rows, {"coefficient", "odds_ratio"})
        for row in serialized:
            row["converged"] = str(row["converged"]).lower()
        return serialized

    def _write_environment(self, path):
        cohort = self.configuration.values["expanded_ground_truth_validation"]
        rows = [
            {"key": "step_id", "value": self.configuration.step_id},
            {"key": "pipeline_version", "value": self.configuration.values["pipeline_version"]},
            {"key": "cohort_id", "value": cohort["cohort_id"]},
            {"key": "ground_truth_definition", "value": "WT=WT_M+WT_F;KO=KO_M"},
            {"key": "class_encoding", "value": "WT=0;KO=1"},
            {"key": "het_female_samples", "value": "excluded_not_loaded"},
            {"key": "validation_method", "value": cohort["method"]},
            {"key": "holdout_unit_field", "value": cohort["holdout_unit_field"]},
            {"key": "probability_threshold", "value": cohort["calling_rule"]["probability_threshold"]},
            {"key": "uncalled_patterns", "value": ",".join(cohort["calling_rule"]["uncalled_pattern_codes"])},
            {"key": "threshold_optimized", "value": "false"},
            {"key": "python_version", "value": platform.python_version()},
            {"key": "numpy_version", "value": np.__version__},
            {"key": "matplotlib_version", "value": matplotlib.__version__},
        ]
        for key in sorted(self.paths):
            rows.append({"key": "{}_path".format(key), "value": self.paths[key]})
            rows.append({"key": "{}_sha256".format(key), "value": sha256_file(self.paths[key])})
        write_tsv(path, ["key", "value"], rows)

    @staticmethod
    def _write_manifest(directory):
        rows = []
        for filename in sorted(os.listdir(directory)):
            path = os.path.join(directory, filename)
            if filename == EXPANDED_MANIFEST_NAME or not os.path.isfile(path):
                continue
            rows.append({"relative_path": filename, "bytes": os.path.getsize(path), "sha256": sha256_file(path)})
        write_tsv(os.path.join(directory, EXPANDED_MANIFEST_NAME), ["relative_path", "bytes", "sha256"], rows)

    @classmethod
    def _verify_manifest(cls, directory):
        manifest_path = os.path.join(directory, EXPANDED_MANIFEST_NAME)
        if not os.path.isfile(manifest_path):
            raise Step05Error("Expanded Step 05 manifest is missing")
        listed = set()
        with open(manifest_path, "r", newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                path = os.path.join(directory, row["relative_path"])
                if row["relative_path"] in listed or not os.path.isfile(path):
                    raise Step05Error("Expanded Step 05 manifest is invalid")
                listed.add(row["relative_path"])
                if os.path.getsize(path) != int(row["bytes"]) or sha256_file(path) != row["sha256"]:
                    raise Step05Error("Expanded Step 05 artifact identity mismatch")
        if listed != cls.EXPECTED_FILES:
            raise Step05Error("Expanded Step 05 manifested file set is unexpected")

    def _verify_existing(self):
        self._verify_manifest(self.output_root)
        with open(os.path.join(self.output_root, EXPANDED_ENVIRONMENT_NAME), "r", newline="") as handle:
            environment = {row["key"]: row["value"] for row in csv.DictReader(handle, delimiter="\t")}
        for key, path in self.paths.items():
            if environment.get("{}_sha256".format(key)) != sha256_file(path):
                raise Step05Error("Expanded Step 05 provenance differs for {}".format(key))

    def publish(
        self,
        predictions,
        confusion_rows,
        per_sample_rows,
        group_rows,
        overall_rows,
        fold_rows,
        full_classifier,
        full_probability_rows,
        male_pattern_rows,
        expanded_pattern_rows,
        comparison_rows,
    ):
        if os.path.exists(self.output_root):
            self._verify_existing()
            print("EXISTING_VALIDATED\t{}".format(self.output_root))
            return
        parent = os.path.dirname(self.output_root)
        staging = tempfile.mkdtemp(prefix=".wt_male_female_vs_ko_male_validation.", dir=parent)
        try:
            write_tsv(
                os.path.join(staging, EXPANDED_PREDICTIONS_NAME),
                EXPANDED_PREDICTION_HEADER,
                self._serialize(predictions, {"predicted_wt_probability", "predicted_ko_probability"}),
            )
            write_tsv(os.path.join(staging, EXPANDED_CONFUSION_NAME), HELD_OUT_CONFUSION_HEADER, confusion_rows)
            write_tsv(
                os.path.join(staging, EXPANDED_PER_SAMPLE_NAME),
                EXPANDED_PER_SAMPLE_HEADER,
                self._serialize(per_sample_rows, {"percent_called", "accuracy_among_called", "ko_sensitivity_where_defined", "ko_specificity_where_defined"}),
            )
            write_tsv(
                os.path.join(staging, EXPANDED_GROUP_METRICS_NAME),
                GROUP_METRICS_HEADER,
                self._serialize(group_rows, {"percent_called", "accuracy_among_called", "ko_sensitivity_where_defined", "ko_specificity_where_defined"}),
            )
            write_tsv(
                os.path.join(staging, EXPANDED_OVERALL_NAME),
                HELD_OUT_OVERALL_HEADER,
                self._overall_rows(overall_rows),
            )
            write_tsv(
                os.path.join(staging, EXPANDED_FOLD_COEFFICIENTS_NAME),
                EXPANDED_FOLD_COEFFICIENT_HEADER,
                self._fold_rows(fold_rows),
            )
            full_coefficients = full_classifier.coefficient_rows()
            write_tsv(
                os.path.join(staging, EXPANDED_FULL_COEFFICIENTS_NAME),
                COEFFICIENT_HEADER,
                self._serialize(full_coefficients, {"coefficient", "odds_ratio"}),
            )
            write_tsv(
                os.path.join(staging, EXPANDED_FULL_PROBABILITIES_NAME),
                PROBABILITY_HEADER,
                self._serialize(full_probability_rows, {"linear_predictor_log_odds_ko", "p_wt_logistic", "p_ko_logistic"}, blank_hard_call=True),
            )
            pattern_float_keys = {"false_ko_rate_among_called", "percent_of_cohort_false_ko_errors"}
            write_tsv(os.path.join(staging, MALE_PATTERN_ERRORS_NAME), PATTERN_ERROR_HEADER, self._serialize(male_pattern_rows, pattern_float_keys))
            write_tsv(os.path.join(staging, EXPANDED_PATTERN_ERRORS_NAME), PATTERN_ERROR_HEADER, self._serialize(expanded_pattern_rows, pattern_float_keys))
            write_tsv(
                os.path.join(staging, COHORT_COMPARISON_NAME),
                COHORT_COMPARISON_HEADER,
                self._serialize(comparison_rows, {"male_only_value", "wt_male_female_ko_male_value", "expanded_minus_male_only"}),
            )
            plot_paths = self.plotter.render_expanded(
                confusion_rows, per_sample_rows, male_pattern_rows, expanded_pattern_rows, comparison_rows, staging
            )
            self.validation.require_equal("expanded_validation_plots_created", sum(os.path.getsize(path) > 0 for path in plot_paths), 5)
            write_tsv(os.path.join(staging, EXPANDED_CHECKS_NAME), VALIDATION_HEADER, self.validation.rows)
            self._write_environment(os.path.join(staging, EXPANDED_ENVIRONMENT_NAME))
            self._write_manifest(staging)
            self._verify_manifest(staging)
            if os.path.exists(self.output_root):
                raise Step05Error("Expanded Step 05 output appeared during publication")
            os.replace(staging, self.output_root)
            staging = None
            print("PUBLISHED\t{}\tfolds=9\tcells={}".format(self.output_root, len(predictions)))
        finally:
            if staging and os.path.exists(staging):
                shutil.rmtree(staging)


class PCDH19LogisticRegressionBaselineStep(object):
    """Coordinate Step 05 using the existing Step 04 classification framework."""

    def __init__(self, args):
        self.configuration = Step05Configuration(args.lock, args.bundle_root)
        self.validation = ValidationRecorder(self.configuration.step_id)
        self.shared_configuration = Step04Configuration(args.step04_lock)
        self.encoder = ProbePatternEncoder(self.configuration.pattern_order)
        self.step03_reader = Step03ClassificationTableReader(
            self.shared_configuration,
            self.encoder,
            self.validation,
            args.step03_root,
        )
        self.step04_reader = Step04EmpiricalModelReader(
            self.configuration, self.validation, args.step04_root
        )
        self.estimator = LogisticRegressionEstimator(self.configuration, self.encoder)
        self.evaluator = LogisticModelEvaluator(self.configuration, self.encoder, self.validation)
        self.plotter = LogisticRegressionPlotter(self.configuration.values["plot_dpi"])
        output_root = os.path.join(os.path.abspath(args.paper3_root), "results", STEP_OUTPUT_DIRECTORY)
        paths = {
            "python_script": __file__,
            "lock": args.lock,
            "requirements": args.requirements,
            "shared_framework_script": self.configuration.framework_script_path,
            "step_04_lock": args.step04_lock,
            "step_03_manifest": self.step03_reader.manifest_path,
            "step_03_table": self.step03_reader.table_path,
            "step_04_manifest": self.step04_reader.manifest_path,
            "step_04_empirical_model": self.step04_reader.model_path,
        }
        self.publisher = Step05OutputPublisher(self.configuration, self.validation, self.plotter, output_root, paths)
        self.base_verifier = ExistingStep05BasePackageVerifier(
            self.configuration, output_root, self.validation
        )
        self.sample_aware_reader = Step03SampleAwarePatternReader(
            self.configuration,
            self.shared_configuration,
            self.encoder,
            self.validation,
            self.step03_reader.table_path,
        )
        self.calling_policy = HeldOutCallingPolicy(
            self.configuration.values["sample_level_validation"]["calling_rule"]
        )
        self.held_out_validator = LeaveOneSampleOutLogisticValidator(
            self.configuration, self.encoder, self.validation, self.calling_policy
        )
        self.held_out_evaluator = HeldOutValidationEvaluator(
            self.configuration, self.validation
        )
        self.held_out_plotter = HeldOutValidationPlotter(
            self.configuration.values["plot_dpi"]
        )
        held_out_root = os.path.join(output_root, HELD_OUT_VALIDATION_DIRECTORY)
        held_out_paths = {
            "python_script": __file__,
            "lock": args.lock,
            "requirements": args.requirements,
            "shared_framework_script": self.configuration.framework_script_path,
            "step_04_lock": args.step04_lock,
            "step_03_manifest": self.step03_reader.manifest_path,
            "step_03_table": self.step03_reader.table_path,
            "step_04_manifest": self.step04_reader.manifest_path,
            "step_04_empirical_model": self.step04_reader.model_path,
            "step_05_base_manifest": os.path.join(output_root, MANIFEST_NAME),
        }
        self.held_out_publisher = HeldOutValidationOutputPublisher(
            self.configuration,
            self.validation,
            self.held_out_plotter,
            held_out_root,
            held_out_paths,
        )
        self.male_only_reader = ExistingMaleOnlyValidationReader(
            self.configuration, self.validation, output_root
        )
        expanded_config = self.configuration.values["expanded_ground_truth_validation"]
        self.expanded_reader = ManifestedProbePatternCohortReader(
            self.configuration,
            self.encoder,
            self.validation,
            args.sample_key,
            args.step02_root,
        )
        self.expanded_calling_policy = HeldOutCallingPolicy(expanded_config["calling_rule"])
        self.expanded_validator = LeaveOneSampleOutLogisticValidator(
            self.configuration,
            self.encoder,
            self.validation,
            self.expanded_calling_policy,
            validation_config=expanded_config,
            check_prefix="expanded_",
        )
        self.expanded_evaluator = HeldOutValidationEvaluator(
            self.configuration,
            self.validation,
            validation_config=expanded_config,
            check_prefix="expanded_",
        )
        self.expanded_group_evaluator = ExpandedCohortEvaluator()
        self.pattern_error_analyzer = PatternErrorAnalyzer()
        self.expanded_full_fit = ExpandedFullFitBuilder(self.configuration, self.encoder)
        self.expanded_plotter = ExpandedValidationPlotter(
            self.configuration.values["plot_dpi"]
        )
        expanded_root = os.path.join(output_root, expanded_config["relative_directory"])
        expanded_paths = {
            "python_script": __file__,
            "lock": args.lock,
            "requirements": args.requirements,
            "shared_framework_script": self.configuration.framework_script_path,
            "step_04_lock": args.step04_lock,
            "sample_key": args.sample_key,
            "step_02a_manifest": self.expanded_reader.manifest_path,
            "step_03_manifest": self.step03_reader.manifest_path,
            "step_03_table": self.step03_reader.table_path,
            "step_04_manifest": self.step04_reader.manifest_path,
            "step_04_empirical_model": self.step04_reader.model_path,
            "step_05_base_manifest": os.path.join(output_root, MANIFEST_NAME),
            "male_only_validation_manifest": self.male_only_reader.manifest_path,
            "male_only_predictions": self.male_only_reader.predictions_path,
        }
        self.expanded_publisher = ExpandedCohortOutputPublisher(
            self.configuration,
            self.validation,
            self.expanded_plotter,
            expanded_root,
            expanded_paths,
        )

    def run(self):
        if os.path.exists(self.publisher.output_root):
            self.base_verifier.verify()
        else:
            recomputed_empirical = EmpiricalPatternEstimator(self.encoder).fit(
                self.step03_reader.iter_observations()
            )
            empirical_rows = self.step04_reader.read_rows()
            classifier = self.estimator.fit(recomputed_empirical)
            probability_rows, comparison_rows, metrics = self.evaluator.evaluate(
                classifier, empirical_rows, recomputed_empirical
            )
            self.publisher.publish(classifier, probability_rows, comparison_rows, metrics)
        if os.path.exists(self.held_out_publisher.output_root):
            male_predictions = self.male_only_reader.read_predictions()
            male_overall = self.male_only_reader.read_overall()
        else:
            records = self.sample_aware_reader.read_records()
            male_predictions, fold_coefficients = self.held_out_validator.run(records)
            confusion_rows, per_sample_rows, overall_rows = self.held_out_evaluator.evaluate(
                male_predictions
            )
            self.held_out_publisher.publish(
                male_predictions,
                confusion_rows,
                per_sample_rows,
                overall_rows,
                fold_coefficients,
            )
            male_overall = {row["metric"]: float(row["value"]) for row in overall_rows}
        if os.path.exists(self.expanded_publisher.output_root):
            self.expanded_publisher.publish(
                None, None, None, None, None, None, None, None, None, None, None
            )
            return
        male_ground_truth_records = self.sample_aware_reader.read_records()
        expanded_records = self.expanded_reader.read_records(male_ground_truth_records)
        expanded_predictions, expanded_fold_rows = self.expanded_validator.run(
            expanded_records
        )
        expanded_confusion, expanded_per_sample, expanded_overall = self.expanded_evaluator.evaluate(
            expanded_predictions
        )
        expanded_group_rows = self.expanded_group_evaluator.group_metrics(
            expanded_predictions
        )
        male_pattern_rows = self.pattern_error_analyzer.analyze_male_only(
            male_predictions
        )
        expanded_pattern_rows = self.pattern_error_analyzer.analyze_expanded(
            expanded_predictions
        )
        comparison_rows = self.expanded_group_evaluator.comparison(
            male_predictions,
            male_overall,
            expanded_predictions,
            expanded_overall,
        )
        full_classifier, full_probability_rows = self.expanded_full_fit.fit(
            expanded_records
        )
        self.validation.require_equal(
            "male_only_wt_false_ko_errors",
            sum(row["predicted_ko_cells_false_ko"] for row in male_pattern_rows),
            8265,
        )
        self.validation.require_equal(
            "expanded_het_female_predictions",
            sum(row["true_genotype"] == "HET" for row in expanded_predictions),
            0,
        )
        self.validation.require_equal(
            "expanded_full_model_converged", full_classifier.fit_diagnostics["converged"], True
        )
        expanded_group_index = {row["evaluation_group"]: row for row in expanded_group_rows}
        self.validation.require_equal(
            "expanded_wt_male_false_ko_errors",
            expanded_group_index["WT_M"]["false_ko_cells"],
            8265,
        )
        self.validation.require_equal(
            "expanded_wt_female_false_ko_errors",
            expanded_group_index["WT_F"]["false_ko_cells"],
            10035,
        )
        self.validation.require_equal(
            "expanded_non000_calls_follow_fixed_threshold",
            sum(
                row["pattern_code"] != "000"
                and row["predicted_genotype"]
                != (
                    "KO" if row["predicted_ko_probability"] > 0.5
                    else ("WT" if row["predicted_wt_probability"] > 0.5 else "uncalled")
                )
                for row in expanded_predictions
            ),
            0,
        )
        self.expanded_publisher.publish(
            expanded_predictions,
            expanded_confusion,
            expanded_per_sample,
            expanded_group_rows,
            expanded_overall,
            expanded_fold_rows,
            full_classifier,
            full_probability_rows,
            male_pattern_rows,
            expanded_pattern_rows,
            comparison_rows,
        )


def build_parser():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--lock", required=True)
    parser.add_argument("--requirements", required=True)
    parser.add_argument("--step04-lock", required=True)
    parser.add_argument("--step03-root", required=True)
    parser.add_argument("--step04-root", required=True)
    parser.add_argument("--step02-root", required=True)
    parser.add_argument("--sample-key", required=True)
    parser.add_argument("--paper3-root", required=True)
    parser.add_argument("--bundle-root", required=True)
    return parser


def main(argv=None):
    args = build_parser().parse_args(argv)
    PCDH19LogisticRegressionBaselineStep(args).run()


if __name__ == "__main__":
    try:
        main()
    except (Step05Error, Step04Error) as error:
        sys.stderr.write("ERROR: {}\n".format(error))
        sys.exit(1)
