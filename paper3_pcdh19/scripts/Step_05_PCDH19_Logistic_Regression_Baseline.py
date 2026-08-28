#!/usr/bin/env python3
"""Fit the Step 05 three-feature PCDH19 logistic-regression baseline.

The model is an unpenalized binary logistic regression with the exact formula

``logit(P(KO)) = intercept + beta_A*A + beta_B*B + beta_C*C``.

Step 03 encoding is retained: WT=0 and KO=1, so positive coefficients increase
KO log odds and negative coefficients decrease KO log odds. The only
predictors are binary ``A_detected``, ``B_detected``, and ``C_detected``. The
fit retains all ``000`` cells; the intercept is their fitted KO log odds.

No UMI counts, transcriptome features, cell types, interactions, nonlinear
terms, penalty, HET cells, hard calls, confidence threshold, train/test split,
or held-out performance estimate is introduced. Step 05 compares fitted
eight-pattern probabilities with the Step 04 empirical probabilities.
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

    def run(self):
        if os.path.exists(self.publisher.output_root):
            self.publisher.publish(None, [], [], {})
            return
        recomputed_empirical = EmpiricalPatternEstimator(self.encoder).fit(
            self.step03_reader.iter_observations()
        )
        empirical_rows = self.step04_reader.read_rows()
        classifier = self.estimator.fit(recomputed_empirical)
        probability_rows, comparison_rows, metrics = self.evaluator.evaluate(
            classifier, empirical_rows, recomputed_empirical
        )
        self.publisher.publish(classifier, probability_rows, comparison_rows, metrics)


def build_parser():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--lock", required=True)
    parser.add_argument("--requirements", required=True)
    parser.add_argument("--step04-lock", required=True)
    parser.add_argument("--step03-root", required=True)
    parser.add_argument("--step04-root", required=True)
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
