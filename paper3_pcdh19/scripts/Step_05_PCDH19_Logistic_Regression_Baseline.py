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
Step 05 also contains a descriptive A-negative raw-count module that asks
whether unnormalized Probe B and C molecule counts provide progressively
stronger genotype evidence. It fits no additional classifier.

A separately manifested paired model comparison retains the original binary
held-out predictions and fits one additional unpenalized model with predictors
``A_detected``, raw ``B_UMI``, and raw ``C_UMI`` on the identical expanded
WT-M+F/KO-M sample folds. These counts are probe-level UMI/ligation evidence,
not transcript-number estimates. Exact ``000`` cells remain uncalled, the
decision threshold remains fixed at 0.5, and no HET cell is loaded.
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
from collections import Counter, defaultdict

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pcdh19_step05_sample_probe_evidence_diagnostics as sample_diagnostic_module

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
from pcdh19_step05_sample_probe_evidence_diagnostics import (
    DiagnosticError,
    SampleProbeEvidenceDiagnosticStep,
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
A_NEGATIVE_RAW_BC_DIRECTORY = "a_negative_raw_bc_umi_evidence"
COUNT_MODEL_COMPARISON_DIRECTORY = "count_informed_vs_binary_validation"

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

A_NEGATIVE_CELL_TABLE_NAME = "step_05_a_negative_raw_bc_cells.tsv"
A_NEGATIVE_COHORT_SUMMARY_NAME = "step_05_a_negative_raw_bc_cohort_summary.tsv"
A_NEGATIVE_DISTRIBUTION_NAME = "step_05_a_negative_raw_umi_distribution.tsv"
A_NEGATIVE_BIN_ENRICHMENT_NAME = "step_05_a_negative_umi_bin_enrichment.tsv"
A_NEGATIVE_JOINT_NAME = "step_05_a_negative_joint_bc_count_combinations.tsv"
A_NEGATIVE_DISTRIBUTION_PLOT_NAME = "step_05_a_negative_raw_b_c_umi_distributions.png"
A_NEGATIVE_ENRICHMENT_PLOT_NAME = "step_05_a_negative_umi_bin_ko_enrichment.png"
A_NEGATIVE_JOINT_PLOT_NAME = "step_05_a_negative_joint_bc_count_combinations.png"
A_NEGATIVE_CHECKS_NAME = "step_05_a_negative_raw_bc_validation_checks.tsv"
A_NEGATIVE_ENVIRONMENT_NAME = "software_environment.tsv"
A_NEGATIVE_MANIFEST_NAME = "output_manifest.tsv"

COUNT_PAIRED_PREDICTIONS_NAME = "step_05_binary_vs_count_informed_held_out_predictions.tsv"
COUNT_CONFUSION_NAME = "step_05_count_informed_confusion_matrix.tsv"
COUNT_PER_SAMPLE_NAME = "step_05_count_informed_per_sample_metrics.tsv"
COUNT_GROUP_METRICS_NAME = "step_05_count_informed_metrics_by_group.tsv"
COUNT_OVERALL_NAME = "step_05_count_informed_overall_metrics.tsv"
COUNT_FOLD_COEFFICIENTS_NAME = "step_05_count_informed_fold_model_coefficients.tsv"
COUNT_FULL_COEFFICIENTS_NAME = "step_05_count_informed_full_fit_coefficients.tsv"
COUNT_OVERALL_COMPARISON_NAME = "step_05_binary_vs_count_informed_overall_comparison.tsv"
COUNT_SAMPLE_COMPARISON_NAME = "step_05_binary_vs_count_informed_per_sample_comparison.tsv"
COUNT_TRANSITIONS_NAME = "step_05_binary_vs_count_informed_correctness_transitions.tsv"
COUNT_METRICS_PLOT_NAME = "step_05_binary_vs_count_informed_called_metrics.png"
COUNT_SAMPLE_PLOT_NAME = "step_05_binary_vs_count_informed_per_sample_accuracy.png"
COUNT_COEFFICIENT_PLOT_NAME = "step_05_count_informed_logistic_odds_ratios.png"
COUNT_CHECKS_NAME = "step_05_count_informed_validation_checks.tsv"
COUNT_ENVIRONMENT_NAME = "software_environment.tsv"
COUNT_MANIFEST_NAME = "output_manifest.tsv"

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

A_NEGATIVE_CELL_HEADER = [
    "biological_sample_id",
    "submitted_sample_name",
    "sex",
    "design_group",
    "true_genotype",
    "cell_barcode",
    "A_UMI",
    "B_UMI",
    "C_UMI",
    "B_plus_C_UMI",
]

A_NEGATIVE_COHORT_SUMMARY_HEADER = [
    "biological_sample_id",
    "submitted_sample_name",
    "sex",
    "design_group",
    "true_genotype",
    "a_negative_cells",
    "b_positive_cells",
    "c_positive_cells",
    "bc_positive_cells",
    "raw_b_umi",
    "raw_c_umi",
    "raw_b_plus_c_umi",
]

A_NEGATIVE_DISTRIBUTION_HEADER = [
    "measure",
    "raw_umi_count",
    "wt_cells",
    "ko_cells",
    "total_cells",
    "p_wt_given_count",
    "p_ko_given_count",
    "p_count_given_wt",
    "p_count_given_ko",
    "ko_to_wt_likelihood_ratio",
]

A_NEGATIVE_BIN_HEADER = [
    "measure",
    "umi_bin",
    "bin_lower_bound",
    "bin_upper_bound",
    "wt_cells",
    "ko_cells",
    "total_cells",
    "p_wt_given_bin",
    "p_ko_given_bin",
    "ko_to_wt_cell_ratio",
    "p_bin_given_wt",
    "p_bin_given_ko",
    "ko_to_wt_likelihood_ratio",
]

A_NEGATIVE_JOINT_HEADER = [
    "B_UMI",
    "C_UMI",
    "B_plus_C_UMI",
    "wt_cells",
    "ko_cells",
    "total_cells",
    "p_wt_given_joint",
    "p_ko_given_joint",
    "p_joint_given_wt",
    "p_joint_given_ko",
    "ko_to_wt_likelihood_ratio",
]

COUNT_PAIRED_PREDICTION_HEADER = [
    "cell_barcode",
    "biological_sample_id",
    "submitted_sample_name",
    "sex",
    "design_group",
    "true_genotype",
    "A_UMI",
    "A_detected",
    "B_UMI",
    "C_UMI",
    "pattern_code",
    "binary_predicted_wt_probability",
    "binary_predicted_ko_probability",
    "binary_predicted_genotype",
    "count_predicted_wt_probability",
    "count_predicted_ko_probability",
    "count_predicted_genotype",
]

COUNT_OVERALL_COMPARISON_HEADER = [
    "metric",
    "denominator",
    "preferred_direction",
    "binary_model_value",
    "count_informed_model_value",
    "count_minus_binary",
]

COUNT_SAMPLE_COMPARISON_HEADER = [
    "biological_sample_id",
    "submitted_sample_name",
    "sex",
    "design_group",
    "true_genotype",
    "called_cells",
    "binary_accuracy_among_called",
    "count_informed_accuracy_among_called",
    "count_minus_binary_accuracy",
]

COUNT_TRANSITION_HEADER = [
    "binary_prediction_correct",
    "count_informed_prediction_correct",
    "cells",
    "interpretation",
]

COUNT_FOLD_COEFFICIENT_HEADER = [
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
        "existing_expanded_ground_truth_package",
        "a_negative_raw_bc_umi_evidence",
        "existing_a_negative_raw_bc_package",
        "count_informed_model_comparison",
        "existing_count_informed_comparison_package",
        "sample_level_probe_evidence_diagnostics",
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
        raw_evidence = self.values["a_negative_raw_bc_umi_evidence"]
        if raw_evidence["a_negative_definition"] != "A_UMI == 0":
            raise Step05Error("Raw B/C evidence must use exact A_UMI == 0")
        if raw_evidence["normalization"] != "none_raw_molecule_counts":
            raise Step05Error("Raw B/C evidence must not normalize counts")
        if raw_evidence["measures"] != ["B_UMI", "C_UMI", "B_plus_C_UMI"]:
            raise Step05Error("Raw B/C evidence measures are not the locked set")
        if raw_evidence["count_bins"] != ["0", "1", "2", "3+"]:
            raise Step05Error("Raw B/C evidence bins must be 0, 1, 2, and 3+")
        prohibited_raw_flags = [
            "het_cells_permitted",
            "classifier_fitting_permitted",
            "threshold_optimization_permitted",
        ]
        if any(raw_evidence[key] is not False for key in prohibited_raw_flags):
            raise Step05Error("Raw B/C evidence scope enables a prohibited operation")
        count_model = self.values["count_informed_model_comparison"]
        if count_model["cohort_id"] != expanded["cohort_id"]:
            raise Step05Error("Count-informed comparison must use the expanded cohort")
        if count_model["predictors"] != ["A_detected", "B_UMI", "C_UMI"]:
            raise Step05Error("Count-informed predictors must be A detection and raw B/C UMIs")
        if count_model["count_interpretation"] != "raw_probe_level_umi_ligation_evidence_not_transcript_numbers":
            raise Step05Error("Count-informed probe-count interpretation is not locked")
        if count_model["normalization"] != "none_raw_molecule_counts":
            raise Step05Error("Count-informed B/C predictors must not be normalized")
        if count_model["model_formula"] != "logit(P(KO)) = intercept + beta_A*A_detected + beta_B*B_UMI + beta_C*C_UMI":
            raise Step05Error("Unexpected count-informed model formula")
        if count_model["fit_method"] != self.values["fit_method"] or count_model["penalty"] != "none":
            raise Step05Error("Count-informed model must use the existing unpenalized fit")
        if count_model["calling_rule"] != expanded["calling_rule"]:
            raise Step05Error("Binary and count-informed calling rules must match")
        prohibited_count_flags = [
            "interaction_terms",
            "nonlinear_transformations",
            "class_weighting",
            "threshold_optimization",
            "het_cells_permitted",
        ]
        if any(count_model[key] is not False for key in prohibited_count_flags):
            raise Step05Error("Count-informed comparison enables a prohibited operation")
        if count_model["expected_total_cells"] != expanded["expected_total_cells"]:
            raise Step05Error("Count-informed cohort total must match the binary comparator")
        diagnostic = self.values["sample_level_probe_evidence_diagnostics"]
        if diagnostic["cohort_id"] != expanded["cohort_id"] or diagnostic["expected_total_cells"] != expanded["expected_total_cells"]:
            raise Step05Error("Sample diagnostic must use the unchanged expanded cohort")
        if diagnostic["normalization"] != "none_raw_molecule_counts" or diagnostic["count_interpretation"] != count_model["count_interpretation"]:
            raise Step05Error("Sample diagnostic raw-count contract is inconsistent")
        if diagnostic["classifier_fitting_permitted"] is not False or diagnostic["prediction_changes_permitted"] is not False or diagnostic["het_cells_permitted"] is not False:
            raise Step05Error("Sample diagnostic enables a prohibited operation")
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


class ANegativeRawCountRecord(object):
    """Represent one known-genotype A-negative cell with raw B/C molecules."""

    __slots__ = (
        "biological_sample_id",
        "submitted_sample_name",
        "sex",
        "design_group",
        "true_genotype",
        "cell_barcode",
        "a_umi",
        "b_umi",
        "c_umi",
    )

    def __init__(self, sample_id, sample_name, sex, design_group, genotype, barcode, a_umi, b_umi, c_umi):
        self.biological_sample_id = sample_id
        self.submitted_sample_name = sample_name
        self.sex = sex
        self.design_group = design_group
        self.true_genotype = genotype
        self.cell_barcode = barcode
        self.a_umi = int(a_umi)
        self.b_umi = int(b_umi)
        self.c_umi = int(c_umi)
        if self.a_umi != 0 or min(self.b_umi, self.c_umi) < 0:
            raise Step05Error("A-negative raw-count record violates its count contract")

    @property
    def bc_umi(self):
        return self.b_umi + self.c_umi

    def as_row(self):
        return {
            "biological_sample_id": self.biological_sample_id,
            "submitted_sample_name": self.submitted_sample_name,
            "sex": self.sex,
            "design_group": self.design_group,
            "true_genotype": self.true_genotype,
            "cell_barcode": self.cell_barcode,
            "A_UMI": self.a_umi,
            "B_UMI": self.b_umi,
            "C_UMI": self.c_umi,
            "B_plus_C_UMI": self.bc_umi,
        }


class ANegativeRawCountCohortReader(object):
    """Reuse Step 03 males and add only manifested WT-female raw counts."""

    def __init__(
        self,
        configuration,
        shared_configuration,
        validation,
        expanded_input_reader,
        step03_table_path,
    ):
        self.configuration = configuration
        self.shared_configuration = shared_configuration
        self.validation = validation
        self.expanded_input_reader = expanded_input_reader
        self.step03_table_path = os.path.abspath(step03_table_path)
        self.scope = configuration.values["a_negative_raw_bc_umi_evidence"]
        self.cohort = configuration.values["expanded_ground_truth_validation"]

    @staticmethod
    def _validated_counts(row, a_key="A_UMI", b_key="B_UMI", c_key="C_UMI", total_key="Pcdh19_total_UMI"):
        try:
            a_umi = int(row[a_key])
            b_umi = int(row[b_key])
            c_umi = int(row[c_key])
            total = int(row[total_key])
        except (KeyError, TypeError, ValueError):
            raise Step05Error("Invalid raw A/B/C molecule count")
        if min(a_umi, b_umi, c_umi, total) < 0 or a_umi + b_umi + c_umi != total:
            raise Step05Error("Raw A/B/C molecule counts are inconsistent")
        return a_umi, b_umi, c_umi

    def _read_step03_males(self, records, seen):
        expected = self.cohort["expected_samples"]
        with open(self.step03_table_path, "r", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if reader.fieldnames != self.shared_configuration.values["required_input_columns"]:
                raise Step05Error("Unexpected Step 03 schema for raw B/C evidence")
            for row in reader:
                sample_id = row["technical_sample_id"]
                definition = expected.get(sample_id)
                if definition is None or definition["sex"] != "M":
                    raise Step05Error("Step 03 raw B/C input contains an unexpected cohort row")
                if row["ground_truth_class"] != definition["genotype"] or row["sex"] != "M":
                    raise Step05Error("Step 03 raw B/C genotype/sex metadata mismatch")
                a_umi, b_umi, c_umi = self._validated_counts(row)
                if a_umi != 0:
                    continue
                key = (sample_id, row["cell_barcode"])
                if key in seen:
                    raise Step05Error("Duplicate A-negative sample/barcode key")
                seen.add(key)
                records.append(
                    ANegativeRawCountRecord(
                        sample_id,
                        row["submitted_sample_name"],
                        row["sex"],
                        row["design_group"],
                        row["ground_truth_class"],
                        row["cell_barcode"],
                        a_umi,
                        b_umi,
                        c_umi,
                    )
                )

    def _read_wt_females(self, records, seen):
        manifest = self.expanded_input_reader._manifest()
        metadata = self.expanded_input_reader._sample_metadata()
        expected = self.cohort["expected_samples"]
        for sample_id, definition in expected.items():
            if definition["design_group"] != "WT_F":
                continue
            relative_path = os.path.join("per_sample", sample_id, "pcdh19_probe_patterns.tsv")
            manifest_row = manifest.get(relative_path)
            path = os.path.join(self.expanded_input_reader.step02_root, relative_path)
            if manifest_row is None or not os.path.isfile(path):
                raise Step05Error("WT-female raw-count input is not manifested")
            if os.path.getsize(path) != int(manifest_row["bytes"]) or sha256_file(path) != manifest_row["sha256"]:
                raise Step05Error("WT-female raw-count input identity mismatch")
            with open(path, "r", newline="") as handle:
                reader = csv.DictReader(handle, delimiter="\t")
                if reader.fieldnames != ManifestedProbePatternCohortReader.REQUIRED_COLUMNS:
                    raise Step05Error("Unexpected WT-female raw-count schema")
                for row in reader:
                    a_umi, b_umi, c_umi = self._validated_counts(row)
                    if a_umi != 0:
                        continue
                    key = (sample_id, row["barcode"])
                    if key in seen:
                        raise Step05Error("Duplicate A-negative sample/barcode key")
                    seen.add(key)
                    records.append(
                        ANegativeRawCountRecord(
                            sample_id,
                            metadata[sample_id]["submitted_sample_name"],
                            definition["sex"],
                            definition["design_group"],
                            definition["genotype"],
                            row["barcode"],
                            a_umi,
                            b_umi,
                            c_umi,
                        )
                    )

    def read_records(self):
        records = []
        seen = set()
        self._read_step03_males(records, seen)
        self._read_wt_females(records, seen)
        class_counts = {genotype: sum(row.true_genotype == genotype for row in records) for genotype in ("WT", "KO")}
        group_counts = {group: sum(row.design_group == group for row in records) for group in ("WT_M", "WT_F", "KO_M")}
        self.validation.require_equal("a_negative_raw_bc_cells", len(records), self.scope["expected_a_negative_cells"])
        self.validation.require_equal("a_negative_unique_sample_barcode_keys", len(seen), len(records))
        for genotype, expected in self.scope["expected_class_cells"].items():
            self.validation.require_equal("a_negative_{}_cells".format(genotype.lower()), class_counts[genotype], expected)
        for group, expected in self.scope["expected_group_cells"].items():
            self.validation.require_equal("a_negative_{}_cells".format(group.lower()), group_counts[group], expected)
        self.validation.require_equal("a_negative_het_female_cells_loaded", sum(row.true_genotype == "HET" for row in records), 0)
        self.validation.require_equal("a_negative_nonzero_a_umi_rows", sum(row.a_umi != 0 for row in records), 0)
        return records


class RawBCUMIEvidenceSummarizer(object):
    """Summarize raw B/C molecule evidence without fitting a classifier."""

    MEASURE_KEYS = ["B_UMI", "C_UMI", "B_plus_C_UMI"]
    BIN_ORDER = ["0", "1", "2", "3+"]

    def __init__(self, configuration, validation):
        self.configuration = configuration
        self.validation = validation
        self.scope = configuration.values["a_negative_raw_bc_umi_evidence"]

    @staticmethod
    def _value(record, measure):
        if measure == "B_UMI":
            return record.b_umi
        if measure == "C_UMI":
            return record.c_umi
        if measure == "B_plus_C_UMI":
            return record.bc_umi
        raise Step05Error("Unknown raw B/C evidence measure")

    @staticmethod
    def _bin(value):
        if value >= 3:
            return "3+"
        return str(value)

    @staticmethod
    def _ratio(numerator, denominator):
        if denominator == 0:
            return math.inf if numerator > 0 else None
        return numerator / denominator

    @classmethod
    def _evidence_fields(cls, wt_cells, ko_cells, wt_total, ko_total, suffix):
        total = wt_cells + ko_cells
        p_wt = wt_cells / total if total else None
        p_ko = ko_cells / total if total else None
        p_given_wt = wt_cells / wt_total
        p_given_ko = ko_cells / ko_total
        return {
            "wt_cells": wt_cells,
            "ko_cells": ko_cells,
            "total_cells": total,
            "p_wt_given_{}".format(suffix): p_wt,
            "p_ko_given_{}".format(suffix): p_ko,
            "p_{}_given_wt".format(suffix): p_given_wt,
            "p_{}_given_ko".format(suffix): p_given_ko,
            "ko_to_wt_likelihood_ratio": cls._ratio(p_given_ko, p_given_wt),
        }

    def cohort_summary(self, records):
        definitions = self.configuration.values["expanded_ground_truth_validation"]["expected_samples"]
        rows = []
        for sample_id, definition in definitions.items():
            subset = [record for record in records if record.biological_sample_id == sample_id]
            rows.append(
                {
                    "biological_sample_id": sample_id,
                    "submitted_sample_name": subset[0].submitted_sample_name,
                    "sex": definition["sex"],
                    "design_group": definition["design_group"],
                    "true_genotype": definition["genotype"],
                    "a_negative_cells": len(subset),
                    "b_positive_cells": sum(record.b_umi > 0 for record in subset),
                    "c_positive_cells": sum(record.c_umi > 0 for record in subset),
                    "bc_positive_cells": sum(record.bc_umi > 0 for record in subset),
                    "raw_b_umi": sum(record.b_umi for record in subset),
                    "raw_c_umi": sum(record.c_umi for record in subset),
                    "raw_b_plus_c_umi": sum(record.bc_umi for record in subset),
                }
            )
        return rows

    def exact_distributions(self, records):
        totals = {genotype: sum(record.true_genotype == genotype for record in records) for genotype in ("WT", "KO")}
        rows = []
        for measure in self.MEASURE_KEYS:
            maximum = max(self._value(record, measure) for record in records)
            self.validation.require_equal(
                "a_negative_{}_maximum_raw_umi".format(measure.lower()),
                maximum,
                self.scope["expected_max_raw_counts"][measure],
            )
            counts = Counter((record.true_genotype, self._value(record, measure)) for record in records)
            for raw_count in range(maximum + 1):
                row = {"measure": measure, "raw_umi_count": raw_count}
                row.update(self._evidence_fields(counts[("WT", raw_count)], counts[("KO", raw_count)], totals["WT"], totals["KO"], "count"))
                rows.append(row)
        return rows

    def binned_enrichment(self, records):
        totals = {genotype: sum(record.true_genotype == genotype for record in records) for genotype in ("WT", "KO")}
        rows = []
        bounds = {"0": (0, 0), "1": (1, 1), "2": (2, 2), "3+": (3, None)}
        for measure in self.MEASURE_KEYS:
            counts = Counter((record.true_genotype, self._bin(self._value(record, measure))) for record in records)
            for bin_name in self.BIN_ORDER:
                lower, upper = bounds[bin_name]
                wt_cells = counts[("WT", bin_name)]
                ko_cells = counts[("KO", bin_name)]
                row = {
                    "measure": measure,
                    "umi_bin": bin_name,
                    "bin_lower_bound": lower,
                    "bin_upper_bound": "Inf" if upper is None else upper,
                    "ko_to_wt_cell_ratio": self._ratio(ko_cells, wt_cells),
                }
                row.update(self._evidence_fields(wt_cells, ko_cells, totals["WT"], totals["KO"], "bin"))
                rows.append(row)
        self.validation.require_equal("a_negative_raw_count_enrichment_rows", len(rows), 12)
        return rows

    def joint_combinations(self, records):
        totals = {genotype: sum(record.true_genotype == genotype for record in records) for genotype in ("WT", "KO")}
        counts = Counter((record.b_umi, record.c_umi, record.true_genotype) for record in records)
        combinations = sorted({(record.b_umi, record.c_umi) for record in records})
        rows = []
        for b_umi, c_umi in combinations:
            row = {"B_UMI": b_umi, "C_UMI": c_umi, "B_plus_C_UMI": b_umi + c_umi}
            row.update(self._evidence_fields(counts[(b_umi, c_umi, "WT")], counts[(b_umi, c_umi, "KO")], totals["WT"], totals["KO"], "joint"))
            rows.append(row)
        self.validation.require_equal("a_negative_joint_observed_combinations", len(rows), 33)
        return rows

    def summarize(self, records):
        return (
            self.cohort_summary(records),
            self.exact_distributions(records),
            self.binned_enrichment(records),
            self.joint_combinations(records),
        )


class CountInformedCellRecord(object):
    """One control cell with raw probe-level UMI/ligation evidence."""

    __slots__ = (
        "cell_barcode",
        "biological_sample_id",
        "submitted_sample_name",
        "sex",
        "design_group",
        "true_genotype",
        "a_umi",
        "a_detected",
        "b_umi",
        "c_umi",
        "pattern_code",
    )

    def __init__(self, barcode, sample_id, sample_name, sex, design_group, genotype, a_umi, b_umi, c_umi, pattern_code):
        self.cell_barcode = barcode
        self.biological_sample_id = sample_id
        self.submitted_sample_name = sample_name
        self.sex = sex
        self.design_group = design_group
        self.true_genotype = genotype
        self.a_umi = int(a_umi)
        self.a_detected = int(self.a_umi > 0)
        self.b_umi = int(b_umi)
        self.c_umi = int(c_umi)
        self.pattern_code = pattern_code
        if min(self.a_umi, self.b_umi, self.c_umi) < 0:
            raise Step05Error("Count-informed record contains a negative raw probe count")

    def features(self):
        return {
            "A_detected": self.a_detected,
            "B_UMI": self.b_umi,
            "C_UMI": self.c_umi,
        }


class CountInformedGroundTruthReader(object):
    """Load the exact binary-comparator cells with raw A/B/C probe counts."""

    def __init__(self, configuration, encoder, validation, expanded_reader, binary_predictions_path):
        self.configuration = configuration
        self.encoder = encoder
        self.validation = validation
        self.expanded_reader = expanded_reader
        self.binary_predictions_path = os.path.abspath(binary_predictions_path)
        self.scope = configuration.values["count_informed_model_comparison"]
        self.cohort = configuration.values["expanded_ground_truth_validation"]

    def _binary_rows(self):
        rows = {}
        with open(self.binary_predictions_path, "r", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for row in reader:
                key = (row["biological_sample_id"], row["cell_barcode"])
                if key in rows:
                    raise Step05Error("Duplicate cell in binary comparator predictions")
                rows[key] = row
        self.validation.require_equal(
            "count_comparison_binary_prediction_rows",
            len(rows),
            self.scope["expected_total_cells"],
        )
        return rows

    def read_records(self):
        manifest = self.expanded_reader._manifest()
        metadata = self.expanded_reader._sample_metadata()
        binary_rows = self._binary_rows()
        records = []
        seen = set()
        for sample_id, definition in self.cohort["expected_samples"].items():
            relative_path = os.path.join("per_sample", sample_id, "pcdh19_probe_patterns.tsv")
            manifest_row = manifest.get(relative_path)
            path = os.path.join(self.expanded_reader.step02_root, relative_path)
            if manifest_row is None or not os.path.isfile(path):
                raise Step05Error("Count-informed input is not present in the Step 02a manifest")
            if os.path.getsize(path) != int(manifest_row["bytes"]) or sha256_file(path) != manifest_row["sha256"]:
                raise Step05Error("Count-informed Step 02a input identity mismatch")
            with open(path, "r", newline="") as handle:
                reader = csv.DictReader(handle, delimiter="\t")
                if reader.fieldnames != ManifestedProbePatternCohortReader.REQUIRED_COLUMNS:
                    raise Step05Error("Unexpected count-informed probe table schema")
                for row in reader:
                    key = (sample_id, row["barcode"])
                    if key in seen:
                        raise Step05Error("Duplicate count-informed sample/barcode key")
                    seen.add(key)
                    try:
                        a_umi, b_umi, c_umi = [int(row[name]) for name in ("A_UMI", "B_UMI", "C_UMI")]
                        total_umi = int(row["Pcdh19_total_UMI"])
                    except (TypeError, ValueError):
                        raise Step05Error("Invalid raw count in count-informed input")
                    if min(a_umi, b_umi, c_umi, total_umi) < 0 or a_umi + b_umi + c_umi != total_umi:
                        raise Step05Error("Count-informed raw probe counts are inconsistent")
                    pattern_code = self.encoder.encode(int(a_umi > 0), int(b_umi > 0), int(c_umi > 0))
                    binary = binary_rows.get(key)
                    if binary is None:
                        raise Step05Error("Count-informed cell is absent from the binary comparator")
                    if (
                        binary["pattern_code"] != pattern_code
                        or binary["true_genotype"] != definition["genotype"]
                        or binary["sex"] != definition["sex"]
                        or binary["design_group"] != definition["design_group"]
                    ):
                        raise Step05Error("Count-informed cell differs from its binary comparator row")
                    records.append(
                        CountInformedCellRecord(
                            row["barcode"],
                            sample_id,
                            metadata[sample_id]["submitted_sample_name"],
                            definition["sex"],
                            definition["design_group"],
                            definition["genotype"],
                            a_umi,
                            b_umi,
                            c_umi,
                            pattern_code,
                        )
                    )
        self.validation.require_equal("count_informed_ground_truth_cells", len(records), self.scope["expected_total_cells"])
        self.validation.require_equal("count_informed_unique_sample_barcode_keys", len(seen), len(records))
        binary_keys = set(binary_rows)
        self.validation.require_equal("count_informed_keys_missing_from_binary", len(seen.difference(binary_keys)), 0)
        self.validation.require_equal("binary_keys_missing_from_count_informed", len(binary_keys.difference(seen)), 0)
        for genotype, expected in self.scope["expected_class_cells"].items():
            self.validation.require_equal(
                "count_informed_{}_cells".format(genotype.lower()),
                sum(record.true_genotype == genotype for record in records),
                expected,
            )
        self.validation.require_equal("count_informed_het_cells_loaded", sum(record.true_genotype == "HET" for record in records), 0)
        self.validation.require_equal("count_informed_max_raw_b_umi", max(record.b_umi for record in records), self.scope["expected_max_raw_counts"]["B_UMI"])
        self.validation.require_equal("count_informed_max_raw_c_umi", max(record.c_umi for record in records), self.scope["expected_max_raw_counts"]["C_UMI"])
        return records, binary_rows


class CountInformedLogisticClassifier(ProbabilisticClassifier):
    """Unpenalized A-detection plus raw-B/C-count probability model."""

    TERM_NAMES = ["intercept", "A_detected", "B_UMI", "C_UMI"]

    def __init__(self, coefficients, diagnostics):
        self.coefficients = np.asarray(coefficients, dtype=float)
        self.fit_diagnostics = dict(diagnostics)
        if self.coefficients.shape != (4,):
            raise Step05Error("Count-informed logistic model requires four coefficients")

    @property
    def model_name(self):
        return "pcdh19_a_detection_raw_b_c_umi_logistic_regression"

    def _feature_vector(self, features):
        a_detected = int(features["A_detected"])
        b_umi = int(features["B_UMI"])
        c_umi = int(features["C_UMI"])
        if a_detected not in (0, 1) or min(b_umi, c_umi) < 0:
            raise Step05Error("Invalid count-informed logistic feature vector")
        return np.asarray([1.0, a_detected, b_umi, c_umi], dtype=float)

    def linear_predictor(self, features):
        return float(np.dot(self._feature_vector(features), self.coefficients))

    def predict_proba(self, features):
        p_ko = LogisticRegressionClassifier._expit(self.linear_predictor(features))
        return {"WT": 1.0 - p_ko, "KO": p_ko}

    def coefficient_rows(self):
        rows = []
        for index, (term, coefficient) in enumerate(zip(self.TERM_NAMES, self.coefficients)):
            if index == 0:
                interpretation = "KO:WT odds when A is absent and raw B=C=0"
            elif term == "A_detected":
                interpretation = "KO:WT odds change for A detected versus absent, holding raw B/C UMIs fixed"
            else:
                interpretation = "KO:WT odds change per one additional raw probe-level UMI/ligation event, holding other predictors fixed"
            rows.append(
                {
                    "term": term,
                    "predictor": term,
                    "outcome_log_odds": "KO_vs_WT; WT=0, KO=1",
                    "coefficient": float(coefficient),
                    "odds_ratio": float(math.exp(coefficient)),
                    "interpretation": interpretation,
                }
            )
        return rows


class CountInformedLogisticEstimator(object):
    """Fit the raw-B/C main-effects model using exact grouped-binomial IRLS."""

    def __init__(self, configuration):
        self.configuration = configuration

    def fit(self, records):
        counts = Counter(
            (record.a_detected, record.b_umi, record.c_umi, record.true_genotype)
            for record in records
        )
        features = sorted({(record.a_detected, record.b_umi, record.c_umi) for record in records})
        design = np.asarray([[1.0, a_value, b_umi, c_umi] for a_value, b_umi, c_umi in features], dtype=float)
        wt_cells = np.asarray([counts[(a, b, c, "WT")] for a, b, c in features], dtype=float)
        ko_cells = np.asarray([counts[(a, b, c, "KO")] for a, b, c in features], dtype=float)
        trials = wt_cells + ko_cells
        coefficients = np.zeros(4, dtype=float)
        coefficients[0] = math.log(float(np.sum(ko_cells)) / float(np.sum(wt_cells)))
        optimizer = self.configuration.values["optimizer"]
        maximum_iterations = int(optimizer["maximum_iterations"])
        tolerance = float(optimizer["coefficient_tolerance"])
        maximum_halvings = int(optimizer["maximum_step_halvings"])
        step_halvings = 0
        converged = False
        for iteration in range(1, maximum_iterations + 1):
            linear = np.dot(design, coefficients)
            probabilities = 1.0 / (1.0 + np.exp(-linear))
            weights = trials * probabilities * (1.0 - probabilities)
            score = np.dot(design.T, ko_cells - trials * probabilities)
            information = np.dot(design.T, weights[:, None] * design)
            try:
                delta = np.linalg.solve(information, score)
            except np.linalg.LinAlgError as error:
                raise Step05Error("Count-informed logistic information matrix is singular: {}".format(error))
            old_likelihood = LogisticRegressionEstimator._log_likelihood(design, trials, ko_cells, coefficients)
            scale = 1.0
            candidate = coefficients + delta
            candidate_likelihood = LogisticRegressionEstimator._log_likelihood(design, trials, ko_cells, candidate)
            halvings = 0
            while candidate_likelihood < old_likelihood and halvings < maximum_halvings:
                scale *= 0.5
                candidate = coefficients + scale * delta
                candidate_likelihood = LogisticRegressionEstimator._log_likelihood(design, trials, ko_cells, candidate)
                halvings += 1
            step_halvings += halvings
            if candidate_likelihood < old_likelihood:
                raise Step05Error("Count-informed IRLS failed to increase likelihood")
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
            "log_likelihood": LogisticRegressionEstimator._log_likelihood(design, trials, ko_cells, coefficients),
            "maximum_absolute_score": float(np.max(np.abs(final_score))),
            "information_matrix_rank": int(np.linalg.matrix_rank(information)),
            "information_matrix_condition_number": float(np.linalg.cond(information)),
            "fit_cells": int(np.sum(trials)),
            "fit_feature_combinations": len(features),
        }
        if not converged or diagnostics["information_matrix_rank"] != 4:
            raise Step05Error("Count-informed logistic regression did not converge at full rank")
        return CountInformedLogisticClassifier(coefficients, diagnostics)


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


class LeaveOneSampleOutCountLogisticValidator(object):
    """Fit the count-informed model on the same registered-sample folds."""

    def __init__(self, configuration, validation, calling_policy):
        self.configuration = configuration
        self.validation = validation
        self.calling_policy = calling_policy
        self.cohort = configuration.values["expanded_ground_truth_validation"]
        self.estimator = CountInformedLogisticEstimator(configuration)

    def run(self, records):
        predictions = []
        coefficient_rows = []
        for held_out_sample, definition in self.cohort["expected_samples"].items():
            training_records = [record for record in records if record.biological_sample_id != held_out_sample]
            test_records = [record for record in records if record.biological_sample_id == held_out_sample]
            training_samples = sorted({record.biological_sample_id for record in training_records})
            self.validation.require_equal("count_{}_training_sample_count".format(held_out_sample), len(training_samples), 8)
            self.validation.require_equal("count_{}_training_excludes_holdout".format(held_out_sample), held_out_sample in training_samples, False)
            self.validation.require_equal("count_{}_held_out_cells".format(held_out_sample), len(test_records), definition["cells"])
            classifier = self.estimator.fit(training_records)
            self.validation.require_equal("count_{}_model_converged".format(held_out_sample), classifier.fit_diagnostics["converged"], True)
            for row in classifier.coefficient_rows():
                coefficient_rows.append(
                    {
                        "held_out_biological_sample_id": held_out_sample,
                        "held_out_true_genotype": definition["genotype"],
                        "held_out_sex": definition["sex"],
                        "held_out_design_group": definition["design_group"],
                        "training_samples": ",".join(training_samples),
                        "training_cells": len(training_records),
                        "term": row["term"],
                        "coefficient": row["coefficient"],
                        "odds_ratio": row["odds_ratio"],
                        "converged": classifier.fit_diagnostics["converged"],
                        "iterations": classifier.fit_diagnostics["iterations"],
                    }
                )
            for record in test_records:
                probabilities = classifier.predict_proba(record.features())
                predictions.append(
                    {
                        "cell_barcode": record.cell_barcode,
                        "biological_sample_id": record.biological_sample_id,
                        "submitted_sample_name": record.submitted_sample_name,
                        "sex": record.sex,
                        "design_group": record.design_group,
                        "true_genotype": record.true_genotype,
                        "A_UMI": record.a_umi,
                        "A_detected": record.a_detected,
                        "B_UMI": record.b_umi,
                        "C_UMI": record.c_umi,
                        "pattern_code": record.pattern_code,
                        "predicted_wt_probability": probabilities["WT"],
                        "predicted_ko_probability": probabilities["KO"],
                        "predicted_genotype": self.calling_policy.call(record.pattern_code, probabilities),
                    }
                )
        self.validation.require_equal("count_held_out_prediction_rows", len(predictions), len(records))
        self.validation.require_equal(
            "count_each_cell_predicted_once",
            len({(row["biological_sample_id"], row["cell_barcode"]) for row in predictions}),
            len(records),
        )
        self.validation.require_equal("count_fold_coefficient_rows", len(coefficient_rows), 36)
        return predictions, coefficient_rows


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


class ExistingExpandedGroundTruthPackageVerifier(object):
    """Protect the published WT-M+F/KO-M Step 05 comparison package."""

    def __init__(self, configuration, validation, base_root):
        self.configuration = configuration
        self.validation = validation
        definition = configuration.values["existing_expanded_ground_truth_package"]
        self.root = os.path.join(os.path.abspath(base_root), definition["relative_directory"])
        self.manifest_path = os.path.join(self.root, EXPANDED_MANIFEST_NAME)
        self.predictions_path = os.path.join(self.root, EXPANDED_PREDICTIONS_NAME)

    def verify(self):
        definition = self.configuration.values["existing_expanded_ground_truth_package"]
        if not os.path.isfile(self.manifest_path):
            raise Step05Error("Expanded ground-truth benchmark package is missing")
        self.validation.require_equal(
            "expanded_ground_truth_manifest_unchanged",
            sha256_file(self.manifest_path),
            definition["output_manifest_sha256"],
        )
        with open(self.manifest_path, "r", newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        self.validation.require_equal("expanded_ground_truth_manifest_files", len(rows), 18)
        for row in rows:
            path = os.path.join(self.root, row["relative_path"])
            if not os.path.isfile(path) or os.path.getsize(path) != int(row["bytes"]):
                raise Step05Error("Expanded ground-truth benchmark is incomplete")
            if sha256_file(path) != row["sha256"]:
                raise Step05Error("Expanded ground-truth benchmark identity mismatch")
        self.validation.require_equal(
            "expanded_ground_truth_predictions_unchanged",
            sha256_file(self.predictions_path),
            definition["held_out_predictions_sha256"],
        )

    def read_predictions(self):
        self.verify()
        rows = []
        with open(self.predictions_path, "r", newline="") as handle:
            rows.extend(csv.DictReader(handle, delimiter="\t"))
        self.validation.require_equal("expanded_binary_comparator_rows", len(rows), 349686)
        return rows


class ExistingRawBCEvidencePackageVerifier(object):
    """Protect the published descriptive raw-count package from source drift."""

    def __init__(self, configuration, validation, base_root):
        self.configuration = configuration
        self.validation = validation
        definition = configuration.values["existing_a_negative_raw_bc_package"]
        self.root = os.path.join(os.path.abspath(base_root), definition["relative_directory"])
        self.manifest_path = os.path.join(self.root, A_NEGATIVE_MANIFEST_NAME)
        self.cell_table_path = os.path.join(self.root, A_NEGATIVE_CELL_TABLE_NAME)

    def verify(self):
        definition = self.configuration.values["existing_a_negative_raw_bc_package"]
        if not os.path.isfile(self.manifest_path):
            raise Step05Error("Published A-negative raw B/C package is missing")
        self.validation.require_equal(
            "a_negative_raw_bc_manifest_unchanged",
            sha256_file(self.manifest_path),
            definition["output_manifest_sha256"],
        )
        with open(self.manifest_path, "r", newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        self.validation.require_equal("a_negative_raw_bc_manifest_files", len(rows), 10)
        for row in rows:
            path = os.path.join(self.root, row["relative_path"])
            if not os.path.isfile(path) or os.path.getsize(path) != int(row["bytes"]):
                raise Step05Error("Published A-negative raw B/C package is incomplete")
            if sha256_file(path) != row["sha256"]:
                raise Step05Error("Published A-negative raw B/C artifact identity mismatch")
        self.validation.require_equal(
            "a_negative_raw_bc_cell_table_unchanged",
            sha256_file(self.cell_table_path),
            definition["cell_table_sha256"],
        )


class ExistingCountModelComparisonPackageVerifier(object):
    """Protect the published paired binary/count comparison package."""

    def __init__(self, configuration, validation, base_root):
        definition = configuration.values["existing_count_informed_comparison_package"]
        self.definition = definition
        self.validation = validation
        self.root = os.path.join(os.path.abspath(base_root), definition["relative_directory"])
        self.manifest_path = os.path.join(self.root, COUNT_MANIFEST_NAME)
        self.predictions_path = os.path.join(self.root, COUNT_PAIRED_PREDICTIONS_NAME)
        self.comparison_path = os.path.join(self.root, COUNT_OVERALL_COMPARISON_NAME)

    def verify(self):
        if not os.path.isfile(self.manifest_path):
            raise Step05Error("Published count-model comparison package is missing")
        self.validation.require_equal("count_comparison_manifest_unchanged", sha256_file(self.manifest_path), self.definition["output_manifest_sha256"])
        with open(self.manifest_path, "r", newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        self.validation.require_equal("count_comparison_manifest_files", len(rows), 15)
        for row in rows:
            path = os.path.join(self.root, row["relative_path"])
            if not os.path.isfile(path) or os.path.getsize(path) != int(row["bytes"]) or sha256_file(path) != row["sha256"]:
                raise Step05Error("Published count-model comparison artifact mismatch")
        self.validation.require_equal("count_comparison_predictions_unchanged", sha256_file(self.predictions_path), self.definition["paired_predictions_sha256"])


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


class PairedBinaryCountModelEvaluator(object):
    """Compare binary and count-informed held-out predictions cell for cell."""

    def __init__(self, configuration, validation):
        self.configuration = configuration
        self.validation = validation
        self.scope = configuration.values["count_informed_model_comparison"]

    @staticmethod
    def _probability_scores(rows):
        probabilities = [float(row["predicted_ko_probability"]) for row in rows]
        outcomes = [1 if row["true_genotype"] == "KO" else 0 for row in rows]
        brier = sum((probability - outcome) ** 2 for probability, outcome in zip(probabilities, outcomes)) / len(rows)
        epsilon = 1e-15
        log_loss = -sum(
            outcome * math.log(max(probability, epsilon))
            + (1 - outcome) * math.log(max(1.0 - probability, epsilon))
            for probability, outcome in zip(probabilities, outcomes)
        ) / len(rows)
        order = sorted(range(len(rows)), key=lambda index: probabilities[index])
        positive_ranks = 0.0
        rank = 1
        index = 0
        while index < len(order):
            end = index + 1
            while end < len(order) and probabilities[order[end]] == probabilities[order[index]]:
                end += 1
            average_rank = (rank + rank + (end - index) - 1) / 2.0
            positive_ranks += average_rank * sum(outcomes[order[position]] for position in range(index, end))
            rank += end - index
            index = end
        positives = sum(outcomes)
        negatives = len(outcomes) - positives
        roc_auc = (positive_ranks - positives * (positives + 1) / 2.0) / (positives * negatives)
        return {"brier_score_all_cells": brier, "log_loss_all_cells": log_loss, "roc_auc_all_cells": roc_auc}

    @staticmethod
    def _called_metrics(rows):
        called = [row for row in rows if row["predicted_genotype"] in ("WT", "KO")]
        counts = Counter((row["true_genotype"], row["predicted_genotype"]) for row in called)
        tn = counts[("WT", "WT")]
        fp = counts[("WT", "KO")]
        fn = counts[("KO", "WT")]
        tp = counts[("KO", "KO")]
        return {
            "percent_called": 100.0 * len(called) / len(rows),
            "accuracy_among_called": (tn + tp) / len(called),
            "ko_sensitivity_among_called": tp / (tp + fn),
            "ko_specificity_among_called": tn / (tn + fp),
        }

    def compare(self, binary_predictions, count_predictions):
        binary_by_key = {(row["biological_sample_id"], row["cell_barcode"]): row for row in binary_predictions}
        count_by_key = {(row["biological_sample_id"], row["cell_barcode"]): row for row in count_predictions}
        self.validation.require_equal("paired_binary_prediction_keys", len(binary_by_key), self.scope["expected_total_cells"])
        self.validation.require_equal("paired_count_prediction_keys", len(count_by_key), self.scope["expected_total_cells"])
        self.validation.require_equal("paired_prediction_key_differences", len(set(binary_by_key).symmetric_difference(count_by_key)), 0)
        paired_rows = []
        for count in count_predictions:
            key = (count["biological_sample_id"], count["cell_barcode"])
            binary = binary_by_key[key]
            for field in ("true_genotype", "sex", "design_group", "pattern_code"):
                if str(binary[field]) != str(count[field]):
                    raise Step05Error("Paired binary/count metadata mismatch for {}".format(field))
            paired_rows.append(
                {
                    "cell_barcode": count["cell_barcode"],
                    "biological_sample_id": count["biological_sample_id"],
                    "submitted_sample_name": count["submitted_sample_name"],
                    "sex": count["sex"],
                    "design_group": count["design_group"],
                    "true_genotype": count["true_genotype"],
                    "A_UMI": count["A_UMI"],
                    "A_detected": count["A_detected"],
                    "B_UMI": count["B_UMI"],
                    "C_UMI": count["C_UMI"],
                    "pattern_code": count["pattern_code"],
                    "binary_predicted_wt_probability": float(binary["predicted_wt_probability"]),
                    "binary_predicted_ko_probability": float(binary["predicted_ko_probability"]),
                    "binary_predicted_genotype": binary["predicted_genotype"],
                    "count_predicted_wt_probability": count["predicted_wt_probability"],
                    "count_predicted_ko_probability": count["predicted_ko_probability"],
                    "count_predicted_genotype": count["predicted_genotype"],
                }
            )
        binary_scores = self._probability_scores(binary_predictions)
        count_scores = self._probability_scores(count_predictions)
        binary_metrics = self._called_metrics(binary_predictions)
        count_metrics = self._called_metrics(count_predictions)
        definitions = [
            ("percent_called", "all held-out cells", "higher", binary_metrics["percent_called"], count_metrics["percent_called"]),
            ("accuracy_among_called", "called held-out cells", "higher", binary_metrics["accuracy_among_called"], count_metrics["accuracy_among_called"]),
            ("ko_sensitivity_among_called", "called true-KO cells", "higher", binary_metrics["ko_sensitivity_among_called"], count_metrics["ko_sensitivity_among_called"]),
            ("ko_specificity_among_called", "called true-WT cells", "higher", binary_metrics["ko_specificity_among_called"], count_metrics["ko_specificity_among_called"]),
            ("brier_score_all_cells", "all held-out cells", "lower", binary_scores["brier_score_all_cells"], count_scores["brier_score_all_cells"]),
            ("log_loss_all_cells", "all held-out cells", "lower", binary_scores["log_loss_all_cells"], count_scores["log_loss_all_cells"]),
            ("roc_auc_all_cells", "pooled predictions from all held-out folds", "higher", binary_scores["roc_auc_all_cells"], count_scores["roc_auc_all_cells"]),
        ]
        overall_rows = [
            {
                "metric": metric,
                "denominator": denominator,
                "preferred_direction": direction,
                "binary_model_value": binary_value,
                "count_informed_model_value": count_value,
                "count_minus_binary": count_value - binary_value,
            }
            for metric, denominator, direction, binary_value, count_value in definitions
        ]
        sample_rows = []
        for sample_id, definition in self.configuration.values["expanded_ground_truth_validation"]["expected_samples"].items():
            binary_sample = [row for row in binary_predictions if row["biological_sample_id"] == sample_id and row["predicted_genotype"] in ("WT", "KO")]
            count_sample = [row for row in count_predictions if row["biological_sample_id"] == sample_id and row["predicted_genotype"] in ("WT", "KO")]
            if len(binary_sample) != len(count_sample):
                raise Step05Error("Binary and count-informed call denominators differ")
            binary_accuracy = sum(row["predicted_genotype"] == row["true_genotype"] for row in binary_sample) / len(binary_sample)
            count_accuracy = sum(row["predicted_genotype"] == row["true_genotype"] for row in count_sample) / len(count_sample)
            sample_rows.append(
                {
                    "biological_sample_id": sample_id,
                    "submitted_sample_name": count_sample[0]["submitted_sample_name"],
                    "sex": definition["sex"],
                    "design_group": definition["design_group"],
                    "true_genotype": definition["genotype"],
                    "called_cells": len(count_sample),
                    "binary_accuracy_among_called": binary_accuracy,
                    "count_informed_accuracy_among_called": count_accuracy,
                    "count_minus_binary_accuracy": count_accuracy - binary_accuracy,
                }
            )
        binary_macro_accuracy = sum(row["binary_accuracy_among_called"] for row in sample_rows) / len(sample_rows)
        count_macro_accuracy = sum(row["count_informed_accuracy_among_called"] for row in sample_rows) / len(sample_rows)
        overall_rows.append(
            {
                "metric": "mean_per_sample_accuracy_among_called",
                "denominator": "nine held-out samples weighted equally",
                "preferred_direction": "higher",
                "binary_model_value": binary_macro_accuracy,
                "count_informed_model_value": count_macro_accuracy,
                "count_minus_binary": count_macro_accuracy - binary_macro_accuracy,
            }
        )
        transition_counts = Counter()
        for pair in paired_rows:
            if pair["binary_predicted_genotype"] not in ("WT", "KO") or pair["count_predicted_genotype"] not in ("WT", "KO"):
                continue
            binary_correct = pair["binary_predicted_genotype"] == pair["true_genotype"]
            count_correct = pair["count_predicted_genotype"] == pair["true_genotype"]
            transition_counts[(binary_correct, count_correct)] += 1
        labels = {
            (False, False): "wrong under both models",
            (False, True): "corrected by count-informed model",
            (True, False): "regressed under count-informed model",
            (True, True): "correct under both models",
        }
        transition_rows = [
            {
                "binary_prediction_correct": str(binary_correct).lower(),
                "count_informed_prediction_correct": str(count_correct).lower(),
                "cells": transition_counts[(binary_correct, count_correct)],
                "interpretation": labels[(binary_correct, count_correct)],
            }
            for binary_correct, count_correct in ((False, False), (False, True), (True, False), (True, True))
        ]
        self.validation.require_equal("paired_comparison_rows", len(paired_rows), self.scope["expected_total_cells"])
        self.validation.require_equal("paired_called_denominator", sum(row["cells"] for row in transition_rows), 65314)
        self.validation.require_equal("count_model_classifier_family", "logistic_regression", "logistic_regression")
        self.validation.require_equal("count_model_threshold_optimized", self.scope["threshold_optimization"], False)
        return paired_rows, overall_rows, sample_rows, transition_rows


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


class RawBCUMIEvidencePlotter(object):
    """Render descriptive raw-count plots without model fitting."""

    WT_COLOR = "#3B6FB6"
    KO_COLOR = "#D95F43"
    MEASURE_COLORS = {"B_UMI": "#4C78A8", "C_UMI": "#F58518", "B_plus_C_UMI": "#8E5BA6"}

    def __init__(self, dpi):
        self.dpi = int(dpi)

    def plot_distributions(self, rows, path):
        measures = ["B_UMI", "C_UMI", "B_plus_C_UMI"]
        titles = ["Raw Probe B UMI", "Raw Probe C UMI", "Raw B+C UMI"]
        figure, axes = plt.subplots(1, 3, figsize=(15, 4.8))
        for axis, measure, title in zip(axes, measures, titles):
            subset = [row for row in rows if row["measure"] == measure]
            x_values = [row["raw_umi_count"] for row in subset]
            axis.plot(x_values, [100.0 * row["p_count_given_wt"] for row in subset], marker="o", color=self.WT_COLOR, label="WT-M+F")
            axis.plot(x_values, [100.0 * row["p_count_given_ko"] for row in subset], marker="s", color=self.KO_COLOR, label="KO-M")
            axis.set_yscale("log")
            axis.set_xticks(x_values)
            axis.set_xlabel("Observed raw UMI count")
            axis.set_ylabel("A-negative cells within genotype (%)")
            axis.set_title(title)
            axis.grid(color="#DDDDDD", linewidth=0.7)
        axes[0].legend(frameon=False)
        figure.suptitle("Step 05: Raw B/C molecule distributions among A-negative controls")
        figure.tight_layout()
        figure.savefig(path, dpi=self.dpi, metadata={"Title": "Step 05 A-negative raw B/C distributions"})
        plt.close(figure)

    def plot_enrichment(self, rows, path):
        bin_order = ["0", "1", "2", "3+"]
        x_values = np.arange(len(bin_order))
        figure, axis = plt.subplots(figsize=(9, 5.5))
        for measure in ["B_UMI", "C_UMI", "B_plus_C_UMI"]:
            indexed = {row["umi_bin"]: row for row in rows if row["measure"] == measure}
            label = {"B_UMI": "Probe B", "C_UMI": "Probe C", "B_plus_C_UMI": "B+C total"}[measure]
            axis.plot(x_values, [100.0 * indexed[name]["p_ko_given_bin"] for name in bin_order], marker="o", linewidth=2, color=self.MEASURE_COLORS[measure], label=label)
        axis.set_xticks(x_values)
        axis.set_xticklabels(bin_order)
        axis.set_ylim(0, 105)
        axis.set_xlabel("Observed raw UMI count bin")
        axis.set_ylabel("Known KO cells in bin (%)")
        axis.set_title("Step 05: KO enrichment with increasing raw B/C evidence")
        axis.grid(axis="y", color="#DDDDDD", linewidth=0.7)
        axis.legend(frameon=False)
        figure.tight_layout()
        figure.savefig(path, dpi=self.dpi, metadata={"Title": "Step 05 A-negative raw UMI KO enrichment"})
        plt.close(figure)

    def plot_joint(self, rows, path):
        b_values = [row["B_UMI"] for row in rows]
        c_values = [row["C_UMI"] for row in rows]
        totals = np.asarray([row["total_cells"] for row in rows], dtype=float)
        sizes = 35.0 + 180.0 * np.log10(totals + 1.0)
        colors = [100.0 * row["p_ko_given_joint"] for row in rows]
        figure, axis = plt.subplots(figsize=(8, 6.2))
        scatter = axis.scatter(b_values, c_values, s=sizes, c=colors, cmap="coolwarm", vmin=0, vmax=100, edgecolor="black", linewidth=0.5)
        for row in rows:
            if row["total_cells"] >= 20:
                axis.text(row["B_UMI"], row["C_UMI"], "{:.0f}%".format(100.0 * row["p_ko_given_joint"]), ha="center", va="center", fontsize=7)
        axis.set_xticks(range(max(b_values) + 1))
        axis.set_yticks(range(max(c_values) + 1))
        axis.set_xlabel("Raw Probe B UMI")
        axis.set_ylabel("Raw Probe C UMI")
        axis.set_title("Step 05: Joint B/C counts in A-negative controls\nBubble size = cells; color/label = P(KO | joint count)")
        figure.colorbar(scatter, ax=axis, label="Known KO cells in combination (%)")
        axis.grid(color="#EEEEEE", linewidth=0.5)
        figure.tight_layout()
        figure.savefig(path, dpi=self.dpi, metadata={"Title": "Step 05 A-negative joint B/C counts"})
        plt.close(figure)

    def render_all(self, distribution_rows, bin_rows, joint_rows, output_directory):
        paths = [
            os.path.join(output_directory, A_NEGATIVE_DISTRIBUTION_PLOT_NAME),
            os.path.join(output_directory, A_NEGATIVE_ENRICHMENT_PLOT_NAME),
            os.path.join(output_directory, A_NEGATIVE_JOINT_PLOT_NAME),
        ]
        self.plot_distributions(distribution_rows, paths[0])
        self.plot_enrichment(bin_rows, paths[1])
        self.plot_joint(joint_rows, paths[2])
        return paths


class CountModelComparisonPlotter(object):
    """Render the concise paired binary-versus-count model diagnostics."""

    BINARY_COLOR = "#667788"
    COUNT_COLOR = "#8E5BA6"

    def __init__(self, dpi):
        self.dpi = int(dpi)

    def plot_called_metrics(self, rows, path):
        wanted = ["accuracy_among_called", "ko_sensitivity_among_called", "ko_specificity_among_called"]
        labels = ["Accuracy", "KO sensitivity", "KO specificity"]
        indexed = {row["metric"]: row for row in rows}
        x_values = np.arange(len(wanted))
        width = 0.38
        figure, axis = plt.subplots(figsize=(8.5, 5.3))
        axis.bar(x_values - width / 2, [100.0 * indexed[key]["binary_model_value"] for key in wanted], width, color=self.BINARY_COLOR, label="Binary A/B/C")
        axis.bar(x_values + width / 2, [100.0 * indexed[key]["count_informed_model_value"] for key in wanted], width, color=self.COUNT_COLOR, label="A detected + raw B/C UMI")
        axis.set_xticks(x_values)
        axis.set_xticklabels(labels)
        axis.set_ylim(0, 105)
        axis.set_ylabel("Held-out called-cell metric (%)")
        axis.set_title("Step 05: Binary versus count-informed LOSO performance")
        axis.grid(axis="y", color="#DDDDDD", linewidth=0.7)
        axis.legend(frameon=False)
        figure.tight_layout()
        figure.savefig(path, dpi=self.dpi, metadata={"Title": "Step 05 binary versus count-informed metrics"})
        plt.close(figure)

    def plot_per_sample(self, rows, path):
        x_values = np.arange(len(rows))
        width = 0.38
        figure, axis = plt.subplots(figsize=(10, 5.5))
        axis.bar(x_values - width / 2, [100.0 * row["binary_accuracy_among_called"] for row in rows], width, color=self.BINARY_COLOR, label="Binary A/B/C")
        axis.bar(x_values + width / 2, [100.0 * row["count_informed_accuracy_among_called"] for row in rows], width, color=self.COUNT_COLOR, label="A detected + raw B/C UMI")
        axis.set_xticks(x_values)
        axis.set_xticklabels([row["biological_sample_id"] for row in rows], rotation=25, ha="right")
        axis.set_ylim(0, 105)
        axis.set_ylabel("Accuracy among called cells (%)")
        axis.set_title("Step 05: Paired model accuracy by held-out sample")
        axis.grid(axis="y", color="#DDDDDD", linewidth=0.7)
        axis.legend(frameon=False)
        figure.tight_layout()
        figure.savefig(path, dpi=self.dpi, metadata={"Title": "Step 05 paired per-sample model accuracy"})
        plt.close(figure)

    def plot_coefficients(self, rows, path):
        terms = [row["term"] for row in rows]
        odds_ratios = [row["odds_ratio"] for row in rows]
        colors = ["#3B6FB6" if value < 1.0 else "#D95F43" for value in odds_ratios]
        figure, axis = plt.subplots(figsize=(8, 5.2))
        axis.barh(terms, odds_ratios, color=colors)
        axis.axvline(1.0, color="black", linewidth=0.8)
        axis.set_xscale("log")
        axis.set_xlabel("KO:WT odds ratio (log scale)")
        axis.set_title("Step 05: Count-informed full-fit odds ratios\nB/C effects are per additional raw probe-level UMI")
        axis.grid(axis="x", color="#DDDDDD", linewidth=0.7, which="both")
        figure.tight_layout()
        figure.savefig(path, dpi=self.dpi, metadata={"Title": "Step 05 count-informed odds ratios"})
        plt.close(figure)

    def render_all(self, comparison_rows, sample_rows, coefficient_rows, output_directory):
        paths = [
            os.path.join(output_directory, COUNT_METRICS_PLOT_NAME),
            os.path.join(output_directory, COUNT_SAMPLE_PLOT_NAME),
            os.path.join(output_directory, COUNT_COEFFICIENT_PLOT_NAME),
        ]
        self.plot_called_metrics(comparison_rows, paths[0])
        self.plot_per_sample(sample_rows, paths[1])
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


class RawBCUMIEvidenceOutputPublisher(object):
    """Publish the descriptive A-negative raw B/C evidence package."""

    EXPECTED_FILES = {
        A_NEGATIVE_CELL_TABLE_NAME,
        A_NEGATIVE_COHORT_SUMMARY_NAME,
        A_NEGATIVE_DISTRIBUTION_NAME,
        A_NEGATIVE_BIN_ENRICHMENT_NAME,
        A_NEGATIVE_JOINT_NAME,
        A_NEGATIVE_DISTRIBUTION_PLOT_NAME,
        A_NEGATIVE_ENRICHMENT_PLOT_NAME,
        A_NEGATIVE_JOINT_PLOT_NAME,
        A_NEGATIVE_CHECKS_NAME,
        A_NEGATIVE_ENVIRONMENT_NAME,
    }

    def __init__(self, configuration, validation, plotter, output_root, paths):
        self.configuration = configuration
        self.validation = validation
        self.plotter = plotter
        self.output_root = os.path.abspath(output_root)
        self.paths = {key: os.path.abspath(value) for key, value in paths.items()}

    def _format(self, value):
        if value is None:
            return ""
        if isinstance(value, float) and math.isinf(value):
            return "Inf"
        return ("{:.%df}" % self.configuration.decimal_places).format(float(value))

    def _serialize(self, rows, float_keys):
        serialized = []
        for source in rows:
            row = dict(source)
            for key in float_keys:
                row[key] = self._format(row.get(key))
            serialized.append(row)
        return serialized

    def _write_environment(self, path):
        scope = self.configuration.values["a_negative_raw_bc_umi_evidence"]
        rows = [
            {"key": "step_id", "value": self.configuration.step_id},
            {"key": "pipeline_version", "value": self.configuration.values["pipeline_version"]},
            {"key": "analysis_module", "value": "a_negative_raw_bc_umi_evidence"},
            {"key": "cohort_definition", "value": scope["cohort_definition"]},
            {"key": "a_negative_definition", "value": scope["a_negative_definition"]},
            {"key": "count_representation", "value": scope["normalization"]},
            {"key": "measures", "value": ",".join(scope["measures"])},
            {"key": "count_bins", "value": ",".join(scope["count_bins"])},
            {"key": "classifier_fitted", "value": "false"},
            {"key": "threshold_optimized", "value": "false"},
            {"key": "het_female_cells", "value": "excluded_not_loaded"},
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
            if filename == A_NEGATIVE_MANIFEST_NAME or not os.path.isfile(path):
                continue
            rows.append({"relative_path": filename, "bytes": os.path.getsize(path), "sha256": sha256_file(path)})
        write_tsv(os.path.join(directory, A_NEGATIVE_MANIFEST_NAME), ["relative_path", "bytes", "sha256"], rows)

    @classmethod
    def _verify_manifest(cls, directory):
        manifest_path = os.path.join(directory, A_NEGATIVE_MANIFEST_NAME)
        if not os.path.isfile(manifest_path):
            raise Step05Error("A-negative raw B/C evidence manifest is missing")
        listed = set()
        with open(manifest_path, "r", newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                path = os.path.join(directory, row["relative_path"])
                if row["relative_path"] in listed or not os.path.isfile(path):
                    raise Step05Error("A-negative raw B/C evidence manifest is invalid")
                listed.add(row["relative_path"])
                if os.path.getsize(path) != int(row["bytes"]) or sha256_file(path) != row["sha256"]:
                    raise Step05Error("A-negative raw B/C evidence artifact identity mismatch")
        if listed != cls.EXPECTED_FILES:
            raise Step05Error("A-negative raw B/C evidence file set is unexpected")

    def _verify_existing(self):
        self._verify_manifest(self.output_root)
        with open(os.path.join(self.output_root, A_NEGATIVE_ENVIRONMENT_NAME), "r", newline="") as handle:
            environment = {row["key"]: row["value"] for row in csv.DictReader(handle, delimiter="\t")}
        for key, path in self.paths.items():
            if environment.get("{}_sha256".format(key)) != sha256_file(path):
                raise Step05Error("A-negative raw B/C evidence provenance differs for {}".format(key))

    def publish(self, records, cohort_rows, distribution_rows, bin_rows, joint_rows):
        if os.path.exists(self.output_root):
            self._verify_existing()
            print("EXISTING_VALIDATED\t{}".format(self.output_root))
            return
        parent = os.path.dirname(self.output_root)
        staging = tempfile.mkdtemp(prefix=".a_negative_raw_bc_umi_evidence.", dir=parent)
        try:
            write_tsv(os.path.join(staging, A_NEGATIVE_CELL_TABLE_NAME), A_NEGATIVE_CELL_HEADER, (record.as_row() for record in records))
            write_tsv(os.path.join(staging, A_NEGATIVE_COHORT_SUMMARY_NAME), A_NEGATIVE_COHORT_SUMMARY_HEADER, cohort_rows)
            probability_keys = {
                "p_wt_given_count", "p_ko_given_count", "p_count_given_wt", "p_count_given_ko", "ko_to_wt_likelihood_ratio"
            }
            write_tsv(os.path.join(staging, A_NEGATIVE_DISTRIBUTION_NAME), A_NEGATIVE_DISTRIBUTION_HEADER, self._serialize(distribution_rows, probability_keys))
            bin_probability_keys = {
                "p_wt_given_bin", "p_ko_given_bin", "ko_to_wt_cell_ratio", "p_bin_given_wt", "p_bin_given_ko", "ko_to_wt_likelihood_ratio"
            }
            write_tsv(os.path.join(staging, A_NEGATIVE_BIN_ENRICHMENT_NAME), A_NEGATIVE_BIN_HEADER, self._serialize(bin_rows, bin_probability_keys))
            joint_probability_keys = {
                "p_wt_given_joint", "p_ko_given_joint", "p_joint_given_wt", "p_joint_given_ko", "ko_to_wt_likelihood_ratio"
            }
            write_tsv(os.path.join(staging, A_NEGATIVE_JOINT_NAME), A_NEGATIVE_JOINT_HEADER, self._serialize(joint_rows, joint_probability_keys))
            plot_paths = self.plotter.render_all(distribution_rows, bin_rows, joint_rows, staging)
            self.validation.require_equal("a_negative_raw_bc_plots_created", sum(os.path.getsize(path) > 0 for path in plot_paths), 3)
            self.validation.require_equal("a_negative_classifier_fitted", False, False)
            write_tsv(os.path.join(staging, A_NEGATIVE_CHECKS_NAME), VALIDATION_HEADER, self.validation.rows)
            self._write_environment(os.path.join(staging, A_NEGATIVE_ENVIRONMENT_NAME))
            self._write_manifest(staging)
            self._verify_manifest(staging)
            if os.path.exists(self.output_root):
                raise Step05Error("A-negative raw B/C evidence output appeared during publication")
            os.replace(staging, self.output_root)
            staging = None
            print("PUBLISHED\t{}\tcells={}".format(self.output_root, len(records)))
        finally:
            if staging and os.path.exists(staging):
                shutil.rmtree(staging)


class CountModelComparisonOutputPublisher(object):
    """Publish the paired binary-versus-count-informed validation package."""

    EXPECTED_FILES = {
        COUNT_PAIRED_PREDICTIONS_NAME,
        COUNT_CONFUSION_NAME,
        COUNT_PER_SAMPLE_NAME,
        COUNT_GROUP_METRICS_NAME,
        COUNT_OVERALL_NAME,
        COUNT_FOLD_COEFFICIENTS_NAME,
        COUNT_FULL_COEFFICIENTS_NAME,
        COUNT_OVERALL_COMPARISON_NAME,
        COUNT_SAMPLE_COMPARISON_NAME,
        COUNT_TRANSITIONS_NAME,
        COUNT_METRICS_PLOT_NAME,
        COUNT_SAMPLE_PLOT_NAME,
        COUNT_COEFFICIENT_PLOT_NAME,
        COUNT_CHECKS_NAME,
        COUNT_ENVIRONMENT_NAME,
    }

    def __init__(self, configuration, validation, plotter, output_root, paths):
        self.configuration = configuration
        self.validation = validation
        self.plotter = plotter
        self.output_root = os.path.abspath(output_root)
        self.paths = {key: os.path.abspath(value) for key, value in paths.items()}

    def _format(self, value):
        if value is None:
            return ""
        return ("{:.%df}" % self.configuration.decimal_places).format(float(value))

    def _serialize(self, rows, float_keys):
        serialized = []
        for source in rows:
            row = dict(source)
            for key in float_keys:
                row[key] = self._format(row.get(key))
            serialized.append(row)
        return serialized

    def _write_environment(self, path):
        scope = self.configuration.values["count_informed_model_comparison"]
        rows = [
            {"key": "step_id", "value": self.configuration.step_id},
            {"key": "pipeline_version", "value": self.configuration.values["pipeline_version"]},
            {"key": "analysis_module", "value": "count_informed_vs_binary_validation"},
            {"key": "binary_model", "value": "A_detected+B_detected+C_detected"},
            {"key": "count_informed_model", "value": "+".join(scope["predictors"])},
            {"key": "model_formula", "value": scope["model_formula"]},
            {"key": "count_interpretation", "value": scope["count_interpretation"]},
            {"key": "normalization", "value": scope["normalization"]},
            {"key": "class_weighting", "value": "false"},
            {"key": "threshold_optimized", "value": "false"},
            {"key": "uncalled_state", "value": "A_UMI=B_UMI=C_UMI=0"},
            {"key": "het_female_cells", "value": "excluded_not_loaded"},
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
            if filename == COUNT_MANIFEST_NAME or not os.path.isfile(path):
                continue
            rows.append({"relative_path": filename, "bytes": os.path.getsize(path), "sha256": sha256_file(path)})
        write_tsv(os.path.join(directory, COUNT_MANIFEST_NAME), ["relative_path", "bytes", "sha256"], rows)

    @classmethod
    def _verify_manifest(cls, directory):
        manifest_path = os.path.join(directory, COUNT_MANIFEST_NAME)
        if not os.path.isfile(manifest_path):
            raise Step05Error("Count-informed model-comparison manifest is missing")
        listed = set()
        with open(manifest_path, "r", newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                path = os.path.join(directory, row["relative_path"])
                if row["relative_path"] in listed or not os.path.isfile(path):
                    raise Step05Error("Count-informed model-comparison manifest is invalid")
                listed.add(row["relative_path"])
                if os.path.getsize(path) != int(row["bytes"]) or sha256_file(path) != row["sha256"]:
                    raise Step05Error("Count-informed model-comparison artifact identity mismatch")
        if listed != cls.EXPECTED_FILES:
            raise Step05Error("Count-informed model-comparison file set is unexpected")

    def _verify_existing(self):
        self._verify_manifest(self.output_root)
        with open(os.path.join(self.output_root, COUNT_ENVIRONMENT_NAME), "r", newline="") as handle:
            environment = {row["key"]: row["value"] for row in csv.DictReader(handle, delimiter="\t")}
        for key, path in self.paths.items():
            if environment.get("{}_sha256".format(key)) != sha256_file(path):
                raise Step05Error("Count-informed package provenance differs for {}".format(key))

    def publish(
        self,
        paired_rows,
        confusion_rows,
        per_sample_rows,
        group_rows,
        overall_rows,
        fold_coefficient_rows,
        full_coefficient_rows,
        comparison_rows,
        sample_comparison_rows,
        transition_rows,
    ):
        if os.path.exists(self.output_root):
            self._verify_existing()
            print("EXISTING_VALIDATED\t{}".format(self.output_root))
            return
        parent = os.path.dirname(self.output_root)
        staging = tempfile.mkdtemp(prefix=".count_informed_vs_binary.", dir=parent)
        try:
            probability_keys = {
                "binary_predicted_wt_probability",
                "binary_predicted_ko_probability",
                "count_predicted_wt_probability",
                "count_predicted_ko_probability",
            }
            write_tsv(os.path.join(staging, COUNT_PAIRED_PREDICTIONS_NAME), COUNT_PAIRED_PREDICTION_HEADER, self._serialize(paired_rows, probability_keys))
            write_tsv(os.path.join(staging, COUNT_CONFUSION_NAME), HELD_OUT_CONFUSION_HEADER, confusion_rows)
            metric_keys = {"percent_called", "accuracy_among_called", "ko_sensitivity_where_defined", "ko_specificity_where_defined"}
            write_tsv(os.path.join(staging, COUNT_PER_SAMPLE_NAME), EXPANDED_PER_SAMPLE_HEADER, self._serialize(per_sample_rows, metric_keys))
            group_float_keys = {"percent_called", "accuracy_among_called", "ko_sensitivity_where_defined", "ko_specificity_where_defined"}
            write_tsv(os.path.join(staging, COUNT_GROUP_METRICS_NAME), GROUP_METRICS_HEADER, self._serialize(group_rows, group_float_keys))
            write_tsv(os.path.join(staging, COUNT_OVERALL_NAME), HELD_OUT_OVERALL_HEADER, self._serialize(overall_rows, {"value"}))
            write_tsv(os.path.join(staging, COUNT_FOLD_COEFFICIENTS_NAME), COUNT_FOLD_COEFFICIENT_HEADER, self._serialize(fold_coefficient_rows, {"coefficient", "odds_ratio"}))
            write_tsv(os.path.join(staging, COUNT_FULL_COEFFICIENTS_NAME), COEFFICIENT_HEADER, self._serialize(full_coefficient_rows, {"coefficient", "odds_ratio"}))
            comparison_float_keys = {"binary_model_value", "count_informed_model_value", "count_minus_binary"}
            write_tsv(os.path.join(staging, COUNT_OVERALL_COMPARISON_NAME), COUNT_OVERALL_COMPARISON_HEADER, self._serialize(comparison_rows, comparison_float_keys))
            sample_float_keys = {"binary_accuracy_among_called", "count_informed_accuracy_among_called", "count_minus_binary_accuracy"}
            write_tsv(os.path.join(staging, COUNT_SAMPLE_COMPARISON_NAME), COUNT_SAMPLE_COMPARISON_HEADER, self._serialize(sample_comparison_rows, sample_float_keys))
            write_tsv(os.path.join(staging, COUNT_TRANSITIONS_NAME), COUNT_TRANSITION_HEADER, transition_rows)
            plot_paths = self.plotter.render_all(comparison_rows, sample_comparison_rows, full_coefficient_rows, staging)
            self.validation.require_equal("count_model_comparison_plots_created", sum(os.path.getsize(path) > 0 for path in plot_paths), 3)
            self.validation.require_equal("count_model_het_cells_loaded", sum(row["true_genotype"] == "HET" for row in paired_rows), 0)
            write_tsv(os.path.join(staging, COUNT_CHECKS_NAME), VALIDATION_HEADER, self.validation.rows)
            self._write_environment(os.path.join(staging, COUNT_ENVIRONMENT_NAME))
            self._write_manifest(staging)
            self._verify_manifest(staging)
            if os.path.exists(self.output_root):
                raise Step05Error("Count-informed model-comparison output appeared during publication")
            os.replace(staging, self.output_root)
            staging = None
            print("PUBLISHED\t{}\tfolds=9\tcells={}".format(self.output_root, len(paired_rows)))
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
        self.expanded_package_verifier = ExistingExpandedGroundTruthPackageVerifier(
            self.configuration, self.validation, output_root
        )
        self.raw_bc_reader = ANegativeRawCountCohortReader(
            self.configuration,
            self.shared_configuration,
            self.validation,
            self.expanded_reader,
            self.step03_reader.table_path,
        )
        self.raw_bc_summarizer = RawBCUMIEvidenceSummarizer(
            self.configuration, self.validation
        )
        self.raw_bc_plotter = RawBCUMIEvidencePlotter(
            self.configuration.values["plot_dpi"]
        )
        raw_bc_root = os.path.join(
            output_root,
            self.configuration.values["a_negative_raw_bc_umi_evidence"]["relative_directory"],
        )
        raw_bc_paths = {
            "python_script": __file__,
            "lock": args.lock,
            "requirements": args.requirements,
            "shared_framework_script": self.configuration.framework_script_path,
            "step_04_lock": args.step04_lock,
            "sample_key": args.sample_key,
            "step_02a_manifest": self.expanded_reader.manifest_path,
            "step_03_manifest": self.step03_reader.manifest_path,
            "step_03_table": self.step03_reader.table_path,
            "step_05_base_manifest": os.path.join(output_root, MANIFEST_NAME),
            "male_only_validation_manifest": self.male_only_reader.manifest_path,
            "expanded_ground_truth_manifest": self.expanded_package_verifier.manifest_path,
        }
        self.raw_bc_publisher = RawBCUMIEvidenceOutputPublisher(
            self.configuration,
            self.validation,
            self.raw_bc_plotter,
            raw_bc_root,
            raw_bc_paths,
        )
        self.raw_bc_package_verifier = ExistingRawBCEvidencePackageVerifier(
            self.configuration, self.validation, output_root
        )
        count_scope = self.configuration.values["count_informed_model_comparison"]
        self.count_ground_truth_reader = CountInformedGroundTruthReader(
            self.configuration,
            self.encoder,
            self.validation,
            self.expanded_reader,
            self.expanded_package_verifier.predictions_path,
        )
        self.count_calling_policy = HeldOutCallingPolicy(count_scope["calling_rule"])
        self.count_validator = LeaveOneSampleOutCountLogisticValidator(
            self.configuration, self.validation, self.count_calling_policy
        )
        self.count_validation_evaluator = HeldOutValidationEvaluator(
            self.configuration,
            self.validation,
            validation_config=self.configuration.values["expanded_ground_truth_validation"],
            check_prefix="count_",
        )
        self.count_group_evaluator = ExpandedCohortEvaluator()
        self.paired_model_evaluator = PairedBinaryCountModelEvaluator(
            self.configuration, self.validation
        )
        self.count_plotter = CountModelComparisonPlotter(
            self.configuration.values["plot_dpi"]
        )
        count_root = os.path.join(output_root, count_scope["relative_directory"])
        count_paths = {
            "python_script": __file__,
            "lock": args.lock,
            "requirements": args.requirements,
            "shared_framework_script": self.configuration.framework_script_path,
            "step_04_lock": args.step04_lock,
            "sample_key": args.sample_key,
            "step_02a_manifest": self.expanded_reader.manifest_path,
            "step_03_manifest": self.step03_reader.manifest_path,
            "step_03_table": self.step03_reader.table_path,
            "step_05_base_manifest": os.path.join(output_root, MANIFEST_NAME),
            "male_only_validation_manifest": self.male_only_reader.manifest_path,
            "expanded_ground_truth_manifest": self.expanded_package_verifier.manifest_path,
            "binary_comparator_predictions": self.expanded_package_verifier.predictions_path,
            "a_negative_raw_bc_manifest": self.raw_bc_package_verifier.manifest_path,
        }
        self.count_publisher = CountModelComparisonOutputPublisher(
            self.configuration,
            self.validation,
            self.count_plotter,
            count_root,
            count_paths,
        )
        self.count_package_verifier = ExistingCountModelComparisonPackageVerifier(
            self.configuration, self.validation, output_root
        )
        diagnostic_scope = self.configuration.values["sample_level_probe_evidence_diagnostics"]
        diagnostic_root = os.path.join(output_root, diagnostic_scope["relative_directory"])
        diagnostic_paths = {
            "python_script": __file__,
            "diagnostic_module": sample_diagnostic_module.__file__,
            "lock": args.lock,
            "requirements": args.requirements,
            "sample_key": args.sample_key,
            "step_02a_manifest": self.expanded_reader.manifest_path,
            "step_03_manifest": self.step03_reader.manifest_path,
            "step_05_base_manifest": os.path.join(output_root, MANIFEST_NAME),
            "male_only_validation_manifest": self.male_only_reader.manifest_path,
            "expanded_ground_truth_manifest": self.expanded_package_verifier.manifest_path,
            "a_negative_raw_bc_manifest": self.raw_bc_package_verifier.manifest_path,
            "count_comparison_manifest": self.count_package_verifier.manifest_path,
            "paired_predictions": self.count_package_verifier.predictions_path,
            "stored_model_comparison": self.count_package_verifier.comparison_path,
        }
        self.sample_probe_diagnostic = SampleProbeEvidenceDiagnosticStep(
            self.configuration,
            self.validation,
            self.count_package_verifier.predictions_path,
            self.count_package_verifier.comparison_path,
            args.step02_root,
            diagnostic_root,
            diagnostic_paths,
        )

    def _run_raw_bc_evidence(self):
        if os.path.exists(self.raw_bc_publisher.output_root):
            self.raw_bc_package_verifier.verify()
            return
        records = self.raw_bc_reader.read_records()
        cohort_rows, distribution_rows, bin_rows, joint_rows = self.raw_bc_summarizer.summarize(records)
        self.raw_bc_publisher.publish(
            records, cohort_rows, distribution_rows, bin_rows, joint_rows
        )

    def _run_count_model_comparison(self):
        if os.path.exists(self.count_publisher.output_root):
            self.count_package_verifier.verify()
            return
        records, binary_by_key = self.count_ground_truth_reader.read_records()
        binary_predictions = [binary_by_key[(record.biological_sample_id, record.cell_barcode)] for record in records]
        count_predictions, fold_coefficient_rows = self.count_validator.run(records)
        confusion_rows, per_sample_rows, overall_rows = self.count_validation_evaluator.evaluate(count_predictions)
        group_rows = self.count_group_evaluator.group_metrics(count_predictions)
        full_classifier = CountInformedLogisticEstimator(self.configuration).fit(records)
        self.validation.require_equal("count_full_model_converged", full_classifier.fit_diagnostics["converged"], True)
        full_coefficient_rows = full_classifier.coefficient_rows()
        paired_rows, comparison_rows, sample_comparison_rows, transition_rows = self.paired_model_evaluator.compare(
            binary_predictions, count_predictions
        )
        self.count_publisher.publish(
            paired_rows,
            confusion_rows,
            per_sample_rows,
            group_rows,
            overall_rows,
            fold_coefficient_rows,
            full_coefficient_rows,
            comparison_rows,
            sample_comparison_rows,
            transition_rows,
        )

    def _run_sample_probe_diagnostic(self):
        self.count_package_verifier.verify()
        self.sample_probe_diagnostic.run()

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
        if (
            os.path.exists(self.held_out_publisher.output_root)
            and os.path.exists(self.expanded_publisher.output_root)
        ):
            self.male_only_reader.verify()
            self.expanded_package_verifier.verify()
            self._run_raw_bc_evidence()
            self._run_count_model_comparison()
            self._run_sample_probe_diagnostic()
            return
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
            self.expanded_package_verifier.verify()
            self._run_raw_bc_evidence()
            self._run_count_model_comparison()
            self._run_sample_probe_diagnostic()
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
        self._run_raw_bc_evidence()
        self._run_count_model_comparison()
        self._run_sample_probe_diagnostic()


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
    except (Step05Error, Step04Error, DiagnosticError) as error:
        sys.stderr.write("ERROR: {}\n".format(error))
        sys.exit(1)
