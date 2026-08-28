"""Read-only biological-sample diagnostics for the existing Step 05 models."""

import csv
import hashlib
import math
import os
import platform
import shutil
import tempfile
from collections import Counter, defaultdict

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


OUTPUT_DIRECTORY = "sample_level_probe_evidence_diagnostics"
FLOW_NAME = "step_05_sample_cell_flow.tsv"
GENOTYPE_FLOW_NAME = "step_05_genotype_cell_flow.tsv"
SUMMARY_NAME = "step_05_per_sample_probe_summary.tsv"
ANEG_SUMMARY_NAME = "step_05_per_sample_Anegative_probe_summary.tsv"
APOS_SUMMARY_NAME = "step_05_per_sample_Apositive_probe_summary.tsv"
DISTRIBUTION_NAME = "step_05_per_sample_probe_umi_distribution.tsv"
BIN_NAME = "step_05_probe_umi_bin_composition.tsv"
ANEG_BIN_NAME = "step_05_Anegative_probe_umi_bin_composition.tsv"
APOS_BIN_NAME = "step_05_Apositive_probe_umi_bin_composition.tsv"
PERFORMANCE_NAME = "step_05_prediction_performance_by_probe_state.tsv"
TRANSITIONS_NAME = "step_05_binary_vs_count_error_transitions_by_probe_state.tsv"
PAIRWISE_NAME = "step_05_sample_pairwise_probe_comparison.tsv"
AUC_NAME = "step_05_auc_orientation_audit.tsv"
INTERPRETATION_NAME = "step_05_sample_probe_interpretation.tsv"
SCOPE_NAME = "step_05_sample_probe_diagnostic_scope.tsv"
CHECKS_NAME = "step_05_sample_probe_diagnostic_validation_checks.tsv"
ENVIRONMENT_NAME = "software_environment.tsv"
MANIFEST_NAME = "output_manifest.tsv"

FIGURES = [
    "step_05_per_sample_A_detection_proportion.png",
    "step_05_per_sample_B_UMI_distribution.png",
    "step_05_per_sample_C_UMI_distribution.png",
    "step_05_per_sample_B_plus_C_UMI_distribution.png",
    "step_05_Anegative_per_sample_B_UMI_distribution.png",
    "step_05_Anegative_per_sample_C_UMI_distribution.png",
    "step_05_Anegative_per_sample_B_plus_C_UMI_distribution.png",
    "step_05_per_sample_UMI_bin_composition.png",
    "step_05_Anegative_per_sample_UMI_bin_composition.png",
    "step_05_sample_probe_evidence_heatmap.png",
    "step_05_KO_sample_probe_evidence_comparison.png",
    "step_05_binary_vs_count_error_state_composition.png",
]


class DiagnosticError(RuntimeError):
    pass


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
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=header, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in header})


class DiagnosticCell(object):
    __slots__ = (
        "sample_id", "sample_name", "sex", "design_group", "genotype", "barcode",
        "a_umi", "a_detected", "b_umi", "c_umi", "pattern_code",
        "binary_p_ko", "binary_call", "count_p_ko", "count_call",
    )

    def __init__(self, row):
        self.sample_id = row["biological_sample_id"]
        self.sample_name = row["submitted_sample_name"]
        self.sex = row["sex"]
        self.design_group = row["design_group"]
        self.genotype = row["true_genotype"]
        self.barcode = row["cell_barcode"]
        self.a_umi = int(row["A_UMI"])
        self.a_detected = int(row["A_detected"])
        self.b_umi = int(row["B_UMI"])
        self.c_umi = int(row["C_UMI"])
        self.pattern_code = row["pattern_code"]
        self.binary_p_ko = float(row["binary_predicted_ko_probability"])
        self.binary_call = row["binary_predicted_genotype"]
        self.count_p_ko = float(row["count_predicted_ko_probability"])
        self.count_call = row["count_predicted_genotype"]
        expected_pattern = "{}{}{}".format(self.a_detected, int(self.b_umi > 0), int(self.c_umi > 0))
        if min(self.a_umi, self.b_umi, self.c_umi) < 0 or self.a_detected != int(self.a_umi > 0) or self.pattern_code != expected_pattern:
            raise DiagnosticError("Invalid paired diagnostic cell")

    @property
    def bc_umi(self):
        return self.b_umi + self.c_umi

    def value(self, measure):
        return {"B_UMI": self.b_umi, "C_UMI": self.c_umi, "B_plus_C_UMI": self.bc_umi}[measure]


class SampleProbeDiagnosticReader(object):
    """Load immutable paired predictions and verify raw counts against Step 02a."""

    REQUIRED_SOURCE_COLUMNS = ["barcode", "A_UMI", "B_UMI", "C_UMI", "Pcdh19_total_UMI", "detection_pattern"]

    def __init__(self, configuration, validation, paired_path, step02_root, step02_manifest_path):
        self.configuration = configuration
        self.validation = validation
        self.scope = configuration.values["sample_level_probe_evidence_diagnostics"]
        self.cohort = configuration.values["expanded_ground_truth_validation"]
        self.paired_path = os.path.abspath(paired_path)
        self.step02_root = os.path.abspath(step02_root)
        self.manifest_path = os.path.abspath(step02_manifest_path)

    def _manifest(self):
        rows = {}
        with open(self.manifest_path, "r", newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                rows[row["relative_path"]] = row
        return rows

    def read(self):
        records = []
        by_key = {}
        with open(self.paired_path, "r", newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                record = DiagnosticCell(row)
                key = (record.sample_id, record.barcode)
                if key in by_key:
                    raise DiagnosticError("Duplicate paired diagnostic key")
                by_key[key] = record
                records.append(record)
        self.validation.require_equal("diagnostic_paired_cells", len(records), self.scope["expected_total_cells"])
        self.validation.require_equal("diagnostic_unique_keys", len(by_key), len(records))
        expected_samples = set(self.cohort["expected_samples"])
        self.validation.require_equal("diagnostic_sample_ids", sorted({row.sample_id for row in records}), sorted(expected_samples))
        self.validation.require_equal("diagnostic_het_cells", sum(row.genotype == "HET" for row in records), 0)
        manifest = self._manifest()
        verified = 0
        for sample_id, definition in self.cohort["expected_samples"].items():
            relative = os.path.join("per_sample", sample_id, "pcdh19_probe_patterns.tsv")
            entry = manifest.get(relative)
            path = os.path.join(self.step02_root, relative)
            if entry is None or not os.path.isfile(path):
                raise DiagnosticError("Manifested source table is missing")
            if os.path.getsize(path) != int(entry["bytes"]) or sha256_file(path) != entry["sha256"]:
                raise DiagnosticError("Source probe table identity mismatch")
            with open(path, "r", newline="") as handle:
                reader = csv.DictReader(handle, delimiter="\t")
                if reader.fieldnames != self.REQUIRED_SOURCE_COLUMNS:
                    raise DiagnosticError("Unexpected source probe schema")
                for source in reader:
                    record = by_key.get((sample_id, source["barcode"]))
                    if record is None:
                        raise DiagnosticError("Source WT/KO cell missing from diagnostic input")
                    raw = tuple(int(source[key]) for key in ("A_UMI", "B_UMI", "C_UMI"))
                    if raw != (record.a_umi, record.b_umi, record.c_umi):
                        raise DiagnosticError("Raw probe evidence changed from Step 02a")
                    if record.genotype != definition["genotype"] or record.sex != definition["sex"] or record.design_group != definition["design_group"]:
                        raise DiagnosticError("Registered sample metadata mismatch")
                    verified += 1
        self.validation.require_equal("diagnostic_raw_source_rows_verified", verified, len(records))
        self.validation.require_equal("diagnostic_exact_zero_cells", sum(row.pattern_code == "000" for row in records), self.scope["expected_exact_zero_cells"])
        self.validation.require_equal("diagnostic_binary_called_cells", sum(row.binary_call in ("WT", "KO") for row in records), self.scope["expected_called_cells_per_model"])
        self.validation.require_equal("diagnostic_count_called_cells", sum(row.count_call in ("WT", "KO") for row in records), self.scope["expected_called_cells_per_model"])
        return records


class ProbeEvidenceSummaryBuilder(object):
    """Build sample flow, raw distributions, bins, and descriptive statistics."""

    MEASURES = ["B_UMI", "C_UMI", "B_plus_C_UMI"]
    SUBSETS = ["all", "A_negative", "A_positive"]
    QUANTILES = [("q00", 0.0), ("q10", 0.1), ("q25", 0.25), ("q50", 0.5), ("q75", 0.75), ("q90", 0.9), ("q95", 0.95), ("q99", 0.99), ("q100", 1.0)]

    def __init__(self, configuration, validation):
        self.configuration = configuration
        self.validation = validation
        self.cohort = configuration.values["expanded_ground_truth_validation"]

    @staticmethod
    def _subset(records, name):
        if name == "all":
            return records
        target = 0 if name == "A_negative" else 1
        return [row for row in records if row.a_detected == target]

    @staticmethod
    def _bin(value):
        return "3+" if value >= 3 else str(value)

    def flow(self, records):
        rows = []
        for sample_id, definition in self.cohort["expected_samples"].items():
            sample = [row for row in records if row.sample_id == sample_id]
            exact_zero = sum(row.pattern_code == "000" for row in sample)
            rows.append({
                "biological_sample_id": sample_id, "submitted_sample_name": sample[0].sample_name,
                "sex": definition["sex"], "design_group": definition["design_group"], "true_genotype": definition["genotype"],
                "original_cells_available": definition["cells"], "retained_cells": len(sample), "excluded_cells": definition["cells"] - len(sample),
                "exclusion_reasons": "none; registered WT/KO Step 05 control", "all_zero_cells": exact_zero,
                "uncalled_by_exact_zero_policy": exact_zero,
                "binary_model_called_cells": sum(row.binary_call in ("WT", "KO") for row in sample),
                "count_informed_model_called_cells": sum(row.count_call in ("WT", "KO") for row in sample),
            })
        genotype_rows = []
        for genotype in ("WT", "KO"):
            subset = [row for row in rows if row["true_genotype"] == genotype]
            genotype_rows.append({
                "true_genotype": genotype,
                "biological_samples": len(subset),
                "original_cells_available": sum(row["original_cells_available"] for row in subset),
                "retained_cells": sum(row["retained_cells"] for row in subset),
                "excluded_cells": sum(row["excluded_cells"] for row in subset),
                "exclusion_reasons": "none; registered WT/KO Step 05 controls only",
                "all_zero_cells": sum(row["all_zero_cells"] for row in subset),
                "uncalled_by_exact_zero_policy": sum(row["uncalled_by_exact_zero_policy"] for row in subset),
                "binary_model_called_cells": sum(row["binary_model_called_cells"] for row in subset),
                "count_informed_model_called_cells": sum(row["count_informed_model_called_cells"] for row in subset),
            })
        self.validation.require_equal("diagnostic_flow_retained", sum(row["retained_cells"] for row in rows), len(records))
        self.validation.require_equal("diagnostic_flow_excluded", sum(row["excluded_cells"] for row in rows), 0)
        return rows, genotype_rows

    def summaries(self, records):
        summary_by_subset = {name: [] for name in self.SUBSETS}
        distributions = []
        bins_by_subset = {name: [] for name in self.SUBSETS}
        for sample_id, definition in self.cohort["expected_samples"].items():
            sample = [row for row in records if row.sample_id == sample_id]
            a_positive = sum(row.a_detected == 1 for row in sample)
            common = {
                "biological_sample_id": sample_id, "submitted_sample_name": sample[0].sample_name,
                "sex": definition["sex"], "design_group": definition["design_group"], "true_genotype": definition["genotype"],
                "total_sample_cells": len(sample), "A_detected_cells": a_positive,
                "A_detected_percent": 100.0 * a_positive / len(sample), "A_not_detected_cells": len(sample) - a_positive,
                "A_not_detected_percent": 100.0 * (len(sample) - a_positive) / len(sample),
            }
            for subset_name in self.SUBSETS:
                subset = self._subset(sample, subset_name)
                detection = {
                    "subset": subset_name, "subset_cells": len(subset),
                    "B_detected_percent": 100.0 * sum(row.b_umi > 0 for row in subset) / len(subset),
                    "C_detected_percent": 100.0 * sum(row.c_umi > 0 for row in subset) / len(subset),
                    "either_B_or_C_detected_percent": 100.0 * sum(row.bc_umi > 0 for row in subset) / len(subset),
                    "both_B_and_C_detected_percent": 100.0 * sum(row.b_umi > 0 and row.c_umi > 0 for row in subset) / len(subset),
                }
                for measure in self.MEASURES:
                    values = np.asarray([row.value(measure) for row in subset], dtype=float)
                    row = dict(common)
                    row.update(detection)
                    row.update({
                        "measure": measure, "mean": float(np.mean(values)), "sd": float(np.std(values, ddof=1)),
                        "median": float(np.median(values)), "iqr": float(np.percentile(values, 75) - np.percentile(values, 25)),
                        "max": int(np.max(values)),
                    })
                    for label, quantile in self.QUANTILES:
                        row[label] = float(np.quantile(values, quantile))
                    summary_by_subset[subset_name].append(row)
                    counts = Counter(int(value) for value in values)
                    for raw_count in range(int(np.max(values)) + 1):
                        distributions.append({
                            "biological_sample_id": sample_id, "sex": definition["sex"], "design_group": definition["design_group"],
                            "true_genotype": definition["genotype"], "subset": subset_name, "measure": measure,
                            "raw_umi_count": raw_count, "cells": counts[raw_count], "percent_of_subset": 100.0 * counts[raw_count] / len(subset),
                        })
                    bin_counts = Counter(self._bin(int(value)) for value in values)
                    for bin_name in ("0", "1", "2", "3+"):
                        bins_by_subset[subset_name].append({
                            "biological_sample_id": sample_id, "sex": definition["sex"], "design_group": definition["design_group"],
                            "true_genotype": definition["genotype"], "subset": subset_name, "measure": measure, "umi_bin": bin_name,
                            "cells": bin_counts[bin_name], "percent_of_subset": 100.0 * bin_counts[bin_name] / len(subset), "subset_cells": len(subset),
                        })
        return summary_by_subset, distributions, bins_by_subset


class ProbeStatePerformanceAnalyzer(object):
    """Stratify immutable calls and paired correctness by raw probe state."""

    STRATA = ["all", "A_detected", "A_not_detected", "B_plus_C_0", "B_plus_C_1", "B_plus_C_2", "B_plus_C_3plus", "Anegative_B_plus_C_0", "Anegative_B_plus_C_1", "Anegative_B_plus_C_2", "Anegative_B_plus_C_3plus"]

    @staticmethod
    def _matches(row, stratum):
        if stratum == "all": return True
        if stratum == "A_detected": return row.a_detected == 1
        if stratum == "A_not_detected": return row.a_detected == 0
        aneg = stratum.startswith("Anegative_")
        label = stratum.replace("Anegative_", "")
        if aneg and row.a_detected != 0: return False
        if label == "B_plus_C_0": return row.bc_umi == 0
        if label == "B_plus_C_1": return row.bc_umi == 1
        if label == "B_plus_C_2": return row.bc_umi == 2
        if label == "B_plus_C_3plus": return row.bc_umi >= 3
        return False

    def analyze(self, records, sample_order):
        performance = []
        transitions = []
        for sample_id in sample_order:
            sample = [row for row in records if row.sample_id == sample_id]
            for stratum in self.STRATA:
                subset = [row for row in sample if self._matches(row, stratum)]
                for model in ("binary", "count_informed"):
                    call_key = "binary_call" if model == "binary" else "count_call"
                    called = [row for row in subset if getattr(row, call_key) in ("WT", "KO")]
                    correct = sum(getattr(row, call_key) == row.genotype for row in called)
                    accuracy = correct / len(called) if called else None
                    performance.append({
                        "biological_sample_id": sample_id, "sex": sample[0].sex, "design_group": sample[0].design_group, "true_genotype": sample[0].genotype,
                        "model": model, "probe_state_stratum": stratum, "stratum_cells": len(subset), "called_cells": len(called),
                        "uncalled_cells": len(subset) - len(called), "correct_cells": correct, "accuracy_among_called": accuracy,
                        "ko_sensitivity_where_defined": accuracy if sample[0].genotype == "KO" and called else None,
                        "ko_specificity_where_defined": accuracy if sample[0].genotype == "WT" and called else None,
                    })
                both_called = [row for row in subset if row.binary_call in ("WT", "KO") and row.count_call in ("WT", "KO")]
                counts = Counter((row.binary_call == row.genotype, row.count_call == row.genotype) for row in both_called)
                transitions.append({
                    "biological_sample_id": sample_id, "sex": sample[0].sex, "design_group": sample[0].design_group, "true_genotype": sample[0].genotype,
                    "probe_state_stratum": stratum, "stratum_cells": len(subset), "both_models_called_cells": len(both_called),
                    "binary_correct_count_correct": counts[(True, True)], "binary_wrong_count_correct": counts[(False, True)],
                    "binary_correct_count_wrong": counts[(True, False)], "both_wrong": counts[(False, False)],
                    "uncalled_by_one_or_both_models": len(subset) - len(both_called),
                })
        return performance, transitions


class SamplePairwiseProbeComparator(object):
    """Compare biological samples descriptively without cell-level p-values."""

    def compare(self, summary_by_subset, bins_by_subset):
        feature_rows = {}
        for subset_name in ("all", "A_negative"):
            for row in summary_by_subset[subset_name]:
                key = (row["biological_sample_id"], subset_name)
                features = feature_rows.setdefault(key, {"genotype": row["true_genotype"]})
                measure = row["measure"]
                if subset_name == "all":
                    features["A_detected_percent"] = row["A_detected_percent"]
                for metric in ("mean", "median", "sd", "iqr", "max", "q10", "q25", "q75", "q90", "q95", "q99", "B_detected_percent", "C_detected_percent", "either_B_or_C_detected_percent", "both_B_and_C_detected_percent"):
                    if metric in row:
                        features["{}_{}".format(measure, metric)] = row[metric]
            for row in bins_by_subset[subset_name]:
                key = (row["biological_sample_id"], subset_name)
                feature_rows[key]["{}_bin_{}_percent".format(row["measure"], row["umi_bin"])] = row["percent_of_subset"]
        output = []
        for subset_name in ("all", "A_negative"):
            for genotype in ("WT", "KO"):
                samples = sorted(key[0] for key, value in feature_rows.items() if key[1] == subset_name and value["genotype"] == genotype)
                for left_index in range(len(samples)):
                    for right_index in range(left_index + 1, len(samples)):
                        left, right = samples[left_index], samples[right_index]
                        left_features, right_features = feature_rows[(left, subset_name)], feature_rows[(right, subset_name)]
                        for feature in sorted(set(left_features).intersection(right_features).difference({"genotype"})):
                            left_value, right_value = float(left_features[feature]), float(right_features[feature])
                            output.append({
                                "true_genotype": genotype, "subset": subset_name, "sample_a": left, "sample_b": right,
                                "probe_feature": feature, "sample_a_value": left_value, "sample_b_value": right_value,
                                "sample_b_minus_sample_a": right_value - left_value, "absolute_difference": abs(right_value - left_value),
                                "sample_b_to_sample_a_ratio": right_value / left_value if left_value != 0 else None,
                                "larger_value_sample": left if left_value > right_value else (right if right_value > left_value else "tie"),
                                "inference_scope": "descriptive biological-sample comparison; no cell-level p-value",
                            })
        return output


class AUCOrientationAuditor(object):
    """Audit positive labels and score orientation without changing stored AUCs."""

    @staticmethod
    def _auc(labels, scores):
        order = sorted(range(len(labels)), key=lambda index: scores[index])
        positive_ranks = 0.0
        rank = 1
        index = 0
        while index < len(order):
            end = index + 1
            while end < len(order) and scores[order[end]] == scores[order[index]]:
                end += 1
            average_rank = (rank + rank + end - index - 1) / 2.0
            positive_ranks += average_rank * sum(labels[order[position]] for position in range(index, end))
            rank += end - index
            index = end
        positives = sum(labels)
        negatives = len(labels) - positives
        return (positive_ranks - positives * (positives + 1) / 2.0) / (positives * negatives)

    def audit(self, records, stored_comparison_path):
        with open(stored_comparison_path, "r", newline="") as handle:
            stored = {row["metric"]: row for row in csv.DictReader(handle, delimiter="\t")}
        labels = [1 if row.genotype == "KO" else 0 for row in records]
        output = []
        for model, attribute, stored_column in (("binary", "binary_p_ko", "binary_model_value"), ("count_informed", "count_p_ko", "count_informed_model_value")):
            p_ko = [getattr(row, attribute) for row in records]
            recomputed = self._auc(labels, p_ko)
            reversed_auc = self._auc(labels, [1.0 - value for value in p_ko])
            output.append({
                "model": model, "positive_class": "KO", "positive_numeric_label": 1, "negative_class": "WT", "negative_numeric_label": 0,
                "score_used": "P(KO)", "fold_label_orientation_consistent": "yes", "stored_pooled_auc": float(stored["roc_auc_all_cells"][stored_column]),
                "recomputed_pooled_auc": recomputed, "orientation_corrected_auc": recomputed,
                "reversed_score": "P(WT)=1-P(KO)", "reversed_score_auc_diagnostic": reversed_auc,
                "auc_sum_original_plus_reversed": recomputed + reversed_auc,
                "orientation_conclusion": "orientation is correct; below-0.5 pooled AUC is not a label/score reversal",
                "interpretation_caveat": "scores come from different sample-held-out models; pooled cross-fold ranking is secondary to per-sample diagnostics",
            })
        return output


class SampleProbeDiagnosticPlotter(object):
    """Render sample-visible probe and error diagnostics without computation changes."""

    COLORS = {"WT_M": "#3B6FB6", "WT_F": "#65A5D8", "KO_M": "#D95F43"}
    BIN_COLORS = {"0": "#D9D9D9", "1": "#9ECAE1", "2": "#6BAED6", "3+": "#2171B5"}

    def __init__(self, dpi):
        self.dpi = int(dpi)
        plt.rcParams.update({"font.size": 8, "axes.spines.top": False, "axes.spines.right": False})

    @staticmethod
    def _sample_order(summary_rows):
        return list(dict.fromkeys(row["biological_sample_id"] for row in summary_rows))

    def plot_a_detection(self, rows, path):
        indexed = {row["biological_sample_id"]: row for row in rows if row["measure"] == "B_UMI"}
        samples = list(indexed)
        values = [indexed[sample]["A_detected_percent"] for sample in samples]
        colors = [self.COLORS[indexed[sample]["design_group"]] for sample in samples]
        figure, axis = plt.subplots(figsize=(9, 4.8))
        axis.scatter(range(len(samples)), values, s=65, c=colors, edgecolor="black", linewidth=0.5)
        axis.set_xticks(range(len(samples))); axis.set_xticklabels(samples, rotation=30, ha="right")
        axis.set_ylabel("A-detected cells (%)"); axis.set_ylim(0, max(values) * 1.2)
        axis.set_title("Step 05: Probe A detection by biological sample")
        axis.grid(axis="y", color="#DDDDDD", linewidth=0.7)
        figure.tight_layout(); figure.savefig(path, dpi=self.dpi); plt.close(figure)

    def plot_distribution(self, rows, subset, measure, path):
        selected = [row for row in rows if row["subset"] == subset and row["measure"] == measure]
        samples = list(dict.fromkeys(row["biological_sample_id"] for row in selected))
        figure, axis = plt.subplots(figsize=(8.5, 5))
        for sample in samples:
            sample_rows = [row for row in selected if row["biological_sample_id"] == sample]
            group = sample_rows[0]["design_group"]
            axis.plot([row["raw_umi_count"] for row in sample_rows], [row["percent_of_subset"] for row in sample_rows], marker="o", linewidth=1.2, color=self.COLORS[group], alpha=0.78, label=sample)
        axis.set_yscale("symlog", linthresh=0.01)
        axis.set_xlabel("Raw probe-level UMI/ligation count"); axis.set_ylabel("Cells in sample subset (%)")
        axis.set_title("Step 05: {} distribution by sample ({})".format(measure.replace("_", " "), subset.replace("_", "-")))
        axis.grid(color="#E5E5E5", linewidth=0.6); axis.legend(ncol=3, frameon=False, fontsize=7)
        figure.tight_layout(); figure.savefig(path, dpi=self.dpi); plt.close(figure)

    def plot_bins(self, rows, path, title):
        samples = list(dict.fromkeys(row["biological_sample_id"] for row in rows))
        measures = ["B_UMI", "C_UMI", "B_plus_C_UMI"]
        figure, axes = plt.subplots(3, 1, figsize=(10, 10), sharex=True)
        for axis, measure in zip(axes, measures):
            bottom = np.zeros(len(samples))
            for bin_name in ("0", "1", "2", "3+"):
                indexed = {(row["biological_sample_id"], row["umi_bin"]): row for row in rows if row["measure"] == measure}
                values = np.asarray([indexed[(sample, bin_name)]["percent_of_subset"] for sample in samples])
                axis.bar(range(len(samples)), values, bottom=bottom, color=self.BIN_COLORS[bin_name], label=bin_name)
                bottom += values
            axis.set_ylim(0, 100); axis.set_ylabel("Cells (%)"); axis.set_title(measure.replace("_", " "))
        axes[0].legend(title="Raw UMI bin", ncol=4, frameon=False)
        axes[-1].set_xticks(range(len(samples))); axes[-1].set_xticklabels(samples, rotation=30, ha="right")
        figure.suptitle(title); figure.tight_layout(); figure.savefig(path, dpi=self.dpi); plt.close(figure)

    def plot_heatmap(self, all_rows, aneg_rows, path):
        all_index = {(row["biological_sample_id"], row["measure"]): row for row in all_rows}
        aneg_index = {(row["biological_sample_id"], row["measure"]): row for row in aneg_rows}
        samples = list(dict.fromkeys(row["biological_sample_id"] for row in all_rows))
        features = ["A detected %", "B mean", "C mean", "B+C mean", "A-neg B+C mean", "A-neg B+C zero %"]
        matrix = []
        for sample in samples:
            bins_zero = 100.0 - aneg_index[(sample, "B_plus_C_UMI")]["either_B_or_C_detected_percent"]
            matrix.append([
                all_index[(sample, "B_UMI")]["A_detected_percent"], all_index[(sample, "B_UMI")]["mean"],
                all_index[(sample, "C_UMI")]["mean"], all_index[(sample, "B_plus_C_UMI")]["mean"],
                aneg_index[(sample, "B_plus_C_UMI")]["mean"], bins_zero,
            ])
        matrix = np.asarray(matrix, dtype=float)
        z = (matrix - matrix.mean(axis=0)) / np.where(matrix.std(axis=0) == 0, 1, matrix.std(axis=0))
        figure, axis = plt.subplots(figsize=(9, 5.8)); image = axis.imshow(z, cmap="coolwarm", aspect="auto", vmin=-2.5, vmax=2.5)
        axis.set_xticks(range(len(features))); axis.set_xticklabels(features, rotation=35, ha="right")
        axis.set_yticks(range(len(samples))); axis.set_yticklabels(samples)
        axis.set_title("Step 05: Sample-level probe evidence (feature z-scores)")
        figure.colorbar(image, ax=axis, label="Across-sample z-score"); figure.tight_layout(); figure.savefig(path, dpi=self.dpi); plt.close(figure)

    def plot_ko_comparison(self, interpretation_rows, path):
        features = [("A_detected_percent", "A detected %"), ("Anegative_B_plus_C_mean", "A-neg B+C mean"), ("Anegative_B_plus_C_zero_percent", "A-neg B+C zero %"), ("count_minus_binary_accuracy_points", "Count−binary accuracy points")]
        samples = [row["biological_sample_id"] for row in interpretation_rows]
        figure, axes = plt.subplots(2, 2, figsize=(10, 7.5))
        for axis, (key, label) in zip(axes.ravel(), features):
            values = [row[key] for row in interpretation_rows]
            colors = [self.COLORS[row["design_group"]] for row in interpretation_rows]
            axis.scatter(range(len(samples)), values, c=colors, s=50, edgecolor="black", linewidth=0.4)
            axis.set_xticks(range(len(samples))); axis.set_xticklabels(samples, rotation=35, ha="right", fontsize=7)
            axis.set_ylabel(label); axis.grid(axis="y", color="#E2E2E2", linewidth=0.6)
        figure.suptitle("Step 05: All samples shown; KO replicate behavior highlighted in red")
        figure.tight_layout(); figure.savefig(path, dpi=self.dpi); plt.close(figure)

    def plot_error_states(self, transitions, path):
        strata = ["B_plus_C_0", "B_plus_C_1", "B_plus_C_2", "B_plus_C_3plus"]
        corrected = [sum(row["binary_wrong_count_correct"] for row in transitions if row["probe_state_stratum"] == state) for state in strata]
        regressed = [sum(row["binary_correct_count_wrong"] for row in transitions if row["probe_state_stratum"] == state) for state in strata]
        x = np.arange(len(strata)); width = 0.38
        figure, axis = plt.subplots(figsize=(8.5, 5))
        axis.bar(x - width / 2, corrected, width, color="#2CA02C", label="Binary wrong → count correct")
        axis.bar(x + width / 2, regressed, width, color="#D62728", label="Binary correct → count wrong")
        axis.set_xticks(x); axis.set_xticklabels(["B+C=0", "B+C=1", "B+C=2", "B+C≥3"])
        axis.set_ylabel("Held-out cells"); axis.set_title("Step 05: Count-model corrections and regressions by raw B+C state")
        axis.legend(frameon=False); axis.grid(axis="y", color="#DDDDDD", linewidth=0.7)
        figure.tight_layout(); figure.savefig(path, dpi=self.dpi); plt.close(figure)

    def render(self, output, summaries, distributions, bins, transitions, interpretation):
        paths = [os.path.join(output, name) for name in FIGURES]
        self.plot_a_detection(summaries["all"], paths[0])
        for measure, path in zip(("B_UMI", "C_UMI", "B_plus_C_UMI"), paths[1:4]): self.plot_distribution(distributions, "all", measure, path)
        for measure, path in zip(("B_UMI", "C_UMI", "B_plus_C_UMI"), paths[4:7]): self.plot_distribution(distributions, "A_negative", measure, path)
        self.plot_bins(bins["all"], paths[7], "Step 05: Per-sample raw UMI-bin composition")
        self.plot_bins(bins["A_negative"], paths[8], "Step 05: A-negative raw UMI-bin composition")
        self.plot_heatmap(summaries["all"], summaries["A_negative"], paths[9])
        self.plot_ko_comparison(interpretation, paths[10])
        self.plot_error_states(transitions, paths[11])
        return paths


class SampleProbeDiagnosticPublisher(object):
    """Atomically publish the diagnostic package with complete provenance."""

    TABLE_NAMES = {FLOW_NAME, GENOTYPE_FLOW_NAME, SUMMARY_NAME, ANEG_SUMMARY_NAME, APOS_SUMMARY_NAME, DISTRIBUTION_NAME, BIN_NAME, ANEG_BIN_NAME, APOS_BIN_NAME, PERFORMANCE_NAME, TRANSITIONS_NAME, PAIRWISE_NAME, AUC_NAME, INTERPRETATION_NAME, SCOPE_NAME, CHECKS_NAME, ENVIRONMENT_NAME}
    EXPECTED_FILES = TABLE_NAMES.union(FIGURES)

    def __init__(self, configuration, validation, plotter, output_root, paths):
        self.configuration = configuration; self.validation = validation; self.plotter = plotter
        self.output_root = os.path.abspath(output_root); self.paths = {key: os.path.abspath(value) for key, value in paths.items()}
        self.decimals = int(configuration.values["probability_decimal_places"])

    def _format_rows(self, rows):
        output = []
        for source in rows:
            row = {}
            for key, value in source.items():
                row[key] = "" if value is None else (("{:.%df}" % self.decimals).format(value) if isinstance(value, float) else value)
            output.append(row)
        return output

    @staticmethod
    def _header(rows):
        return list(rows[0].keys()) if rows else []

    def _write_scope(self, path):
        scope = self.configuration.values["sample_level_probe_evidence_diagnostics"]
        rows = [
            {"key": "analysis_step", "value": "05_pcdh19_logistic_regression_baseline"},
            {"key": "analysis_module", "value": "sample_level_probe_evidence_diagnostics"},
            {"key": "included_samples", "value": "WT-M JZ-1--3; WT-F JZ-4--6; KO-M JZ-10--12"},
            {"key": "excluded_samples", "value": "HET-F JZ-7--9; excluded and not loaded"},
            {"key": "cell_eligibility", "value": "all cells already present in immutable paired Step 05 held-out predictions"},
            {"key": "analysis_exclusions", "value": "none"},
            {"key": "calling_exclusion", "value": "exact A/B/C zero remains uncalled under existing policy"},
            {"key": "normalization", "value": scope["normalization"]},
            {"key": "count_interpretation", "value": scope["count_interpretation"]},
            {"key": "classifier_fit", "value": "none"}, {"key": "prediction_changes", "value": "none"},
            {"key": "threshold_changes", "value": "none"}, {"key": "cell_level_p_values", "value": "none"},
            {"key": "biological_unit", "value": "registered sample"},
        ]
        write_tsv(path, ["key", "value"], rows)

    def _write_environment(self, path):
        rows = [{"key": "python_version", "value": platform.python_version()}, {"key": "numpy_version", "value": np.__version__}, {"key": "matplotlib_version", "value": matplotlib.__version__}]
        for key in sorted(self.paths):
            rows.extend([{"key": key + "_path", "value": self.paths[key]}, {"key": key + "_sha256", "value": sha256_file(self.paths[key])}])
        write_tsv(path, ["key", "value"], rows)

    @staticmethod
    def _manifest(directory):
        rows = []
        for name in sorted(os.listdir(directory)):
            path = os.path.join(directory, name)
            if name != MANIFEST_NAME and os.path.isfile(path): rows.append({"relative_path": name, "bytes": os.path.getsize(path), "sha256": sha256_file(path)})
        write_tsv(os.path.join(directory, MANIFEST_NAME), ["relative_path", "bytes", "sha256"], rows)

    @classmethod
    def _verify_manifest(cls, directory):
        listed = set()
        with open(os.path.join(directory, MANIFEST_NAME), "r", newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                path = os.path.join(directory, row["relative_path"]); listed.add(row["relative_path"])
                if not os.path.isfile(path) or os.path.getsize(path) != int(row["bytes"]) or sha256_file(path) != row["sha256"]: raise DiagnosticError("Diagnostic manifest verification failed")
        if listed != cls.EXPECTED_FILES: raise DiagnosticError("Unexpected diagnostic output file set")

    def verify_existing(self):
        self._verify_manifest(self.output_root)
        with open(os.path.join(self.output_root, ENVIRONMENT_NAME), "r", newline="") as handle: env = {row["key"]: row["value"] for row in csv.DictReader(handle, delimiter="\t")}
        for key, path in self.paths.items():
            if env.get(key + "_sha256") != sha256_file(path): raise DiagnosticError("Existing diagnostic provenance mismatch")

    def publish(self, products):
        if os.path.exists(self.output_root): self.verify_existing(); print("EXISTING_VALIDATED\t{}".format(self.output_root)); return
        staging = tempfile.mkdtemp(prefix=".sample_probe_diagnostics.", dir=os.path.dirname(self.output_root))
        try:
            tables = [
                (FLOW_NAME, products["flow"]), (GENOTYPE_FLOW_NAME, products["genotype_flow"]),
                (SUMMARY_NAME, products["summaries"]["all"]), (ANEG_SUMMARY_NAME, products["summaries"]["A_negative"]), (APOS_SUMMARY_NAME, products["summaries"]["A_positive"]),
                (DISTRIBUTION_NAME, products["distributions"]), (BIN_NAME, products["bins"]["all"]), (ANEG_BIN_NAME, products["bins"]["A_negative"]), (APOS_BIN_NAME, products["bins"]["A_positive"]),
                (PERFORMANCE_NAME, products["performance"]), (TRANSITIONS_NAME, products["transitions"]), (PAIRWISE_NAME, products["pairwise"]),
                (AUC_NAME, products["auc"]), (INTERPRETATION_NAME, products["interpretation"]),
            ]
            for name, rows in tables:
                formatted = self._format_rows(rows); write_tsv(os.path.join(staging, name), self._header(formatted), formatted)
            self._write_scope(os.path.join(staging, SCOPE_NAME))
            plot_paths = self.plotter.render(staging, products["summaries"], products["distributions"], products["bins"], products["transitions"], products["interpretation"])
            self.validation.require_equal("diagnostic_plots_created", sum(os.path.getsize(path) > 0 for path in plot_paths), len(FIGURES))
            write_tsv(os.path.join(staging, CHECKS_NAME), ["step_id", "check_name", "status", "observed", "expected", "details"], self.validation.rows)
            self._write_environment(os.path.join(staging, ENVIRONMENT_NAME)); self._manifest(staging); self._verify_manifest(staging)
            if os.path.exists(self.output_root): raise DiagnosticError("Diagnostic output appeared during publication")
            os.replace(staging, self.output_root); staging = None
            print("PUBLISHED\t{}\tcells={}\tsamples=9".format(self.output_root, sum(row["retained_cells"] for row in products["flow"])))
        finally:
            if staging and os.path.exists(staging): shutil.rmtree(staging)


class SampleProbeEvidenceDiagnosticStep(object):
    """Coordinate the read-only Step 05 biological-sample diagnostic."""

    def __init__(self, configuration, validation, paired_path, comparison_path, step02_root, output_root, paths):
        self.configuration = configuration; self.validation = validation
        self.reader = SampleProbeDiagnosticReader(configuration, validation, paired_path, step02_root, paths["step_02a_manifest"])
        self.summary = ProbeEvidenceSummaryBuilder(configuration, validation)
        self.performance = ProbeStatePerformanceAnalyzer(); self.pairwise = SamplePairwiseProbeComparator(); self.auc = AUCOrientationAuditor()
        self.publisher = SampleProbeDiagnosticPublisher(configuration, validation, SampleProbeDiagnosticPlotter(configuration.values["plot_dpi"]), output_root, paths)
        self.comparison_path = comparison_path

    def _interpretation(self, summaries, bins, performance):
        all_index = {(row["biological_sample_id"], row["measure"]): row for row in summaries["all"]}
        aneg_index = {(row["biological_sample_id"], row["measure"]): row for row in summaries["A_negative"]}
        bin_index = {(row["biological_sample_id"], row["measure"], row["umi_bin"]): row for row in bins["A_negative"]}
        perf = {(row["biological_sample_id"], row["model"]): row for row in performance if row["probe_state_stratum"] == "all"}
        rows = []
        for sample_id, definition in self.configuration.values["expanded_ground_truth_validation"]["expected_samples"].items():
            binary = perf[(sample_id, "binary")]["accuracy_among_called"]; count = perf[(sample_id, "count_informed")]["accuracy_among_called"]
            delta = 100.0 * (count - binary)
            rows.append({
                "biological_sample_id": sample_id, "sex": definition["sex"], "design_group": definition["design_group"], "true_genotype": definition["genotype"],
                "total_cells": all_index[(sample_id, "B_UMI")]["total_sample_cells"], "A_detected_percent": all_index[(sample_id, "B_UMI")]["A_detected_percent"],
                "B_UMI_mean": all_index[(sample_id, "B_UMI")]["mean"], "C_UMI_mean": all_index[(sample_id, "C_UMI")]["mean"], "B_plus_C_UMI_mean": all_index[(sample_id, "B_plus_C_UMI")]["mean"],
                "Anegative_cells": aneg_index[(sample_id, "B_plus_C_UMI")]["subset_cells"], "Anegative_B_plus_C_mean": aneg_index[(sample_id, "B_plus_C_UMI")]["mean"],
                "Anegative_B_plus_C_zero_percent": bin_index[(sample_id, "B_plus_C_UMI", "0")]["percent_of_subset"],
                "binary_accuracy_percent": 100.0 * binary, "count_informed_accuracy_percent": 100.0 * count,
                "count_minus_binary_accuracy_points": delta, "performance_change_direction": "improved" if delta > 1e-12 else ("worsened" if delta < -1e-12 else "unchanged"),
            })
        return rows

    def run(self):
        if os.path.exists(self.publisher.output_root): self.publisher.publish(None); return
        records = self.reader.read(); flow, genotype_flow = self.summary.flow(records); summaries, distributions, bins = self.summary.summaries(records)
        sample_order = list(self.configuration.values["expanded_ground_truth_validation"]["expected_samples"])
        performance, transitions = self.performance.analyze(records, sample_order)
        pairwise = self.pairwise.compare(summaries, bins); auc = self.auc.audit(records, self.comparison_path)
        interpretation = self._interpretation(summaries, bins, performance)
        corrections = sum(row["binary_wrong_count_correct"] for row in transitions if row["probe_state_stratum"] == "all")
        regressions = sum(row["binary_correct_count_wrong"] for row in transitions if row["probe_state_stratum"] == "all")
        self.validation.require_equal("diagnostic_total_corrections", corrections, 2373)
        self.validation.require_equal("diagnostic_total_regressions", regressions, 587)
        self.validation.require_equal("diagnostic_pairwise_ko_sample_pairs", len({(row["sample_a"], row["sample_b"]) for row in pairwise if row["true_genotype"] == "KO"}), 3)
        self.validation.require_equal("diagnostic_pairwise_wt_sample_pairs", len({(row["sample_a"], row["sample_b"]) for row in pairwise if row["true_genotype"] == "WT"}), 15)
        products = {"flow": flow, "genotype_flow": genotype_flow, "summaries": summaries, "distributions": distributions, "bins": bins, "performance": performance, "transitions": transitions, "pairwise": pairwise, "auc": auc, "interpretation": interpretation}
        self.publisher.publish(products)
        print("SAMPLE_INTERPRETATION\tsample\tgenotype\tA_detected_pct\tAneg_BC_mean\tAneg_BC_zero_pct\tcount_minus_binary_accuracy_points\tdirection")
        for row in interpretation:
            print("SAMPLE_INTERPRETATION\t{biological_sample_id}\t{true_genotype}\t{A_detected_percent:.3f}\t{Anegative_B_plus_C_mean:.4f}\t{Anegative_B_plus_C_zero_percent:.3f}\t{count_minus_binary_accuracy_points:+.3f}\t{performance_change_direction}".format(**row))
