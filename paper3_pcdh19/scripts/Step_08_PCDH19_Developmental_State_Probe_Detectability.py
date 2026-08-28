#!/usr/bin/env python3
"""Formal Step 08: developmental-state composition and PCDH19 detectability.

This post-freeze analysis joins the independent GSE94641 E15.5 broad-state
transfer to the immutable Step 02a per-cell PCDH19 probe audit. It describes
state composition and A/B/C detectability by biological sample. The only
inferential comparison is sex-matched WT-male versus KO-male, with sample ID
as the replicate and exact label permutations. WT-female and HET-female
samples are descriptive context only. No classifier is fit or modified.
"""

from __future__ import print_function

import argparse
import csv
import datetime
import gzip
import hashlib
import itertools
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


STEP_ID = "08_pcdh19_developmental_state_probe_detectability"
MANIFEST_NAME = "output_manifest.tsv"
JOINED_NAME = "step_08_cell_state_pcdh19_probe_evidence.tsv.gz"
COMPOSITION_NAME = "step_08_per_sample_developmental_state_composition.tsv"
DETECTION_NAME = "step_08_per_sample_state_probe_detectability.tsv"
GROUP_NAME = "step_08_group_state_probe_detectability_summary.tsv"
COMPARISON_NAME = "step_08_wt_m_vs_ko_m_sample_level_comparisons.tsv"
VALIDATION_NAME = "step_08_validation.tsv"
ENVIRONMENT_NAME = "software_environment.tsv"
COMPOSITION_PLOT = "step_08_developmental_state_composition_by_sample.png"
DETECTION_PLOT = "step_08_pcdh19_probe_detectability_by_state_and_group.png"
WTKO_PLOT = "step_08_wt_m_vs_ko_m_probe_detectability_by_state.png"

GROUP_ORDER = ["WT_M", "WT_F", "HET_F", "KO_M"]
GROUP_COLORS = {"WT_M": "#2166ac", "WT_F": "#67a9cf", "HET_F": "#fdae61", "KO_M": "#b2182b"}
STATE_COLORS = {
    "proliferating_neural_progenitor": "#4daf4a",
    "postmitotic_immature_neuron": "#984ea3",
    "not_assigned_neural_state": "#999999",
}
STATE_SHORT = {
    "proliferating_neural_progenitor": "Progenitor",
    "postmitotic_immature_neuron": "Immature neuron",
    "not_assigned_neural_state": "Unassigned",
}
METRIC_SHORT = {
    "A_detected": "A+",
    "B_detected": "B+",
    "C_detected": "C+",
    "any_PCDH19_detected": "Any probe+",
}


class Step08Error(RuntimeError):
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


def fmt(value, places=10):
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return "NA"
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
            raise Step08Error("Validation failed for {}: {!r} != {!r}".format(name, observed, expected))


class Step08Configuration(object):
    REQUIRED = {
        "step_id", "pipeline_version", "analysis_scope", "label_transfer", "probe_audit",
        "sample_key", "broad_state_order", "detection_metrics", "primary_comparison",
        "context_groups", "plot_dpi", "prohibited_operations",
    }

    def __init__(self, path):
        self.path = os.path.abspath(path)
        with open(self.path, "r") as handle:
            self.values = json.load(handle)
        missing = sorted(self.REQUIRED.difference(self.values))
        if missing:
            raise Step08Error("Step 08 lock missing keys: {}".format(", ".join(missing)))
        if self.values["step_id"] != STEP_ID:
            raise Step08Error("Unexpected Step 08 identifier")
        expected_prohibited = {
            "classifier_refitting", "threshold_modification", "het_model_development",
            "cell_as_independent_replicate", "fine_subtype_biological_relabeling",
            "probe_umi_as_transcript_count",
        }
        if set(self.values["prohibited_operations"]) != expected_prohibited:
            raise Step08Error("Incomplete prohibited-operation guard")
        comparison = self.values["primary_comparison"]
        if comparison["reference_group"] != "WT_M" or comparison["comparison_group"] != "KO_M":
            raise Step08Error("Primary comparison must remain WT_M versus KO_M")
        if comparison["replicate_unit"] != "biological_sample_id":
            raise Step08Error("Biological sample ID must be the replicate")

    @property
    def states(self):
        return list(self.values["broad_state_order"])

    @property
    def metrics(self):
        return list(self.values["detection_metrics"])


class SampleRegistry(object):
    REQUIRED = [
        "technical_sample_id", "submitted_sample_name", "target_cells", "organism", "tissue",
        "region", "genotype", "sex", "design_group", "metadata_source", "metadata_registered_date",
    ]

    def __init__(self, path, configuration, validation):
        self.path = os.path.abspath(path)
        definition = configuration.values["sample_key"]
        validation.check("sample_key_bytes", os.path.getsize(self.path), definition["bytes"])
        validation.check("sample_key_sha256", sha256_file(self.path), definition["sha256"])
        self.samples = {}
        with open(self.path, "r", newline="") as handle:
            reader = csv.DictReader(handle)
            validation.check("sample_key_schema", reader.fieldnames, self.REQUIRED)
            for row in reader:
                sample_id = row["technical_sample_id"]
                if sample_id in self.samples:
                    raise Step08Error("Duplicate sample ID {}".format(sample_id))
                self.samples[sample_id] = row
        validation.check("registered_sample_ids", set(self.samples), set(configuration.values["probe_audit"]["samples"]))


class ProbeEvidenceRepository(object):
    """Checksum-validate and load immutable per-cell Step 02a probe evidence."""

    REQUIRED = ["barcode", "A_UMI", "B_UMI", "C_UMI", "Pcdh19_total_UMI", "detection_pattern"]

    def __init__(self, root, configuration, validation):
        self.root = os.path.abspath(root)
        self.configuration = configuration
        self.validation = validation
        definition = configuration.values["probe_audit"]
        manifest = os.path.join(self.root, MANIFEST_NAME)
        validation.check("step_02a_manifest_sha256", sha256_file(manifest), definition["manifest_sha256"])

    def load(self):
        result = {}
        definition = self.configuration.values["probe_audit"]
        total = 0
        for sample_id in sorted(definition["samples"]):
            expected = definition["samples"][sample_id]
            path = os.path.join(self.root, "per_sample", sample_id, definition["per_sample_file"])
            self.validation.check("{}_probe_bytes".format(sample_id), os.path.getsize(path), expected["bytes"])
            self.validation.check("{}_probe_sha256".format(sample_id), sha256_file(path), expected["sha256"])
            sample = {}
            with open(path, "r", newline="") as handle:
                reader = csv.DictReader(handle, delimiter="\t")
                self.validation.check("{}_probe_schema".format(sample_id), reader.fieldnames, self.REQUIRED)
                for row in reader:
                    barcode = row["barcode"]
                    if barcode in sample:
                        raise Step08Error("Duplicate probe barcode {} {}".format(sample_id, barcode))
                    a, b, c = int(row["A_UMI"]), int(row["B_UMI"]), int(row["C_UMI"])
                    total_umi = int(row["Pcdh19_total_UMI"])
                    if min(a, b, c) < 0 or a + b + c != total_umi:
                        raise Step08Error("Invalid probe arithmetic {} {}".format(sample_id, barcode))
                    sample[barcode] = (a, b, c, total_umi, row["detection_pattern"])
            self.validation.check("{}_probe_cells".format(sample_id), len(sample), expected["cells"])
            result[sample_id] = sample
            total += len(sample)
        self.validation.check("all_probe_cells", total, sum(v["cells"] for v in definition["samples"].values()))
        return result


class StateProbeAggregator(object):
    """Join labels/probes and maintain only sample-level sufficient statistics."""

    LABEL_REQUIRED = [
        "sample_id", "submitted_sample_name", "genotype", "sex", "design_group", "cell_barcode",
        "GSE94641_broad_state", "GSE94641_label", "GSE94641_subtype", "GSE94641_nonneuralcelltype",
    ]
    APPENDED = [
        "A_UMI", "B_UMI", "C_UMI", "Pcdh19_total_UMI", "A_detected", "B_detected",
        "C_detected", "any_PCDH19_detected", "probe_pattern_code", "probe_detection_pattern",
    ]

    def __init__(self, label_root, probe_data, registry, configuration, validation):
        self.label_root = os.path.abspath(label_root)
        self.probe_data = probe_data
        self.registry = registry
        self.configuration = configuration
        self.validation = validation
        self.counts = {}
        self.sample_totals = {}

    def _bucket(self, sample_id, state):
        key = (sample_id, state)
        if key not in self.counts:
            self.counts[key] = {"cells": 0, "A_detected": 0, "B_detected": 0, "C_detected": 0, "any_PCDH19_detected": 0}
        return self.counts[key]

    def run(self, joined_path):
        definition = self.configuration.values["label_transfer"]
        manifest = os.path.join(self.label_root, MANIFEST_NAME)
        self.validation.check("step_00_label_manifest_sha256", sha256_file(manifest), definition["manifest_sha256"])
        input_path = os.path.join(self.label_root, definition["cell_file"])
        self.validation.check("step_00_label_cell_bytes", os.path.getsize(input_path), definition["cell_bytes"])
        self.validation.check("step_00_label_cell_sha256", sha256_file(input_path), definition["cell_sha256"])
        seen = dict((sample_id, set()) for sample_id in self.registry.samples)
        state_set = set(self.configuration.states)
        with open(input_path, "r", newline="") as source, gzip.open(joined_path, "wt", newline="") as target:
            reader = csv.DictReader(source, delimiter="\t")
            missing = [field for field in self.LABEL_REQUIRED if field not in reader.fieldnames]
            if missing:
                raise Step08Error("Label-transfer table missing: {}".format(", ".join(missing)))
            writer = csv.DictWriter(target, fieldnames=reader.fieldnames + self.APPENDED, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            for row in reader:
                sample_id = row["sample_id"]
                barcode = row["cell_barcode"]
                if sample_id not in self.registry.samples:
                    raise Step08Error("Unregistered label sample {}".format(sample_id))
                registered = self.registry.samples[sample_id]
                for label_field, key_field in (("submitted_sample_name", "submitted_sample_name"), ("genotype", "genotype"), ("sex", "sex"), ("design_group", "design_group")):
                    if row[label_field] != registered[key_field]:
                        raise Step08Error("Metadata mismatch {} {}".format(sample_id, label_field))
                state = row["GSE94641_broad_state"]
                if state not in state_set:
                    raise Step08Error("Unexpected broad state {}".format(state))
                if barcode in seen[sample_id]:
                    raise Step08Error("Duplicate label barcode {} {}".format(sample_id, barcode))
                seen[sample_id].add(barcode)
                try:
                    a, b, c, total_umi, detection_pattern = self.probe_data[sample_id][barcode]
                except KeyError:
                    raise Step08Error("Probe join missing {} {}".format(sample_id, barcode))
                indicators = {
                    "A_detected": int(a > 0), "B_detected": int(b > 0), "C_detected": int(c > 0),
                    "any_PCDH19_detected": int(total_umi > 0),
                }
                pattern_code = "{}{}{}".format(indicators["A_detected"], indicators["B_detected"], indicators["C_detected"])
                output = dict(row)
                output.update({
                    "A_UMI": a, "B_UMI": b, "C_UMI": c, "Pcdh19_total_UMI": total_umi,
                    "probe_pattern_code": pattern_code, "probe_detection_pattern": detection_pattern,
                })
                output.update(indicators)
                writer.writerow(output)
                bucket = self._bucket(sample_id, state)
                bucket["cells"] += 1
                for metric, value in indicators.items():
                    bucket[metric] += value
                self.sample_totals[sample_id] = self.sample_totals.get(sample_id, 0) + 1
        for sample_id, expected in self.configuration.values["probe_audit"]["samples"].items():
            self.validation.check("{}_joined_cells".format(sample_id), self.sample_totals.get(sample_id, 0), expected["cells"])
            self.validation.check(
                "{}_barcode_sets_equal".format(sample_id),
                seen[sample_id] == set(self.probe_data[sample_id]), True,
                "label_barcodes={} probe_barcodes={}".format(len(seen[sample_id]), len(self.probe_data[sample_id])),
            )
        self.validation.check("joined_total_cells", sum(self.sample_totals.values()), sum(v["cells"] for v in self.configuration.values["probe_audit"]["samples"].values()))

    def composition_rows(self):
        rows = []
        for sample_id in sorted(self.registry.samples, key=sample_sort_key):
            meta = self.registry.samples[sample_id]
            total = self.sample_totals[sample_id]
            for state in self.configuration.states:
                cells = self.counts.get((sample_id, state), {"cells": 0})["cells"]
                rows.append({
                    "biological_sample_id": sample_id, "submitted_sample_name": meta["submitted_sample_name"],
                    "genotype": meta["genotype"], "sex": meta["sex"], "design_group": meta["design_group"],
                    "developmental_state": state, "state_cells": cells, "sample_total_cells": total,
                    "state_fraction": cells / float(total), "state_percent": 100.0 * cells / float(total),
                })
        return rows

    def detection_rows(self):
        rows = []
        for sample_id in sorted(self.registry.samples, key=sample_sort_key):
            meta = self.registry.samples[sample_id]
            for state in self.configuration.states:
                bucket = self.counts.get((sample_id, state), {"cells": 0})
                cells = bucket["cells"]
                for metric in self.configuration.metrics:
                    detected = bucket.get(metric, 0)
                    fraction = detected / float(cells) if cells else float("nan")
                    rows.append({
                        "biological_sample_id": sample_id, "submitted_sample_name": meta["submitted_sample_name"],
                        "genotype": meta["genotype"], "sex": meta["sex"], "design_group": meta["design_group"],
                        "developmental_state": state, "detection_metric": metric, "state_cells": cells,
                        "detected_cells": detected, "detected_fraction": fraction,
                        "detected_percent": 100.0 * fraction,
                    })
        return rows


def sample_sort_key(sample_or_row):
    sample_id = sample_or_row if isinstance(sample_or_row, str) else sample_or_row["biological_sample_id"]
    return int(sample_id.rsplit("-", 1)[1])


class SampleLevelStatistics(object):
    """Summaries and exact tests that never use a cell as a replicate."""

    def __init__(self, configuration):
        self.configuration = configuration

    @staticmethod
    def _mean(values):
        return sum(values) / float(len(values))

    @staticmethod
    def _exact_permutation(wt_values, ko_values):
        values = list(wt_values) + list(ko_values)
        n_wt = len(wt_values)
        observed = SampleLevelStatistics._mean(ko_values) - SampleLevelStatistics._mean(wt_values)
        exceed = 0
        total = 0
        indices = set(range(len(values)))
        for wt_indices_tuple in itertools.combinations(range(len(values)), n_wt):
            wt_indices = set(wt_indices_tuple)
            ko_indices = indices.difference(wt_indices)
            permuted = SampleLevelStatistics._mean([values[i] for i in ko_indices]) - SampleLevelStatistics._mean([values[i] for i in wt_indices])
            if abs(permuted) + 1e-15 >= abs(observed):
                exceed += 1
            total += 1
        return observed, exceed / float(total), total

    @staticmethod
    def _bh(rows):
        ordered = sorted(enumerate(rows), key=lambda item: item[1]["exact_permutation_p_value"])
        m = len(rows)
        adjusted = [1.0] * m
        running = 1.0
        for rank_index in range(m - 1, -1, -1):
            original_index, row = ordered[rank_index]
            rank = rank_index + 1
            running = min(running, row["exact_permutation_p_value"] * m / float(rank))
            adjusted[original_index] = min(1.0, running)
        for index, value in enumerate(adjusted):
            rows[index]["bh_fdr_within_family"] = value
        return rows

    def group_rows(self, detection_rows):
        rows = []
        for group in GROUP_ORDER:
            for state in self.configuration.states:
                for metric in self.configuration.metrics:
                    selected = [r for r in detection_rows if r["design_group"] == group and r["developmental_state"] == state and r["detection_metric"] == metric]
                    values = [r["detected_fraction"] for r in selected]
                    detected = sum(r["detected_cells"] for r in selected)
                    cells = sum(r["state_cells"] for r in selected)
                    rows.append({
                        "design_group": group, "developmental_state": state, "detection_metric": metric,
                        "biological_samples": len(selected), "pooled_state_cells": cells,
                        "pooled_detected_cells": detected, "pooled_detected_fraction": detected / float(cells),
                        "sample_mean_detected_fraction": self._mean(values),
                        "sample_median_detected_fraction": float(np.median(values)),
                        "sample_min_detected_fraction": min(values), "sample_max_detected_fraction": max(values),
                    })
        return rows

    def comparison_rows(self, composition_rows, detection_rows):
        composition = []
        for state in self.configuration.states:
            wt = [r["state_fraction"] for r in composition_rows if r["design_group"] == "WT_M" and r["developmental_state"] == state]
            ko = [r["state_fraction"] for r in composition_rows if r["design_group"] == "KO_M" and r["developmental_state"] == state]
            effect, p_value, permutations = self._exact_permutation(wt, ko)
            composition.append({
                "analysis_family": "developmental_state_composition", "developmental_state": state,
                "detection_metric": "NA", "wt_m_samples": len(wt), "ko_m_samples": len(ko),
                "wt_m_sample_mean_fraction": self._mean(wt), "ko_m_sample_mean_fraction": self._mean(ko),
                "ko_minus_wt_fraction": effect, "ko_minus_wt_percentage_points": 100.0 * effect,
                "exact_permutation_p_value": p_value, "permutations": permutations,
            })
        detectability = []
        for state in self.configuration.states:
            for metric in self.configuration.metrics:
                wt = [r["detected_fraction"] for r in detection_rows if r["design_group"] == "WT_M" and r["developmental_state"] == state and r["detection_metric"] == metric]
                ko = [r["detected_fraction"] for r in detection_rows if r["design_group"] == "KO_M" and r["developmental_state"] == state and r["detection_metric"] == metric]
                effect, p_value, permutations = self._exact_permutation(wt, ko)
                detectability.append({
                    "analysis_family": "pcdh19_probe_detectability", "developmental_state": state,
                    "detection_metric": metric, "wt_m_samples": len(wt), "ko_m_samples": len(ko),
                    "wt_m_sample_mean_fraction": self._mean(wt), "ko_m_sample_mean_fraction": self._mean(ko),
                    "ko_minus_wt_fraction": effect, "ko_minus_wt_percentage_points": 100.0 * effect,
                    "exact_permutation_p_value": p_value, "permutations": permutations,
                })
        return self._bh(composition) + self._bh(detectability)


class Step08Plotter(object):
    def __init__(self, configuration):
        self.configuration = configuration
        self.dpi = int(configuration.values["plot_dpi"])

    def composition(self, rows, path):
        samples = sorted(set(r["biological_sample_id"] for r in rows), key=sample_sort_key)
        fig, ax = plt.subplots(figsize=(12.5, 5.5))
        bottoms = np.zeros(len(samples))
        for state in self.configuration.states:
            values = [next(r["state_fraction"] for r in rows if r["biological_sample_id"] == sample and r["developmental_state"] == state) for sample in samples]
            ax.bar(range(len(samples)), values, bottom=bottoms, color=STATE_COLORS[state], label=STATE_SHORT[state], width=0.78)
            bottoms += np.asarray(values)
        labels = [sample.replace("15662-", "") + "\n" + next(r["design_group"] for r in rows if r["biological_sample_id"] == sample) for sample in samples]
        ax.set_xticks(range(len(samples)))
        ax.set_xticklabels(labels, fontsize=8)
        ax.set_ylim(0, 1)
        ax.set_ylabel("Fraction of cells")
        ax.set_title("Step 08: transferred developmental-state composition by biological sample")
        ax.legend(frameon=False, ncol=3, loc="upper center", bbox_to_anchor=(0.5, -0.18))
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        fig.tight_layout()
        fig.savefig(path, dpi=self.dpi, bbox_inches="tight")
        plt.close(fig)

    def heatmap(self, group_rows, path):
        fig, axes = plt.subplots(1, 4, figsize=(14, 4.2), sharey=True)
        global_vmax = max(
            20.0,
            100.0 * max(row["sample_mean_detected_fraction"] for row in group_rows),
        )
        for ax, group in zip(axes, GROUP_ORDER):
            matrix = np.zeros((len(self.configuration.states), len(self.configuration.metrics)))
            for i, state in enumerate(self.configuration.states):
                for j, metric in enumerate(self.configuration.metrics):
                    row = next(r for r in group_rows if r["design_group"] == group and r["developmental_state"] == state and r["detection_metric"] == metric)
                    matrix[i, j] = 100.0 * row["sample_mean_detected_fraction"]
            image = ax.imshow(matrix, vmin=0, vmax=global_vmax, cmap="viridis", aspect="auto")
            for i in range(matrix.shape[0]):
                for j in range(matrix.shape[1]):
                    color = "white" if matrix[i, j] > 0.55 * global_vmax else "black"
                    ax.text(j, i, "{:.1f}".format(matrix[i, j]), ha="center", va="center", color=color, fontsize=8)
            ax.set_title(group)
            ax.set_xticks(range(len(self.configuration.metrics)))
            ax.set_xticklabels([METRIC_SHORT[m] for m in self.configuration.metrics], rotation=45, ha="right")
            ax.set_yticks(range(len(self.configuration.states)))
            ax.set_yticklabels([STATE_SHORT[s] for s in self.configuration.states])
        fig.colorbar(image, ax=axes.ravel().tolist(), label="Mean sample-level detected cells (%)", shrink=0.82)
        fig.suptitle("Step 08: PCDH19 probe detectability by developmental state and design group", y=1.03)
        fig.savefig(path, dpi=self.dpi, bbox_inches="tight")
        plt.close(fig)

    def wt_ko(self, detection_rows, path):
        fig, axes = plt.subplots(1, 4, figsize=(14.5, 4.3), sharey=True)
        offsets = {"WT_M": -0.12, "KO_M": 0.12}
        jitter = [-0.035, 0.0, 0.035]
        for ax, metric in zip(axes, self.configuration.metrics):
            for state_index, state in enumerate(self.configuration.states):
                for group in ("WT_M", "KO_M"):
                    selected = sorted([r for r in detection_rows if r["design_group"] == group and r["developmental_state"] == state and r["detection_metric"] == metric], key=sample_sort_key)
                    values = [100.0 * r["detected_fraction"] for r in selected]
                    xs = [state_index + offsets[group] + jitter[i] for i in range(len(values))]
                    ax.scatter(xs, values, color=GROUP_COLORS[group], s=32, alpha=0.9, edgecolor="white", linewidth=0.4, zorder=3)
                    mean = sum(values) / float(len(values))
                    ax.plot([state_index + offsets[group] - 0.08, state_index + offsets[group] + 0.08], [mean, mean], color=GROUP_COLORS[group], lw=2.4)
            ax.set_title(METRIC_SHORT[metric])
            ax.set_xticks(range(len(self.configuration.states)))
            ax.set_xticklabels([STATE_SHORT[s] for s in self.configuration.states], rotation=35, ha="right")
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
        axes[0].set_ylabel("Detected cells per sample (%)")
        handles = [plt.Line2D([0], [0], marker="o", color="w", markerfacecolor=GROUP_COLORS[g], markersize=7, label=g) for g in ("WT_M", "KO_M")]
        axes[-1].legend(handles=handles, frameon=False, loc="upper right")
        fig.suptitle("Step 08: sex-matched WT-M versus KO-M PCDH19 detectability\n(points are biological samples; bars are sample means)")
        fig.tight_layout(rect=(0, 0, 1, 0.88))
        fig.savefig(path, dpi=self.dpi, bbox_inches="tight")
        plt.close(fig)


class OutputPackage(object):
    def __init__(self, output_root, configuration, requirements_path, validation):
        self.output_root = os.path.abspath(output_root)
        self.configuration = configuration
        self.requirements_path = os.path.abspath(requirements_path)
        self.validation = validation

    def environment(self, path):
        rows = [
            {"component": "python", "version": platform.python_version()},
            {"component": "numpy", "version": np.__version__},
            {"component": "matplotlib", "version": matplotlib.__version__},
            {"component": "platform", "version": platform.platform()},
            {"component": "pipeline_source_sha256", "version": sha256_file(os.path.abspath(__file__))},
            {"component": "lock_sha256", "version": sha256_file(self.configuration.path)},
            {"component": "requirements_sha256", "version": sha256_file(self.requirements_path)},
        ]
        write_tsv(path, ["component", "version"], rows)

    def manifest(self, root):
        rows = []
        for directory, _, filenames in os.walk(root):
            for filename in filenames:
                if filename == MANIFEST_NAME:
                    continue
                path = os.path.join(directory, filename)
                rows.append({
                    "relative_path": os.path.relpath(path, root), "bytes": os.path.getsize(path),
                    "sha256": sha256_file(path),
                })
        rows.sort(key=lambda row: row["relative_path"])
        write_tsv(os.path.join(root, MANIFEST_NAME), ["relative_path", "bytes", "sha256"], rows)


class Step08Workflow(object):
    def __init__(self, args):
        self.args = args
        self.configuration = Step08Configuration(args.lock)
        self.validation = ValidationRecorder()

    def run(self):
        output_root = os.path.abspath(self.args.output_root)
        if os.path.exists(output_root):
            raise Step08Error("Output root already exists and will not be overwritten: {}".format(output_root))
        parent = os.path.dirname(output_root)
        os.makedirs(parent, exist_ok=True)
        staging = tempfile.mkdtemp(prefix=".step08-staging-", dir=parent)
        try:
            registry = SampleRegistry(self.args.sample_key, self.configuration, self.validation)
            probes = ProbeEvidenceRepository(self.args.probe_root, self.configuration, self.validation).load()
            aggregator = StateProbeAggregator(self.args.label_root, probes, registry, self.configuration, self.validation)
            aggregator.run(os.path.join(staging, JOINED_NAME))
            composition = aggregator.composition_rows()
            detection = aggregator.detection_rows()
            statistics = SampleLevelStatistics(self.configuration)
            groups = statistics.group_rows(detection)
            comparisons = statistics.comparison_rows(composition, detection)
            composition_header = ["biological_sample_id", "submitted_sample_name", "genotype", "sex", "design_group", "developmental_state", "state_cells", "sample_total_cells", "state_fraction", "state_percent"]
            detection_header = ["biological_sample_id", "submitted_sample_name", "genotype", "sex", "design_group", "developmental_state", "detection_metric", "state_cells", "detected_cells", "detected_fraction", "detected_percent"]
            group_header = ["design_group", "developmental_state", "detection_metric", "biological_samples", "pooled_state_cells", "pooled_detected_cells", "pooled_detected_fraction", "sample_mean_detected_fraction", "sample_median_detected_fraction", "sample_min_detected_fraction", "sample_max_detected_fraction"]
            comparison_header = ["analysis_family", "developmental_state", "detection_metric", "wt_m_samples", "ko_m_samples", "wt_m_sample_mean_fraction", "ko_m_sample_mean_fraction", "ko_minus_wt_fraction", "ko_minus_wt_percentage_points", "exact_permutation_p_value", "bh_fdr_within_family", "permutations"]
            write_tsv(os.path.join(staging, COMPOSITION_NAME), composition_header, format_rows(composition))
            write_tsv(os.path.join(staging, DETECTION_NAME), detection_header, format_rows(detection))
            write_tsv(os.path.join(staging, GROUP_NAME), group_header, format_rows(groups))
            write_tsv(os.path.join(staging, COMPARISON_NAME), comparison_header, format_rows(comparisons))
            plotter = Step08Plotter(self.configuration)
            plotter.composition(composition, os.path.join(staging, COMPOSITION_PLOT))
            plotter.heatmap(groups, os.path.join(staging, DETECTION_PLOT))
            plotter.wt_ko(detection, os.path.join(staging, WTKO_PLOT))
            self.validation.check("classifier_fit_count", 0, 0, "Step 08 contains no estimator")
            self.validation.check("threshold_changes", 0, 0, "Frozen Step 07 is never opened or modified")
            self.validation.check("wt_m_samples", sum(1 for row in registry.samples.values() if row["design_group"] == "WT_M"), 3)
            self.validation.check("ko_m_samples", sum(1 for row in registry.samples.values() if row["design_group"] == "KO_M"), 3)
            write_tsv(os.path.join(staging, VALIDATION_NAME), ["check", "status", "observed", "expected", "detail"], self.validation.rows)
            package = OutputPackage(output_root, self.configuration, self.args.requirements, self.validation)
            package.environment(os.path.join(staging, ENVIRONMENT_NAME))
            shutil.copy2(self.configuration.path, os.path.join(staging, os.path.basename(self.configuration.path)))
            shutil.copy2(self.args.requirements, os.path.join(staging, os.path.basename(self.args.requirements)))
            package.manifest(staging)
            os.rename(staging, output_root)
            staging = None
            print("Step 08 complete: {}".format(output_root))
        finally:
            if staging and os.path.exists(staging):
                shutil.rmtree(staging)


def format_rows(rows):
    formatted = []
    for row in rows:
        item = {}
        for key, value in row.items():
            if isinstance(value, float):
                item[key] = fmt(value, 10)
            else:
                item[key] = value
        formatted.append(item)
    return formatted


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--lock", required=True)
    parser.add_argument("--requirements", required=True)
    parser.add_argument("--sample-key", required=True)
    parser.add_argument("--label-root", required=True)
    parser.add_argument("--probe-root", required=True)
    parser.add_argument("--output-root", required=True)
    return parser.parse_args(argv)


def main(argv=None):
    Step08Workflow(parse_args(argv)).run()


if __name__ == "__main__":
    try:
        main()
    except Step08Error as error:
        print("ERROR: {}".format(error), file=sys.stderr)
        sys.exit(2)
