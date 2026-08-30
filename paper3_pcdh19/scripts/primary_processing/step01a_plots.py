"""Visualize every Step 01a per-sample MAD boundary and candidate summary."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

from .step01a_models import Step01aSettings
from .step01a_outliers import Step01aResult


class Step01aPlotPublisher:
    """Publish per-sample boundary panels and cross-sample flag summaries."""

    COLORS = {3.0: "#d73027", 4.0: "#fc8d59", 5.0: "#4575b4"}
    PATTERN_COLORS = {
        "counts_only": "#4c78a8",
        "genes_only": "#59a14f",
        "mt_only": "#e15759",
        "counts+genes": "#76b7b2",
        "counts+mt": "#b07aa1",
        "genes+mt": "#f28e2b",
        "all_three": "#7f3c8d",
    }

    def __init__(self, settings: Step01aSettings, figures_dir: Path):
        """Store rendering settings and create the run-scoped figure directory."""

        self.settings = settings
        self.figures_dir = figures_dir
        self.per_sample_dir = figures_dir / "per_sample"
        self.per_sample_dir.mkdir(parents=True, exist_ok=True)

    def publish_all(self, obs: pd.DataFrame, result: Step01aResult) -> pd.DataFrame:
        """Render all requested visual sensitivity diagnostics in PNG and PDF."""

        records: list[dict[str, str]] = []
        for sample, sample_obs in obs.groupby(self.settings.sample_field, sort=False, observed=True):
            sample_boundaries = result.boundaries.loc[result.boundaries["sample"] == str(sample)]
            figure = self._sample_boundary_figure(str(sample), sample_obs, sample_boundaries)
            records.extend(self._save(figure, self.per_sample_dir / f"{sample}_mad_boundaries", str(sample), "Per-sample distributions with 3/4/5 scaled-MAD candidate boundaries."))
        figure = self._flag_percentage_heatmap(result.flag_summary)
        records.extend(self._save(figure, self.figures_dir / "candidate_flag_percentages_by_sample", "sample_summary", "Individual and joint candidate percentages for every sample and MAD stringency."))
        figure = self._overlap_figure(result.overlap_summary)
        records.extend(self._save(figure, self.figures_dir / "candidate_flag_overlap_by_sample", "sample_overlap", "Exact overlap composition of candidate flags by sample and MAD stringency."))
        figure = self._boundary_comparison(result.flag_summary)
        records.extend(self._save(figure, self.figures_dir / "candidate_boundaries_by_sample", "boundary_summary", "Original-scale candidate boundaries across samples and stringencies."))
        return pd.DataFrame(records)

    def _sample_boundary_figure(self, sample: str, obs: pd.DataFrame, boundaries: pd.DataFrame) -> plt.Figure:
        """Overlay all candidate boundaries on one sample's existing distributions."""

        figure, axes = plt.subplots(1, 3, figsize=(17, 5.2))
        panels = (
            ("total_counts", np.log10(obs["total_counts"].to_numpy(dtype=float) + 1), "log10(total counts + 1)", "#4c78a8"),
            ("n_genes_by_counts", obs["n_genes_by_counts"].to_numpy(dtype=float), "Detected genes", "#59a14f"),
            ("pct_counts_mt", obs["pct_counts_mt"].to_numpy(dtype=float), "Mitochondrial counts (%)", "#e15759"),
        )
        for axis, (metric, values, label, color) in zip(axes, panels, strict=True):
            axis.hist(values, bins=self.settings.histogram_bins, color=color, alpha=0.78)
            metric_rows = boundaries.loc[boundaries["metric"] == metric].sort_values("stringency_mad")
            for row in metric_rows.itertuples(index=False):
                boundary = float(row.boundary_original_scale)
                plotted_boundary = np.log10(boundary + 1) if metric == "total_counts" else boundary
                axis.axvline(plotted_boundary, color=self.COLORS[float(row.stringency_mad)], linewidth=2.0, linestyle="--")
            axis.set_xlabel(label)
            axis.set_ylabel("Cells")
            axis.set_title(metric.replace("_", " "))
        legend = [Line2D([0], [0], color=self.COLORS[value], linestyle="--", linewidth=2, label=f"{value:g} MAD") for value in self.settings.stringencies]
        figure.legend(handles=legend, loc="upper center", ncol=len(legend), frameon=False, bbox_to_anchor=(0.5, 0.95))
        figure.suptitle(f"{sample}: per-sample technical-outlier sensitivity only; no filtering", fontsize=15)
        figure.tight_layout(rect=(0, 0, 1, 0.88))
        return figure

    def _flag_percentage_heatmap(self, summary: pd.DataFrame) -> plt.Figure:
        """Display individual and joint candidate percentages for all samples."""

        metrics = [
            ("low_total_counts_pct", "Low counts"),
            ("low_n_genes_pct", "Low genes"),
            ("high_pct_counts_mt_pct", "High mt %"),
            ("any_candidate_pct", "Any criterion"),
            ("all_three_candidate_pct", "All three"),
        ]
        samples = summary["technical_sample_id"].drop_duplicates().tolist()
        maximum = max(float(summary[column].max()) for column, _ in metrics)
        figure, axes = plt.subplots(1, len(self.settings.stringencies), figsize=(18, 7), sharey=True)
        for axis, stringency in zip(axes, self.settings.stringencies, strict=True):
            selected = summary.loc[summary["stringency_mad"] == stringency].set_index("technical_sample_id").loc[samples]
            matrix = selected[[column for column, _ in metrics]].to_numpy(dtype=float)
            image = axis.imshow(matrix, aspect="auto", cmap="magma", vmin=0, vmax=maximum)
            axis.set_title(f"{stringency:g} MAD")
            axis.set_xticks(range(len(metrics)), [label for _, label in metrics], rotation=35, ha="right")
            axis.set_yticks(range(len(samples)), samples)
            for row in range(matrix.shape[0]):
                for column in range(matrix.shape[1]):
                    axis.text(column, row, f"{matrix[row, column]:.2f}", ha="center", va="center", fontsize=7, color="black" if matrix[row, column] > maximum * 0.55 else "white")
        figure.colorbar(image, ax=axes, label="Cells flagged (%)", fraction=0.02, pad=0.02)
        figure.suptitle("Per-sample MAD candidate flags — descriptive sensitivity only; no filtering", fontsize=15)
        figure.subplots_adjust(left=0.10, right=0.94, bottom=0.18, top=0.88, wspace=0.28)
        return figure

    def _overlap_figure(self, overlap: pd.DataFrame) -> plt.Figure:
        """Plot exact non-null flag-pattern percentages as stacked bars."""

        patterns = list(self.PATTERN_COLORS)
        samples = overlap["technical_sample_id"].drop_duplicates().tolist()
        figure, axes = plt.subplots(len(self.settings.stringencies), 1, figsize=(17, 12), sharex=True)
        for axis, stringency in zip(axes, self.settings.stringencies, strict=True):
            selected = overlap.loc[(overlap["stringency_mad"] == stringency) & (overlap["candidate_pattern"] != "none")]
            pivot = selected.pivot(index="technical_sample_id", columns="candidate_pattern", values="pct_cells").reindex(samples).fillna(0)
            bottom = np.zeros(len(samples))
            for pattern in patterns:
                values = pivot[pattern].to_numpy(dtype=float) if pattern in pivot else np.zeros(len(samples))
                axis.bar(samples, values, bottom=bottom, label=pattern, color=self.PATTERN_COLORS[pattern], width=0.78)
                bottom += values
            axis.set_ylabel("Cells (%)")
            axis.set_title(f"{stringency:g} MAD: exact candidate overlap")
            axis.grid(axis="y", alpha=0.2)
        axes[-1].tick_params(axis="x", rotation=35)
        handles, labels = axes[0].get_legend_handles_labels()
        figure.legend(handles, labels, loc="upper center", ncol=7, frameon=False)
        figure.suptitle("Joint candidate-flag composition by technical sample; no cells removed", fontsize=15)
        figure.tight_layout(rect=(0, 0, 1, 0.93))
        return figure

    def _boundary_comparison(self, summary: pd.DataFrame) -> plt.Figure:
        """Compare original-scale boundary values across samples and stringencies."""

        panels = (
            ("low_total_counts_boundary", "Low total-count boundary"),
            ("low_n_genes_boundary", "Low detected-gene boundary"),
            ("high_pct_counts_mt_boundary", "High mitochondrial boundary (%)"),
        )
        samples = summary["technical_sample_id"].drop_duplicates().tolist()
        x = np.arange(len(samples))
        figure, axes = plt.subplots(1, 3, figsize=(19, 5.5))
        for axis, (field, title) in zip(axes, panels, strict=True):
            for stringency in self.settings.stringencies:
                selected = summary.loc[summary["stringency_mad"] == stringency].set_index("technical_sample_id").loc[samples]
                axis.plot(x, selected[field], marker="o", color=self.COLORS[stringency], label=f"{stringency:g} MAD")
            axis.set_xticks(x, samples, rotation=40, ha="right")
            axis.set_title(title)
            axis.grid(axis="y", alpha=0.25)
        axes[0].set_ylabel("Boundary on original metric scale")
        handles, labels = axes[0].get_legend_handles_labels()
        figure.legend(handles, labels, loc="upper center", ncol=3, frameon=False)
        figure.suptitle("Per-sample robust candidate boundaries; no upper count/gene boundary", fontsize=15)
        figure.tight_layout(rect=(0, 0, 1, 0.90))
        return figure

    def _save(self, figure: plt.Figure, stem: Path, scope: str, description: str) -> list[dict[str, str]]:
        """Save one figure as review PNG and vector PDF and return manifest rows."""

        records = []
        for extension in ("png", "pdf"):
            path = stem.with_suffix(f".{extension}")
            figure.savefig(path, dpi=self.settings.plot_dpi if extension == "png" else None, bbox_inches="tight")
            records.append({"relative_path": str(path.relative_to(self.figures_dir.parent)), "scope": scope, "format": extension, "description": description})
        plt.close(figure)
        return records
