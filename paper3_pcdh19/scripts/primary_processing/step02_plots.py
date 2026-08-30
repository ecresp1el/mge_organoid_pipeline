"""Review plots for approved Step 02 QC filtering."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .step02_models import Step02Settings


class Step02PlotPublisher:
    """Publish before/after and exclusion-reason review figures."""

    def __init__(self, settings: Step02Settings, figures_dir: Path):
        """Store rendering settings and the run-scoped figure directory."""

        self.settings = settings
        self.figures_dir = figures_dir

    def publish_all(self, by_sample: pd.DataFrame, by_design: pd.DataFrame) -> pd.DataFrame:
        """Render all Step 02 review figures in PNG and PDF formats."""

        records = []
        for figure, stem, scope, description in (
            (self._before_after(by_sample, "Technical sample"), "cells_before_after_by_sample", "sample", "Before/after cell counts and removed percentage by technical sample."),
            (self._before_after(by_design, "Design group"), "cells_before_after_by_design_group", "design_group", "Before/after cell counts and removed percentage by design group; filters were not defined by design group."),
            (self._reason_plot(by_sample), "excluded_cells_by_reason_and_sample", "sample_reason", "Selected 5-MAD exclusion criteria by sample; exact overlaps retained."),
        ):
            records.extend(self._save(figure, stem, scope, description))
        return pd.DataFrame(records)

    @staticmethod
    def _before_after(summary: pd.DataFrame, label: str) -> plt.Figure:
        """Plot retained cells and removed percentages for one grouping table."""

        labels = summary["group"].astype(str).tolist()
        x = np.arange(len(labels))
        figure, axes = plt.subplots(1, 2, figsize=(16, 5.5), gridspec_kw={"width_ratios": [1.7, 1]})
        axes[0].bar(x, summary["cells_before"], color="#bdbdbd", label="Before")
        axes[0].bar(x, summary["cells_after"], color="#4c78a8", label="After")
        axes[0].set_xticks(x, labels, rotation=35, ha="right")
        axes[0].set_ylabel("Cells")
        axes[0].legend(frameon=False)
        axes[1].bar(x, summary["pct_removed"], color="#e15759")
        axes[1].set_xticks(x, labels, rotation=35, ha="right")
        axes[1].set_ylabel("Removed (%)")
        for index, value in enumerate(summary["pct_removed"]):
            axes[1].text(index, value, f"{value:.2f}%", ha="center", va="bottom", fontsize=8)
        figure.suptitle(f"Step 02 cells before/after by {label.lower()} — approved per-sample 5-MAD union", fontsize=15)
        figure.tight_layout(rect=(0, 0, 1, 0.93))
        return figure

    @staticmethod
    def _reason_plot(summary: pd.DataFrame) -> plt.Figure:
        """Show criterion-specific removal counts for every technical sample."""

        labels = summary["group"].astype(str).tolist()
        x = np.arange(len(labels))
        figure, axis = plt.subplots(figsize=(16, 6))
        fields = (
            ("low_total_counts_n", "Low counts", "#4c78a8"),
            ("low_n_genes_n", "Low genes", "#59a14f"),
            ("high_pct_counts_mt_n", "High mt %", "#e15759"),
        )
        bottom = np.zeros(len(labels))
        for field, label, color in fields:
            values = summary[field].to_numpy(dtype=float)
            axis.bar(x, values, bottom=bottom, label=label, color=color)
            bottom += values
        axis.set_xticks(x, labels, rotation=35, ha="right")
        axis.set_ylabel("Criterion flags (cells)")
        axis.legend(frameon=False, ncol=3)
        axis.set_title("Step 02 exclusion evidence by sample; no doublet criterion")
        figure.tight_layout()
        return figure

    def _save(self, figure: plt.Figure, stem: str, scope: str, description: str) -> list[dict[str, str]]:
        """Save one figure in review and vector formats and return manifest rows."""

        rows = []
        for extension in ("png", "pdf"):
            path = self.figures_dir / f"{stem}.{extension}"
            figure.savefig(path, dpi=self.settings.plot_dpi if extension == "png" else None, bbox_inches="tight")
            rows.append({"relative_path": f"figures/{path.name}", "scope": scope, "format": extension, "description": description})
        plt.close(figure)
        return rows
