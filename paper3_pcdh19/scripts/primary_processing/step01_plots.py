"""Review-oriented QC plots for pooled, per-sample, and design-group data."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .step01_models import Step01Settings


class QCPlotPublisher:
    """Generate diagnostic figures without defining thresholds or exclusions."""

    DESIGN_COLORS = {
        "WT_M": "#4C78A8",
        "WT_F": "#72B7B2",
        "HET_F": "#F58518",
        "KO_M": "#E45756",
    }

    def __init__(self, settings: Step01Settings, figures_dir: Path):
        """Store rendering settings and create the staged figure directory."""

        self.settings = settings
        self.figures_dir = figures_dir
        self.figures_dir.mkdir(parents=True, exist_ok=True)
        (self.figures_dir / "per_sample").mkdir(exist_ok=True)
        self.records: list[dict[str, str]] = []

    def publish_all(self, obs: pd.DataFrame) -> pd.DataFrame:
        """Render pooled, 12 sample-specific, and grouped comparison figures."""

        self._diagnostic_panel(obs, "Pooled: all 12 samples", "pooled_qc_diagnostics", "pooled")
        for sample_id, frame in obs.groupby("technical_sample_id", observed=True, sort=False):
            genotype = frame["genotype"].iloc[0]
            sex = frame["sex"].iloc[0]
            design = frame["design_group"].iloc[0]
            title = f"{sample_id}: {genotype} {sex} ({design})"
            stem = f"per_sample/{sample_id}_qc_diagnostics"
            self._diagnostic_panel(frame, title, stem, str(sample_id))
        sample_order = list(dict.fromkeys(obs["technical_sample_id"].astype(str)))
        design_order = [value for value in ("WT_M", "WT_F", "HET_F", "KO_M") if value in set(obs["design_group"])]
        self._comparison_boxplots(obs, "technical_sample_id", sample_order, "QC distributions by sample", "qc_distributions_by_sample", "sample_comparison")
        self._comparison_boxplots(obs, "design_group", design_order, "QC distributions by design group", "qc_distributions_by_design_group", "design_group_comparison")
        return pd.DataFrame(self.records)

    def _diagnostic_panel(
        self,
        obs: pd.DataFrame,
        title: str,
        stem: str,
        scope: str,
    ) -> None:
        """Plot three distributions and the three requested pairwise relationships."""

        total = obs["total_counts"].to_numpy(dtype=float)
        genes = obs["n_genes_by_counts"].to_numpy(dtype=float)
        pct_mt = obs["pct_counts_mt"].to_numpy(dtype=float)
        log_total = np.log10(total + 1.0)
        figure, axes = plt.subplots(2, 3, figsize=(15, 9), constrained_layout=True)
        axes[0, 0].hist(log_total, bins=self.settings.histogram_bins, color="#4C78A8", alpha=0.9)
        axes[0, 0].set(xlabel="log10(total counts + 1)", ylabel="Cells", title="Total counts")
        axes[0, 1].hist(genes, bins=self.settings.histogram_bins, color="#59A14F", alpha=0.9)
        axes[0, 1].set(xlabel="Detected genes", ylabel="Cells", title="Detected genes")
        axes[0, 2].hist(pct_mt, bins=self.settings.histogram_bins, color="#E15759", alpha=0.9)
        axes[0, 2].set(xlabel="Mitochondrial counts (%)", ylabel="Cells", title="Mitochondrial fraction")
        self._hexbin(axes[1, 0], log_total, genes, "log10(total counts + 1)", "Detected genes")
        self._hexbin(axes[1, 1], log_total, pct_mt, "log10(total counts + 1)", "Mitochondrial counts (%)")
        self._hexbin(axes[1, 2], genes, pct_mt, "Detected genes", "Mitochondrial counts (%)")
        figure.suptitle(f"{title} — metrics only; no thresholds", fontsize=14)
        self._save(figure, stem, scope, "QC distributions and requested pairwise relationships; no threshold lines.")

    def _hexbin(
        self,
        axis: plt.Axes,
        x: np.ndarray,
        y: np.ndarray,
        xlabel: str,
        ylabel: str,
    ) -> None:
        """Render a density-aware relationship plot for large cell collections."""

        artist = axis.hexbin(x, y, gridsize=70, bins="log", mincnt=1, cmap="viridis", rasterized=True)
        axis.set(xlabel=xlabel, ylabel=ylabel)
        axis.figure.colorbar(artist, ax=axis, label="log10(cell count)")

    def _comparison_boxplots(
        self,
        obs: pd.DataFrame,
        group_column: str,
        order: list[str],
        title: str,
        stem: str,
        scope: str,
    ) -> None:
        """Compare full metric distributions without implying filtering cutoffs."""

        metrics = (
            ("total_counts", "log10(total counts + 1)", True),
            ("n_genes_by_counts", "Detected genes", False),
            ("pct_counts_mt", "Mitochondrial counts (%)", False),
        )
        figure, axes = plt.subplots(1, 3, figsize=(18, 6), constrained_layout=True)
        for axis, (metric, ylabel, log_transform) in zip(axes, metrics):
            values = []
            for group in order:
                group_values = obs.loc[obs[group_column].astype(str) == group, metric].to_numpy(dtype=float)
                values.append(np.log10(group_values + 1.0) if log_transform else group_values)
            artists = axis.boxplot(values, labels=order, whis=(5, 95), showfliers=False, patch_artist=True)
            for patch, group in zip(artists["boxes"], order):
                design = obs.loc[obs[group_column].astype(str) == group, "design_group"].iloc[0] if group_column == "technical_sample_id" else group
                patch.set_facecolor(self.DESIGN_COLORS.get(str(design), "#9D9D9D"))
                patch.set_alpha(0.8)
            axis.set(ylabel=ylabel, xlabel=group_column.replace("_", " "))
            axis.tick_params(axis="x", rotation=45 if len(order) > 4 else 0)
            axis.set_title("Box: IQR; whiskers: 5th–95th percentile")
        figure.suptitle(f"{title} — descriptive only; no thresholds", fontsize=14)
        self._save(figure, stem, scope, "Distribution comparison using IQR and 5th–95th percentile whiskers.")

    def _save(self, figure: plt.Figure, stem: str, scope: str, description: str) -> None:
        """Save both PNG and PDF copies and append plot-manifest rows."""

        for extension in ("png", "pdf"):
            relative = Path(f"{stem}.{extension}")
            output = self.figures_dir / relative
            output.parent.mkdir(parents=True, exist_ok=True)
            figure.savefig(output, dpi=self.settings.plot_dpi if extension == "png" else None, bbox_inches="tight")
            self.records.append({"relative_path": str(Path("figures") / relative), "scope": scope, "format": extension, "description": description})
        plt.close(figure)
