"""Generate review plots for Step 03 calls, scores, composition, and PCA."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .step03_models import Step03Settings


class Step03PlotPublisher:
    """Publish deterministic primary, composition, reproducibility, and PCA plots."""

    def __init__(self, settings: Step03Settings, figures_dir: Path):
        """Store rendering settings and the staged figure directory."""

        self.settings = settings
        self.figures_dir = figures_dir

    def publish_all(self, frame: pd.DataFrame, pca: np.ndarray) -> pd.DataFrame:
        """Render every required diagnostic in PNG and vector PDF forms."""

        records: list[dict[str, str]] = []
        specifications = (
            (self._score_distribution(frame), "primary_score_distribution", "score", "Primary scDblFinder score distribution and empirical call boundary."),
            (self._group_score_distributions(frame, self.settings.sample_field, "Technical sample"), "score_distributions_by_technical_sample", "sample_score", "Primary score distributions for each retained probe-barcode technical sample."),
            (self._group_score_distributions(frame, self.settings.design_field, "Design group"), "score_distributions_by_design_group", "design_score", "Primary score distributions by design group for reporting only."),
            (self._called_fraction(frame, self.settings.sample_field, "Technical sample"), "called_fraction_by_technical_sample", "sample", "Called-doublet fraction for each retained probe-barcode technical sample."),
            (self._called_fraction(frame, self.settings.design_field, "Design group"), "called_fraction_by_design_group", "design_group", "Called-doublet fraction by design group for reporting only."),
            (self._called_fraction(frame, "primary_cluster", "Generated cluster"), "called_fraction_by_generated_cluster", "generated_cluster", "Called-doublet fraction by the cluster generated internally because clusters=TRUE."),
            (self._reproducibility(frame), "primary_vs_replicate_reproducibility", "reproducibility", "Primary-versus-repeat score density and class concordance under a second seed."),
            (self._pca_score(frame, pca), "internal_pca_mean_score", "expression_space", "Mean primary score across the exact internal scDblFinder PCA expression space."),
            (self._pca_called_fraction(frame, pca), "internal_pca_called_fraction", "expression_space", "Called-doublet fraction across the exact internal scDblFinder PCA expression space."),
        )
        for figure, stem, scope, description in specifications:
            records.extend(self._save(figure, stem, scope, description))
        return pd.DataFrame(records)

    @staticmethod
    def _score_distribution(frame: pd.DataFrame) -> plt.Figure:
        """Plot the complete primary score distribution split by final call."""

        figure, axis = plt.subplots(figsize=(10, 6))
        bins = np.linspace(0, 1, 151)
        for label, color in (("singlet", "#4c78a8"), ("doublet", "#e15759")):
            values = frame.loc[frame["primary_class"] == label, "primary_score"].to_numpy(float)
            axis.hist(values, bins=bins, histtype="stepfilled", alpha=0.55, color=color, label=f"{label} (n={len(values):,})")
        doublet_scores = frame.loc[frame["primary_class"] == "doublet", "primary_score"]
        if len(doublet_scores):
            axis.axvline(float(doublet_scores.min()), color="#8b0000", linestyle="--", linewidth=1.5, label="Lowest called-doublet score")
        axis.set_yscale("log")
        axis.set_xlabel("scDblFinder primary score")
        axis.set_ylabel("Cells (log scale)")
        axis.set_title("Step 03 primary scDblFinder score distribution — one GEX_1 capture")
        axis.legend(frameon=False)
        figure.tight_layout()
        return figure

    @staticmethod
    def _group_score_distributions(frame: pd.DataFrame, field: str, title: str) -> plt.Figure:
        """Plot complete score histograms in small multiples for one grouping."""

        groups = sorted(frame[field].astype(str).unique())
        columns = 4 if len(groups) > 4 else 2
        rows = int(np.ceil(len(groups) / columns))
        figure, axes = plt.subplots(rows, columns, figsize=(4.2 * columns, 3.1 * rows), sharex=True, sharey=True, squeeze=False)
        bins = np.linspace(0, 1, 101)
        for axis, group in zip(axes.ravel(), groups):
            values = frame.loc[frame[field].astype(str) == group, "primary_score"].to_numpy(float)
            axis.hist(values, bins=bins, color="#4c78a8", alpha=0.8)
            axis.set_yscale("log")
            axis.set_title(f"{group} (n={len(values):,})", fontsize=10)
        for axis in axes.ravel()[len(groups) :]:
            axis.set_visible(False)
        figure.supxlabel("scDblFinder primary score")
        figure.supylabel("Cells (log scale)")
        figure.suptitle(f"Step 03 primary score distributions by {title.lower()}", fontsize=15)
        figure.tight_layout(rect=(0.02, 0.02, 1, 0.96))
        return figure

    @staticmethod
    def _called_fraction(frame: pd.DataFrame, field: str, title: str) -> plt.Figure:
        """Plot called fractions for one reporting-only metadata grouping."""

        summary = frame.groupby(field, observed=True)["primary_class"].agg(cells="size", doublets=lambda x: int((x == "doublet").sum())).reset_index()
        summary["fraction"] = summary["doublets"] / summary["cells"]
        summary = summary.sort_values(field)
        figure, axis = plt.subplots(figsize=(max(9, min(24, 0.55 * len(summary))), 6))
        x = np.arange(len(summary))
        axis.bar(x, 100 * summary["fraction"], color="#e15759")
        axis.set_xticks(x, summary[field].astype(str), rotation=35, ha="right")
        axis.set_ylabel("Called doublets (%)")
        axis.set_title(f"Step 03 called-doublet fraction by {title.lower()}")
        for index, row in summary.reset_index(drop=True).iterrows():
            axis.text(index, 100 * row["fraction"], f"{100*row['fraction']:.2f}%\n({row['doublets']:,})", ha="center", va="bottom", fontsize=8)
        figure.tight_layout()
        return figure

    @staticmethod
    def _reproducibility(frame: pd.DataFrame) -> plt.Figure:
        """Show score agreement and the two-run classification confusion matrix."""

        figure, axes = plt.subplots(1, 2, figsize=(14, 6))
        axes[0].hexbin(frame["primary_score"], frame["replicate_score"], gridsize=100, bins="log", mincnt=1, cmap="viridis")
        axes[0].plot([0, 1], [0, 1], color="white", linewidth=1, linestyle="--")
        axes[0].set(xlabel="Primary score", ylabel="Second-seed score", title="Score reproducibility (hexbin density)")
        confusion = pd.crosstab(frame["primary_class"], frame["replicate_class"]).reindex(index=["singlet", "doublet"], columns=["singlet", "doublet"], fill_value=0)
        image = axes[1].imshow(confusion.to_numpy(), cmap="Blues")
        axes[1].set_xticks([0, 1], confusion.columns)
        axes[1].set_yticks([0, 1], confusion.index)
        axes[1].set_xlabel("Second-seed call")
        axes[1].set_ylabel("Primary call")
        axes[1].set_title("Call concordance")
        for row in range(2):
            for column in range(2):
                axes[1].text(column, row, f"{confusion.iloc[row, column]:,}", ha="center", va="center", color="black")
        figure.colorbar(image, ax=axes[1], fraction=0.046)
        figure.tight_layout()
        return figure

    @staticmethod
    def _binned_pca(pca: np.ndarray, values: np.ndarray, reducer: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Aggregate all cells into stable PCA-plane bins for interpretable plots."""

        x, y = pca[:, 0], pca[:, 1]
        x_limits = np.quantile(x, [0.001, 0.999])
        y_limits = np.quantile(y, [0.001, 0.999])
        x_edges = np.linspace(*x_limits, 181)
        y_edges = np.linspace(*y_limits, 181)
        x = np.clip(x, np.nextafter(x_limits[0], x_limits[1]), np.nextafter(x_limits[1], x_limits[0]))
        y = np.clip(y, np.nextafter(y_limits[0], y_limits[1]), np.nextafter(y_limits[1], y_limits[0]))
        counts, _, _ = np.histogram2d(x, y, bins=(x_edges, y_edges))
        sums, _, _ = np.histogram2d(x, y, bins=(x_edges, y_edges), weights=values)
        with np.errstate(divide="ignore", invalid="ignore"):
            aggregate = sums / counts if reducer == "mean" else counts
        aggregate[counts == 0] = np.nan
        return x_edges, y_edges, aggregate.T

    def _pca_score(self, frame: pd.DataFrame, pca: np.ndarray) -> plt.Figure:
        """Plot local mean primary score on the internal PCA plane."""

        x_edges, y_edges, values = self._binned_pca(pca, frame["primary_score"].to_numpy(float), "mean")
        return self._pca_heatmap(x_edges, y_edges, values, "Mean scDblFinder score", "Step 03 internal PCA — local mean primary score", "magma", 0, 1)

    def _pca_called_fraction(self, frame: pd.DataFrame, pca: np.ndarray) -> plt.Figure:
        """Plot local called-doublet fraction on the internal PCA plane."""

        calls = (frame["primary_class"].to_numpy(str) == "doublet").astype(float)
        x_edges, y_edges, values = self._binned_pca(pca, calls, "mean")
        return self._pca_heatmap(x_edges, y_edges, values, "Called-doublet fraction", "Step 03 internal PCA — local called-doublet fraction", "viridis", 0, 1)

    @staticmethod
    def _pca_heatmap(x_edges: np.ndarray, y_edges: np.ndarray, values: np.ndarray, label: str, title: str, cmap: str, vmin: float, vmax: float) -> plt.Figure:
        """Render one complete-cell binned PCA diagnostic."""

        figure, axis = plt.subplots(figsize=(9, 7))
        image = axis.pcolormesh(x_edges, y_edges, values, cmap=cmap, shading="auto", vmin=vmin, vmax=vmax)
        figure.colorbar(image, ax=axis, label=label)
        axis.set(xlabel="scDblFinder internal PC1", ylabel="scDblFinder internal PC2", title=title)
        figure.tight_layout()
        return figure

    def _save(self, figure: plt.Figure, stem: str, scope: str, description: str) -> list[dict[str, str]]:
        """Save one diagnostic as review PNG and vector PDF with manifest rows."""

        records = []
        for extension in ("png", "pdf"):
            path = self.figures_dir / f"{stem}.{extension}"
            figure.savefig(path, dpi=self.settings.plot_dpi if extension == "png" else None, bbox_inches="tight")
            records.append({"relative_path": f"figures/{path.name}", "scope": scope, "format": extension, "description": description})
        plt.close(figure)
        return records
