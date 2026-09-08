"""Render the four-page A-L report from already-computed diagnostics.

``Step06ReportPlotter.publish`` reuses one balanced display subset and shared
PCA/UMAP coordinates. Pages A-D show composition/QC/PCA; E-G show UMAP views;
H-I show provisional states/composition; J-L show mixing/similarity/flags.
No fitting, normalization, clustering, or annotation is performed here.

Numerical companion tables describe all cells. Scatter plots and the panel J
boxplot use the rendering subset; panel J's printed medians use all cells.
Colors, grids, marker sizes and font sizes are code-level presentation knobs.
The E-G grid assumes at most 12 samples; it is not a general layout engine.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages

from .step06_models import Step06Artifacts, Step06Results, Step06Settings
from .step06_progress import NullStep06ProgressTracker


class Step06Palette:
    """Provide fixed colors for samples and biological design variables.

    Notes
    -----
    Own deterministic sample colors and fixed genotype/sex/design mappings.
    Construct once per report and reuse across pages. It controls presentation only
    and does not modify labels or group membership.
    """

    def __init__(self, samples: list[str]):
        """Create stable mappings from sorted labels.

        Parameters
        ----------
        samples : list[str]
            Stable ordered sample IDs used across all report pages.

        Notes
        -----
        Store sample colors in supplied order and fixed WT/HET/KO, F/M and design-group
        colors. Caller supplies sorted samples. Colors affect display only.
        """

        colors = sns.color_palette("husl", len(samples))
        self.sample = dict(zip(samples, colors, strict=True))
        self.genotype = {"WT": "#4c78a8", "HET": "#f58518", "KO": "#e45756"}
        self.sex = {"F": "#b279a2", "M": "#59a14f"}
        self.design = {
            "WT_M": "#4c78a8",
            "WT_F": "#72b7b2",
            "HET_F": "#f58518",
            "KO_M": "#e45756",
        }


class Step06ReportPlotter:
    """Assemble one primary multi-page A-L report from shared coordinates.

    Notes
    -----
    Own PDF/PNG assembly from existing coordinates and summary tables. publish()
    selects saved display indices, calls four page builders, writes assets and
    closes figures. The plotting class never recomputes PCA, neighbors, UMAP or
    cell-state labels.
    """

    def __init__(self, settings: Step06Settings, output_dir: Path, progress=None):
        """Store plotting settings, output location, and progress publisher.

        Parameters
        ----------
        settings : Step06Settings
            Resolved scientific/rendering controls; see the settings class and tuning
            guide.
        output_dir : pathlib.Path
            Directory for rendered PDF and PNG assets; created if absent.
        progress : Step06ProgressTracker or None
            Optional event publisher; None selects the no-I/O tracker.

        Notes
        -----
        Create output_dir and store rendering controls/progress adapter. Does not fit or
        sample cells.
        """

        self.settings = settings
        self.output_dir = output_dir
        self.progress = progress or NullStep06ProgressTracker()
        self.output_dir.mkdir(parents=True, exist_ok=True)

    def publish(
        self,
        obs: pd.DataFrame,
        artifacts: Step06Artifacts,
        results: Step06Results,
    ) -> pd.DataFrame:
        """Write the primary PDF and page previews with a manifest.

        Parameters
        ----------
        obs : pandas.DataFrame
            All retained cells in artifact row order, with sample/design and raw QC
            columns.
        artifacts : Step06Artifacts
            Full-cell coordinates, graphs, scores and HVG/pseudobulk results in input
            order.
        results : Step06Results
            Companion tables plus the automatic, unreviewed diagnostic outcome.

        Returns
        -------
        pandas.DataFrame
            asset/role manifest for one four-page PDF and four PNG previews.

        Notes
        -----
        Slice the saved rendering positions once, build each page from shared coordinates and
        full-data tables, save figures, then close them. plot_dpi affects PNGs; PDF scatter
        artists are rasterized. Only rendering uses the subset.
        """

        sample = self.settings.sample_field
        samples = sorted(obs[sample].astype(str).unique())
        palette = Step06Palette(samples)
        # These indices affect display only; do not recompute embeddings on this subset.
        positions = artifacts.rendering_cells["row_position"].to_numpy(int)
        plot_obs = obs.iloc[positions].copy()
        plot_annotations = artifacts.annotations.iloc[positions].copy()
        pca = artifacts.pca[positions]
        umap = artifacts.umap[positions]
        pdf_path = self.output_dir / "step06_primary_diagnostic_report.pdf"
        page_builders = [
            (
                "page_1_panels_A_D",
                "Step06ReportPlotter._page_a_d",
                lambda: self._page_a_d(
                    obs, plot_obs, pca, artifacts, results, palette
                ),
            ),
            (
                "page_2_panels_E_G",
                "Step06ReportPlotter._page_e_g",
                lambda: self._page_e_g(plot_obs, umap, palette),
            ),
            (
                "page_3_panels_H_I",
                "Step06ReportPlotter._page_h_i",
                lambda: self._page_h_i(
                    plot_obs, plot_annotations, umap, artifacts, results
                ),
            ),
            (
                "page_4_panels_J_L",
                "Step06ReportPlotter._page_j_l",
                lambda: self._page_j_l(
                    plot_obs, plot_annotations, results, palette
                ),
            ),
        ]
        pages = []
        for name, function, builder in page_builders:
            with self.progress.track(
                "output.report.page",
                function,
                {
                    "page": name,
                    "all_numerical_cells": len(obs),
                    "rendering_cells": len(plot_obs),
                    "rendering_only_downsampling": True,
                    "coordinate_sources": ["X_pca", "X_umap"],
                },
            ) as event:
                figure = builder()
                pages.append((name, figure))
                event.outputs["figure_axes"] = len(figure.axes)
        with PdfPages(pdf_path) as pdf:
            for _, figure in pages:
                pdf.savefig(figure, bbox_inches="tight")
        rows = [{"asset": pdf_path.name, "role": "primary A-L multi-page review report"}]
        for name, figure in pages:
            path = self.output_dir / f"{name}.png"
            figure.savefig(path, dpi=self.settings.plot_dpi, bbox_inches="tight")
            plt.close(figure)
            rows.append({"asset": path.name, "role": "page preview"})
        return pd.DataFrame(rows)

    def _page_a_d(self, obs, plot_obs, pca, artifacts, results, palette):
        """Create panels A-D for composition, QC, PCA, and centroids.

        Parameters
        ----------
        obs : pandas.DataFrame
            All retained cells in artifact row order, with sample/design and raw QC
            columns.
        plot_obs : pandas.DataFrame
            Rendering-subset metadata aligned with the supplied coordinate rows.
        pca : numpy.ndarray
            Rendering-subset PCA coordinates; first two components are displayed.
        artifacts : Step06Artifacts
            Full-cell coordinates, graphs, scores and HVG/pseudobulk results in input
            order.
        results : Step06Results
            Companion tables plus the automatic, unreviewed diagnostic outcome.
        palette : Step06Palette
            Stable mappings from sample/design labels to display colors.

        Returns
        -------
        matplotlib.figure.Figure
            Unsaved figure; publish() writes and closes it.

        Notes
        -----
        A/B counts and QC summaries use all cells. C displays PCs 1/2 of the saved subset by
        sample/genotype/sex. D plots all-cell sample centers on PCs 1/2, although numeric
        distances use centroid_components PCs. Figure dimensions, legends and fonts are
        code-level knobs.
        """

        figure = plt.figure(figsize=(21, 15))
        grid = figure.add_gridspec(2, 3)
        axes = np.asarray([
            [figure.add_subplot(grid[0, 0]), figure.add_subplot(grid[0, 1]), figure.add_subplot(grid[0, 2])],
            [figure.add_subplot(grid[1, 0]), figure.add_subplot(grid[1, 1]), figure.add_subplot(grid[1, 2])],
        ])
        sample = self.settings.sample_field
        counts = obs[sample].astype(str).value_counts().sort_index()
        axes[0, 0].barh(counts.index, counts.values, color=[palette.sample[x] for x in counts.index])
        axes[0, 0].set_xlabel("Retained cells")
        axes[0, 0].set_title("A  Sample composition (all cells)")

        qc = results.tables["qc_summary_by_sample.tsv"].set_index(sample).reindex(counts.index)
        x = np.arange(len(qc))
        axes[0, 1].plot(x, qc["total_counts_median"], "o-", label="median UMIs")
        axes[0, 1].plot(x, qc["n_genes_median"], "s-", label="median genes")
        axes[0, 1].set_xticks(x, qc.index, rotation=75)
        axes[0, 1].set_ylabel("Median")
        twin = axes[0, 1].twinx()
        twin.plot(x, qc["pct_mt_median"], "^-", color="#e45756", label="median % mt")
        twin.set_ylabel("Median mitochondrial %")
        axes[0, 1].set_title("B  Per-sample QC summaries (all cells)")

        variance = artifacts.variance_ratio
        self._pca_scatter(
            axes[0, 2], pca, plot_obs[sample], palette.sample,
            "C  Unintegrated PCA by biological sample", variance,
        )
        self._pca_scatter(
            axes[1, 0], pca, plot_obs[self.settings.genotype_field], palette.genotype,
            "C  Same PCA by genotype", variance,
        )
        self._pca_scatter(
            axes[1, 1], pca, plot_obs[self.settings.sex_field], palette.sex,
            "C  Same PCA by sex", variance,
        )

        centroids = results.tables["pca_centroids.tsv"]
        for row in centroids.itertuples(index=False):
            value = getattr(row, sample)
            axes[1, 2].scatter(row.PC1, row.PC2, s=90, color=palette.sample[value])
            axes[1, 2].text(row.PC1, row.PC2, value, fontsize=7)
        axes[1, 2].set_xlabel("PC1 centroid")
        axes[1, 2].set_ylabel("PC2 centroid")
        axes[1, 2].set_title("D  Sample centroids in unintegrated PCA")
        figure.suptitle("Step 06: sample composition, QC, and unintegrated PCA", fontsize=18)
        figure.tight_layout()
        return figure

    def _page_e_g(self, plot_obs, umap, palette):
        """Create panels E-G from one shared global UMAP.

        Parameters
        ----------
        plot_obs : pandas.DataFrame
            Rendering-subset metadata aligned with the supplied coordinate rows.
        umap : numpy.ndarray
            Rendering-subset coordinates from the one global UMAP, shape (n_display, 2).
        palette : Step06Palette
            Stable mappings from sample/design labels to display colors.

        Returns
        -------
        matplotlib.figure.Figure
            Unsaved figure; publish() writes and closes it.

        Notes
        -----
        E/G recolor identical global UMAP coordinates. F highlights one sample at a time
        against the same rendering subset. No sample-specific UMAP is recomputed. The fixed
        4x4 grid allocates twelve sample facets; layout must change for more samples.
        """

        figure = plt.figure(figsize=(20, 18))
        grid = figure.add_gridspec(4, 4)
        sample = self.settings.sample_field
        ax_e = figure.add_subplot(grid[0, :2])
        self._categorical_scatter(ax_e, umap, plot_obs[sample], palette.sample, "E  Unintegrated UMAP by sample")
        design_grid = grid[0, 2:].subgridspec(1, 3)
        child = figure.add_subplot(design_grid[0, 0])
        self._categorical_scatter(child, umap, plot_obs[self.settings.genotype_field], palette.genotype, "G  Genotype")
        child = figure.add_subplot(design_grid[0, 1])
        self._categorical_scatter(child, umap, plot_obs[self.settings.sex_field], palette.sex, "G  Sex")
        child = figure.add_subplot(design_grid[0, 2])
        self._categorical_scatter(child, umap, plot_obs[self.settings.design_field], palette.design, "G  Design group")
        samples = sorted(plot_obs[sample].astype(str).unique())
        for index, value in enumerate(samples):
            ax = figure.add_subplot(grid[1 + index // 4, index % 4])
            ax.scatter(umap[:, 0], umap[:, 1], s=0.4, color="#d9d9d9", rasterized=True)
            mask = plot_obs[sample].astype(str).to_numpy() == value
            ax.scatter(umap[mask, 0], umap[mask, 1], s=0.8, color=palette.sample[value], rasterized=True)
            ax.set_title(value, fontsize=9)
            ax.set_xticks([])
            ax.set_yticks([])
        figure.suptitle(
            "Step 06: shared unintegrated UMAP views\n"
            "F: identical global coordinates; gray = full rendering subset, color = selected sample",
            fontsize=18,
        )
        figure.tight_layout(rect=(0, 0, 1, 0.95))
        return figure

    def _page_h_i(self, plot_obs, plot_annotations, umap, artifacts, results):
        """Create panels H-I for provisional biology and sample contributions.

        Parameters
        ----------
        plot_obs : pandas.DataFrame
            Rendering-subset metadata aligned with the supplied coordinate rows.
        plot_annotations : pandas.DataFrame
            Rendering-subset program scores and provisional states in display row order.
        umap : numpy.ndarray
            Rendering-subset coordinates from the one global UMAP, shape (n_display, 2).
        artifacts : Step06Artifacts
            Full-cell coordinates, graphs, scores and HVG/pseudobulk results in input
            order.
        results : Step06Results
            Companion tables plus the automatic, unreviewed diagnostic outcome.

        Returns
        -------
        matplotlib.figure.Figure
            Unsaved figure; publish() writes and closes it.

        Notes
        -----
        H displays automatic cluster-program winner labels and mean standardized scores. I
        uses all-cell counts and percentages from companion tables. State labels are
        provisional, not validated identities. The scatter uses only the rendering subset.
        """

        figure = plt.figure(figsize=(22, 16))
        grid = figure.add_gridspec(2, 3)
        ax_state = figure.add_subplot(grid[0, :2])
        ax_program = figure.add_subplot(grid[0, 2])
        ax_counts = figure.add_subplot(grid[1, 0])
        ax_sample_pct = figure.add_subplot(grid[1, 1])
        ax_state_pct = figure.add_subplot(grid[1, 2])
        states = sorted(plot_annotations["step06_provisional_state"].unique())
        colors = dict(zip(states, sns.color_palette("tab10", len(states)), strict=True))
        self._categorical_scatter(ax_state, umap, plot_annotations["step06_provisional_state"], colors, "H  Preliminary broad states")

        cluster = artifacts.cluster_program_scores.set_index("leiden")
        matrix = cluster.drop(columns="provisional_state")
        sns.heatmap(matrix, cmap="vlag", center=0, ax=ax_program, cbar_kws={"label": "Mean standardized score"})
        ax_program.set_title("H  Marker programs by descriptive Leiden cluster")
        ax_program.set_xlabel("Program")

        sample = self.settings.sample_field
        counts = results.tables["provisional_state_counts.tsv"].pivot(index=sample, columns="provisional_state", values="cells")
        pct = results.tables["provisional_state_within_sample_percentages.tsv"].pivot(index=sample, columns="provisional_state", values="within_sample_pct")
        state_pct = results.tables["provisional_state_within_state_percentages.tsv"].pivot(index=sample, columns="provisional_state", values="within_state_pct")
        sns.heatmap(counts, cmap="mako", ax=ax_counts, cbar_kws={"label": "Cells"})
        ax_counts.set_title("I  Absolute sample × state counts")
        sns.heatmap(pct, cmap="rocket", ax=ax_sample_pct, cbar_kws={"label": "Within-sample %"})
        ax_sample_pct.set_title("I  Within-sample state percentages")
        sns.heatmap(state_pct, cmap="crest", ax=ax_state_pct, cbar_kws={"label": "Within-state %"})
        ax_state_pct.set_title("I  Within-state sample contributions")
        figure.suptitle("Step 06: descriptive biological structure (not final annotation)", fontsize=18)
        figure.tight_layout()
        return figure

    def _page_j_l(self, plot_obs, plot_annotations, results, palette):
        """Create panels J-L for mixing, similarity, and decision review.

        Parameters
        ----------
        plot_obs : pandas.DataFrame
            Rendering-subset metadata aligned with the supplied coordinate rows.
        plot_annotations : pandas.DataFrame
            Rendering-subset program scores and provisional states in display row order.
        results : Step06Results
            Companion tables plus the automatic, unreviewed diagnostic outcome.
        palette : Step06Palette
            Stable mappings from sample/design labels to display colors.

        Returns
        -------
        matplotlib.figure.Figure
            Unsaved figure; publish() writes and closes it.

        Notes
        -----
        J boxplots use the rendering subset; printed entropy/effective-sample medians use ALL
        cells. K plots HVG pseudobulk Pearson correlations on a fixed 0..1 color scale
        (negative values would saturate). L copies heuristic sample flags and headline without
        making a new decision.
        """

        figure = plt.figure(figsize=(20, 15))
        grid = figure.add_gridspec(2, 2, height_ratios=[1, 1.2])
        ax_j = figure.add_subplot(grid[0, 0])
        diversity = results.tables["local_sample_diversity_per_cell.tsv.gz"]
        render_ids = set(plot_obs.index.astype(str))
        # Panel J boxplot is display-subset data; numeric medians immediately below use all cells.
        shown = diversity.loc[diversity["cell_id"].astype(str).isin(render_ids)]
        sns.boxplot(data=shown, x="provisional_state", y="local_effective_samples", color="#72b7b2", showfliers=False, ax=ax_j)
        ax_j.tick_params(axis="x", rotation=70)
        ax_j.set_title("J  Local effective sample diversity by provisional state")
        global_median = diversity["local_effective_samples"].median()
        global_entropy = diversity["local_sample_entropy_normalized"].median()
        ax_j.text(
            0.02, 0.98,
            f"All-cell median effective samples: {global_median:.2f}\nAll-cell median normalized entropy: {global_entropy:.3f}",
            transform=ax_j.transAxes, va="top", fontsize=8,
            bbox={"facecolor": "white", "alpha": 0.8, "edgecolor": "none"},
        )

        ax_k = figure.add_subplot(grid[0, 1])
        corr = results.tables["pseudobulk_correlations.tsv"].pivot(index="sample_1", columns="sample_2", values="pseudobulk_pearson_r")
        metadata = results.tables["sample_metadata_and_cells.tsv"].set_index(self.settings.sample_field)
        labels = [
            f"{value} | {metadata.loc[value, 'genotype']}/{metadata.loc[value, 'sex']}"
            for value in corr.index
        ]
        corr.index = labels
        corr.columns = labels
        sns.heatmap(corr, vmin=0, vmax=1, cmap="vlag", square=True, ax=ax_k)
        ax_k.set_title("K  Sample pseudobulk correlations across HVGs")

        ax_l = figure.add_subplot(grid[1, :])
        ax_l.axis("off")
        decision = results.tables["preintegration_sample_decisions.tsv"].copy()
        columns = [
            self.settings.sample_field,
            "cells",
            "median_total_counts",
            "median_n_genes",
            "median_pct_mt",
            "major_pca_displacement",
            "major_umap_displacement",
            "broad_states_represented",
            "evidence_of_missing_states",
            "unusual_sample_specific_structure",
            "interpretation",
        ]
        table = ax_l.table(cellText=decision[columns].astype(str).values, colLabels=columns, loc="upper center", cellLoc="left")
        table.auto_set_font_size(False)
        table.set_fontsize(6)
        table.scale(1, 1.35)
        ax_l.set_title(f"L  Provisional outcome {results.outcome_code}: {results.outcome_label}\n{results.outcome_rationale}\nNo integration decision is automatic.", fontsize=13, pad=20)
        figure.suptitle("Step 06: quantitative mixing, similarity, and decision summary", fontsize=18)
        figure.tight_layout()
        return figure

    @staticmethod
    def _categorical_scatter(ax, coordinates, labels, colors, title):
        """Draw one rasterized categorical scatter with stable colors.

        Parameters
        ----------
        ax : matplotlib.axes.Axes
            Existing axes to mutate with a scatter plot.
        coordinates : numpy.ndarray
            Coordinates aligned row-for-row with labels; only columns 0 and 1 are drawn.
        labels : pandas.Series
            One categorical label per coordinate row.
        colors : dict
            Label-to-color mapping; unknown labels use gray.
        title : str
            Title displayed on the axes.

        Returns
        -------
        None
            Effects are described below.

        Notes
        -----
        Draw coordinate columns 0/1 per label, with fixed size=0.7, alpha=0.45 and
        rasterization. Missing color keys become gray; hide ticks and add legend. These are
        presentation-only code literals.
        """

        values = labels.astype(str).to_numpy()
        for value in sorted(np.unique(values)):
            mask = values == value
            ax.scatter(coordinates[mask, 0], coordinates[mask, 1], s=0.7, alpha=0.45, rasterized=True, color=colors.get(value, "#777777"), label=value)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_title(title)
        ax.legend(markerscale=6, fontsize=6, ncol=2)

    @staticmethod
    def _pca_scatter(ax, coordinates, labels, colors, title, variance_ratio):
        """Draw one view of the exact shared unintegrated PCA coordinates.

        Parameters
        ----------
        ax : matplotlib.axes.Axes
            Existing axes to mutate with a scatter plot.
        coordinates : numpy.ndarray
            Coordinates aligned row-for-row with labels; only columns 0 and 1 are drawn.
        labels : pandas.Series
            One categorical label per coordinate row.
        colors : dict
            Label-to-color mapping; unknown labels use gray.
        title : str
            Title displayed on the axes.
        variance_ratio : numpy.ndarray
            PCA variance fractions; first two values appear as axis percentages.

        Returns
        -------
        None
            Effects are described below.

        Notes
        -----
        Delegate categorical rendering and label axes with the first two variance percentages.
        Reuses existing PCA; no fitting occurs.
        """

        Step06ReportPlotter._categorical_scatter(ax, coordinates, labels, colors, title)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_xlabel(f"PC1 ({100 * variance_ratio[0]:.1f}%)")
        ax.set_ylabel(f"PC2 ({100 * variance_ratio[1]:.1f}%)")
