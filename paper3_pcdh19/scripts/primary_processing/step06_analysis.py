"""Full-cell unintegrated preprocessing and descriptive state analysis."""

from __future__ import annotations

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse

from .step06_models import Step06Artifacts, Step06Settings
from .step06_progress import NullStep06ProgressTracker


def marker_programs() -> dict[str, tuple[str, ...]]:
    """Return the fixed broad marker programs used for descriptive states."""

    return {
        "Cycling progenitor": ("Mki67", "Top2a", "Pcna", "Cdk1", "Ccnb1"),
        "RG/progenitor-like": ("Sox2", "Hes1", "Hes5", "Fabp7", "Nes", "Vim"),
        "IPC/neuroblast": ("Eomes", "Neurog2", "Neurod1", "Dcx"),
        "Immature inhibitory neuron": (
            "Dcx", "Tubb3", "Gad1", "Gad2", "Dlx1", "Dlx2", "Dlx5", "Dlx6", "Slc32a1"
        ),
        "MGE identity": ("Nkx2-1", "Lhx6", "Lhx8", "Ascl1"),
        "Excitatory-like": ("Slc17a6", "Slc17a7", "Neurog2", "Eomes"),
        "Gliogenic/oligodendroglial": ("Aldoc", "Aqp4", "Pdgfra", "Olig1", "Olig2", "Mbp"),
        "Immune/microglial": ("Ptprc", "C1qa", "C1qb", "C1qc"),
        "Endothelial/mesenchymal": ("Pecam1", "Kdr", "Col1a1", "Col1a2"),
        "Ependymal/choroid-like": ("Foxj1", "Ttr"),
    }


class PseudobulkBuilder:
    """Aggregate all raw counts by sample and return log-CPM expression."""

    def build(self, adata, sample_field: str, gene_mask: np.ndarray) -> pd.DataFrame:
        """Calculate sample pseudobulks using every retained cell."""

        samples = pd.Index(adata.obs[sample_field].astype(str).unique()).sort_values()
        codes = pd.Categorical(adata.obs[sample_field].astype(str), categories=samples).codes
        design = sparse.csr_matrix(
            (np.ones(adata.n_obs, dtype=np.float32), (codes, np.arange(adata.n_obs))),
            shape=(len(samples), adata.n_obs),
        )
        counts = design @ adata.X[:, gene_mask]
        totals = np.asarray(counts.sum(axis=1)).ravel()
        if np.any(totals <= 0):
            raise ValueError("A sample pseudobulk has zero total counts")
        logcpm = np.log1p(counts.toarray() * (1_000_000.0 / totals[:, None]))
        return pd.DataFrame(logcpm, index=samples, columns=adata.var_names[gene_mask])


class ProgramScorer:
    """Calculate transparent mean-log-expression programs and cluster labels."""

    def score(self, adata, clusters: pd.Series) -> tuple[pd.DataFrame, pd.DataFrame]:
        """Return per-cell scores and cluster-level standardized program means."""

        symbols = pd.Series(np.arange(adata.n_vars), index=adata.var["gene_symbol"].astype(str))
        scores: dict[str, np.ndarray] = {}
        for name, genes in marker_programs().items():
            present = [gene for gene in genes if gene in symbols.index]
            if not present:
                raise ValueError(f"No genes available for marker program: {name}")
            positions = symbols.loc[present].to_numpy(int)
            scores[name] = np.asarray(adata.X[:, positions].mean(axis=1)).ravel()
        frame = pd.DataFrame(scores, index=adata.obs_names)
        standard = (frame - frame.mean()) / frame.std(ddof=0).replace(0, 1)
        grouped = standard.assign(leiden=clusters.astype(str).to_numpy()).groupby(
            "leiden", observed=True
        ).mean()
        grouped["provisional_state"] = grouped.idxmax(axis=1)
        safe_names = {
            name: "program_" + "".join(
                character if character.isalnum() else "_" for character in name
            ).strip("_")
            for name in frame.columns
        }
        return frame.rename(columns=safe_names), grouped.reset_index()


class RenderingSampler:
    """Choose a deterministic balanced subset exclusively for plotting."""

    def __init__(self, settings: Step06Settings):
        """Store the fixed cap and random seed."""

        self.settings = settings

    def select(self, obs: pd.DataFrame) -> pd.DataFrame:
        """Return selected cell IDs with their full-data row positions."""

        rng = np.random.default_rng(self.settings.random_seed)
        groups = obs.groupby(self.settings.sample_field, observed=True).indices
        per_sample = max(1, self.settings.render_max_cells // len(groups))
        selected: list[int] = []
        for sample in sorted(groups):
            positions = np.asarray(groups[sample], dtype=int)
            take = min(per_sample, len(positions))
            selected.extend(rng.choice(positions, size=take, replace=False).tolist())
        selected = sorted(selected)
        return pd.DataFrame(
            {"cell_id": obs.index[selected].astype(str), "row_position": selected}
        )


class Step06Analyzer:
    """Run the documented diagnostic workflow without batch correction."""

    def __init__(self, settings: Step06Settings, progress=None):
        """Store immutable settings and the run-scoped progress publisher."""

        self.settings = settings
        self.progress = progress or NullStep06ProgressTracker()

    def run(self, adata) -> Step06Artifacts:
        """Compute HVGs, PCA, neighbors, UMAP, Leiden, states, and pseudobulk."""

        self.progress.note(
            "analysis",
            "Step06Analyzer.run",
            {
                "object_type": type(adata).__name__,
                "shape": list(adata.shape),
                "matrix_type": type(adata.X).__name__,
                "matrix_dtype": str(adata.X.dtype),
                "sample_field": self.settings.sample_field,
                "sample_count": int(adata.obs[self.settings.sample_field].nunique()),
            },
            {"analysis_mode": "unintegrated", "cells_removed": 0},
        )
        sc.settings.n_jobs = self.settings.n_jobs
        with self.progress.track(
            "analysis.hvg",
            "scanpy.pp.highly_variable_genes",
            {
                "input_representation": "approved Step 02 raw integer adata.X",
                "input_shape": list(adata.shape),
                "flavor": self.settings.hvg_flavor,
                "n_top_genes": self.settings.n_top_genes,
                "batch_key": self.settings.sample_field,
                "subset": False,
            },
        ) as event:
            sc.pp.highly_variable_genes(
                adata,
                flavor=self.settings.hvg_flavor,
                n_top_genes=self.settings.n_top_genes,
                batch_key=self.settings.sample_field,
                subset=False,
            )
            event.outputs["highly_variable_genes"] = int(
                adata.var["highly_variable"].sum()
            )
        hvg = adata.var["highly_variable"].to_numpy(bool)
        hvg_columns = [
            column
            for column in (
                "gene_symbol",
                "highly_variable",
                "highly_variable_rank",
                "highly_variable_nbatches",
                "means",
                "variances",
                "variances_norm",
            )
            if column in adata.var
        ]
        hvg_table = adata.var[hvg_columns].copy()
        hvg_table.insert(0, "gene_id", adata.var_names.astype(str))
        with self.progress.track(
            "analysis.pseudobulk",
            "PseudobulkBuilder.build",
            {
                "input_representation": "approved Step 02 raw integer adata.X",
                "cells": adata.n_obs,
                "genes": int(hvg.sum()),
                "sample_field": self.settings.sample_field,
                "aggregation": "sum raw counts by sample; log1p CPM",
            },
        ) as event:
            pseudobulk = PseudobulkBuilder().build(
                adata, self.settings.sample_field, hvg
            )
            event.outputs.update(
                {"samples": pseudobulk.shape[0], "genes": pseudobulk.shape[1]}
            )

        with self.progress.track(
            "analysis.normalization",
            "scanpy.pp.normalize_total",
            {
                "input_representation": "approved Step 02 raw integer adata.X",
                "input_shape": list(adata.shape),
                "target_sum": self.settings.target_sum,
                "exclude_highly_expressed": False,
                "inplace": True,
            },
        ) as event:
            sc.pp.normalize_total(adata, target_sum=self.settings.target_sum)
            event.outputs["representation"] = "library-size-normalized adata.X"
        with self.progress.track(
            "analysis.normalization",
            "scanpy.pp.log1p",
            {
                "input_representation": "library-size-normalized adata.X",
                "input_shape": list(adata.shape),
                "base": "natural logarithm",
            },
        ) as event:
            sc.pp.log1p(adata)
            event.outputs["representation"] = "log1p-normalized adata.X"
        with self.progress.track(
            "analysis.scaling",
            "scanpy.pp.scale",
            {
                "input_representation": "log1p-normalized HVG subset",
                "input_shape": [adata.n_obs, int(hvg.sum())],
                "zero_center": False,
                "max_value": 10,
            },
        ) as event:
            diagnostic = adata[:, hvg].copy()
            sc.pp.scale(diagnostic, zero_center=False, max_value=10)
            event.outputs.update(
                {
                    "output_shape": list(diagnostic.shape),
                    "output_matrix_type": type(diagnostic.X).__name__,
                }
            )
        with self.progress.track(
            "analysis.pca",
            "scanpy.tl.pca",
            {
                "input_representation": "scaled log1p HVGs",
                "input_shape": list(diagnostic.shape),
                "n_comps": self.settings.pca_components,
                "zero_center": True,
                "svd_solver": "arpack",
                "random_state": self.settings.random_seed,
            },
        ) as event:
            sc.tl.pca(
                diagnostic,
                n_comps=self.settings.pca_components,
                zero_center=True,
                svd_solver="arpack",
                random_state=self.settings.random_seed,
            )
            event.outputs["X_pca_shape"] = list(diagnostic.obsm["X_pca"].shape)
        with self.progress.track(
            "analysis.neighbors",
            "scanpy.pp.neighbors",
            {
                "input_representation": "X_pca",
                "input_shape": list(diagnostic.obsm["X_pca"].shape),
                "n_neighbors": self.settings.n_neighbors,
                "n_pcs": self.settings.pca_components,
                "metric": self.settings.neighbor_metric,
                "random_state": self.settings.random_seed,
                "method": "umap",
                "n_jobs": self.settings.n_jobs,
                "batch_correction": "none",
            },
        ) as event:
            sc.pp.neighbors(
                diagnostic,
                n_neighbors=self.settings.n_neighbors,
                n_pcs=self.settings.pca_components,
                metric=self.settings.neighbor_metric,
                random_state=self.settings.random_seed,
                method="umap",
            )
            event.outputs.update(
                {
                    "distances_shape": list(diagnostic.obsp["distances"].shape),
                    "distances_nnz": diagnostic.obsp["distances"].nnz,
                    "connectivities_nnz": diagnostic.obsp["connectivities"].nnz,
                }
            )
        with self.progress.track(
            "analysis.umap",
            "scanpy.tl.umap",
            {
                "input_representation": "unintegrated neighbors graph",
                "cells": diagnostic.n_obs,
                "min_dist": self.settings.umap_min_dist,
                "spread": self.settings.umap_spread,
                "random_state": self.settings.random_seed,
            },
        ) as event:
            sc.tl.umap(
                diagnostic,
                min_dist=self.settings.umap_min_dist,
                spread=self.settings.umap_spread,
                random_state=self.settings.random_seed,
            )
            event.outputs["X_umap_shape"] = list(diagnostic.obsm["X_umap"].shape)
        with self.progress.track(
            "analysis.clustering",
            "scanpy.tl.leiden",
            {
                "input_representation": "unintegrated neighbors graph",
                "cells": diagnostic.n_obs,
                "resolution": self.settings.leiden_resolution,
                "random_state": self.settings.random_seed,
                "flavor": "igraph",
                "directed": False,
                "n_iterations": 2,
                "key_added": "step06_leiden",
            },
        ) as event:
            sc.tl.leiden(
                diagnostic,
                resolution=self.settings.leiden_resolution,
                random_state=self.settings.random_seed,
                flavor="igraph",
                directed=False,
                n_iterations=2,
                key_added="step06_leiden",
            )
            event.outputs["clusters"] = int(
                diagnostic.obs["step06_leiden"].nunique()
            )

        with self.progress.track(
            "analysis.marker_programs",
            "ProgramScorer.score",
            {
                "input_representation": "full log1p-normalized adata.X",
                "input_shape": list(adata.shape),
                "cluster_field": "step06_leiden",
                "marker_programs": {
                    name: list(genes) for name, genes in marker_programs().items()
                },
                "purpose": "descriptive provisional states only",
            },
        ) as event:
            scores, cluster_scores = ProgramScorer().score(
                adata, diagnostic.obs["step06_leiden"]
            )
            event.outputs.update(
                {
                    "cell_score_columns": list(scores.columns),
                    "cluster_score_rows": len(cluster_scores),
                }
            )
        state_map = cluster_scores.set_index("leiden")["provisional_state"]
        annotations = scores.copy()
        annotations["step06_leiden"] = diagnostic.obs["step06_leiden"].astype(str)
        annotations["step06_provisional_state"] = (
            annotations["step06_leiden"].map(state_map).astype(str)
        )
        with self.progress.track(
            "analysis.rendering_subset",
            "RenderingSampler.select",
            {
                "input_cells": adata.n_obs,
                "sample_field": self.settings.sample_field,
                "sample_count": int(adata.obs[self.settings.sample_field].nunique()),
                "render_max_cells": self.settings.render_max_cells,
                "random_seed": self.settings.random_seed,
                "numerical_analysis_downsampled": False,
            },
        ) as event:
            render = RenderingSampler(self.settings).select(adata.obs)
            event.outputs["rendering_cells"] = len(render)
        return Step06Artifacts(
            annotations=annotations,
            pca=np.asarray(diagnostic.obsm["X_pca"], dtype=np.float32),
            umap=np.asarray(diagnostic.obsm["X_umap"], dtype=np.float32),
            connectivities=diagnostic.obsp["connectivities"].tocsr(),
            distances=diagnostic.obsp["distances"].tocsr(),
            hvg_table=hvg_table,
            variance_ratio=np.asarray(diagnostic.uns["pca"]["variance_ratio"]),
            cluster_program_scores=cluster_scores,
            pseudobulk_logcpm=pseudobulk,
            rendering_cells=render,
        )
