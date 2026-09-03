"""Full-cell unintegrated preprocessing and descriptive state analysis."""

from __future__ import annotations

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse

from .step06_models import Step06Artifacts, Step06Settings


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

    def __init__(self, settings: Step06Settings):
        """Store immutable analysis settings."""

        self.settings = settings

    def run(self, adata) -> Step06Artifacts:
        """Compute HVGs, PCA, neighbors, UMAP, Leiden, states, and pseudobulk."""

        sc.settings.n_jobs = self.settings.n_jobs
        sc.pp.highly_variable_genes(
            adata,
            flavor=self.settings.hvg_flavor,
            n_top_genes=self.settings.n_top_genes,
            batch_key=self.settings.sample_field,
            subset=False,
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
        pseudobulk = PseudobulkBuilder().build(adata, self.settings.sample_field, hvg)

        sc.pp.normalize_total(adata, target_sum=self.settings.target_sum)
        sc.pp.log1p(adata)
        diagnostic = adata[:, hvg].copy()
        sc.pp.scale(diagnostic, zero_center=False, max_value=10)
        sc.tl.pca(
            diagnostic,
            n_comps=self.settings.pca_components,
            zero_center=True,
            svd_solver="arpack",
            random_state=self.settings.random_seed,
        )
        sc.pp.neighbors(
            diagnostic,
            n_neighbors=self.settings.n_neighbors,
            n_pcs=self.settings.pca_components,
            metric=self.settings.neighbor_metric,
            random_state=self.settings.random_seed,
            method="umap",
        )
        sc.tl.umap(
            diagnostic,
            min_dist=self.settings.umap_min_dist,
            spread=self.settings.umap_spread,
            random_state=self.settings.random_seed,
        )
        sc.tl.leiden(
            diagnostic,
            resolution=self.settings.leiden_resolution,
            random_state=self.settings.random_seed,
            flavor="igraph",
            directed=False,
            n_iterations=2,
            key_added="step06_leiden",
        )

        scores, cluster_scores = ProgramScorer().score(
            adata, diagnostic.obs["step06_leiden"]
        )
        state_map = cluster_scores.set_index("leiden")["provisional_state"]
        annotations = scores.copy()
        annotations["step06_leiden"] = diagnostic.obs["step06_leiden"].astype(str)
        annotations["step06_provisional_state"] = (
            annotations["step06_leiden"].map(state_map).astype(str)
        )
        render = RenderingSampler(self.settings).select(adata.obs)
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
