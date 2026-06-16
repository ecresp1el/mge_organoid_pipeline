#!/usr/bin/env python3
"""Assign DIV90 cells to the neuron-S9 seven-tip URD grouping."""

from __future__ import annotations

import argparse
import os
from pathlib import Path

import anndata as ad
import matplotlib

matplotlib.use("Agg")

import pandas as pd
import scanpy as sc


PROJECT_ROOT_DEFAULT = "/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder"
NEURON_PARENT_CLUSTERS = ["0", "1", "2", "3", "5", "8", "11"]
NEURON_KEY = "neuron_leiden_r0_6"

NEURON_TIP_BY_CLUSTER = {
    "1": "tip_lhx6_nfia_epha5_mef2c_cortical",
    "2": "tip_lhx6_nfia_epha5_mef2c_cortical",
    "3": "tip_crabp1_angpt2_fetal_precursor",
    "4": "tip_lhx8_isl1_state1",
    "5": "tip_nr2f1_nr2f2",
    "6": "tip_nr2f1_nr2f2",
    "7": "tip_lhx8_isl1_state2",
    "8": "tip_nr2f1_nr2f2",
}

TIP_BY_ORIGINAL_CLUSTER = {
    "4": "tip_astrocytes",
    "10": "tip_astrocytes",
    "9": "tip_opc",
}

TIP_DISPLAY = {
    "tip_lhx8_isl1_state1": "LHX8/ISL1-like ventral MGE branch state 1 (neuron 4)",
    "tip_lhx8_isl1_state2": "LHX8/ISL1-like ventral MGE branch state 2 (neuron 7)",
    "tip_crabp1_angpt2_fetal_precursor": "CRABP1/ANGPT2-like fetal precursor branch (neuron 3)",
    "tip_lhx6_nfia_epha5_mef2c_cortical": "LHX6/NFIA + EPHA5/MEF2C cortical interneuron branch (neurons 1+2)",
    "tip_nr2f1_nr2f2": "NR2F1/NR2F2-like branch (neurons 5+6+8)",
    "tip_astrocytes": "combined astrocyte endpoint (original clusters 4+10)",
    "tip_opc": "OPC endpoint (original cluster 9)",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", default=os.environ.get("PROJECT_ROOT", PROJECT_ROOT_DEFAULT))
    parser.add_argument("--h5ad", default=None)
    parser.add_argument("--out", required=True)
    parser.add_argument("--n-top-genes", type=int, default=2500)
    parser.add_argument("--n-pcs", type=int, default=40)
    parser.add_argument("--n-neighbors", type=int, default=30)
    return parser.parse_args()


def run_neuron_reclustering(adata: ad.AnnData, args: argparse.Namespace) -> pd.Series:
    neuron = adata[adata.obs["cluster_id"].astype(str).isin(NEURON_PARENT_CLUSTERS)].copy()
    work = neuron.copy()
    sc.pp.highly_variable_genes(work, n_top_genes=min(args.n_top_genes, work.n_vars), flavor="seurat")
    work = work[:, work.var["highly_variable"].to_numpy()].copy()
    sc.pp.scale(work, max_value=10)
    sc.tl.pca(work, n_comps=min(args.n_pcs, work.n_obs - 2, work.n_vars - 1), svd_solver="arpack")
    sc.pp.neighbors(work, n_neighbors=min(args.n_neighbors, work.n_obs - 1), n_pcs=min(args.n_pcs, work.obsm["X_pca"].shape[1]))
    sc.tl.leiden(work, resolution=0.6, key_added=NEURON_KEY, random_state=0)
    return pd.Series(work.obs[NEURON_KEY].astype(str).to_numpy(), index=work.obs["cell_id"].astype(str).to_numpy(), name=NEURON_KEY)


def main() -> int:
    args = parse_args()
    project_root = Path(args.project_root).expanduser().resolve()
    h5ad = Path(args.h5ad).expanduser().resolve() if args.h5ad else project_root / "results/python_anndata/varela_div90.h5ad"
    out = Path(args.out).expanduser().resolve()
    out.parent.mkdir(parents=True, exist_ok=True)

    adata = sc.read_h5ad(h5ad)
    adata.obs["cluster_id"] = adata.obs["cluster_id"].astype(str)
    neuron_labels = run_neuron_reclustering(adata, args)

    df = adata.obs[["cell_id", "cluster_id", "cluster_number_name"]].copy()
    df["cell_id"] = df["cell_id"].astype(str)
    df["cluster_id"] = df["cluster_id"].astype(str)
    df[NEURON_KEY] = df["cell_id"].map(neuron_labels)
    df["div90_jia_tip_group"] = pd.NA
    df["div90_jia_tip_group"] = df[NEURON_KEY].map(NEURON_TIP_BY_CLUSTER)
    df.loc[df["cluster_id"].isin(TIP_BY_ORIGINAL_CLUSTER), "div90_jia_tip_group"] = df.loc[
        df["cluster_id"].isin(TIP_BY_ORIGINAL_CLUSTER), "cluster_id"
    ].map(TIP_BY_ORIGINAL_CLUSTER)
    df["div90_jia_tip_display"] = df["div90_jia_tip_group"].map(TIP_DISPLAY)
    df["div90_jia_urd_role"] = "retained_unassigned_or_excluded"
    df.loc[df["cluster_id"] == "12", "div90_jia_urd_role"] = "root_pool_cluster12_dividing_cells"
    df.loc[df["div90_jia_tip_group"].notna(), "div90_jia_urd_role"] = df.loc[df["div90_jia_tip_group"].notna(), "div90_jia_tip_group"]
    df.loc[df[NEURON_KEY] == "0", "div90_jia_urd_role"] = "retained_upstream_lhx8_isl1_neuron0"
    df.loc[df["cluster_id"].isin(["6", "7"]), "div90_jia_urd_role"] = "excluded_stressed"

    df.to_csv(out, sep="\t", index=False)
    summary = (
        df.groupby(["div90_jia_urd_role", "div90_jia_tip_group", "div90_jia_tip_display"], dropna=False)
        .size()
        .reset_index(name="n_cells")
    )
    summary.to_csv(out.with_name(out.stem + "_summary.tsv"), sep="\t", index=False)
    print(f"Wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
