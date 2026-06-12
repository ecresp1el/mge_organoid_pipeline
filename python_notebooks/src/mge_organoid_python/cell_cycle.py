"""Cell-cycle scoring helpers for Notebook 01.

The S and G2M gene lists are the standard human Regev/Tirosh lists used by the
Scanpy and Seurat cell-cycle scoring examples. Notebook 01 filters these lists
to genes present in the current AnnData object before scoring.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

import anndata as ad
import pandas as pd
import scanpy as sc


CELL_CYCLE_GENE_SOURCE = "regev_lab_cell_cycle_genes_from_scanpy_seurat_tutorial"
CCDIFFERENCE_KEY = "CCDifference"

S_PHASE_GENES = (
    "MCM5",
    "PCNA",
    "TYMS",
    "FEN1",
    "MCM2",
    "MCM4",
    "RRM1",
    "UNG",
    "GINS2",
    "MCM6",
    "CDCA7",
    "DTL",
    "PRIM1",
    "UHRF1",
    "MLF1IP",
    "HELLS",
    "RFC2",
    "RPA2",
    "NASP",
    "RAD51AP1",
    "GMNN",
    "WDR76",
    "SLBP",
    "CCNE2",
    "UBR7",
    "POLD3",
    "MSH2",
    "ATAD2",
    "RAD51",
    "RRM2",
    "CDC45",
    "CDC6",
    "EXO1",
    "TIPIN",
    "DSCC1",
    "BLM",
    "CASP8AP2",
    "USP1",
    "CLSPN",
    "POLA1",
    "CHAF1B",
    "BRIP1",
    "E2F8",
)

G2M_PHASE_GENES = (
    "HMGB2",
    "CDK1",
    "NUSAP1",
    "UBE2C",
    "BIRC5",
    "TPX2",
    "TOP2A",
    "NDC80",
    "CKS2",
    "NUF2",
    "CKS1B",
    "MKI67",
    "TMPO",
    "CENPF",
    "TACC3",
    "FAM64A",
    "SMC4",
    "CCNB2",
    "CKAP2L",
    "CKAP2",
    "AURKB",
    "BUB1",
    "KIF11",
    "ANP32E",
    "TUBB4B",
    "GTSE1",
    "KIF20B",
    "HJURP",
    "CDCA3",
    "HN1",
    "CDC20",
    "TTK",
    "CDC25C",
    "KIF2C",
    "RANGAP1",
    "NCAPD2",
    "DLGAP5",
    "CDCA2",
    "CDCA8",
    "ECT2",
    "KIF23",
    "HMMR",
    "AURKA",
    "PSRC1",
    "ANLN",
    "LBR",
    "CKAP5",
    "CENPE",
    "CTCF",
    "NEK2",
    "G2E3",
    "GAS2L3",
    "CBX5",
    "CENPA",
)


@dataclass(frozen=True)
class CellCycleGeneSelection:
    """Cell-cycle gene lists filtered against an AnnData object's genes."""

    s_genes: tuple[str, ...]
    g2m_genes: tuple[str, ...]
    cell_cycle_genes: tuple[str, ...]
    gene_table: pd.DataFrame
    summary_table: pd.DataFrame


def _present_genes(source_genes: Sequence[str], var_names: set[str]) -> tuple[str, ...]:
    """Return source genes present in `var_names`, preserving source order."""
    return tuple(gene for gene in source_genes if gene in var_names)


def cell_cycle_gene_table(
    var_names: Sequence[str],
    *,
    s_genes: Sequence[str] = S_PHASE_GENES,
    g2m_genes: Sequence[str] = G2M_PHASE_GENES,
    source: str = CELL_CYCLE_GENE_SOURCE,
) -> pd.DataFrame:
    """Return one row per source cell-cycle gene with present/missing status."""
    var_name_set = set(map(str, var_names))
    records = []
    for phase, genes in [("S", s_genes), ("G2M", g2m_genes)]:
        for order, gene in enumerate(genes):
            records.append(
                {
                    "phase": phase,
                    "gene": gene,
                    "source_order": order,
                    "present_in_adata": gene in var_name_set,
                    "source": source,
                }
            )
    return pd.DataFrame(records)


def cell_cycle_gene_summary(gene_table: pd.DataFrame) -> pd.DataFrame:
    """Summarize source/present/missing cell-cycle genes by phase and total."""
    records = []
    for phase, frame in gene_table.groupby("phase", sort=False):
        n_source = int(len(frame))
        n_present = int(frame["present_in_adata"].sum())
        records.append(
            {
                "phase": phase,
                "n_source_genes": n_source,
                "n_present_genes": n_present,
                "n_missing_genes": n_source - n_present,
            }
        )
    n_source_total = int(len(gene_table))
    n_present_total = int(gene_table["present_in_adata"].sum())
    records.append(
        {
            "phase": "all",
            "n_source_genes": n_source_total,
            "n_present_genes": n_present_total,
            "n_missing_genes": n_source_total - n_present_total,
        }
    )
    return pd.DataFrame(records)


def select_cell_cycle_genes(var_names: Sequence[str]) -> CellCycleGeneSelection:
    """Filter source S/G2M genes to those present in `var_names`."""
    var_name_set = set(map(str, var_names))
    s_present = _present_genes(S_PHASE_GENES, var_name_set)
    g2m_present = _present_genes(G2M_PHASE_GENES, var_name_set)
    gene_table = cell_cycle_gene_table(var_names)
    summary_table = cell_cycle_gene_summary(gene_table)
    return CellCycleGeneSelection(
        s_genes=s_present,
        g2m_genes=g2m_present,
        cell_cycle_genes=(*s_present, *g2m_present),
        gene_table=gene_table,
        summary_table=summary_table,
    )


def score_cell_cycle_and_ccdifference(
    adata: ad.AnnData,
    *,
    selection: CellCycleGeneSelection | None = None,
    ccdifference_key: str = CCDIFFERENCE_KEY,
) -> CellCycleGeneSelection:
    """Score S/G2M phase and create `CCDifference = S_score - G2M_score`.

    The input `adata.X` should contain normalized/log1p expression. The scoring
    function mutates `adata.obs` by adding `S_score`, `G2M_score`, `phase`, and
    the requested CCDifference column.
    """
    selected = selection or select_cell_cycle_genes(adata.var_names)
    if not selected.s_genes:
        raise ValueError("No S-phase genes are present in adata.var_names.")
    if not selected.g2m_genes:
        raise ValueError("No G2M-phase genes are present in adata.var_names.")

    sc.tl.score_genes_cell_cycle(
        adata,
        s_genes=list(selected.s_genes),
        g2m_genes=list(selected.g2m_genes),
    )
    adata.obs[ccdifference_key] = adata.obs["S_score"] - adata.obs["G2M_score"]
    adata.uns["notebook01_cell_cycle"] = {
        "source": CELL_CYCLE_GENE_SOURCE,
        "s_genes": list(selected.s_genes),
        "g2m_genes": list(selected.g2m_genes),
        "cell_cycle_genes": list(selected.cell_cycle_genes),
        "ccdifference_key": ccdifference_key,
        "ccdifference_formula": "S_score - G2M_score",
    }
    return selected


def cell_cycle_score_summary(
    adata: ad.AnnData,
    *,
    scope: str,
    run_sample_id: object | None = None,
    ccdifference_key: str = CCDIFFERENCE_KEY,
) -> pd.DataFrame:
    """Summarize phase counts and score distributions for one scored object."""
    records = []
    score_columns = ["S_score", "G2M_score", ccdifference_key]
    for phase, frame in adata.obs.groupby("phase", observed=True, sort=False):
        record = {
            "scope": scope,
            "run_sample_id": "" if run_sample_id is None else str(run_sample_id),
            "phase": str(phase),
            "n_cells": int(len(frame)),
        }
        for column in score_columns:
            record[f"{column}_mean"] = float(frame[column].mean())
            record[f"{column}_median"] = float(frame[column].median())
            record[f"{column}_std"] = float(frame[column].std())
            record[f"{column}_min"] = float(frame[column].min())
            record[f"{column}_max"] = float(frame[column].max())
        records.append(record)
    return pd.DataFrame(records)
