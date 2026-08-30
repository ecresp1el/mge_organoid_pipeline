"""Apply the approved 5-MAD union while preserving exact cell dispositions."""

from __future__ import annotations

from dataclasses import dataclass

import anndata as ad
import numpy as np
import pandas as pd

from .step02_models import Step02Settings


@dataclass(frozen=True)
class Step02FilterResult:
    """Collect the filtered object and complete cell-level filtering evidence."""

    filtered: ad.AnnData
    disposition: pd.DataFrame
    by_sample: pd.DataFrame
    by_design: pd.DataFrame
    reason_summary: pd.DataFrame


class ApprovedMADFilter:
    """Filter only the approved per-sample 5-MAD union of three criteria."""

    def __init__(self, settings: Step02Settings):
        """Store the selected stringency and grouping-field contract."""

        self.settings = settings
        suffix = f"{settings.selected_stringency}mad"
        self.flag_columns = (
            f"low_total_counts_{suffix}",
            f"low_n_genes_by_counts_{suffix}",
            f"high_pct_counts_mt_{suffix}",
        )

    def apply(self, adata: ad.AnnData, flags: pd.DataFrame) -> Step02FilterResult:
        """Return a raw-count cell subset and full-universe exclusion ledger."""

        aligned = flags.set_index("cell_id").loc[adata.obs_names]
        selected = self._boolean_flags(aligned)
        excluded = selected.any(axis=1)
        reasons = selected.apply(self._reason_string, axis=1)
        disposition = pd.DataFrame({
            "cell_id": adata.obs_names.astype(str),
            self.settings.sample_field: adata.obs[self.settings.sample_field].astype(str).to_numpy(),
            self.settings.design_field: adata.obs[self.settings.design_field].astype(str).to_numpy(),
            "total_counts": adata.obs["total_counts"].to_numpy(),
            "n_genes_by_counts": adata.obs["n_genes_by_counts"].to_numpy(),
            "pct_counts_mt": adata.obs["pct_counts_mt"].to_numpy(),
            self.flag_columns[0]: selected.iloc[:, 0].to_numpy(),
            self.flag_columns[1]: selected.iloc[:, 1].to_numpy(),
            self.flag_columns[2]: selected.iloc[:, 2].to_numpy(),
            "step02_excluded": excluded.to_numpy(),
            "step02_disposition": np.where(excluded.to_numpy(), "excluded", "retained"),
            "step02_exclusion_reasons": reasons.to_numpy(),
        })
        keep = ~excluded.to_numpy()
        filtered = adata[keep, :].copy()
        filtered.obs["step02_qc_pass"] = True
        filtered.obs["step02_qc_exclusion_reasons"] = ""
        primary = dict(filtered.uns.get("primary_processing", {}))
        primary["latest_step"] = "02_qc_filtering"
        primary["step02_qc_filtering"] = {
            "source": "approved Step 01a per-sample MAD candidate flags",
            "stringency_mad": self.settings.selected_stringency,
            "combination": "low_total_counts OR low_n_genes_by_counts OR high_pct_counts_mt",
            "flag_columns": list(self.flag_columns),
            "cells_before": adata.n_obs,
            "cells_removed": int(excluded.sum()),
            "cells_after": filtered.n_obs,
            "genes_removed": 0,
            "normalization": "none",
            "doublet_filtering": "none",
        }
        filtered.uns["primary_processing"] = primary
        return Step02FilterResult(
            filtered=filtered,
            disposition=disposition,
            by_sample=self._group_summary(disposition, self.settings.sample_field),
            by_design=self._group_summary(disposition, self.settings.design_field),
            reason_summary=self._reason_summary(disposition),
        )

    def _group_summary(self, disposition: pd.DataFrame, field: str) -> pd.DataFrame:
        """Report before, removed, and after counts for one registered grouping."""

        rows = []
        for group, frame in disposition.groupby(field, sort=False, observed=True):
            before = len(frame)
            removed = int(frame["step02_excluded"].sum())
            rows.append({
                "grouping": field,
                "group": group,
                "cells_before": before,
                "cells_removed": removed,
                "cells_after": before - removed,
                "pct_removed": 100.0 * removed / before,
                "low_total_counts_n": int(frame[self.flag_columns[0]].sum()),
                "low_n_genes_n": int(frame[self.flag_columns[1]].sum()),
                "high_pct_counts_mt_n": int(frame[self.flag_columns[2]].sum()),
                "multiple_reasons_n": int((frame[list(self.flag_columns)].sum(axis=1) > 1).sum()),
            })
        return pd.DataFrame(rows)

    def _reason_summary(self, disposition: pd.DataFrame) -> pd.DataFrame:
        """Count every exact disposition/reason combination globally and by sample."""

        rows = []
        for scope, frame in [("pooled", disposition), *[(sample, group) for sample, group in disposition.groupby(self.settings.sample_field, sort=False, observed=True)]]:
            counts = frame["step02_exclusion_reasons"].replace("", "retained").value_counts()
            for reason, count in counts.items():
                rows.append({"scope": scope, "reason_combination": reason, "n_cells": int(count), "pct_cells": 100.0 * int(count) / len(frame)})
        return pd.DataFrame(rows)

    def _reason_string(self, row: pd.Series) -> str:
        """Join every true selected criterion without discarding overlaps."""

        return ";".join(column for column in self.flag_columns if bool(row[column]))

    def _boolean_flags(self, aligned: pd.DataFrame) -> pd.DataFrame:
        """Parse only canonical boolean encodings and reject ambiguous flags."""

        selected = aligned.loc[:, list(self.flag_columns)].copy()
        for column in self.flag_columns:
            if pd.api.types.is_bool_dtype(selected[column]):
                continue
            normalized = selected[column].astype(str).str.strip().str.lower()
            unexpected = sorted(set(normalized) - {"true", "false"})
            if unexpected:
                raise ValueError(f"Unexpected boolean values in {column}: {unexpected}")
            selected[column] = normalized.map({"true": True, "false": False}).astype(bool)
        return selected
