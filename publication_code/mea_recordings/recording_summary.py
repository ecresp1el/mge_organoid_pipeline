"""Recording-level endpoint summaries."""
from __future__ import annotations
import pandas as pd


def summarize_recordings(unit_metrics: pd.DataFrame, metric_columns: list[str],
                         channel_column="best_channel_id") -> pd.DataFrame:
    """Compute channel-averaged class summaries within each recording.

    Inputs: unit metrics and numerical endpoint column names.
    Units are averaged within best channel/class first, then channels are averaged
    within recording/class. Returns endpoint means, unit counts, and channel counts.
    """
    keys = ["recording", "rs_fs_class"]
    channels = unit_metrics.groupby(keys + [channel_column], sort=False, dropna=False)[metric_columns].mean()
    means = channels.groupby(level=keys, sort=False, dropna=False).mean()
    means["unit_count"] = unit_metrics.groupby(keys, sort=False, dropna=False).size()
    means["best_channel_count"] = channels.groupby(level=keys, sort=False, dropna=False).size()
    return means.reset_index()
