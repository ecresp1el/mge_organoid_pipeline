# Ventral MEA activity

- `spike_activity.py`: spike-train burst and activity metrics.
- `recording_summary.py`: recording-level endpoint summaries.
- `figures/`: prepared MEA PNG.

Dependencies: Python 3.10 or newer, NumPy, and pandas.

## Inputs

Spike times are supplied in seconds, with recording duration in seconds.
Endpoint tables contain recording and unit IDs, FS/RS labels, best-channel IDs,
and numerical activity metrics.

## Usage

```python
from spike_activity import BurstDetectionParameters, summarize_unit_activity
from recording_summary import summarize_recordings

summary, bursts = summarize_unit_activity(
    spike_times_s, recording_duration_s, BurstDetectionParameters()
)
recording_metrics = summarize_recordings(unit_metrics, metric_columns)
```

Activity endpoints include firing rate, burst rate, firing rate within bursts,
burst duration, interburst interval, and spikes per burst. Recording summaries
average unit endpoints within channels and then across channels for each
recording/class. Missing endpoint values are omitted from the relevant mean.
