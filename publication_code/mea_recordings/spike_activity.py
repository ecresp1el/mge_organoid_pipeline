"""Spike-train burst detection and activity summaries."""
from __future__ import annotations
from dataclasses import dataclass, asdict
import numpy as np

@dataclass(frozen=True)
class BurstDetectionParameters:
    """Parameters for contiguous-ISI burst detection."""

    max_isi_ms: float = 100.0
    min_spikes: int = 3
    min_duration_ms: float = 100.0

    def validate(self) -> None:
        """Check that ISI, minimum spike count, and duration thresholds are valid."""
        if not np.isfinite(self.max_isi_ms) or self.max_isi_ms <= 0:
            raise ValueError("max_isi_ms must be finite and > 0")
        if self.min_spikes < 2:
            raise ValueError("min_spikes must be >= 2")
        if not np.isfinite(self.min_duration_ms) or self.min_duration_ms < 0:
            raise ValueError("min_duration_ms must be finite and >= 0")

@dataclass(frozen=True)
class BurstEvent:
    """One accepted burst in one sorted unit."""

    burst_index: int
    first_spike_index: int
    last_spike_index: int
    start_time_s: float
    end_time_s: float
    duration_ms: float
    spike_count: int
    firing_rate_hz: float
    interval_rate_hz: float

    def to_dict(self) -> dict[str, int | float]:
        """Return burst boundaries, spike counts, duration, and rate as plain fields."""
        return asdict(self)



def detect_bursts(
    spike_times_s: np.ndarray,
    parameters: BurstDetectionParameters,
) -> list[BurstEvent]:
    """Detect bursts as contiguous runs whose adjacent ISIs do not exceed a limit."""

    parameters.validate()
    spikes = _validated_spike_times(spike_times_s)
    if spikes.size < parameters.min_spikes:
        return []

    max_isi_s = parameters.max_isi_ms / 1000.0
    min_duration_s = parameters.min_duration_ms / 1000.0
    # A small numerical tolerance keeps decimal interval boundaries stable.
    tolerance_s = max(np.finfo(float).eps * 16.0, 1e-12)
    split_after = np.flatnonzero(np.diff(spikes) > max_isi_s + tolerance_s)
    starts = np.concatenate(([0], split_after + 1))
    stops = np.concatenate((split_after + 1, [spikes.size]))

    events: list[BurstEvent] = []
    for start, stop in zip(starts, stops, strict=True):
        spike_count = int(stop - start)
        if spike_count < parameters.min_spikes:
            continue
        duration_s = float(spikes[stop - 1] - spikes[start])
        if duration_s + tolerance_s < min_duration_s:
            continue
        firing_rate_hz = float(spike_count / duration_s) if duration_s > 0 else np.nan
        interval_rate_hz = float((spike_count - 1) / duration_s) if duration_s > 0 else np.nan
        events.append(
            BurstEvent(
                burst_index=len(events),
                first_spike_index=int(start),
                last_spike_index=int(stop - 1),
                start_time_s=float(spikes[start]),
                end_time_s=float(spikes[stop - 1]),
                duration_ms=duration_s * 1000.0,
                spike_count=spike_count,
                firing_rate_hz=firing_rate_hz,
                interval_rate_hz=interval_rate_hz,
            )
        )
    return events


def summarize_unit_activity(
    spike_times_s: np.ndarray,
    recording_duration_s: float,
    parameters: BurstDetectionParameters,
) -> tuple[dict[str, int | float | bool], list[BurstEvent]]:
    """Return one unit summary and its accepted burst events."""

    spikes = _validated_spike_times(spike_times_s)
    if not np.isfinite(recording_duration_s) or recording_duration_s <= 0:
        raise ValueError("recording_duration_s must be finite and > 0")
    events = detect_bursts(spikes, parameters)
    burst_spike_count = int(sum(event.spike_count for event in events))
    if burst_spike_count > spikes.size:
        raise AssertionError("burst spike assignment exceeded total unit spikes")

    ibis_s = np.asarray(
        [events[index].start_time_s - events[index - 1].end_time_s for index in range(1, len(events))],
        dtype=float,
    )
    durations_ms = np.asarray([event.duration_ms for event in events], dtype=float)
    spikes_per_burst = np.asarray([event.spike_count for event in events], dtype=float)
    burst_rates_hz = np.asarray([event.firing_rate_hz for event in events], dtype=float)
    interval_rates_hz = np.asarray([event.interval_rate_hz for event in events], dtype=float)
    spike_count = int(spikes.size)

    summary: dict[str, int | float | bool] = {
        "spike_count_recomputed": spike_count,
        "recording_duration_s_recomputed": float(recording_duration_s),
        "firing_rate_hz_recomputed": float(spike_count / recording_duration_s),
        "burst_count": int(len(events)),
        "burst_positive": bool(events),
        "burst_rate_per_min": float(len(events) / (recording_duration_s / 60.0)),
        "burst_spike_count": burst_spike_count,
        "fraction_spikes_in_bursts": float(burst_spike_count / spike_count) if spike_count else 0.0,
        "mean_firing_rate_within_bursts_hz": _mean_or_nan(burst_rates_hz),
        "mean_interval_rate_within_bursts_hz": _mean_or_nan(interval_rates_hz),
        "mean_burst_duration_ms": _mean_or_nan(durations_ms),
        "median_burst_duration_ms": _median_or_nan(durations_ms),
        "mean_spikes_per_burst": _mean_or_nan(spikes_per_burst),
        "median_spikes_per_burst": _median_or_nan(spikes_per_burst),
        "max_spikes_per_burst": float(np.max(spikes_per_burst)) if spikes_per_burst.size else np.nan,
        "interburst_interval_count": int(ibis_s.size),
        "mean_interburst_interval_s": _mean_or_nan(ibis_s),
        "median_interburst_interval_s": _median_or_nan(ibis_s),
    }
    return summary, events


def _validated_spike_times(spike_times_s: np.ndarray) -> np.ndarray:
    """Validate finite, sorted, one-dimensional spike times.
    
    Callers apply endpoint-specific recording-duration checks."""
    spikes = np.asarray(spike_times_s, dtype=float)
    if spikes.ndim != 1:
        raise ValueError("spike_times_s must be one-dimensional")
    if not np.isfinite(spikes).all():
        raise ValueError("spike_times_s contains non-finite values")
    if spikes.size > 1 and np.any(np.diff(spikes) < 0):
        raise ValueError("spike_times_s must be sorted")
    return spikes


def _mean_or_nan(values: np.ndarray) -> float:
    """Return the arithmetic mean, or NaN for an empty array."""
    return float(np.mean(values)) if values.size else np.nan


def _median_or_nan(values: np.ndarray) -> float:
    """Return the median, or NaN for an empty array."""
    return float(np.median(values)) if values.size else np.nan
