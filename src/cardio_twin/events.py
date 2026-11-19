"""Hypotension event construction with explicit temporal rules."""
from __future__ import annotations

import numpy as np


def _validate_series(time: np.ndarray, values: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    time_array = np.asarray(time, dtype=float)
    value_array = np.asarray(values, dtype=float)
    if time_array.ndim != 1 or value_array.ndim != 1:
        raise ValueError("time and values must be one-dimensional")
    if time_array.size == 0 or time_array.shape != value_array.shape:
        raise ValueError("time and values must be non-empty and equal length")
    if not np.isfinite(time_array).all():
        raise ValueError("time must contain only finite values")
    if time_array.size > 1 and np.any(np.diff(time_array) <= 0):
        raise ValueError("time must be strictly increasing")
    return time_array, value_array


def _sample_interval_seconds(time: np.ndarray) -> float:
    if time.size < 2:
        return 0.0
    return float(np.median(np.diff(time)))


def threshold_episodes(
    time: np.ndarray,
    values: np.ndarray,
    threshold: float,
    minimum_duration: float,
    merge_gap: float = 0.0,
) -> list[tuple[int, int]]:
    """Return half-open intervals ``[start, end)`` below ``threshold``.

    Durations and merge gaps use the same time unit as ``time``. Non-finite
    values are not hypotensive and interrupt an episode. Candidate runs are
    merged when the elapsed gap between runs is no greater than ``merge_gap``.
    The minimum-duration rule is applied after merging.
    """
    time_array, value_array = _validate_series(time, values)
    if not np.isfinite(threshold):
        raise ValueError("threshold must be finite")
    if minimum_duration < 0 or not np.isfinite(minimum_duration):
        raise ValueError("minimum_duration must be finite and non-negative")
    if merge_gap < 0 or not np.isfinite(merge_gap):
        raise ValueError("merge_gap must be finite and non-negative")

    below = np.isfinite(value_array) & (value_array < threshold)
    changes = np.diff(np.r_[False, below, False].astype(int))
    starts = np.flatnonzero(changes == 1)
    ends = np.flatnonzero(changes == -1)
    if starts.size == 0:
        return []

    merged: list[list[int]] = [[int(starts[0]), int(ends[0])]]
    sample_interval = _sample_interval_seconds(time_array)
    for start, end in zip(starts[1:], ends[1:]):
        previous_end = merged[-1][1]
        gap = float(time_array[int(start)] - time_array[previous_end - 1] - sample_interval)
        if gap <= merge_gap + np.finfo(float).eps:
            merged[-1][1] = int(end)
        else:
            merged.append([int(start), int(end)])

    episodes: list[tuple[int, int]] = []
    for start, end in merged:
        duration = float(time_array[end - 1] - time_array[start] + sample_interval)
        if duration >= minimum_duration:
            episodes.append((start, end))
    return episodes


def prediction_intervals(
    event_start_seconds: float,
    predictor_window_seconds: float,
    exclusion_gap_seconds: float,
    prediction_horizon_seconds: float,
) -> dict[str, tuple[float, float]]:
    """Construct leakage-safe predictor, gap, and outcome intervals.

    The outcome interval ends at the event onset and spans the prespecified
    prediction horizon. The predictor interval ends before the exclusion gap.
    All intervals are half-open.
    """
    values = np.asarray(
        [event_start_seconds, predictor_window_seconds, exclusion_gap_seconds,
         prediction_horizon_seconds],
        dtype=float,
    )
    if not np.isfinite(values).all():
        raise ValueError("all timing values must be finite")
    if event_start_seconds < 0:
        raise ValueError("event_start_seconds must be non-negative")
    if predictor_window_seconds <= 0:
        raise ValueError("predictor_window_seconds must be positive")
    if exclusion_gap_seconds < 0:
        raise ValueError("exclusion_gap_seconds must be non-negative")
    if prediction_horizon_seconds <= 0:
        raise ValueError("prediction_horizon_seconds must be positive")

    outcome_start = event_start_seconds - prediction_horizon_seconds
    predictor_end = outcome_start - exclusion_gap_seconds
    predictor_start = predictor_end - predictor_window_seconds
    if predictor_start < 0:
        raise ValueError("insufficient history before the event")
    return {
        "predictor": (predictor_start, predictor_end),
        "exclusion_gap": (predictor_end, outcome_start),
        "outcome_horizon": (outcome_start, event_start_seconds),
    }
