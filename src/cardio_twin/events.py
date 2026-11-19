"""Event construction with explicit threshold and minimum duration."""
from __future__ import annotations
import numpy as np


def threshold_episodes(time: np.ndarray, values: np.ndarray, threshold: float, minimum_duration: float):
    """Return half-open index intervals [start, end) below threshold.

    Threshold and duration are inputs because they must be pre-specified and justified.
    """
    time = np.asarray(time, dtype=float)
    values = np.asarray(values, dtype=float)
    if time.shape != values.shape or time.ndim != 1:
        raise ValueError("time and values must be equal one-dimensional arrays")
    below = np.isfinite(values) & (values < threshold)
    changes = np.diff(np.r_[False, below, False].astype(int))
    starts = np.flatnonzero(changes == 1)
    ends = np.flatnonzero(changes == -1)
    episodes = []
    for start, end in zip(starts, ends):
        last_idx = min(end, len(time) - 1)
        duration = time[last_idx] - time[start]
        if duration >= minimum_duration:
            episodes.append((int(start), int(end)))
    return episodes
