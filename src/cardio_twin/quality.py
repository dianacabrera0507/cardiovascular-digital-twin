"""Signal-quality placeholders. Thresholds must be protocol-defined."""
from __future__ import annotations
import numpy as np


def summarize_signal(values: np.ndarray) -> dict[str, float]:
    values = np.asarray(values, dtype=float)
    finite = np.isfinite(values)
    return {
        "n_samples": int(values.size),
        "finite_fraction": float(finite.mean()) if values.size else 0.0,
        "mean": float(np.nanmean(values)) if finite.any() else float("nan"),
        "std": float(np.nanstd(values)) if finite.any() else float("nan"),
    }
