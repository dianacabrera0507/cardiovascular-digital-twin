"""Synthetic signals for method development and parameter-recovery tests."""
from __future__ import annotations
import numpy as np
from .models.windkessel import simulate_two_element


def pulsatile_inflow(time: np.ndarray, heart_rate_bpm: float, amplitude: float) -> np.ndarray:
    if heart_rate_bpm <= 0 or amplitude < 0:
        raise ValueError("heart_rate_bpm must be positive and amplitude nonnegative")
    phase = (time * heart_rate_bpm / 60.0) % 1.0
    pulse = np.where(phase < 0.30, np.sin(np.pi * phase / 0.30), 0.0)
    return amplitude * pulse


def make_synthetic_record(config: dict, seed: int = 2026) -> dict[str, np.ndarray]:
    duration = float(config["duration_seconds"])
    fs = float(config["sampling_hz"])
    time = np.arange(0.0, duration, 1.0 / fs)
    inflow = pulsatile_inflow(time, float(config["heart_rate_bpm"]), float(config["inflow_amplitude"]))
    clean = simulate_two_element(
        time,
        inflow,
        float(config["resistance"]),
        float(config["compliance"]),
        float(config["venous_pressure"]),
    )
    rng = np.random.default_rng(seed)
    observed = clean + rng.normal(0.0, float(config["noise_std"]), size=clean.shape)
    missing_fraction = float(config.get("missing_fraction", 0.0))
    if not 0.0 <= missing_fraction < 1.0:
        raise ValueError("missing_fraction must be in [0, 1)")
    mask = rng.random(len(time)) < missing_fraction
    observed[mask] = np.nan
    return {"time": time, "inflow": inflow, "pressure_clean": clean, "pressure_observed": observed}
