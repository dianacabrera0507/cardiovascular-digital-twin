"""Minimal two-element Windkessel utilities for synthetic experiments."""
from __future__ import annotations
import numpy as np


def simulate_two_element(
    time: np.ndarray,
    inflow: np.ndarray,
    resistance: float,
    compliance: float,
    venous_pressure: float = 5.0,
    initial_pressure: float = 80.0,
) -> np.ndarray:
    """Integrate C dP/dt = Qin - (P-Pv)/R with explicit Euler.

    This implementation is a development baseline, not a clinical model.
    """
    time = np.asarray(time, dtype=float)
    inflow = np.asarray(inflow, dtype=float)
    if time.ndim != 1 or inflow.shape != time.shape:
        raise ValueError("time and inflow must be one-dimensional arrays of equal size")
    if len(time) < 2 or np.any(np.diff(time) <= 0):
        raise ValueError("time must be strictly increasing with at least two samples")
    if resistance <= 0 or compliance <= 0:
        raise ValueError("resistance and compliance must be positive")

    pressure = np.empty_like(time)
    pressure[0] = initial_pressure
    for idx in range(1, len(time)):
        dt = time[idx] - time[idx - 1]
        derivative = (
            inflow[idx - 1]
            - (pressure[idx - 1] - venous_pressure) / resistance
        ) / compliance
        pressure[idx] = pressure[idx - 1] + dt * derivative
    return pressure


def diastolic_time_constant(resistance: float, compliance: float) -> float:
    """Return tau = R*C; does not imply separate identifiability of R and C."""
    if resistance <= 0 or compliance <= 0:
        raise ValueError("resistance and compliance must be positive")
    return resistance * compliance
