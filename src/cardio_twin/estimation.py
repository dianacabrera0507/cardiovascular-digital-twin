"""Parameter estimation baselines."""
from __future__ import annotations
import numpy as np
from scipy.optimize import least_squares
from .models.windkessel import simulate_two_element


def estimate_resistance_compliance(
    time: np.ndarray,
    inflow: np.ndarray,
    pressure: np.ndarray,
    venous_pressure: float = 5.0,
    initial_pressure: float | None = None,
    initial_guess: tuple[float, float] = (1.0, 1.0),
) -> dict[str, float]:
    """Estimate positive R and C on observed finite samples.

    This baseline assumes inflow is observed; real-data use must justify that assumption.
    """
    pressure = np.asarray(pressure, dtype=float)
    finite = np.isfinite(pressure)
    if finite.sum() < 3:
        raise ValueError("at least three finite pressure samples are required")
    p0 = float(pressure[np.flatnonzero(finite)[0]] if initial_pressure is None else initial_pressure)

    def residual(log_params: np.ndarray) -> np.ndarray:
        resistance, compliance = np.exp(log_params)
        predicted = simulate_two_element(time, inflow, resistance, compliance, venous_pressure, p0)
        return predicted[finite] - pressure[finite]

    result = least_squares(residual, np.log(np.asarray(initial_guess, dtype=float)))
    resistance, compliance = np.exp(result.x)
    return {
        "resistance": float(resistance),
        "compliance": float(compliance),
        "time_constant": float(resistance * compliance),
        "cost": float(result.cost),
        "success": bool(result.success),
    }
