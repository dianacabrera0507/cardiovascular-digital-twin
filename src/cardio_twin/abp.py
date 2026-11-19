"""Herramientas reutilizables para analizar presión arterial invasiva (ABP).

Estas funciones permiten resumir la calidad de una ventana, detectar picos
sistólicos, caracterizar latidos y estimar una constante de tiempo diastólica
efectiva. No identifican resistencia y compliancia por separado desde ABP.
"""

from __future__ import annotations

import numpy as np
from scipy.optimize import least_squares
from scipy.signal import find_peaks, savgol_filter


def _as_1d(values: np.ndarray, name: str) -> np.ndarray:
    array = np.asarray(values, dtype=float)

    if array.ndim != 1 or array.size == 0:
        raise ValueError(
            f"{name} must be a non-empty one-dimensional array"
        )

    return array


def summarize_abp_window(
    values: np.ndarray,
) -> dict[str, float]:
    """Devuelve descriptores robustos de una ventana ABP en mmHg."""

    values = _as_1d(values, "values")
    finite = values[np.isfinite(values)]

    if finite.size == 0:
        return {
            "n_samples": int(values.size),
            "finite_fraction": 0.0,
            "median_mmhg": float("nan"),
            "p05_mmhg": float("nan"),
            "p95_mmhg": float("nan"),
            "robust_amplitude_mmhg": float("nan"),
        }

    p05, median, p95 = np.percentile(
        finite,
        [5, 50, 95],
    )

    return {
        "n_samples": int(values.size),
        "finite_fraction": float(
            finite.size / values.size
        ),
        "median_mmhg": float(median),
        "p05_mmhg": float(p05),
        "p95_mmhg": float(p95),
        "robust_amplitude_mmhg": float(p95 - p05),
    }


def detect_systolic_peaks(
    values: np.ndarray,
    sampling_hz: float,
    maximum_heart_rate_bpm: float = 180.0,
    minimum_prominence_mmhg: float = 5.0,
) -> np.ndarray:
    """Localiza picos sistólicos candidatos en una ventana ABP."""

    values = _as_1d(values, "values")

    if not np.isfinite(values).all():
        raise ValueError(
            "values must be finite for peak detection"
        )

    if sampling_hz <= 0:
        raise ValueError(
            "sampling_hz must be positive"
        )

    if maximum_heart_rate_bpm <= 0:
        raise ValueError(
            "maximum_heart_rate_bpm must be positive"
        )

    smooth_samples = max(
        5,
        int(round(0.11 * sampling_hz)),
    )

    if smooth_samples % 2 == 0:
        smooth_samples += 1

    if smooth_samples >= values.size:
        smooth_samples = (
            values.size - 1
            if values.size % 2 == 0
            else values.size
        )

    if smooth_samples < 5:
        raise ValueError(
            "values contain too few samples "
            "for peak detection"
        )

    smoothed = savgol_filter(
        values,
        window_length=smooth_samples,
        polyorder=2,
    )

    robust_amplitude = float(
        np.percentile(values, 95)
        - np.percentile(values, 5)
    )

    prominence = max(
        minimum_prominence_mmhg,
        0.15 * robust_amplitude,
    )

    minimum_distance = max(
        1,
        int(
            round(
                60.0
                * sampling_hz
                / maximum_heart_rate_bpm
            )
        ),
    )

    peaks, _ = find_peaks(
        smoothed,
        distance=minimum_distance,
        prominence=prominence,
    )

    return peaks.astype(int)


def characterize_beats(
    values: np.ndarray,
    systolic_peaks: np.ndarray,
    sampling_hz: float,
) -> list[dict[str, float | int | bool]]:
    """Caracteriza intervalos entre picos sistólicos consecutivos."""

    values = _as_1d(values, "values")
    peaks = np.asarray(systolic_peaks, dtype=int)

    if peaks.ndim != 1:
        raise ValueError(
            "systolic_peaks must be one-dimensional"
        )

    if peaks.size > 1 and np.any(np.diff(peaks) <= 0):
        raise ValueError(
            "systolic_peaks must be strictly increasing"
        )

    if sampling_hz <= 0:
        raise ValueError(
            "sampling_hz must be positive"
        )

    if peaks.size and (
        peaks[0] < 0
        or peaks[-1] >= values.size
    ):
        raise ValueError(
            "systolic_peaks contain an out-of-range index"
        )

    rows: list[
        dict[str, float | int | bool]
    ] = []

    for beat_number, (start, stop) in enumerate(
        zip(peaks[:-1], peaks[1:]),
        start=1,
    ):
        segment = values[start:stop]

        if (
            segment.size < 3
            or not np.isfinite(segment).all()
        ):
            continue

        diastolic_index = int(
            start + np.argmin(segment)
        )

        duration = float(
            (stop - start) / sampling_hz
        )

        systolic = float(values[start])
        diastolic = float(
            values[diastolic_index]
        )

        pulse_pressure = systolic - diastolic

        plausible = bool(
            0.35 <= duration <= 1.50
            and 40.0 <= systolic <= 250.0
            and 20.0 <= diastolic <= 160.0
            and 10.0 <= pulse_pressure <= 120.0
        )

        rows.append(
            {
                "beat": beat_number,
                "start_peak": int(start),
                "diastolic_index": diastolic_index,
                "duration_s": duration,
                "heart_rate_bpm": float(
                    60.0 / duration
                ),
                "systolic_mmhg": systolic,
                "diastolic_mmhg": diastolic,
                "pulse_pressure_mmhg": float(
                    pulse_pressure
                ),
                "plausible": plausible,
            }
        )

    return rows


def estimate_effective_diastolic_tau(
    values: np.ndarray,
    sampling_hz: float,
    start_peak: int,
    diastolic_index: int,
    venous_pressure_mmhg: float = 5.0,
    fit_start_fraction: float = 0.35,
) -> dict[str, float | int | bool]:
    """Ajusta una constante de tiempo diastólica efectiva.

    Modelo:
        P(t) = Pv + (P0 - Pv) * exp(-t / tau)

    El resultado depende del modelo y de sus supuestos.
    No representa R y C identificados por separado.
    """

    values = _as_1d(values, "values")

    if sampling_hz <= 0:
        raise ValueError(
            "sampling_hz must be positive"
        )

    if not 0.0 <= fit_start_fraction < 1.0:
        raise ValueError(
            "fit_start_fraction must be in [0, 1)"
        )

    if not (
        0 <= start_peak
        < diastolic_index
        < values.size
    ):
        raise ValueError(
            "invalid beat indices"
        )

    beat_length = (
        diastolic_index - start_peak
    )

    fit_start = start_peak + max(
        2,
        int(
            round(
                fit_start_fraction
                * beat_length
            )
        ),
    )

    pressure = values[
        fit_start:diastolic_index + 1
    ]

    finite = (
        np.isfinite(pressure)
        & (pressure > venous_pressure_mmhg)
    )

    pressure = pressure[finite]

    time = (
        np.arange(
            diastolic_index + 1 - fit_start,
            dtype=float,
        )[finite]
        / sampling_hz
    )

    if pressure.size < 8:
        raise ValueError(
            "at least eight valid diastolic "
            "samples are required"
        )

    initial_pressure = float(pressure[0])

    def residual(
        log_tau: np.ndarray,
    ) -> np.ndarray:
        tau = float(np.exp(log_tau[0]))

        predicted = (
            venous_pressure_mmhg
            + (
                initial_pressure
                - venous_pressure_mmhg
            )
            * np.exp(-time / tau)
        )

        return predicted - pressure

    result = least_squares(
        residual,
        x0=np.log([1.5]),
        bounds=(
            np.log([0.10]),
            np.log([10.0]),
        ),
    )

    tau = float(np.exp(result.x[0]))

    predicted = (
        venous_pressure_mmhg
        + (
            initial_pressure
            - venous_pressure_mmhg
        )
        * np.exp(-time / tau)
    )

    errors = pressure - predicted

    residual_sum_squares = float(
        np.sum(errors ** 2)
    )

    total_sum_squares = float(
        np.sum(
            (
                pressure
                - np.mean(pressure)
            ) ** 2
        )
    )

    r_squared = (
        float(
            1.0
            - residual_sum_squares
            / total_sum_squares
        )
        if total_sum_squares > 0
        else float("nan")
    )

    rmse = float(
        np.sqrt(
            np.mean(errors ** 2)
        )
    )

    acceptable = bool(
        result.success
        and 0.20 <= tau <= 5.00
        and r_squared >= 0.90
        and rmse <= 2.5
    )

    return {
        "tau_seconds": tau,
        "r_squared": r_squared,
        "rmse_mmhg": rmse,
        "samples_used": int(pressure.size),
        "success": bool(result.success),
        "acceptable_fit": acceptable,
    }