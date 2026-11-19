import numpy as np
import pytest

from cardio_twin.abp import (
    characterize_beats,
    detect_systolic_peaks,
    estimate_effective_diastolic_tau,
    summarize_abp_window,
)


def test_summarize_abp_window_handles_missing_values():
    values = np.array(
        [60.0, 80.0, 100.0, np.nan]
    )

    summary = summarize_abp_window(values)

    assert summary["n_samples"] == 4
    assert summary["finite_fraction"] == pytest.approx(
        0.75
    )
    assert summary["median_mmhg"] == pytest.approx(
        80.0
    )
    assert summary["robust_amplitude_mmhg"] > 0


def test_detect_and_characterize_periodic_abp():
    sampling_hz = 100.0

    time = np.arange(
        0.0,
        10.0,
        1.0 / sampling_hz,
    )

    phase = np.mod(time, 1.0)

    signal = (
        70.0
        + 40.0
        * np.exp(-phase / 0.25)
    )

    peaks = detect_systolic_peaks(
        signal,
        sampling_hz,
    )

    beats = characterize_beats(
        signal,
        peaks,
        sampling_hz,
    )

    assert 8 <= len(peaks) <= 11
    assert len(beats) >= 7

    assert all(
        row["plausible"]
        for row in beats
    )


def test_tau_recovery_on_exact_exponential_decay():
    sampling_hz = 100.0
    tau_true = 2.0
    venous_pressure = 5.0

    time = (
        np.arange(101)
        / sampling_hz
    )

    values = (
        venous_pressure
        + (
            100.0
            - venous_pressure
        )
        * np.exp(
            -time / tau_true
        )
    )

    result = estimate_effective_diastolic_tau(
        values,
        sampling_hz,
        start_peak=0,
        diastolic_index=100,
        venous_pressure_mmhg=venous_pressure,
        fit_start_fraction=0.35,
    )

    assert result["success"]
    assert result["acceptable_fit"]

    assert result[
        "tau_seconds"
    ] == pytest.approx(
        tau_true,
        abs=0.02,
    )

    assert result["r_squared"] > 0.999


def test_peak_detection_rejects_missing_values():
    values = np.array(
        [70.0, np.nan, 90.0, 80.0, 70.0]
    )

    with pytest.raises(ValueError):
        detect_systolic_peaks(
            values,
            100.0,
        )