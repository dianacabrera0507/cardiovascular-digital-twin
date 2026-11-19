import numpy as np
import pytest

from cardio_twin.events import prediction_intervals, threshold_episodes


def test_threshold_episode_exact_minimum_duration():
    time = np.arange(0.0, 121.0)
    values = np.full(time.shape, 80.0)
    values[20:80] = 60.0
    assert threshold_episodes(time, values, 65.0, 60.0) == [(20, 80)]


def test_short_episode_is_rejected():
    time = np.arange(0.0, 121.0)
    values = np.full(time.shape, 80.0)
    values[20:79] = 60.0
    assert threshold_episodes(time, values, 65.0, 60.0) == []


def test_episodes_merge_across_prespecified_gap():
    time = np.arange(0.0, 181.0)
    values = np.full(time.shape, 80.0)
    values[20:50] = 60.0
    values[70:100] = 60.0
    assert threshold_episodes(time, values, 65.0, 60.0, merge_gap=20.0) == [(20, 100)]


def test_nonfinite_value_interrupts_episode():
    time = np.arange(0.0, 121.0)
    values = np.full(time.shape, 60.0)
    values[60] = np.nan
    assert threshold_episodes(time, values, 65.0, 60.0) == [
        (0, 60),
        (61, 121),
    ]


def test_nonfinite_interruption_prevents_short_runs_from_passing():
    time = np.arange(0.0, 120.0)
    values = np.full(time.shape, 60.0)
    values[59] = np.nan
    assert threshold_episodes(time, values, 65.0, 60.0) == [(60, 120)]


def test_prediction_intervals_are_leakage_safe():
    intervals = prediction_intervals(
        event_start_seconds=3600.0,
        predictor_window_seconds=300.0,
        exclusion_gap_seconds=300.0,
        prediction_horizon_seconds=900.0,
    )
    assert intervals == {
        "predictor": (2100.0, 2400.0),
        "exclusion_gap": (2400.0, 2700.0),
        "outcome_horizon": (2700.0, 3600.0),
    }
    assert intervals["predictor"][1] <= intervals["outcome_horizon"][0]


def test_prediction_intervals_require_history():
    with pytest.raises(ValueError, match="insufficient history"):
        prediction_intervals(1000.0, 300.0, 300.0, 900.0)


def test_time_must_be_strictly_increasing():
    with pytest.raises(ValueError, match="strictly increasing"):
        threshold_episodes(
            np.array([0.0, 1.0, 1.0]),
            np.array([80.0, 60.0, 60.0]),
            65.0,
            1.0,
        )
