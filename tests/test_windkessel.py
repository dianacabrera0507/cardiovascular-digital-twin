import numpy as np
import pytest
from cardio_twin.models.windkessel import simulate_two_element, diastolic_time_constant


def test_positive_time_constant():
    assert diastolic_time_constant(2.0, 3.0) == 6.0


def test_simulation_shape_and_finiteness():
    time = np.linspace(0, 2, 201)
    inflow = np.ones_like(time) * 10
    pressure = simulate_two_element(time, inflow, 1.0, 1.0)
    assert pressure.shape == time.shape
    assert np.isfinite(pressure).all()


def test_rejects_nonpositive_parameters():
    time = np.linspace(0, 1, 10)
    with pytest.raises(ValueError):
        simulate_two_element(time, np.ones(10), 0.0, 1.0)
