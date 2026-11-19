import numpy as np
from cardio_twin.models.windkessel import simulate_two_element
from cardio_twin.estimation import estimate_resistance_compliance


def test_parameter_recovery_noise_free():
    time = np.linspace(0, 20, 2001)
    inflow = 15 + 10 * np.maximum(np.sin(2 * np.pi * 1.2 * time), 0)
    pressure = simulate_two_element(time, inflow, 1.3, 1.7)
    result = estimate_resistance_compliance(time, inflow, pressure)
    assert result["success"]
    assert abs(result["resistance"] - 1.3) < 0.05
    assert abs(result["compliance"] - 1.7) < 0.05
