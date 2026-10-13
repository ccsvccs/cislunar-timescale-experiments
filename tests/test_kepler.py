"""Exact synthetic Sun-Earth ellipse: mechanism and quadrature only."""
import numpy as np
import pytest

from tt_tdb.integration import SegmentedIntegrator
from tt_tdb.validation import KeplerReference


@pytest.mark.parametrize("method,order", [("romberg", 10), ("gauss_legendre", 24)])
def test_exact_kepler_integral_annual_amplitude(method, order):
    orbit = KeplerReference()
    times = np.linspace(-orbit.period, orbit.period, 65)
    integration = SegmentedIntegrator(orbit.rate, anchor=0, channels=1,
                                       method=method, order=order,
                                       fingerprint="explicit-synthetic-Kepler-test")
    expected = orbit.integral(0, times)
    actual = integration.integral(times)[:, 0]
    np.testing.assert_allclose(actual, expected, atol=1e-15, rtol=0)
    assert 0.0016 < orbit.amplitude < 0.0017
    assert 365 < orbit.period / 86400 < 366


def test_kepler_rate_matches_analytic_time_derivative():
    orbit = KeplerReference()
    time = np.array([0, orbit.period / 4, orbit.period / 2])
    anomaly = orbit.eccentric_anomaly(time)
    analytic_derivative = (orbit.amplitude * np.cos(anomaly)
                           * orbit.n / (1 - orbit.e * np.cos(anomaly)))
    np.testing.assert_allclose(orbit.rate(time)[:, 0], analytic_derivative, rtol=2e-13, atol=1e-23)
