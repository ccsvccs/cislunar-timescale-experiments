"""Analytic numerical checks; these are deliberately not DE440 validation."""
import math

import mpmath as mp
import numpy as np
import pytest

from tt_tdb.integration import SegmentedIntegrator, integrate_rate


def polynomial_rate(t):
    return np.column_stack((t ** 4 - 2 * t + 1, 3e-16 * np.ones_like(t)))


def polynomial_integral(a, b):
    def primitive(t):
        return t ** 5 / 5 - t ** 2 + t
    return np.array([primitive(b) - primitive(a), 3e-16 * (b - a)])


@pytest.mark.parametrize("method,order", [("romberg", 10), ("gauss_legendre", 16)])
@pytest.mark.parametrize("a,b", [(-2, 3), (3, -2), (1.25, 1.25)])
def test_independent_quadratures_analytic_polynomial(method, order, a, b):
    result = integrate_rate(polynomial_rate, a, b, method=method, order=order)
    np.testing.assert_allclose(result.value, polynomial_integral(a, b), rtol=5e-15, atol=1e-28)
    assert np.all(result.error >= 0)
    assert result.evaluations > 0
    if method == "romberg":
        assert result.evaluations == 513
        assert result.levels == 10


def test_romberg_order_convergence_to_mpmath():
    def rate(t):
        return np.column_stack((np.exp(t), 1e-16 * np.cos(t)))
    with mp.workdps(60):
        expected = np.array([float(mp.e - 1), float(mp.mpf("1e-16") * mp.sin(1))])
    errors = [abs(integrate_rate(rate, 0, 1, order=n).value[0] - expected[0])
              for n in (2, 4, 6, 10)]
    assert errors[2] < errors[1] < errors[0]
    assert errors[-1] < 1e-15
    np.testing.assert_allclose(integrate_rate(rate, 0, 1, method="gauss_legendre", order=16).value,
                               expected, rtol=3e-15, atol=1e-30)


@pytest.mark.parametrize("method", ["romberg", "gauss_legendre"])
def test_forward_backward_anchored_cache_and_shared_queries(method):
    cache = SegmentedIntegrator(polynomial_rate, anchor=0.75, step_seconds=0.5,
                                 method=method, fingerprint="analytic-polynomial-only")
    queries = np.array([[0.75, -1.23], [3.17, -0.75]])
    result = cache.integral(queries)
    expected = np.stack([polynomial_integral(0.75, t) for t in queries.ravel()]).reshape(2, 2, 2)
    np.testing.assert_allclose(result, expected, rtol=5e-15, atol=1e-27)
    evaluations = cache.stats["evaluations"]
    np.testing.assert_array_equal(cache.integral(queries), result)
    assert cache.stats["evaluations"] == evaluations
    assert cache.stats["query_intervals"] == 2
    assert cache.integral(0.75).shape == (2,)
    assert cache.integral(np.array([])).shape == (0, 2)


def test_cache_resume_and_fingerprint_validation(tmp_path):
    cache = SegmentedIntegrator(polynomial_rate, anchor=0, step_seconds=0.5,
                                 fingerprint={"model": "synthetic", "sha256": "test"})
    cache.build(-2, 2)
    path = tmp_path / "checkpoint.npz"
    cache.save(path)
    resumed = SegmentedIntegrator.load(path, polynomial_rate,
                                       fingerprint={"sha256": "test", "model": "synthetic"})
    original_evaluations = resumed.stats["evaluations"]
    resumed.build(3)
    assert resumed.stats["evaluations"] - original_evaluations == 2 * 513
    np.testing.assert_allclose(resumed.integral([-2, 3]),
                               np.stack([polynomial_integral(0, -2), polynomial_integral(0, 3)]),
                               rtol=5e-15, atol=1e-29)
    np.testing.assert_array_equal(resumed.segment_errors[:len(cache.segment_errors)],
                                  cache.segment_errors)
    with pytest.raises(ValueError, match="fingerprint"):
        SegmentedIntegrator.load(path, polynomial_rate, fingerprint="changed-model")
    with pytest.raises(ValueError, match="fingerprint"):
        SegmentedIntegrator.load(path, polynomial_rate, fingerprint=cache.model_fingerprint, order=8)
    with pytest.raises(ValueError, match="fingerprint"):
        SegmentedIntegrator.load(path, polynomial_rate, fingerprint=cache.model_fingerprint, anchor=1)


def test_corrupt_cache_rejected(tmp_path):
    cache = SegmentedIntegrator(polynomial_rate, anchor=0, step_seconds=1,
                                 fingerprint="synthetic-only")
    cache.build(1)
    path = tmp_path / "checkpoint.npz"
    cache.save(path)
    with np.load(path, allow_pickle=False) as saved:
        arrays = {name: saved[name].copy() for name in saved.files}
    arrays["totals"][1, 0] += 1
    np.savez_compressed(path, **arrays)
    with pytest.raises(ValueError, match="checksum"):
        SegmentedIntegrator.load(path, polynomial_rate, fingerprint="synthetic-only")


def test_compensated_accumulation_keeps_small_channel():
    # Exact piecewise-constant rate with enormous cancellation; GL nodes stay
    # within each segment. The expected sum is independently math.fsum'd.
    increments = np.array([1e16, 1, -1e16] * 100)
    def rate(t):
        indices = np.minimum(np.floor(t).astype(int), len(increments) - 1)
        return np.column_stack((increments[indices], 1e-16 * np.ones_like(t)))
    cache = SegmentedIntegrator(rate, anchor=0, step_seconds=1,
                                 method="gauss_legendre", order=2,
                                 fingerprint="synthetic-cancellation-test")
    result = cache.integral(len(increments))
    assert result[0] == math.fsum(increments)
    assert abs(result[1] - len(increments) * 1e-16) < 1e-28


def test_vectorization_does_not_call_rate_once_per_node():
    calls = []
    def rate(t):
        calls.append(t.shape)
        return np.column_stack((np.sin(t / 100), np.ones_like(t) * 1e-16))
    cache = SegmentedIntegrator(rate, anchor=0, step_seconds=1, order=10,
                                 batch_segments=32, fingerprint="synthetic-batching")
    cache.build(64)
    assert len(calls) == 20
    assert cache.stats["evaluations"] == 64 * 513
    assert max(shape[0] for shape in calls) == 32 * 256


@pytest.mark.parametrize("step", [86400, 43200, 21600])
def test_annual_periodic_synthetic_mechanism(step):
    # Leading-order Kepler test, amplitude 2 e GM_sun/(a c^2 n) ~1.66ms.
    # This verifies mechanism and integration only, not real ephemeris accuracy.
    year = 365.25 * 86400
    angular_rate = 2 * np.pi / year
    amplitude = 2 * 0.0167 * 1.32712440018e20 / (149597870700 * 299792458 ** 2 * angular_rate)
    def rate(t):
        return np.column_stack((amplitude * angular_rate * np.cos(angular_rate * t),
                                np.zeros_like(t)))
    cache = SegmentedIntegrator(rate, anchor=0, step_seconds=step, order=8,
                                 fingerprint="synthetic-leading-order-kepler")
    times = np.linspace(0, year, 33)
    np.testing.assert_allclose(cache.integral(times)[:, 0], amplitude * np.sin(angular_rate * times),
                               rtol=0, atol=2e-17)
    assert 0.0016 < amplitude < 0.0017


def test_rate_validation():
    with pytest.raises(ValueError, match="rate must return"):
        integrate_rate(lambda t: np.zeros_like(t), 0, 1)
    with pytest.raises(ValueError, match="nonfinite"):
        integrate_rate(lambda t: np.full((len(t), 2), np.nan), 0, 1)
    with pytest.raises(ValueError, match="fingerprint"):
        SegmentedIntegrator(polynomial_rate, anchor=0, fingerprint={})
