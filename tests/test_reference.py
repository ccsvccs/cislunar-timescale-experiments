import numpy as np
from tt_tdb.reference import difference_metrics, sofa_geocentric


def test_geocentric_reference_and_metrics():
    # Independent published ERFA/SOFA test value at jd=2448939.5+0.123,
    # but the published case is topocentric; here check geocentric magnitude
    # and independently evaluate a known linear diagnostic.
    d = sofa_geocentric(2451545.0, np.arange(370, dtype=float))
    assert np.max(np.abs(d)) < .002
    assert np.ptp(d) > .003
    m = difference_metrics(np.arange(10)*86400, 2e-9+np.arange(10)*3e-12)
    assert abs(m['linear_slope_seconds_per_day']-3e-12) < 1e-24
    assert abs(m['mean_seconds']-(2e-9+4.5*3e-12)) < 1e-24
