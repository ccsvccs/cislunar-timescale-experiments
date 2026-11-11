from pathlib import Path
import numpy as np
import pytest
import spiceypy as spice
from tt_tdb.constants import T0_TDB, TDB0
from tt_tdb.lte_reference import LTE440Reference


def test_missing_and_wrong_companion_are_rejected(tmp_path):
    with pytest.raises(FileNotFoundError):
        LTE440Reference(tmp_path/'missing.bsp')
    if not Path('data/lte440.bsp').is_file():
        pytest.skip('Fixed LTE440 reference not installed')
    bad = tmp_path/'wrong.tpc'
    bad.write_text('KPL/PCK\n\\begindata\nBODY1000000005_RATE = (0)\n')
    with pytest.raises(ValueError, match='SHA256'):
        LTE440Reference('data/lte440.bsp', bad)


@pytest.mark.integration
def test_original_author_samples_anchor_endpoints_and_spice():
    if not Path('data/lte440.bsp').is_file():
        pytest.skip('Fixed LTE440 reference not installed')
    with LTE440Reference() as r:
        assert r.rate == 6.798355238724282e-10
        anchor = r.evaluate_jd(*T0_TDB.jd_parts())[0]
        assert abs(anchor+TDB0) < 1e-16
        values, rates = r.evaluate_jd(np.array([2305445., 2451545., 2524595.]))
        np.testing.assert_allclose(values, [-8.0882739246330395, .49330749643254862, 4.7842468668519906], rtol=0, atol=2e-15)
        # Actual record boundaries, both coverage endpoints, and J2000.
        start, interval, coefficients = r.segment.load_array()
        dates = np.array([r.start_jd, start+interval, 2451545., r.end_jd])
        values, rates = r.evaluate_jd(dates)
        spice.furnsh(str(r.path))
        try:
            for i, date in enumerate(dates):
                t = (date-2451545.)*86400.
                pv, _ = spice.spkez(1000000005, float(t), 'J2000', 'NONE', 1000000000)
                linear = r.rate*(np.longdouble(t)-np.longdouble(T0_TDB.seconds_since_j2000()))
                assert abs(values[i]-(np.longdouble(pv[0])+linear)) < 5e-15
                assert abs(rates[i]-(pv[3]+r.rate)) < 1e-21
        finally:
            spice.unload(str(r.path))
        with pytest.raises(ValueError, match='coverage'):
            r.evaluate_jd(r.end_jd, 1e-8)
