"""Numerical reference units and sign verified against independent SPICE."""
from pathlib import Path
import os
import numpy as np
import pytest
import spiceypy as spice
from tt_tdb.jpl_reference import JPLTimeEphemeris


@pytest.fixture
def time_path():
    path = Path(os.environ.get('JPL_TIME_BSP', 'data/de440t.bsp'))
    if not path.is_file():
        pytest.skip('Genuine JPL numerical TT−TDB reference not downloaded')
    return path


def test_missing_reference_is_not_sofa(tmp_path):
    with pytest.raises(FileNotFoundError):
        JPLTimeEphemeris(tmp_path/'absent.bsp')


@pytest.mark.integration
def test_time_reference_spice_sign_and_derivative(time_path):
    spice.furnsh(str(time_path))
    try:
        with JPLTimeEphemeris(time_path) as ref:
            times = np.array([-725803167.8160655, 0., 365.25*86400., 25*365.25*86400.])
            delta, rate = ref.evaluate(times)
            for i, time in enumerate(times):
                state, _ = spice.spkezr('1000000001', float(time), 'J2000', 'NONE', '1000000000')
                assert delta[i] == pytest.approx(-state[0], rel=0, abs=2e-16)
                assert rate[i] == pytest.approx(-state[3], rel=0, abs=3e-22)
            assert ref.evaluate(0.)[0] < 0.
            # Differentiate a real time polynomial independently by symmetric
            # finite differences; the hour interval makes rounding negligible.
            h = 3600.
            finite = (ref.evaluate(h)[0] - ref.evaluate(-h)[0])/(2*h)
            assert finite == pytest.approx(ref.evaluate(0.)[1], rel=2e-8)
            with pytest.raises(ValueError):
                ref.evaluate_jd(ref.start_jd-1.)
    finally:
        spice.unload(str(time_path))
