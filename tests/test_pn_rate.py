import numpy as np
import mpmath as mp
import pytest
from tt_tdb.pn_rate import rate_c2,rate_c4,relativistic_rate,PNModel
from tt_tdb.ephemeris import DE440Ephemeris


def test_formula_against_high_precision_scalar_and_scaling():
    with mp.workdps(60):
        v2,u,c=map(mp.mpf,('9e8','8.9e8','299792458'))
        expected2=(v2/2+u)/c**2
        expected4=(v2**2/8+3*v2*u/2-u**2/2)/c**4
        assert abs(mp.mpf(float(rate_c2(float(v2),float(u))))-expected2)<mp.mpf('3e-24')
        assert abs(mp.mpf(float(rate_c4(float(v2),float(u))))-expected4)<mp.mpf('4e-32')
    v=np.array([7e8,9e8]);u=np.array([8e8,10e8])
    assert np.allclose(relativistic_rate(v,u,include_c4=False),rate_c2(v,u),rtol=0,atol=0)
    # Consistent spatial and GM scaling leaves GM/r, velocities and F invariant.
    assert np.all(rate_c4(v,u)>0)
    assert np.all(rate_c4(v,u)/rate_c2(v,u)<1e-7)


@pytest.mark.de440
def test_genuine_geocentric_rates(de440_paths):
    with DE440Ephemeris(*de440_paths) as eph:
        d=PNModel(eph).diagnostics(0.)
        assert 8e8<d['sun_potential']<1e9
        assert 8e8<d['velocity_squared']<1e9
        assert 1.4e-8<d['F2']<1.6e-8
        assert 8e-17<d['F4']<2e-16
        assert d['observer_position_term_seconds']==0.
        assert 'earth' not in d['potential_by_body']
        assert len(d['potential_by_body'])==10
