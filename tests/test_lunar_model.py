from pathlib import Path
import numpy as np
import pytest
import spiceypy as spice
import erfa
from tt_tdb.constants import L_B, TDB0, T0_TDB
from tt_tdb.lunar_model import (vondrak_poles, earth_pole_park, external_j2_potential,
                               MoonObserverEphemeris, LunarModel, LunarSmallBodyPotential)


def test_long_term_poles_against_independent_erfa_reference():
    years = np.array([1550., 1977., 2000., 2026., 2650.])
    seconds = (years-2000.)*365.25*86400.
    n, k = vondrak_poles(seconds)
    np.testing.assert_allclose(n, erfa.ltpequ(years), rtol=0, atol=3e-16)
    np.testing.assert_allclose(k, erfa.ltpecl(years), rtol=0, atol=3e-16)
    pole = earth_pole_park(seconds)
    np.testing.assert_allclose(np.linalg.norm(pole, axis=-1), 1., rtol=0, atol=5e-16)
    # Precession is physical, not a fixed J2000 z-axis placeholder.
    assert pole[0, 0] < -.04
    assert pole[-1, 0] > .06
    assert np.max(np.linalg.norm(pole-n, axis=-1)) < 1e-4


def test_earth_j2_sign_at_equator_and_poles():
    u = external_j2_potential(np.array([[10., 0., 0.], [0., 0., 10.]]),
                              1., 1., .01, [0., 0., 1.])
    np.testing.assert_allclose(u, [5e-6, -1e-5], rtol=1e-15)


def test_declared_frame_offset_units_order_and_auxiliary_exclusion():
    seconds = 820497600.
    fx, fy = .005343991604414805, -.01712882484031753
    # ERFA supplies an independent inertial->date matrix. CSPICE supplies
    # the passive elementary rotations with Park's specified order.
    date = erfa.ltp(2000.+seconds/(365.25*86400.))
    zero = earth_pole_park(seconds)
    expected = date.T @ spice.rotate(-fx*np.pi/(180*3600), 1) @ spice.rotate(-fy*np.pi/(180*3600), 2) @ date @ zero
    actual = earth_pole_park(seconds, x_offset_arcsec=fx, y_offset_arcsec=fy)
    np.testing.assert_allclose(actual, expected, rtol=0, atol=3e-16)
    model = LunarModel.__new__(LunarModel)
    model.centering_rate = 1.4825e-8
    integrals = np.zeros((1, 11))
    integrals[:, 10] = 1e6  # diagnostic channel must never enter physics.
    assert model.restore([T0_TDB.seconds_since_j2000()], integrals)[0] == -np.longdouble(TDB0)


def test_exact_anchor_and_centering_rate_is_not_a_physical_parameter():
    anchor = T0_TDB.seconds_since_j2000()
    t = np.array([anchor, anchor+86400., anchor-4e9])
    f = np.longdouble('1.4825362e-8')
    elapsed = t.astype(np.longdouble)-np.longdouble(anchor)
    expected = -np.longdouble(TDB0)+(np.longdouble(L_B)-f)*elapsed/(1-np.longdouble(L_B))
    for center in (1.4825e-8, 1.48e-8):
        model = LunarModel.__new__(LunarModel)
        model.centering_rate = center
        integrals = ((f-np.longdouble(center))*elapsed)[:, None]
        q = model.restore(t, integrals)
        np.testing.assert_allclose(q, expected, rtol=0, atol=1e-18)
        assert q[0] == -np.longdouble(TDB0)


@pytest.mark.integration
def test_moon_observer_states_and_external_body_identity(de440_paths):
    spice.furnsh(str(de440_paths[0]))
    try:
        with MoonObserverEphemeris(*de440_paths) as eph:
            t = np.array([0., 757339200., 820497600.])
            state = eph.states(t)
            assert 'moon' not in state.body_positions
            assert 'earth_external' in state.body_positions
            assert state.earth_gm == eph.all_gm[301]
            assert state.gm['earth_external'] == eph.all_gm[399]
            for i, seconds in enumerate(t):
                moon, _ = spice.spkezr('301', float(seconds), 'J2000', 'NONE', '0')
                np.testing.assert_allclose(state.earth_position[i], moon[:3]*1000., rtol=0, atol=5e-5)
                np.testing.assert_allclose(state.earth_velocity[i], moon[3:]*1000., rtol=0, atol=1e-10)
            model = LunarModel(eph)
            channels = model(t)
            assert channels.shape == (3, 11)
            assert np.all(np.isfinite(channels))
            assert model.fingerprint()['reference_used_in_integrand'] is False
            if Path('data/sb441-n373s.bsp').is_file():
                with LunarSmallBodyPotential(eph) as bodies:
                    lunar_rates = bodies(t)
                    assert np.all(lunar_rates > 0)
                    expected = ((L_B-model.centering_rate)-channels[:, :10].sum(-1)-lunar_rates.sum(-1))/(1-L_B)
                    np.testing.assert_allclose(model.rate(t, lunar_rates), expected, rtol=0, atol=2e-25)
    finally:
        spice.unload(str(de440_paths[0]))
