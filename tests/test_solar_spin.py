"""Geometric, invariance and fixed-anchor checks of the separate spin diagnostic."""
import numpy as np
import pytest
from tt_tdb.constants import L_B
from tt_tdb.solar_spin import spin_rate_components, solar_spin_parameter, solar_axis


def test_prograde_retrograde_radial_and_spin_reversal():
    # Unit analytic fixture: prograde angular momentum parallel +z solar spin.
    r = np.array([2., 0., 0.]); ve = np.array([0., 3., 0.]); vs = np.zeros(3)
    gs = np.array([0., 0., 5.])
    pro = spin_rate_components(r, ve, vs, gs, c=1.)
    assert pro['scalar'] == 0
    assert pro['vector'] == pro['combined'] == pro['relative'] == -7.5
    retro = spin_rate_components(r, -ve, vs, gs, c=1.)
    assert retro['combined'] == 7.5
    assert spin_rate_components(r, ve, vs, -gs, c=1.)['combined'] == 7.5
    assert spin_rate_components(r, r, vs, gs, c=1.)['combined'] == 0


def test_scalar_vector_closure_and_common_velocity_cancellation():
    # A common barycentric boost changes individual scalar/vector pieces;
    # their complete contribution must depend only on relative velocity.
    r=np.array([[2., 1., -3.], [4., -2., 1.]])
    ve=np.array([[2., 3., 4.], [-1., 5., 2.]])
    vs=np.array([[1., -2., 3.], [2., 4., 1.]])
    gs=np.array([2., -4., 5.]); boost=np.array([9., -7., 11.])
    first=spin_rate_components(r, ve, vs, gs, c=1.)
    shifted=spin_rate_components(r, ve+boost, vs+boost, gs, c=1.)
    np.testing.assert_allclose(first['combined'], first['relative'], rtol=2e-15, atol=1e-15)
    np.testing.assert_allclose(shifted['combined'], first['combined'], rtol=3e-15, atol=2e-15)
    assert not np.allclose(first['scalar'], shifted['scalar'])
    same=spin_rate_components(r, ve, ve, gs, c=1.)
    np.testing.assert_allclose(same['combined'], np.zeros(2), rtol=0, atol=1e-15)


def test_consistent_coordinate_spin_scaling_and_rotation_invariance():
    r=np.array([3., -4., 5.]); ve=np.array([7., 2., 1.]); vs=np.array([1., -3., 2.])
    gs=np.array([2., 3., 4.]); k=.7
    a=spin_rate_components(r, ve, vs, gs)
    b=spin_rate_components(k*r, ve, vs, k*k*gs)
    q=np.array([[0., 1., 0.], [0., 0., 1.], [1., 0., 0.]])
    rotated=spin_rate_components(q@r, q@ve, q@vs, q@gs)
    for name in a:
        # This fixture deliberately has cancellation between scalar/vector.
        assert b[name] == pytest.approx(a[name], rel=1e-14, abs=0)
        assert rotated[name] == pytest.approx(a[name], rel=1e-14, abs=0)
    np.testing.assert_allclose(solar_spin_parameter(), (1-L_B)**2*solar_spin_parameter(tdb_compatible=False), rtol=2e-16)
    assert np.linalg.norm(solar_axis()) == pytest.approx(1., abs=2e-16)


def test_invalid_geometry_rejected():
    for r in (np.zeros(3), np.array([np.nan, 0, 1]), np.zeros(2)):
        with pytest.raises(ValueError):
            spin_rate_components(r, np.ones(3), np.zeros(3), np.ones(3))


@pytest.mark.de440
def test_real_states_anchor_and_combined_quadrature(de440_paths):
    from tt_tdb.extended_ephemeris import ExtendedDE440Ephemeris
    from tt_tdb.solar_spin import SolarSpinDiagnostic
    from tt_tdb.integration import SegmentedIntegrator
    from tt_tdb.constants import T0_TDB
    with ExtendedDE440Ephemeris(*de440_paths) as e:
        diagnostic=SolarSpinDiagnostic(e)
        anchor=T0_TDB.seconds_since_j2000()
        t=np.array([anchor, anchor+1e6, 8e8])
        rates=diagnostic(t)
        np.testing.assert_allclose(rates[:, 2], rates[:, 3], rtol=5e-16, atol=1e-36)
        assert np.all(rates[:, 2] < 0)
        integ=SegmentedIntegrator(diagnostic, anchor=anchor, method='gauss_legendre', order=12,
                                 channels=5, fingerprint=diagnostic.fingerprint())
        value=integ.integral([anchor, anchor+86400])
        np.testing.assert_array_equal(value[0], np.zeros(5))
        assert value[1, 2] < 0
        assert value[1, 2] == pytest.approx(value[1, 0]+value[1, 1], rel=5e-16)
