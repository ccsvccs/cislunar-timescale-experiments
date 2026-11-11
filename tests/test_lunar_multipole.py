"""Analytic gravity conventions and independent NAIF PA orientation checks."""
from pathlib import Path
import numpy as np
import pytest
import spiceypy as spice
from tt_tdb.constants import C_LIGHT, T0_TDB
from tt_tdb.ephemeris import DE440Ephemeris
from tt_tdb.jpl_model import read_de440_constants
from tt_tdb.lunar_multipole import (degree2_potential, LunarPAOrientation,
                                     LunarMultipolePotential)


def test_degree2_equator_poles_and_principal_axes():
    gm, radius, distance, j2, c22 = 5., 2., 20., 0.01, 0.003
    r = distance*np.array([[1., 0., 0.], [0., 1., 0.], [0., 0., 1.], [0., 0., -1.]])
    u = degree2_potential(r, gm, radius, j2=j2, c22=c22)
    scale = gm*radius**2/distance**3
    np.testing.assert_allclose(u[:, 0], scale*j2*np.array([.5, .5, -1., -1.]), rtol=2e-15)
    np.testing.assert_allclose(u[:, 3], scale*3*c22*np.array([1., -1., 0., 0.]), rtol=2e-15)
    np.testing.assert_array_equal(u[:, [1, 2, 4]], 0.)


def test_unnormalized_tesseral_sign_and_factors():
    # Unit directions isolate P21=3xz, P22*sin(2lambda)=6xy.
    r = 10*np.array([[1., 0., 1.], [0., 1., 1.], [1., 1., 0.]])/np.sqrt(2)
    u = degree2_potential(r, 1., 1., j2=0., c21=1., s21=1., s22=1.)
    np.testing.assert_allclose(u[:, 1], [1.5e-3, 0., 0.], atol=0, rtol=2e-15)
    np.testing.assert_allclose(u[:, 2], [0., 1.5e-3, 0.], atol=0, rtol=2e-15)
    np.testing.assert_allclose(u[:, 4], [0., 0., 3.e-3], atol=0, rtol=2e-15)
    np.testing.assert_allclose(degree2_potential(-r, 1., 1., j2=0., c21=1., s21=1., s22=1.), u)
    with pytest.raises(ValueError, match='external-potential'):
        degree2_potential([0., 0., .5], 1., 1., j2=.1)


def test_header_j2_c22_consistent_with_undistorted_inertia():
    h = read_de440_constants('docs/sources/de440_header.440')
    j, b, g = h['J2M'], h['LBET'], h['LGAM']
    # Park (2021) Eqs. (47)--(51), divided by m_M R_M².
    a = 2*(1-b*g)*j/(2*b-g+b*g)
    b_moment = 2*(1+g)*j/(2*b-g+b*g)
    assert h['AM'] == 1738.
    assert (b_moment-a)/4 == pytest.approx(h['C22M'], rel=4e-13)
    assert h['S22M'] == 0.


def test_missing_orientation_has_no_toy_fallback(tmp_path):
    with pytest.raises(FileNotFoundError, match='authentic DE440'):
        LunarPAOrientation(tmp_path/'missing.bpc')


@pytest.mark.integration
def test_real_orientation_against_naif_and_published_direction(de440_paths):
    pck = Path('data/moon_pa_de440_200625.bpc')
    fk = Path('data/moon_de440_250416.tf')
    if not pck.is_file() or not fk.is_file():
        pytest.skip('Official DE440 lunar binary PCK/FK unavailable')
    for path in (de440_paths[0], pck, fk):
        spice.furnsh(str(path))
    try:
        with DE440Ephemeris(*de440_paths) as eph, LunarMultipolePotential(eph) as model:
            times = np.array([T0_TDB.seconds_since_j2000(), 0., 717768000., 820497600.])
            r = model.orientation.matrices(times)
            np.testing.assert_allclose(r@np.swapaxes(r, -1, -2), np.broadcast_to(np.eye(3), r.shape), atol=8e-16)
            np.testing.assert_allclose(np.linalg.det(r), 1., atol=8e-16)
            independent = np.stack([spice.pxform('J2000', 'MOON_PA_DE440', float(t)) for t in times])
            np.testing.assert_allclose(r, independent, rtol=0, atol=5e-12)
            positions = np.stack([spice.spkpos('399', float(t), 'MOON_PA_DE440', 'NONE', '301')[0]
                                  for t in times])*1000.
            np.testing.assert_allclose(model.displacement_pa(times), positions, rtol=0, atol=.003)
            # Official FK example: 2022-09-30 TDB, km rounded to 1 m.
            np.testing.assert_allclose(positions[2]/1000, [373997.028, -23558.987, 10284.057], rtol=0, atol=.0006)
            n = positions/np.linalg.norm(positions, axis=-1)[:, None]
            assert np.all(n[:, 0] > .98)  # Earth near long +X PA axis, with actual libration.
            assert np.max(abs(n[:, 1])) > .01
            expected = degree2_potential(positions, model.gm, model.radius, **model.coefficients)/C_LIGHT**2
            np.testing.assert_allclose(model(times), expected, rtol=1e-9, atol=1e-33)
            with pytest.raises(ValueError, match='coverage'):
                model.orientation.matrices(model.orientation.end+1.)
            # Independent orientation decoding on both sides of the binary
            # PCK split; this is a format check, not a long-time physics test.
            boundary = model.orientation.segments[0].final_second
            split_times = np.array([boundary-1., boundary, boundary+1.])
            split = model.orientation.matrices(split_times)
            split_reference = np.stack([spice.pxform('J2000', 'MOON_PA_DE440', float(t)) for t in split_times])
            np.testing.assert_allclose(split, split_reference, rtol=0, atol=1e-10)
    finally:
        for path in (fk, pck, de440_paths[0]):
            spice.unload(str(path))
