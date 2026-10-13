"""Independent SPICE checks of genuine DE440 interpolation and units."""
import numpy as np
import pytest
import spiceypy as spice

from tt_tdb.ephemeris import DE440Ephemeris, BODY_IDS, read_gm_kernel


def test_gm_parser_ignores_comment_and_converts_si(tmp_path):
    from pathlib import Path
    source = Path('data/gm_de440.tpc')
    if not source.exists():
        pytest.skip('Real GM PCK not downloaded')
    gm = read_gm_kernel(source)
    assert gm[10] == pytest.approx(1.3271244004127942e20, rel=1e-15)
    assert gm[3] == pytest.approx(gm[399] + gm[301], rel=2e-15)


def test_missing_data_never_synthesized(tmp_path):
    with pytest.raises(FileNotFoundError):
        DE440Ephemeris(tmp_path/'absent.bsp', tmp_path/'absent.tpc')


@pytest.mark.integration
def test_states_against_spice(de440_paths):
    bsp, gm = de440_paths
    spice.furnsh(str(bsp))
    try:
        with DE440Ephemeris(bsp, gm) as eph:
            epochs = np.array([-725803167.8160655, 0., 30*86400., 25*365.25*86400.])
            result = eph.states(epochs)
            assert result.earth_position.shape == (4, 3)
            for index, epoch in enumerate(epochs):
                reference, _ = spice.spkezr('399', float(epoch), 'J2000', 'NONE', '0')
                np.testing.assert_allclose(result.earth_position[index], reference[:3]*1000, rtol=0, atol=2e-4)
                np.testing.assert_allclose(result.earth_velocity[index], reference[3:]*1000, rtol=0, atol=1e-10)
                for name, body in BODY_IDS.items():
                    expected, _ = spice.spkpos(str(body), float(epoch), 'J2000', 'NONE', '0')
                    np.testing.assert_allclose(result.body_positions[name][index], expected*1000, rtol=0, atol=2e-3)
    finally:
        spice.unload(str(bsp))


@pytest.mark.integration
def test_shapes_physical_units_coverage(de440_paths):
    with DE440Ephemeris(*de440_paths) as eph:
        scalar = eph.states(0.)
        array = eph.states(np.zeros((2, 3)))
        assert scalar.earth_position.shape == (3,)
        assert array.earth_position.shape == (2, 3, 3)
        assert 1.4e11 < np.linalg.norm(scalar.earth_position) < 1.6e11
        assert 2.8e4 < np.linalg.norm(scalar.earth_velocity) < 3.2e4
        assert 3.5e8 < np.linalg.norm(scalar.body_positions['moon'] - scalar.earth_position) < 4.2e8
        assert 'earth' not in scalar.gm and 'emb' not in scalar.gm
        np.testing.assert_array_equal(eph.states_jd(2451545., 0.).earth_position, scalar.earth_position)
        with pytest.raises(ValueError, match='coverage'):
            eph.states_jd(eph.start_jd - 1.)
        with pytest.raises(ValueError, match='coverage'):
            eph.states(float('nan'))
