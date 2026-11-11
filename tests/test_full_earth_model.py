"""Composition and reference-independence checks for the explicit Earth model."""
from pathlib import Path
import numpy as np
import pytest
from tt_tdb.constants import A_SCALE
from tt_tdb.full_earth_model import FullEarthModel, anchored_prefix


def test_fixed_anchor_split_and_bidirectional_compensated_accumulation():
    # Integral of f(t)=t from a physical anchor t=1.25 to grid 0,1,2,3.
    increments=np.array([[.5],[1.5],[2.5]])
    value=anchored_prefix(increments,1,np.array([(1.-1.25**2)/2]))
    np.testing.assert_array_equal(value[:,0],(np.arange(4.)**2-1.25**2)/2)
    # A tiny increment between cancelling large values must survive summation.
    result=anchored_prefix(np.array([[1e16],[1.],[-1e16]]),0,np.array([0.]))
    assert result[-1,0]==1.


def test_configuration_excludes_spin_and_literal_lunar_coefficients():
    assert FullEarthModel.axis_index not in FullEarthModel.primary_indices
    assert len(FullEarthModel.channel_names)==14
    assert not any('spin' in n or 'C21' in n or 'S21' in n or 'S22' in n for n in FullEarthModel.channel_names)
    source=Path('src/tt_tdb/full_earth_model.py').read_text()
    assert 'jpl_reference' not in source and 'solar_spin import' not in source


@pytest.mark.de440
def test_exact_component_composition_and_modern_archived_rates(de440_paths):
    for path in ('data/sb441-n373s.bsp','data/moon_pa_de440_200625.bpc','data/moon_de440_250416.tf'):
        if not Path(path).exists():pytest.skip('Actual independent component data unavailable')
    times=np.array([757339200.,788961600.,820497600.])
    with FullEarthModel() as model:
        values=model(times)
        lunar=model.lunar(times)
        np.testing.assert_array_equal(values[:,9],lunar[:,0])
        np.testing.assert_array_equal(values[:,10],lunar[:,3])
        np.testing.assert_array_equal(values[:,12:],model.small(times))
        park=A_SCALE*values[:,:9].sum(axis=-1)
        # A separate call to the unchanged Park core reconstructs the same
        # model with expected floating-order differences in leading subtraction.
        np.testing.assert_allclose(park,A_SCALE*model.park(times).sum(axis=-1),rtol=0,atol=2e-23)
        old=np.genfromtxt('results/tables/phase2_sb441_time_series.csv',delimiter=',',names=True)
        indices=np.searchsorted(old['tdb_j2000_seconds'],times)
        actual=A_SCALE*values[:,list(model.primary_indices)].sum(axis=-1)
        expected=old['sb441_rate'][indices]+A_SCALE*(lunar[:,0]+lunar[:,3])
        np.testing.assert_allclose(actual,expected,rtol=0,atol=2e-23)
        assert model.start_jd==2287184.5 and model.end_jd==2688976.5
