"""Evidence-based target/GM mapping and independent SB441 SPICE checks."""
import json
from pathlib import Path
import numpy as np
import pytest
import spiceypy as spice
from tt_tdb.small_body import build_mapping, SmallBodyPotential
from tt_tdb.ephemeris import DE440Ephemeris
from tt_tdb.constants import C_LIGHT


def test_all_373_iom_targets_match_de440_masses_without_adjustment():
    m=build_mapping();rows=m['entries']
    assert len(rows)==373
    assert len({r['target'] for r in rows})==373
    assert len({r['header_name'] for r in rows})==373
    ast=[r for r in rows if r['group']=='asteroids'];kbo=[r for r in rows if r['group']=='kbos']
    assert len(ast)==343 and len(kbo)==30
    original=json.loads(Path('data/small_body_gm.json').read_text())['groups']
    for group, subset in [('asteroids',ast),('kbos',kbo)]:
        expected=original['asteroids' if group=='asteroids' else 'individual_kbos']['gm_m3_s2']
        for row in subset: assert row['gm_si'] == pytest.approx(expected[row['header_name']],rel=2e-15)
    assert next(r for r in kbo if r['target']==2136199)['header_name']=='MA8001'  # Eris, not MA8001==NAIF target!
    assert next(r for r in ast if r['target']==2000001)['header_name']=='MA0001'


def test_mapping_rejects_wrong_or_duplicated_inputs(tmp_path):
    text=Path('docs/sources/sb441_perturbers.txt').read_text()
    bad=tmp_path/'wrong.txt';bad.write_text(text.replace('1.3964518123081070e-13','1.4964518123081070e-13'))
    with pytest.raises(ValueError,match='disagree'):build_mapping(bad)
    bad.write_text(text+'\nCeres 1 2000001 1.3964518123081070e-13\n')
    with pytest.raises(ValueError,match='373 unique'):build_mapping(bad)


def test_missing_kernel_has_no_synthetic_fallback(tmp_path):
    with pytest.raises(FileNotFoundError,match='actual SB441'):
        SmallBodyPotential(None,tmp_path/'missing.bsp')


@pytest.mark.integration
def test_sb441_potentials_against_independent_spice(de440_paths):
    path=Path('data/sb441-n373s.bsp')
    if not path.is_file():pytest.skip('Official SB441373s kernel unavailable')
    spice.furnsh(str(de440_paths[0]));spice.furnsh(str(path))
    try:
        with DE440Ephemeris(*de440_paths) as e,SmallBodyPotential(e,path) as model:
            times=np.array([0.,(2460310.5-2451545.)*86400,(2300000.5-2451545.)*86400,(2399474.5-2451545.)*86400])
            d=model.diagnostics(times,by_body=True)
            assert model(times).shape==(4,2)
            for target in (2000001,2000004,2136199,2001036):
                row=next(r for r in model.mapping['entries'] if r['target']==target)
                for index,t in enumerate(times):
                    relative,_=spice.spkpos(str(target),float(t),'J2000','NONE','399')
                    expected=row['gm_si']/np.linalg.norm(relative*1000)
                    assert d['potential_by_body'][row['header_name']][index] == pytest.approx(expected,rel=3e-15)
            assert 0.4 < d['potential_asteroids'][0] < 1.
            assert .1 < d['potential_kbos'][0] < 1.
            np.testing.assert_allclose(model(times).sum(-1),np.sum(np.stack(list(d['potential_by_body'].values())),axis=0)/C_LIGHT**2,rtol=2e-15,atol=0)
            assert model.start_jd < 2443144.5 and model.end_jd > 2461041.5
    finally:
        spice.unload(str(path));spice.unload(str(de440_paths[0]))
