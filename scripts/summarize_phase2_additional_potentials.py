"""Combine separately verified candidate increments; never fit the reference."""
import json
from pathlib import Path
import sys
import _common
import numpy as np
from tt_tdb.solar_spin import sha256


ROOT = Path('results/tables')


def raw_metrics(values):
    return {'max_abs_seconds': float(np.max(np.abs(values))),
            'rms_seconds': float(np.sqrt(np.mean(values**2))),
            'first_seconds': float(values[0]), 'last_seconds': float(values[-1])}


def verify_sources(mapping, base=Path('.')):
    for name, digest in mapping.items():
        path = base/name
        if path.is_absolute():
            # Producer records __file__ absolutely; preserve clone portability.
            path = Path('scripts')/path.name
        if sha256(path) != digest:
            raise ValueError(f'Producer source or input changed: {path}')


def main():
    files = {n: ROOT/n for n in (
        'phase2_solar_spin_diagnostic.json', 'phase2_solar_spin_overlay.csv',
        'phase2_solar_spin_time_series.csv', 'phase2_lunar_multipole_diagnostic.json',
        'phase2_lunar_multipole_modern.csv')}
    spin = json.loads(files['phase2_solar_spin_diagnostic.json'].read_text())
    moon = json.loads(files['phase2_lunar_multipole_diagnostic.json'].read_text())
    ss, ms = spin['model_signature'], moon['signature']
    verify_sources(ss['sources'], Path('src/tt_tdb'))
    verify_sources(ms['source_sha256'], Path('src/tt_tdb'))
    verify_sources(ms['experiment_source_sha256'])
    verify_sources(spin['overlay']['verified_input_sha256'], ROOT)
    verify_sources({'scripts/compare_phase2_solar_spin.py': ss['script_sha256'],
                    'tests/test_solar_spin.py': ss['tests_sha256']})
    verify_sources({'data/de440.bsp': ss['bsp_sha256'],
                    'data/gm_de440.tpc': ss['gm_sha256'],
                    'docs/sources/de440_header.440': ms['header_sha256'],
                    'data/moon_pa_de440_200625.bpc': ms['pck_sha256'],
                    'data/moon_de440_250416.tf': ms['fk_sha256']})
    if ss['bsp_sha256'] != ms['earth_moon_bsp_sha256']:
        raise ValueError('DE440 geometry identities differ')
    if ss['gm_sha256'] != ms['gm_kernel_sha256'] or ss['header_sha256'] != ms['header_sha256']:
        raise ValueError('GM/header identities differ')
    if spin['anchor_seconds'] != moon['anchor_tdb_j2000_seconds']:
        raise ValueError('Original integration anchors differ')
    overlay = np.genfromtxt(files['phase2_solar_spin_overlay.csv'], delimiter=',', names=True)
    spin_table = np.genfromtxt(files['phase2_solar_spin_time_series.csv'], delimiter=',', names=True)
    lunar = np.genfromtxt(files['phase2_lunar_multipole_modern.csv'], delimiter=',', names=True)
    times = overlay['tdb_j2000_seconds']
    for table in (spin_table, lunar):
        if not np.array_equal(times, table['tdb_j2000_seconds']):
            raise ValueError('Diagnostic epoch grids differ')
    delta_spin = spin_table['combined_delta_s']
    delta_moon = lunar['primary_J2_C22_delta_s']
    columns = [times, delta_spin, delta_moon]
    names = ['tdb_j2000_seconds', 'solar_spin_delta_s', 'lunar_static_J2_C22_delta_s']
    results = {}
    for axis, label in (('solar_equator', 'solar_equator_J2'), ('ecliptic', 'literal_ecliptic_J2')):
        before = overlay[axis+'_before_s']
        variants = {'before': before, 'solar_spin_only': before+delta_spin,
                    'lunar_static_J2_C22_only': before+delta_moon,
                    'both_candidates': before+delta_spin+delta_moon}
        results[label] = {name: raw_metrics(value) for name, value in variants.items()}
        for name, value in variants.items():
            names.append(axis+'_'+name+'_residual_s')
            columns.append(value)
    result = {'command': ' '.join([sys.executable, *sys.argv]),
              'script_sha256': sha256(__file__),
              'classification': 'Mixed SB441 residual plus independent sourced candidate increments; no acceptance',
              'input_sha256': {n: sha256(p) for n, p in files.items()},
              'producer_sources_verified': True,
              'shared_ephemeris_GM_header_anchor_and_epochs_verified': True,
              'no_fitting_or_reanchoring': True, 'sample_count': len(times),
              'scope': 'Solar direct spin and static lunar J2+C22 only; no default model change',
              'residual_metrics': results}
    _common.write_json(ROOT/'phase2_additional_potentials_overlay.json', result)
    _common.csv_table(ROOT/'phase2_additional_potentials_overlay.csv', names, np.column_stack(columns))
    print(json.dumps(results, indent=2))


if __name__ == '__main__':
    main()
