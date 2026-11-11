"""Independent static lunar multipole diagnostic, 1977 anchor to 2026.

Run: .venv/bin/python scripts/compare_phase2_lunar_multipole.py
The official DE440 PCK/FK must already exist; no toy orientation fallback.
"""
import datetime
import json
import time
from pathlib import Path
import _common
import numpy as np
import spiceypy as spice
from tt_tdb.constants import A_SCALE, T0_TDB
from tt_tdb.ephemeris import DE440Ephemeris
from tt_tdb.integration import SegmentedIntegrator
from tt_tdb.lunar_multipole import LunarMultipolePotential
from tt_tdb.reference import difference_metrics
from tt_tdb.small_body import file_sha256


def seconds_at_jan1(year):
    return (datetime.datetime(year, 1, 1)-datetime.datetime(2000, 1, 1, 12)).total_seconds()


def verify_orientation(model):
    paths = (model.ephemeris.path, model.orientation.path, model.fk_path)
    for path in paths:
        spice.furnsh(str(path))
    try:
        times = np.unique(np.r_[np.linspace(T0_TDB.seconds_since_j2000(), seconds_at_jan1(2026), 129), 717768000.])
        matrices = model.orientation.matrices(times)
        independent = np.stack([spice.pxform('J2000', 'MOON_PA_DE440', float(t)) for t in times])
        vector = model.displacement_pa(times)
        independent_vector = np.stack([spice.spkpos('399', float(t), 'MOON_PA_DE440', 'NONE', '301')[0]
                                       for t in times])*1000.
        norms = np.linalg.norm(vector, axis=-1)
        n = vector/norms[:, None]
        result = {'sample_count': len(times), 'frame_id': spice.namfrm('MOON_PA_DE440'),
                  'max_matrix_difference_from_cspice': float(np.max(abs(matrices-independent))),
                  'max_orthogonality_error': float(np.max(abs(matrices@np.swapaxes(matrices, -1, -2)-np.eye(3)))),
                  'max_abs_determinant_minus_one': float(np.max(abs(np.linalg.det(matrices)-1))),
                  'max_vector_difference_from_cspice_m': float(np.max(np.linalg.norm(vector-independent_vector, axis=-1))),
                  'earth_unit_direction_in_PA_min': n.min(axis=0).tolist(),
                  'earth_unit_direction_in_PA_max': n.max(axis=0).tolist(),
                  'earth_near_positive_long_axis_verified': bool(np.all(n[:, 0] > .97))}
        # The FK example is explicitly evaluated, independent of sample ordering.
        result['kernel_FK_example_position_at_717768000_tdb_seconds_km'] = (model.displacement_pa(717768000.)/1000.).tolist()
        if result['frame_id'] != 31008 or result['max_matrix_difference_from_cspice'] > 1e-10 or result['max_vector_difference_from_cspice_m'] > .02 or not result['earth_near_positive_long_axis_verified']:
            raise ValueError('Independent CSPICE PA orientation validation failed')
        return result
    finally:
        for path in reversed(paths):
            spice.unload(str(path))


def main():
    began = time.perf_counter()
    anchor = T0_TDB.seconds_since_j2000()
    modern = seconds_at_jan1(2024)+np.arange(732)*86400.
    history = np.r_[anchor, [seconds_at_jan1(y) for y in range(1978, 2027)]]
    times = np.unique(np.r_[history, modern])
    with DE440Ephemeris() as eph, LunarMultipolePotential(eph) as model:
        identity = model.fingerprint()
        identity['experiment_source_sha256'] = {str(p): file_sha256(p) for p in
            (Path(__file__), Path('tests/test_lunar_multipole.py'),
             Path('docs/sources/park_2021_de440_de441.pdf'),
             Path('docs/sources/moon_pa_de440_200625.cmt'))}
        orientation = verify_orientation(model)
        values, statistics = {}, {}
        for label, method, order, step in (
                ('gl24_halfday', 'gauss_legendre', 24, 43200.),
                ('romberg8_halfday', 'romberg', 8, 43200.),
                ('gl24_day', 'gauss_legendre', 24, 86400.)):
            integrator = SegmentedIntegrator(model, anchor=anchor, method=method, order=order,
                                            step_seconds=step, channels=5, fingerprint=identity,
                                            batch_segments=128)
            started = time.perf_counter()
            values[label] = A_SCALE*integrator.integral(times)
            statistics[label] = dict(integrator.stats, wall_seconds=time.perf_counter()-started,
                                     method=method, order=order, step_seconds=step)
            print(label, json.dumps(statistics[label]), flush=True)
        rate = A_SCALE*model(times)
        directions = model.displacement_pa(times)
        directions /= np.linalg.norm(directions, axis=-1)[:, None]
    selected = np.searchsorted(times, modern)
    baseline = values['gl24_halfday']
    primary = baseline[:, 0]+baseline[:, 3]
    primary_rate = rate[:, 0]+rate[:, 3]
    comparisons = {}
    for name in ('romberg8_halfday', 'gl24_day'):
        difference = baseline-values[name]
        comparisons[name] = {'max_abs_each_channel_seconds': abs(difference).max(axis=0).tolist(),
                             'max_abs_primary_sum_seconds': float(abs(difference[:, 0]+difference[:, 3]).max())}
    result = {'purpose': 'independent static lunar external U2/c2 diagnostic; default physical model unchanged',
              'command': '.venv/bin/python scripts/compare_phase2_lunar_multipole.py',
              'generated_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
              'anchor_tdb_j2000_seconds': anchor, 'anchor_definition': 'T0 + TDB0, IAU fixed original boundary',
              'modern_sample_count': len(modern), 'history_sample_count': len(history),
              'modern_TDB_dates': ['2024-01-01', '2026-01-01'],
              'history_TDB_dates': ['1977 anchor', '1978 through 2026 Jan 1 annually'],
              'signature': identity, 'A_SCALE': A_SCALE,
              'normalization_equation': 'U2=mu R²/r³[-J2(3nz²-1)/2+3C21 nx nz+3S21 ny nz+3C22(nx²-ny²)+6S22 nx ny]',
              'output_sign': 'add A_SCALE integral(U2/c²)dTDB to existing TDB-TT; no fitted correction',
              'summary_helper_linear_fit_fields_diagnostic_only': True,
              'primary_static_J2_C22_metrics_modern': difference_metrics(modern, primary[selected]),
              'primary_static_J2_C22_first_last_seconds': [float(primary[selected[0]]), float(primary[selected[-1]])],
              'primary_max_abs_rate_modern': float(abs(primary_rate[selected]).max()),
              'channels_modern': {name: {'time_metrics': difference_metrics(modern, baseline[selected, i]),
                                          'max_abs_rate': float(abs(rate[selected, i]).max())}
                                  for i, name in enumerate(LunarMultipolePotential.channel_names)},
              'literal_static_header_all_five_modern_metrics': difference_metrics(modern, baseline[selected].sum(axis=-1)),
              'orientation_validation': orientation, 'quadrature_comparisons_all_samples': comparisons,
              'integration_statistics': statistics, 'versions': _common.versions(),
              'reference_time_data_used_in_integrand': False,
              'reference_overlay_performed': False,
              'not_used_in_primary_acceptance': True,
              'model_limitation': 'J2/C22 are the static undistorted component; time-variable tide/rotation corrections omitted. C21/S21 header role unresolved, retained only as separate literal-static sensitivities. Not established as JPL time-ephemeris terms.',
              'wall_elapsed_seconds': time.perf_counter()-began}
    output = Path('results/tables/phase2_lunar_multipole_diagnostic.json')
    _common.write_json(output, result)
    names = ['tdb_j2000_seconds']+[name+'_delta_s' for name in LunarMultipolePotential.channel_names]
    names += ['primary_J2_C22_delta_s', 'literal_all_delta_s']
    names += [name+'_rate' for name in LunarMultipolePotential.channel_names]+['primary_rate', 'nx_pa', 'ny_pa', 'nz_pa']
    table = np.column_stack((times, baseline, primary, baseline.sum(axis=-1), rate, primary_rate, directions))
    _common.csv_table('results/tables/phase2_lunar_multipole_all_epochs.csv', names, table)
    _common.csv_table('results/tables/phase2_lunar_multipole_modern.csv', names, table[selected])
    directory = _common.make_plots_directory()
    import matplotlib.pyplot as plt
    from matplotlib import font_manager
    font_path=Path('/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc')
    if font_path.is_file():
        font_manager.fontManager.addfont(str(font_path))
        plt.rcParams['font.family']=font_manager.FontProperties(fname=str(font_path)).get_name()
    plt.rcParams['axes.unicode_minus']=False
    fig, axes = plt.subplots(2, 1, figsize=(8, 6), constrained_layout=True)
    years = 2000.+(times/86400.-.5)/365.25
    for i, label in ((0, '静态 J2'), (3, '静态 C22')):
        axes[0].plot(years, baseline[:, i]*1e12, label=label)
    axes[0].plot(years, primary*1e12, label='J2 + C22')
    axes[0].set(xlabel='TDB 年份（近似）', ylabel='对 TDB−TT 的锚定增量（ps）',
                title='月球静态二阶外势：DE440 系数与真实 PA 姿态')
    axes[0].legend()
    axes[1].plot((modern-modern[0])/86400, primary_rate[selected]/1e-22)
    axes[1].set(xlabel='自 2024-01-01 TDB 起的天数', ylabel='J2＋C22 速率（10⁻²²）',
                title='原始速率诊断（未施加拟合修正）')
    fig.savefig(directory/'phase2_lunar_multipole_diagnostic.png', dpi=160)
    print(json.dumps({'result': str(output), 'modern_primary_ps': (primary[selected[[0, -1]]]*1e12).tolist(),
                      'wall_seconds': result['wall_elapsed_seconds'], 'quadrature': comparisons}, indent=2))


if __name__ == '__main__':
    main()
