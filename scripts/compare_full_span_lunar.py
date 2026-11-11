"""Original-anchor moon-center TCL−TDB vs fixed published LTE440 product.

Run --modern-only first to check convergence/performance, then the default
full common coverage. Physical scalar/PN and SB441 leading channels are
independently integrated; the reference is read only after integration.
"""
import argparse
import gzip
import inspect
import json
import time
from pathlib import Path
import _common
import numpy as np
from tt_tdb.constants import L_B, T0_TDB
from tt_tdb.integration import _quadrature_intervals
from tt_tdb.lunar_model import LunarModel, MoonObserverEphemeris, LunarSmallBodyPotential
from tt_tdb.lte_reference import LTE440Reference
from tt_tdb.small_body import file_sha256


def raw_metrics(t, residual, *, unit='seconds'):
    x = np.asarray(residual, dtype=np.longdouble)
    return {'count': int(x.size), 'max_abs_'+unit: float(abs(x).max()),
            'rms_'+unit: float(np.sqrt(np.mean(x*x))), 'mean_'+unit: float(np.mean(x)),
            'min_'+unit: float(x.min()), 'max_'+unit: float(x.max()),
            'max_epoch_tdb_seconds': float(np.asarray(t)[np.argmax(abs(x))])}


def anchored_integral(rate, times, anchor, signature, name, method, order, channels,
                      *, subdivisions=1, use_cache=True):
    """Batch actual quadrature increments, then accumulate from the anchor.

    No physical rate interpolation. No left-end total is subtracted to form
    an anchor value. Long-double cumsum is done separately in both directions.
    """
    import hashlib
    identity = dict(model=signature, method=method, order=order, subdivisions=subdivisions,
                    times_sha256=hashlib.sha256(times.tobytes()).hexdigest(),
                    integration_driver_sha256=hashlib.sha256(inspect.getsource(anchored_integral).encode()).hexdigest(),
                    accumulation='anchor_outward_longdouble_cumsum_v1')
    encoded = json.dumps(identity, sort_keys=True)
    cache = Path('results/cache')/(name+'.npz')
    began = time.perf_counter()
    cached = False
    if use_cache and cache.is_file():
        with np.load(cache, allow_pickle=False) as z:
            if str(z['identity']) == encoded:
                increments = z['increments'].copy()
                stats = json.loads(str(z['statistics']))
                cached = True
    if not cached:
        increments = np.zeros((len(times)-1, channels))
        error_sum = np.zeros(channels)
        evaluations = 0
        batch = 256
        for lo in range(0, len(times)-1, batch):
            hi = min(lo+batch, len(times)-1)
            starts, ends = times[lo:hi], times[lo+1:hi+1]
            for split in range(subdivisions):
                a = starts+(ends-starts)*split/subdivisions
                b = starts+(ends-starts)*(split+1)/subdivisions
                answer = _quadrature_intervals(rate, a, b, method, order, channels)
                increments[lo:hi] += answer.value
                error_sum += answer.error.sum(axis=0)
                evaluations += answer.evaluations
            if lo and lo % (batch*200) == 0:
                print(f'{name}: {hi}/{len(times)-1} intervals, {time.perf_counter()-began:.1f}s', flush=True)
        stats = {'method': method, 'order': order, 'subdivisions_per_output_interval': subdivisions,
                 'interval_count': len(times)-1, 'evaluations': evaluations,
                 'sum_local_order_difference_indicators_seconds': error_sum.tolist(),
                 'integration_wall_seconds': time.perf_counter()-began}
        cache.parent.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(cache, identity=encoded, increments=increments,
                            statistics=json.dumps(stats))
    index = np.flatnonzero(times == anchor)
    if len(index) != 1:
        raise ValueError('Exactly one actual original anchor must be in the integration grid')
    i = int(index[0])
    values = np.zeros((len(times), channels), dtype=np.longdouble)
    values[i+1:] = np.cumsum(increments[i:].astype(np.longdouble), axis=0, dtype=np.longdouble)
    values[:i] = -np.cumsum(increments[:i][::-1].astype(np.longdouble), axis=0, dtype=np.longdouble)[::-1]
    stats = dict(stats, loaded_cache=cached, call_wall_seconds=time.perf_counter()-began,
                 anchor_index=i, cumulative_directions='independent forward and backward from original anchor',
                 cache_identity_sha256=hashlib.sha256(encoded.encode()).hexdigest(),
                 longdouble_mantissa_bits=int(np.finfo(np.longdouble).nmant))
    print(name, json.dumps(stats), flush=True)
    return values, stats


def convergence(model, bodies):
    modern = 757339200.+np.arange(733)*86400.
    a, b = modern[:-1], modern[1:]
    answer = {}
    for name, f, channels, configs in (
            ('major', model, 11, [('gauss_legendre', 8, 1), ('romberg', 5, 1), ('gauss_legendre', 12, 2)]),
            ('SB441', bodies, 2, [('gauss_legendre', 4, 1), ('romberg', 4, 1), ('gauss_legendre', 12, 2)])):
        vals, wall = [], []
        for method, order, n in configs:
            began = time.perf_counter()
            inc = np.zeros((len(a), channels))
            for split in range(n):
                inc += _quadrature_intervals(f, a+(b-a)*split/n,
                                            a+(b-a)*(split+1)/n,
                                            method, order, channels).value
            vals.append(np.cumsum(inc.astype(np.longdouble), axis=0))
            wall.append(time.perf_counter()-began)
        answer[name] = {'scope': 'modern 732 intervals, computational convergence anchor 2024-01-01 only; not physical origin',
                        'configs': configs, 'wall_seconds': wall,
                        'max_primary_minus_highorder_seconds': float(abs((vals[0]-vals[2]).sum(axis=-1)).max()),
                        'max_independent_minus_highorder_seconds': float(abs((vals[1]-vals[2]).sum(axis=-1)).max())}
    return answer


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--modern-only', action='store_true')
    parser.add_argument('--no-cache', action='store_true')
    args = parser.parse_args()
    start_wall = time.perf_counter()
    anchor = T0_TDB.seconds_since_j2000()
    with MoonObserverEphemeris() as eph, LunarSmallBodyPotential(eph) as bodies:
        model = LunarModel(eph)
        signature, sb_signature = model.fingerprint(), bodies.fingerprint()
        conv = convergence(model, bodies)
        print('convergence', json.dumps(conv), flush=True)
        # Product metadata fixes the comparison domain, not the integrand.
        with LTE440Reference() as reference:
            first_jd = max(eph.start_jd, bodies.start_jd, reference.start_jd)
            last_jd = min(eph.end_jd, bodies.end_jd, reference.end_jd)
        first, last = (first_jd-2451545.)*86400., (last_jd-2451545.)*86400.
        modern = 757339200.+np.arange(732)*86400.
        if args.modern_only:
            # History must still be fully integrated from the true 1977 origin.
            grid = np.arange(np.floor(anchor/86400.)*86400., modern[-1]+1., 86400.)
        else:
            grid = first+np.arange(int(round((last-first)/86400.))+1)*86400.
        times = np.unique(np.r_[grid, anchor, modern])
        mode = 'modern' if args.modern_only else 'full'
        major, major_stats = anchored_integral(model, times, anchor, signature,
                f'full_span_lunar_{mode}_major_gl8', 'gauss_legendre', 8, 11, use_cache=not args.no_cache)
        small, small_stats = anchored_integral(bodies, times, anchor, sb_signature,
                f'full_span_lunar_{mode}_SB441_gl4', 'gauss_legendre', 4, 2, use_cache=not args.no_cache)
        independent_major, independent_major_stats = anchored_integral(model, times, anchor, signature,
                f'full_span_lunar_{mode}_major_romberg5', 'romberg', 5, 11, use_cache=not args.no_cache)
        independent_small, independent_small_stats = anchored_integral(bodies, times, anchor, sb_signature,
                f'full_span_lunar_{mode}_SB441_romberg4', 'romberg', 4, 2, use_cache=not args.no_cache)
        value = model.restore(times, major, small)
        independent = model.restore(times, independent_major, independent_small)
        rate_parts = np.empty((len(times), 11))
        sb_rate = np.empty((len(times), 2))
        for lo in range(0, len(times), 4096):
            hi = min(lo+4096, len(times))
            rate_parts[lo:hi] = model(times[lo:hi])
            sb_rate[lo:hi] = bodies(times[lo:hi])
        rate = ((L_B-model.centering_rate)-rate_parts[:, :10].sum(-1)-sb_rate.sum(-1))/(1-L_B)
        # Independent product values are read only after physical integration.
        with LTE440Reference() as reference:
            ref, ref_rate = reference.evaluate(times)
            ref_signature = reference.fingerprint()
        residual, rate_residual = value-ref, rate-ref_rate
        modern_indices = np.searchsorted(times, modern)
        ai = int(np.flatnonzero(times == anchor)[0])
        result = {'mode': mode, 'quantity': 'Moon-center TCL-TDB, seconds; no TL/proper-time mapping',
                  'command': '.venv/bin/python scripts/compare_full_span_lunar.py'+(' --modern-only' if args.modern_only else ''),
                  'sample_count': len(times), 'daily_grid_sample_count': len(grid),
                  'common_coverage_jd': [first_jd, last_jd], 'actual_sample_tdb_seconds': [float(times[0]), float(times[-1])],
                  'anchor_tdb_seconds': anchor, 'model_anchor_seconds': float(value[ai]),
                  'reference_anchor_seconds': float(ref[ai]),
                  'model_signature': signature, 'SB441_signature': sb_signature,
                  'reference_signature': ref_signature,
                  'script_sha256': file_sha256(__file__),
                  'test_sha256': {str(p): file_sha256(p) for p in map(Path, ['tests/test_lunar_model.py', 'tests/test_lte_reference.py'])},
                  'time_residual_raw': raw_metrics(times, residual),
                  'rate_residual_raw': raw_metrics(times, rate_residual, unit='dimensionless'),
                  'modern_time_residual_raw': raw_metrics(modern, residual[modern_indices]),
                  'modern_rate_residual_raw': raw_metrics(modern, rate_residual[modern_indices], unit='dimensionless'),
                  'independent_quadrature_difference': raw_metrics(times, value-independent),
                  'independent_major_integral_sum_max_seconds': float(abs((major-independent_major)[:, :10].sum(-1)).max()),
                  'independent_SB441_integral_sum_max_seconds': float(abs((small-independent_small).sum(-1)).max()),
                  'Earth_J2_zero_offset_minus_primary_max_seconds': float(abs(major[:, 10]/(1-L_B)).max()),
                  'modern_convergence': conv,
                  'methods': {'major_primary': major_stats, 'SB441_primary': small_stats,
                              'major_independent': independent_major_stats, 'SB441_independent': independent_small_stats},
                  'linear_or_zero_fit_applied': False, 'reference_used_in_integrand': False,
                  'acceptance': {'threshold_seconds': 1e-12, 'passed_raw_sampled_time': bool(abs(residual).max() < 1e-12),
                                 'scope': 'daily samples plus exact anchor and endpoints, not continuous-domain guarantee'},
                  'limitations': signature['excluded']+['SB441 orbits are separately fitted with DE441; original DE440 small-body states not obtained',
                                                       'LTE440 release predates paper v1; production source/input manifest unavailable'],
                  'versions': _common.versions(), 'wall_elapsed_seconds_before_plot': time.perf_counter()-start_wall}
        prefix = f'results/tables/full_span_lunar_{mode}'
        _common.write_json(prefix+'_diagnostic.json', result)
        output = prefix+'_daily.csv.gz'
        names = ['tdb_j2000_seconds', 'model_TCL_minus_TDB_s', 'LTE440_TCL_minus_TDB_s',
                 'raw_time_residual_s', 'model_rate', 'LTE440_rate', 'raw_rate_residual',
                 'independent_quadrature_difference_s', 'Earth_J2_contribution_s', 'SB441_contribution_s',
                 'Earth_J2_zero_offset_minus_primary_s']
        table = np.column_stack((times, value, ref, residual, rate, ref_rate, rate_residual,
                                 value-independent, -major[:, 9]/(1-L_B), -small.sum(-1)/(1-L_B),
                                 -major[:, 10]/(1-L_B)))
        with gzip.open(output, 'wt') as stream:
            np.savetxt(stream, table, delimiter=',', header=','.join(names), comments='', fmt='%.20e')
        sample = np.unique(np.r_[0, len(times)-1, ai, modern_indices, np.arange(0, len(times), 3652)]).astype(int)
        _common.csv_table(prefix+'_selected_epochs.csv', names, table[sample])
        result['daily_table_sha256'] = file_sha256(output)
        _common.write_json(prefix+'_diagnostic.json', result)
        plot(times, residual, rate_residual, value-independent, mode)
        print(json.dumps({k: result[k] for k in ('sample_count', 'time_residual_raw', 'modern_time_residual_raw',
                                                'independent_quadrature_difference', 'wall_elapsed_seconds_before_plot')}, indent=2), flush=True)


def plot(times, residual, rate_residual, quadrature, mode):
    directory = _common.make_plots_directory()
    import matplotlib.pyplot as plt
    from matplotlib import font_manager
    font = Path('/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc')
    if not font.is_file():
        raise FileNotFoundError('Chinese plotting font unavailable')
    font_manager.fontManager.addfont(str(font))
    plt.rcParams['font.family'] = font_manager.FontProperties(fname=str(font)).get_name()
    plt.rcParams['axes.unicode_minus'] = False
    years = 2000.+times/(365.25*86400.)
    fig, axes = plt.subplots(3, 1, figsize=(11, 8), sharex=True, constrained_layout=True)
    axes[0].plot(years, rate_residual/1e-21, lw=.55)
    axes[0].set_ylabel(r'原始速率残差 / $10^{-21}$')
    axes[1].plot(years, np.asarray(residual, float)*1e12, lw=.55)
    axes[1].set_ylabel('TCL−TDB 原始残差 / ps')
    axes[2].plot(years, np.asarray(quadrature, float)*1e12, lw=.55)
    axes[2].set_ylabel('独立求积差 / ps')
    axes[2].set_xlabel('TDB 年份（近似）')
    fig.suptitle('月心 TCL−TDB 对 LTE440 发布产品：1977 固定锚点，未拟合零点或趋势')
    for ax in axes:
        ax.grid(alpha=.25)
    fig.savefig(directory/f'full_span_lunar_{mode}_raw_residuals.png', dpi=170)
    plt.close(fig)


if __name__ == '__main__':
    main()
