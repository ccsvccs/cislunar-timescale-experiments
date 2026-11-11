"""Integrate sourced direct solar spin from the fixed 1977 anchor, without fitting."""
import json
from pathlib import Path
import sys
import time
import _common
import numpy as np
import spiceypy as spice
from tt_tdb.constants import A_SCALE, C_LIGHT, L_B, T0_TDB
from tt_tdb.extended_ephemeris import ExtendedDE440Ephemeris
from tt_tdb.integration import SegmentedIntegrator
from tt_tdb.jpl_model import JPLModel, read_de440_constants, solar_j2_potential
from tt_tdb.jpl_reference import JPLTimeEphemeris
from tt_tdb.solar_spin import SolarSpinDiagnostic, solar_axis, sha256, G_STANDARD_UNCERTAINTY, G_CODATA2022


PREFIX='results/tables/phase2_solar_spin_'


def raw_metrics(values):
    """No fitted zero point, trend or fitted summary statistics."""
    values=np.asarray(values, float)
    return {'count':int(values.size), 'min':float(values.min()), 'max':float(values.max()),
            'max_abs':float(abs(values).max()), 'mean':float(values.mean()),
            'rms':float(np.sqrt(np.mean(values**2))), 'peak_to_peak':float(np.ptp(values))}


def checked_existing_overlay(eph, times, delta_spin, rate_spin, axis_integral):
    """Validate archived inputs before adding a diagnostic to existing residuals."""
    root=Path('results/tables')
    paths={n:root/n for n in ('phase2_sb441_diagnostic.json','phase2_sb441_time_series.csv',
                             'phase2_solar_axis_diagnostic.json','phase2_solar_axis_diagnostic.csv',
                             'phase2_jpl_comparison.json','phase2_time_series.csv')}
    docs={n:json.loads(p.read_text()) for n,p in paths.items() if p.suffix=='.json'}
    sb=docs['phase2_sb441_diagnostic.json']; primary=docs['phase2_jpl_comparison.json']
    axis=docs['phase2_solar_axis_diagnostic.json']
    if sb['primary_model_signature'] != primary['model_signature']:
        raise ValueError('SB441 and primary model signatures differ')
    if json.loads(JPLModel(eph).fingerprint()) != primary['model_signature']['model']:
        raise ValueError('Archived primary model signature differs from current complete Park signature')
    if sha256(paths['phase2_time_series.csv']) != sb['primary_result_sha256']:
        raise ValueError('Primary CSV SHA differs from SB441 provenance')
    sbsig=sb['model_signature']
    for name, digest in sbsig['sources'].items():
        if sha256(Path('src/tt_tdb')/name) != digest:
            raise ValueError('Archived SB441 source differs: '+name)
    if sha256('data/sb441-n373s.bsp') != sbsig['bsp_sha256']:
        raise ValueError('Archived SB441 BSP identity differs')
    if sha256(eph.path) != sbsig['earth_bsp_sha256'] or sha256(eph.path) != axis['signature']['bsp_sha256']:
        raise ValueError('Archived diagnostics differ from DE440 geometry')
    header_sha=sha256('docs/sources/de440_header.440')
    if header_sha != sbsig['mapping']['header_sha256'] or header_sha != axis['signature']['header_sha256']:
        raise ValueError('Archived diagnostics differ from DE440 header')
    if axis['anchor_seconds'] != T0_TDB.seconds_since_j2000() or primary['anchor_j2000_seconds'] != axis['anchor_seconds']:
        raise ValueError('Archived anchors differ')
    if sb['anchor_jd_parts'] != list(T0_TDB.jd_parts()):
        raise ValueError('Archived two-part SB441 anchor differs')
    np.testing.assert_array_equal(axis['signature']['axes'][0], solar_axis())
    np.testing.assert_array_equal(axis['signature']['axes'][1], spice.pxform('ECLIPJ2000','J2000',0.)[:,2])
    tables={n:np.genfromtxt(p,delimiter=',',names=True) for n,p in paths.items() if p.suffix=='.csv'}
    for table in tables.values():
        np.testing.assert_array_equal(table['tdb_j2000_seconds'], times)
    sbtable=tables['phase2_sb441_time_series.csv']; ptable=tables['phase2_time_series.csv']
    atable=tables['phase2_solar_axis_diagnostic.csv']
    np.testing.assert_array_equal(sbtable['primary_delta_s'], ptable['extended_delta_s'])
    np.testing.assert_array_equal(sbtable['jpl_delta_s'], ptable['jpl_delta_s'])
    np.testing.assert_allclose(sbtable['sb441_delta_s'], sbtable['primary_delta_s']+
                               sbtable['asteroids_delta_s']+sbtable['kbos_delta_s'],rtol=0,atol=3e-19)
    axis_difference=atable['ecliptic_J2_delta_s']-atable['spin_J2_delta_s']
    axis_closure=float(abs(axis_difference-axis_integral).max())
    if axis_closure > 1e-25:
        raise ValueError('Fresh anchor-based axis integral differs from archived axis CSV')
    with JPLTimeEphemeris() as reference:
        if sha256(reference.path) != primary['reference']['sha256']:
            raise ValueError('Reference BSP identity differs')
        ref_delta,ref_rate=reference.evaluate(times)
    np.testing.assert_array_equal(ref_delta, sbtable['jpl_delta_s'])
    np.testing.assert_array_equal(ref_rate, sbtable['jpl_rate'])
    original=sbtable['sb441_delta_s']-ref_delta
    original_rate=sbtable['sb441_rate']-ref_rate
    ecliptic=original+axis_difference
    ecliptic_rate=original_rate+atable['ecliptic_J2_rate']-atable['spin_J2_rate']
    outputs={}
    for name, residual, rr in (('solar_equator_J2',original,original_rate),('literal_ecliptic_J2',ecliptic,ecliptic_rate)):
        outputs[name]={'before_time_seconds':raw_metrics(residual),
                       'after_time_seconds':raw_metrics(residual+delta_spin),
                       'before_rate_dimensionless':raw_metrics(rr),
                       'after_rate_dimensionless':raw_metrics(rr+rate_spin)}
    values=np.column_stack((times,original,original+delta_spin,ecliptic,ecliptic+delta_spin,
                            original_rate,original_rate+rate_spin,ecliptic_rate,ecliptic_rate+rate_spin))
    _common.csv_table(PREFIX+'overlay.csv', ['tdb_j2000_seconds','solar_equator_before_s','solar_equator_plus_spin_s',
                      'ecliptic_before_s','ecliptic_plus_spin_s','solar_equator_before_rate',
                      'solar_equator_plus_spin_rate','ecliptic_before_rate','ecliptic_plus_spin_rate'],values)
    return {'classification':'Previously measured mixed SB441 diagnostic plus independently integrated solar spin; no model acceptance',
            'verified_input_sha256':{n:sha256(p) for n,p in paths.items()},
            'fresh_axis_integral_closure_max_seconds':axis_closure,
            'current_Park_signature_matches_archived_complete_signature':True,
            'reference_rechecked_at_all_epochs':True, 'residuals':outputs}, values


def main():
    started=time.perf_counter()
    times=757339200.+np.arange(732)*86400.
    anchor=T0_TDB.seconds_since_j2000()
    # Exact anchor, 1980/1990/2000/2010/2020-01-01 00:00 TDB.
    historical_jd=np.array([2444239.5,2447892.5,2451544.5,2455197.5,2458849.5])
    history=np.r_[anchor,(historical_jd-2451545.)*86400.]
    queries=np.r_[history,times]
    with ExtendedDE440Ephemeris() as eph:
        diagnostic=SolarSpinDiagnostic(eph)
        header=read_de440_constants('docs/sources/de440_header.440')
        axes=(solar_axis(),spice.pxform('ECLIPJ2000','J2000',0.)[:,2])
        def rates(t):
            spin=diagnostic(t)
            state=eph.states(t)
            j2=[solar_j2_potential(state.earth_position,state.body_positions['sun'],state.gm['sun'],
                    header['J2SUN'],header['ASUN']*1000,axis)/C_LIGHT**2 for axis in axes]
            return np.column_stack((spin,j2[1]-j2[0]))
        signature=diagnostic.fingerprint()
        signature.update(script_sha256=sha256(__file__), tests_sha256=sha256('tests/test_solar_spin.py'),
                         header_sha256=sha256('docs/sources/de440_header.440'),
                         axis_integral_verification_source_sha256=sha256('src/tt_tdb/jpl_model.py'),
                         sixth_channel='Independent ecliptic-minus-equator J2 integral solely to validate archived overlay input')
        values={};methods={}
        for label,method,order,step in (('GL24_halfday','gauss_legendre',24,43200.),
                                        ('Romberg6_halfday','romberg',6,43200.),
                                        ('GL16_day','gauss_legendre',16,86400.)):
            began=time.perf_counter()
            integ=SegmentedIntegrator(rates, anchor=anchor, method=method, order=order,
                    step_seconds=step,channels=6,fingerprint=signature,batch_segments=256)
            values[label]=A_SCALE*integ.integral(queries)
            methods[label]={'method':method,'order':order,'step_seconds':step,
                            'elapsed_wall_seconds':time.perf_counter()-began, 'statistics':integ.stats,
                            'fingerprint':integ.fingerprint,
                            'sum_segment_error_indicators_seconds':(A_SCALE*integ.cumulative_error_diagnostic).tolist()}
            print(label,methods[label]['elapsed_wall_seconds'],flush=True)
        rate_values=A_SCALE*rates(queries)
        main=values['GL24_halfday']; modern=main[len(history):]; mr=rate_values[len(history):]
        overlay,overlay_values=checked_existing_overlay(eph,times,modern[:,2],mr[:,2],modern[:,5])
        # Sample the complete anchor-to-end interval daily for an observed rate
        # envelope. This is not a continuous supremum nor a physical error bound.
        daily=np.r_[anchor+np.arange(int((times[-1]-anchor)//86400)+1)*86400,times[-1]]
        daily_rates=A_SCALE*diagnostic(daily)
    names=diagnostic.channel_names
    fields=['tdb_j2000_seconds']+[n+'_delta_s' for n in names]+[n+'_rate' for n in names]
    _common.csv_table(PREFIX+'time_series.csv', fields, np.column_stack((times,modern[:,:5],mr[:,:5])))
    _common.csv_table(PREFIX+'history.csv', fields, np.column_stack((history,main[:len(history),:5],rate_values[:len(history),:5])))
    channels={n:{'modern_time_seconds':raw_metrics(modern[:,i]),'modern_rate_dimensionless':raw_metrics(mr[:,i]),
                 'anchor_to_first_seconds':float(modern[0,i]),'anchor_to_last_seconds':float(modern[-1,i]),
                 'daily_anchor_to_end_rate_dimensionless':raw_metrics(daily_rates[:,i])} for i,n in enumerate(names)}
    result={'command':' '.join([sys.executable,*sys.argv]),'versions':_common.versions(),'model_signature':signature,
            'sign':'additive correction to Delta=TDB-TT, seconds; rates dDelta/dTDB',
            'source_parameters_chosen_before_reference_comparison':True, 'primary_models_unchanged':True,
            'no_fitting_or_reanchoring':True, 'DE440_time_model_spin_inclusion':'unknown; Park Eq3-7 do not include direct spin',
            'anchor_jd_parts':T0_TDB.jd_parts(),'anchor_seconds':anchor,'anchor_correction_seconds':main[0,:5].tolist(),
            'modern_epoch_count':len(times),'modern_start_jd':2460310.5,'modern_end_jd':2461041.5,
            'historical_jd':historical_jd.tolist(),'daily_rate_sample_count':len(daily),'methods':methods,
            'channels':channels,'quadrature_max_abs_difference_each_channel_seconds':
                {n:abs(main-v).max(axis=0).tolist() for n,v in values.items() if n!='GL24_halfday'},
            'closure':{'instantaneous_combined_minus_relative_max_abs':float(abs(rate_values[:,2]-rate_values[:,3]).max()),
                       'integrated_combined_minus_vector_plus_scalar_max_seconds':float(abs(main[:,2]-main[:,0]-main[:,1]).max()),
                       'integrated_combined_minus_relative_max_seconds':float(abs(main[:,2]-main[:,3]).max())},
            'GS_scale_sensitivity_max_seconds':float(abs(main[:,4]-main[:,2]).max()),
            'G_standard_uncertainty_only_max_seconds':float(abs(main[:,2]).max()*G_STANDARD_UNCERTAINTY/G_CODATA2022),
            'uncertainty_limit':'G-only sensitivity is not a total uncertainty; source S, fixed-axis and omitted-physics uncertainties not supplied',
            'overlay':overlay,'elapsed_wall_seconds':time.perf_counter()-started}
    _common.write_json(PREFIX+'diagnostic.json',result)
    make_plots(history,main,times,modern,mr,overlay_values)
    print(json.dumps({k:result[k] for k in ('channels','closure','quadrature_max_abs_difference_each_channel_seconds',
                                           'GS_scale_sensitivity_max_seconds','elapsed_wall_seconds')},indent=2))


def make_plots(history,main,times,modern,rates,overlay):
    directory=_common.make_plots_directory()
    import matplotlib.pyplot as plt
    from matplotlib import font_manager
    font_path=Path('/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc')
    if font_path.is_file():
        font_manager.fontManager.addfont(str(font_path))
        plt.rcParams['font.family']=font_manager.FontProperties(fname=str(font_path)).get_name()
    plt.rcParams['axes.unicode_minus']=False
    years=2000.+np.r_[history,times]/(365.25*86400)
    fig,ax=plt.subplots(figsize=(10,5))
    for i,label in enumerate(('向量自旋项','标量自旋项','合并自旋项')):
        ax.plot(years,main[:,i]*1e12,label=label)
    ax.set(xlabel='TDB 年份（近似）',ylabel='从 1977 固定锚点累计的修正（ps）',
           title='太阳直接自旋势的独立诊断（主模型不含此项）')
    ax.grid(alpha=.3);ax.legend();fig.tight_layout()
    fig.savefig(directory/'phase2_solar_spin_history.png',dpi=170);plt.close(fig)
    days=(times-times[0])/86400
    fig,axs=plt.subplots(3,1,figsize=(10,10),sharex=True)
    for i,label in enumerate(('向量项','标量项','合并项')):
        axs[0].plot(days,modern[:,i]*1e12,label=label)
        axs[1].plot(days,rates[:,i]*1e21,label=label)
    for i,label in ((1,'太阳赤道 J2 原结果'),(2,'太阳赤道 J2＋自旋'),
                    (3,'黄道 J2 原结果'),(4,'黄道 J2＋自旋')):
        axs[2].plot(days,overlay[:,i]*1e12,label=label)
    axs[0].set(ylabel='锚定累计修正（ps）',title='太阳自旋的现代区间诊断（不进入主模型）')
    axs[1].set_ylabel('速率修正（10$^{-21}$）')
    axs[2].set(ylabel='相对 DE440t 原始残差（ps）',xlabel='自 2024-01-01 TDB 起的天数')
    for ax in axs:ax.grid(alpha=.3);ax.legend(fontsize=8)
    fig.tight_layout();fig.savefig(directory/'phase2_solar_spin_modern.png',dpi=170);plt.close(fig)


if __name__=='__main__':
    main()
