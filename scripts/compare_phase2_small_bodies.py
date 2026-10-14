"""Separately labelled SB441 direct scalar diagnostic; no result fitting.

Read the measured Park candidate baseline after independent small-body
quadrature. Do not relabel the separately fitted SB441 orbits as DE440 states.
"""
import argparse
import hashlib
import json
from pathlib import Path
import time
import _common
import numpy as np
from tt_tdb.small_body import SmallBodyPotential
from tt_tdb.ephemeris import DE440Ephemeris
from tt_tdb.extended_ephemeris import ExtendedDE440Ephemeris
from tt_tdb.jpl_model import JPLModel
from tt_tdb.integration import SegmentedIntegrator
from tt_tdb.constants import A_SCALE,TDB0,T0_TDB,C_LIGHT
from tt_tdb.reference import difference_metrics
from tt_tdb.jpl_reference import JPLTimeEphemeris


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--no-cache',action='store_true');args=parser.parse_args()
    times=(2460310.5-2451545.)*86400+np.arange(732)*86400
    anchor=T0_TDB.seconds_since_j2000();started=time.perf_counter()
    with DE440Ephemeris() as e,SmallBodyPotential(e) as bodies:
        signature=bodies.fingerprint();methods={};integrals={}
        modern_rates=A_SCALE*bodies(times)
        for method,order in [('gauss_legendre',24),('romberg',6)]:
            path=Path(f'results/cache/phase2_sb441_{method}{order}.npz')
            config=dict(anchor=anchor,step_seconds=43200.,method=method,order=order,channels=2,fingerprint=signature,batch_segments=256)
            if path.is_file() and not args.no_cache:
                try:integ=SegmentedIntegrator.load(path,bodies,**config)
                except ValueError:integ=SegmentedIntegrator(bodies,**config)
            else:integ=SegmentedIntegrator(bodies,**config)
            began=time.perf_counter();value=A_SCALE*integ.integral(times);integ.save(path)
            integrals[method]=value
            methods[method]={'order':order,'step_seconds':43200.,'elapsed_wall_seconds':time.perf_counter()-began,
                     'statistics':integ.stats,'cache_fingerprint':integ.fingerprint,
                     'sum_error_indicators_seconds':integ.cumulative_error_diagnostic.tolist()}
            print(f'{method} done {methods[method]["elapsed_wall_seconds"]:.1f}s',flush=True)
        # A 1-day high-order integration is an additional separate convergence
        # diagnostic rather than an adjustment to the half-day result.
        gl1=SegmentedIntegrator(bodies,anchor=anchor,step_seconds=86400.,method='gauss_legendre',order=16,channels=2,fingerprint=signature,batch_segments=256)
        day_integral=A_SCALE*gl1.integral(times)
        main_integrals=integrals['gauss_legendre']
        raw=Path('results/tables/phase2_time_series.csv')
        primary=json.loads(Path('results/tables/phase2_jpl_comparison.json').read_text())
        current_sha=hashlib.sha256(Path('src/tt_tdb/jpl_model.py').read_bytes()).hexdigest()
        if primary['model_signature']['model']['sources']['jpl_model.py'] != current_sha:
            raise ValueError('Primary measured comparison does not match current Park model source')
        table=np.genfromtxt(raw,names=True,delimiter=',')
        if not np.array_equal(table['tdb_j2000_seconds'],times):raise ValueError('Primary and diagnostic epochs disagree')
        base=table['extended_delta_s'];result=base+main_integrals.sum(-1)
        with ExtendedDE440Ephemeris() as extended:
            park=JPLModel(extended)
            rate=A_SCALE*park(times).sum(-1)+modern_rates.sum(-1)
            # Model-specific conditional sensitivity of baseline mixed/squared
            # c^-4 terms if the missing direct U is inserted. Not used to
            # change the reported diagnostic model or acceptance result.
            d=park.diagnostics(times);added_u=modern_rates.sum(-1)/A_SCALE*C_LIGHT**2
            mixed_sensitivity=(1.5*d['velocity_squared']*added_u-d['potential']*added_u-.5*added_u**2)/C_LIGHT**4
        with JPLTimeEphemeris() as reference:
            delta_ref,rate_ref=reference.evaluate(times)
        residual=result-delta_ref;rate_residual=rate-rate_ref
        contribution=main_integrals
        document={'input_classification':'SB441 diagnostic: separate DE441-refitted heliocentric states, DE440-compatible header GMs and DE440 Earth-Sun geometry',
              'no_reference_fitting':True,'model_signature':signature,'versions':_common.versions(),
              'primary_result_sha256':hashlib.sha256(raw.read_bytes()).hexdigest(),
              'primary_model_signature':primary['model_signature'],
              'epoch_count':732,'start_jd':'2460310.5','end_jd':'2461041.5',
              'anchor_jd_parts':T0_TDB.jd_parts(),'anchor_delta_seconds':TDB0,
              'methods':methods,'time_residual':difference_metrics(times,residual),
              'rate_residual':difference_metrics(times,rate_residual,value_unit='dimensionless'),
              'acceptance_diagnostic':{'max_abs_raw_seconds':float(abs(residual).max()),'threshold_seconds':1e-12,
                    'passed_measured_threshold':bool(abs(residual).max()<1e-12),
                    'original_DE440_small_body_states_obtained':False,
                    'full_Park_all_body_PN_implemented':False},
              'quadrature_agreement':difference_metrics(times,(integrals['gauss_legendre']-integrals['romberg']).sum(-1)),
              'halfday_GL24_vs_day_GL16':difference_metrics(times,(main_integrals-day_integral).sum(-1)),
              'term_attribution':{name:{'max_abs_rate':float(abs(modern_rates[:,i]).max()),
                        'anchor_to_first_seconds':float(contribution[0,i]),'anchor_to_last_seconds':float(contribution[-1,i]),
                        'modern_peak_to_peak_seconds':float(np.ptp(contribution[:,i])),
                        'raw_residual_without_term':difference_metrics(times,residual-contribution[:,i])}
                         for i,name in enumerate(bodies.channel_names)},
              'omitted_mixed_squared_c4_sensitivity':{'max_abs_rate':float(abs(mixed_sensitivity).max()),
                         'scope':'Instantaneous 2024-2026 conditional expansion only; not an integrated PN error bound'},
              'elapsed_wall_seconds':time.perf_counter()-started}
    _common.write_json('results/tables/phase2_sb441_diagnostic.json',document)
    _common.csv_table('results/tables/phase2_sb441_time_series.csv',
       ['tdb_j2000_seconds','primary_delta_s','sb441_delta_s','jpl_delta_s','asteroids_delta_s','kbos_delta_s',
        'sb441_rate','jpl_rate','asteroids_rate','kbos_rate'],
       np.column_stack((times,base,result,delta_ref,contribution,rate,rate_ref,modern_rates)))
    import matplotlib.pyplot as plt
    fig,axs=plt.subplots(2,1,figsize=(10,7),sharex=True);days=(times-times[0])/86400
    axs[0].plot(days,residual*1e12,label='SB441 direct-potential diagnostic')
    axs[0].axhline(1,color='black',linestyle=':',linewidth=.7);axs[0].axhline(-1,color='black',linestyle=':',linewidth=.7)
    axs[0].set_ylabel('Raw time residual (ps)')
    axs[1].plot(days,rate_residual*1e22)
    axs[1].set_ylabel('Raw rate residual (10⁻²²)');axs[1].set_xlabel('Days since 2024-01-01 TDB')
    for ax in axs:ax.grid(alpha=.3)
    axs[0].legend();fig.tight_layout();fig.savefig(_common.make_plots_directory()/'phase2_sb441_diagnostic.png',dpi=170);plt.close(fig)
    print(json.dumps({k:document[k] for k in ('acceptance_diagnostic','time_residual','rate_residual','quadrature_agreement')},indent=2),flush=True)

if __name__=='__main__':main()
