"""Reproduce Eq. (21) from its 1977 anchor with two independent quadratures."""
import argparse
from pathlib import Path
import hashlib
import json
import subprocess
import sys
import time
import _common
import numpy as np
from tt_tdb.constants import A_SCALE,C_EXACT,K_LINEAR,L_C_PAPER,TDB0,T0,T0_TDB,decimal_constants
from tt_tdb.ephemeris import DE440Ephemeris
from tt_tdb.pn_rate import PNModel
from tt_tdb.integration import SegmentedIntegrator
from tt_tdb.transforms import TTTDBConverter
from tt_tdb.timebase import TwoPartJD,time_precision_experiment
from tt_tdb.reference import sofa_geocentric,difference_metrics
from tt_tdb.jpl_reference import JPLTimeEphemeris


def get_integrator(rate,signature,method,order,channels,anchor,args,label):
    path=Path('results/cache')/(label+'.npz')
    parameters=dict(fingerprint=signature,anchor=anchor,method=method,order=order,
                    channels=channels,step_seconds=43200,batch_segments=64)
    if path.exists() and not args.no_cache:
        try:
            integrator=SegmentedIntegrator.load(path,rate,**parameters)
            print(f'Validated checkpoint cache: {label}',flush=True)
            return integrator,path
        except ValueError as error:
            print(f'Rebuilding incompatible cache {label}: {error}',flush=True)
    return SegmentedIntegrator(rate,**parameters),path


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--start-jd',default='2460310.5',help='Default 2024-01-01 TDB')
    parser.add_argument('--end-jd',default='2461041.5',help='Default 2026-01-01 TDB')
    parser.add_argument('--no-cache',action='store_true')
    parser.add_argument('--skip-short-validation',action='store_true')
    args=parser.parse_args()
    began=time.perf_counter()
    if not args.skip_short_validation:
        print('First perform genuine short-interval convergence study',flush=True)
        subprocess.run([sys.executable,'scripts/convergence_study.py'],check=True)
    start=TwoPartJD.from_decimal_jd(args.start_jd).seconds_since_j2000()
    end=TwoPartJD.from_decimal_jd(args.end_jd).seconds_since_j2000()
    if end<=start:
        raise ValueError('Output end must follow start')
    times=np.arange(start,end+1,86400.)
    if times[-1]!=end:
        times=np.append(times,end)
    anchor=T0_TDB.seconds_since_j2000()
    with DE440Ephemeris() as eph:
        model=PNModel(eph)
        signature=json.dumps({'model':json.loads(model.fingerprint()),
                              'integrator_source_sha256':hashlib.sha256(Path('src/tt_tdb/integration.py').read_bytes()).hexdigest(),
                              'versions':_common.versions()},sort_keys=True)
        romberg,rpath=get_integrator(model,signature,'romberg',10,2,anchor,args,'romberg10_halfday')
        gauss,gpath=get_integrator(model,signature,'gauss_legendre',16,2,anchor,args,'gauss16_halfday')
        query=np.concatenate(([anchor,0.],times))
        print(f'Full integration: anchor {anchor:.10f} s → end {end:.1f} s',flush=True)
        channels=romberg.integral(query)
        romberg.save(rpath)
        print(f'Romberg complete: {romberg.stats}',flush=True)
        independent=gauss.integral(query)
        gauss.save(gpath)
        print(f'Independent Gauss complete: {gauss.stats}',flush=True)
        full=TDB0+A_SCALE*np.sum(channels,axis=-1)
        c2=TDB0+A_SCALE*channels[:,0]
        full_gauss=TDB0+A_SCALE*np.sum(independent,axis=-1)
        diagnostics=model.diagnostics(times)
        reference=sofa_geocentric(2451545.,times/86400.)
        try:
            with JPLTimeEphemeris() as jpl:
                jpl_delta,jpl_rate=jpl.evaluate(times)
                jpl_origin,jpl_origin_rate=jpl.evaluate(anchor)
                jpl_j2000,jpl_j2000_rate=jpl.evaluate(0.)
                jpl_metadata={'available':True,'path':str(jpl.path),'segment':[1000000000,1000000001],
                              'anchor_delta_seconds':jpl_origin,'anchor_model_minus_jpl_seconds':TDB0-jpl_origin,
                              'j2000_delta_seconds':jpl_j2000,'j2000_rate':jpl_j2000_rate}
        except FileNotFoundError as error:
            jpl_delta=None
            jpl_metadata={'available':False,'reason':str(error)}
        modern_full,modern_c2,modern_gauss=full[2:],c2[2:],full_gauss[2:]
        model_rate=A_SCALE*diagnostics['residual_rate']
        headers=['tdb_j2000_seconds','delta_full_s','delta_c2_s','delta_gauss_s','sofa_s','F2','F4','model_rate']
        values=[times,modern_full,modern_c2,modern_gauss,reference,diagnostics['F2'],diagnostics['F4'],model_rate]
        if jpl_delta is not None:
            headers+=['jpl_s','jpl_rate']
            values+=[jpl_delta,jpl_rate]
        _common.csv_table('results/tables/time_series.csv',headers,np.column_stack(values))
        roundtrips=[]
        converter=TTTDBConverter(romberg)
        for second in np.linspace(start,end,13):
            epoch=TwoPartJD.from_seconds_since_j2000(second).add_seconds(.123456789)
            tt=converter.tdb_to_tt(epoch)
            back=converter.tt_to_tdb(tt)
            diag1=converter.iteration_diagnostics()
            tdb=converter.tt_to_tdb(epoch)
            diag2=converter.iteration_diagnostics()
            back_tt=converter.tdb_to_tt(tdb)
            roundtrips.append({'tdb_seconds':float(second),
                              'TDB_TT_TDB_error_seconds':back.difference_seconds(epoch),
                              'TT_TDB_TT_error_seconds':back_tt.difference_seconds(epoch),
                              'TDB_TT_TDB_solver':diag1,'TT_TDB_TT_solver':diag2})
        _common.write_json('results/tables/roundtrips.json',roundtrips)
        print('Measuring raw / constant-subtracted / fitted-average floating-point forms',flush=True)
        fitted_mean=float(np.mean(diagnostics['F2']+diagnostics['F4']))
        def stability_rate(t):
            d=model.diagnostics(t)
            f=d['F2']+d['F4']
            return np.stack((f,f-L_C_PAPER,f-fitted_mean),axis=-1)
        stability_signature=json.dumps({'base':signature,'purpose':'stability','fitted_mean':fitted_mean},sort_keys=True)
        stability,spath=get_integrator(stability_rate,stability_signature,'gauss_legendre',16,3,anchor,args,'stability_gauss16')
        stability_values=stability.integral(times)
        stability.save(spath)
        duration=times-anchor
        raw=A_SCALE*(TDB0+stability_values[:,0])-K_LINEAR*(times-T0.seconds_since_j2000())
        paper=TDB0+A_SCALE*(stability_values[:,1]+(L_C_PAPER-C_EXACT)*duration)
        fitted=TDB0+A_SCALE*(stability_values[:,2]+(fitted_mean-C_EXACT)*duration)
        wrong_paper=TDB0+A_SCALE*stability_values[:,1]
        stability_results={'reference':'independent Gauss16 residual channels, no detrending',
                           'fitted_mean_rate':fitted_mean,'raw_vs_residual':difference_metrics(times,raw-modern_gauss),
                           'paper_mean_restored_vs_residual':difference_metrics(times,paper-modern_gauss),
                           'fitted_mean_restored_vs_residual':difference_metrics(times,fitted-modern_gauss),
                           'paper_mean_without_restoration_vs_correct':difference_metrics(times,wrong_paper-modern_gauss),
                           'stats':stability.stats}
        _common.write_json('results/tables/stability.json',stability_results)
        summary={'command':' '.join(sys.argv),'versions':_common.versions(),'model_signature':json.loads(signature),
                 'anchor_jd_parts':T0_TDB.jd_parts(),'anchor_j2000_seconds':anchor,'anchor_delta_seconds':float(full[0]),
                 'output_start_jd':args.start_jd,'output_end_jd':args.end_jd,'sample_count':len(times),
                 'j2000_delta_seconds':float(full[1]),'j2000_c4_seconds':float(full[1]-c2[1]),
                 'delta_min_seconds':float(np.min(modern_full)),'delta_max_seconds':float(np.max(modern_full)),
                 'half_peak_to_peak_seconds':float(np.ptp(modern_full)/2),
                 'sofa_comparison':difference_metrics(times,modern_full-reference),
                 'long_romberg_vs_gauss':difference_metrics(times,modern_full-modern_gauss),
                 'long_channel_max_difference_seconds':np.max(np.abs(channels-independent),axis=0).tolist(),
                 'F2_range': [float(np.min(diagnostics['F2'])),float(np.max(diagnostics['F2']))],
                 'F4_range': [float(np.min(diagnostics['F4'])),float(np.max(diagnostics['F4']))],
                 'F4_over_F2_range':[float(np.min(diagnostics['F4_over_F2'])),float(np.max(diagnostics['F4_over_F2']))],
                 'c4_accumulated_range_seconds':[float(np.min(modern_full-modern_c2)),float(np.max(modern_full-modern_c2))],
                 'potential_at_j2000':{k:float(v) for k,v in model.diagnostics(0.)['potential_by_body'].items()},
                 'romberg_stats':romberg.stats,'gauss_stats':gauss.stats,
                 'sum_segment_error_indicators_seconds':romberg.cumulative_error_diagnostic.tolist(),
                 'max_partial_error_indicators_seconds':np.max(romberg.query_errors,axis=0).tolist(),
                 'jpl_reference':jpl_metadata,'elapsed_wall_seconds':time.perf_counter()-began}
        if jpl_delta is not None:
            summary['jpl_comparison']=difference_metrics(times,modern_full-jpl_delta)
            summary['jpl_rate_difference']=difference_metrics(times,model_rate-jpl_rate,value_unit='dimensionless')
        _common.write_json('results/tables/experiment_summary.json',summary)
        _common.write_json('results/tables/constants.json',{k:str(v) for k,v in decimal_constants().items()})
        _common.write_json('results/tables/time_precision.json',time_precision_experiment())
        output=_common.make_plots_directory()
        import matplotlib.pyplot as plt
        days=(times-times[0])/86400.
        plt.figure(figsize=(10,4))
        plt.plot(days,modern_full*1e3,label='Independent Eq. (21), DE440')
        plt.plot(days,reference*1e3,'--',alpha=.6,label='ERFA comparison')
        plt.xlabel('Days since first sample (TDB)');plt.ylabel('TDB − TT (ms)')
        plt.grid(alpha=.3);plt.legend();plt.tight_layout();plt.savefig(output/'delta_years.png',dpi=170);plt.close()
        fig,axes=plt.subplots(3,1,figsize=(10,8),sharex=True)
        axes[0].plot(days,modern_c2*1e3,label='c⁻² only');axes[0].plot(days,modern_full*1e3,'--',label='c⁻² + c⁻⁴')
        axes[0].set_ylabel('TDB − TT (ms)');axes[0].legend();axes[0].grid(alpha=.3)
        axes[1].plot(days,(modern_full-modern_c2)*1e9)
        axes[1].set_ylabel('Accumulated c⁻⁴ (ns)');axes[1].grid(alpha=.3)
        axes[2].plot(days,diagnostics['F4']*1e16)
        axes[2].set_ylabel('F4 rate (10⁻¹⁶)');axes[2].set_xlabel('Days since first sample (TDB)');axes[2].grid(alpha=.3)
        fig.tight_layout();fig.savefig(output/'pn_contribution.png',dpi=170);plt.close(fig)
        print(json.dumps({k:v for k,v in summary.items() if k in ('delta_min_seconds','delta_max_seconds','long_romberg_vs_gauss','sofa_comparison','jpl_comparison','elapsed_wall_seconds')},indent=2),flush=True)
    subprocess.run([sys.executable,'scripts/compare_sofa.py'],check=True)
    subprocess.run([sys.executable,'scripts/compare_jpl.py'],check=True)


if __name__=='__main__':
    main()
