"""Compare independent Eq. (21) with verified JPL de440t numerical time data."""
import argparse
import _common
import numpy as np
from tt_tdb.jpl_reference import JPLTimeEphemeris
from tt_tdb.reference import difference_metrics


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--table',default='results/tables/time_series.csv')
    parser.add_argument('--reference',default=None)
    args=parser.parse_args()
    table=np.genfromtxt(args.table,delimiter=',',names=True)
    times=table['tdb_j2000_seconds']
    try:
        with JPLTimeEphemeris(args.reference) as reference:
            delta,rate=reference.evaluate(times)
            path=str(reference.path)
    except FileNotFoundError as error:
        status={'available':False,'reason':str(error)}
        _common.write_json('results/tables/jpl_comparison.json',status)
        print(status)
        return
    difference=table['delta_full_s']-delta
    rate_difference=table['model_rate']-rate
    output={'available':True,'reference':path,'sign':'TDB-TT in SI seconds',
            'time_difference':difference_metrics(times,difference),
            'derivative_difference':difference_metrics(times,rate_difference,value_unit='dimensionless'),
            'first_tdb_j2000_seconds':float(times[0]),'last_tdb_j2000_seconds':float(times[-1])}
    _common.write_json('results/tables/jpl_comparison.json',output)
    figures=_common.make_plots_directory()
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(2,1,figsize=(10,6),sharex=True)
    days=(times-times[0])/86400.
    axes[0].plot(days,difference*1e9)
    axes[0].set_ylabel('Model − JPL (ns)');axes[0].grid(alpha=.3)
    axes[1].plot(days,rate_difference*1e16)
    axes[1].set_ylabel('Rate difference (10⁻¹⁶)');axes[1].grid(alpha=.3)
    axes[1].set_xlabel('Days since first sample (TDB)')
    fig.tight_layout();fig.savefig(figures/'jpl_difference.png',dpi=170);plt.close(fig)
    print(output)


if __name__=='__main__':
    main()
