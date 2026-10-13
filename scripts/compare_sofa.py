"""Compare saved independent integral results with the geocentric ERFA series."""
import argparse
import _common
import numpy as np
from tt_tdb.reference import sofa_geocentric, difference_metrics


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--table', default='results/tables/time_series.csv')
    args = parser.parse_args()
    table = np.genfromtxt(args.table, delimiter=',', names=True)
    t = table['tdb_j2000_seconds']
    reference = sofa_geocentric(2451545., t/86400.)
    difference = table['delta_full_s'] - reference
    metrics = difference_metrics(t, difference)
    _common.write_json('results/tables/sofa_comparison.json', metrics)
    import matplotlib.pyplot as plt
    output = _common.make_plots_directory()
    plt.figure(figsize=(10,4))
    plt.plot((t-t[0])/86400., difference*1e6)
    plt.xlabel('Days since first sample (TDB)')
    plt.ylabel('Model − ERFA (microseconds)')
    plt.grid(alpha=.3)
    plt.tight_layout()
    plt.savefig(output/'sofa_difference.png', dpi=170)
    plt.close()
    print(metrics)


if __name__ == '__main__':
    main()
