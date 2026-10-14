"""Separate the measured Kuiper-ring ns contribution from smaller ps terms."""
import _common
import numpy as np


def main():
    table = np.genfromtxt('results/tables/phase2_time_series.csv', delimiter=',', names=True)
    columns = table.dtype.names
    terms = [name[:-8] for name in columns if name.endswith('_delta_s')
             and name not in ('phase1_delta_s', 'extended_delta_s', 'gauss_delta_s', 'jpl_delta_s')]
    days = (table['tdb_j2000_seconds']-table['tdb_j2000_seconds'][0])/86400.
    figures = _common.make_plots_directory()
    import matplotlib.pyplot as plt
    from matplotlib.ticker import FormatStrFormatter
    fig, axes = plt.subplots(4, 1, figsize=(11, 11), sharex=True)
    ring = 'kuiper_ring_c2'
    if ring in terms:
        axes[0].plot(days, table[ring+'_delta_s']*1e9, label=ring)
        axes[2].plot(days, table[ring+'_rate']*1e17, label=ring)
    for name in terms:
        if name == ring:
            continue
        axes[1].plot(days, table[name+'_delta_s']*1e12, label=name)
        axes[3].plot(days, table[name+'_rate']*1e19, label=name)
    for ax, label in zip(axes, ('Ring integral from 1977 (ns)',
                               'Other integrals from 1977 (ps)',
                               'Ring rate (10⁻¹⁷)', 'Other rates (10⁻¹⁹)')):
        ax.set_ylabel(label)
        ax.grid(alpha=.3)
        ax.legend(fontsize=8, ncol=2)
    axes[2].yaxis.set_major_formatter(FormatStrFormatter('%.7f'))
    axes[-1].set_xlabel('Days since first TDB epoch')
    fig.suptitle('Actual model contributions; original anchor, no fitted offsets')
    fig.tight_layout()
    fig.savefig(figures/'phase2_term_contributions_detail.png', dpi=170)
    plt.close(fig)


if __name__ == '__main__':
    main()
