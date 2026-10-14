"""A few fixed-anchor historical raw checks; not long-span acceptance."""
from datetime import datetime
import json
from pathlib import Path
import _common
import numpy as np
import spiceypy as spice
from tt_tdb.constants import A_SCALE, C_LIGHT, TDB0, T0_TDB
from tt_tdb.ephemeris import DE440Ephemeris
from tt_tdb.extended_ephemeris import ExtendedDE440Ephemeris
from tt_tdb.integration import SegmentedIntegrator
from tt_tdb.jpl_model import JPLModel, solar_j2_potential
from tt_tdb.jpl_reference import JPLTimeEphemeris
from tt_tdb.pn_rate import PNModel
from tt_tdb.small_body import SmallBodyPotential


def main():
    anchor = T0_TDB.seconds_since_j2000()
    labels = ['1977 defining anchor', '1980-01-01', '1990-01-01', '2000-01-01',
              '2010-01-01', '2020-01-01', '2024-01-01', '2026-01-01']
    years = [1980, 1990, 2000, 2010, 2020, 2024, 2026]
    origin = datetime(2000, 1, 1, 12)
    times = np.array([anchor]+[(datetime(year, 1, 1)-origin).total_seconds() for year in years])
    primary = json.loads(Path('results/tables/phase2_jpl_comparison.json').read_text())
    with DE440Ephemeris() as e, ExtendedDE440Ephemeris() as extended, SmallBodyPotential(e) as small:
        baseline, model = PNModel(e), JPLModel(extended)
        if json.loads(model.fingerprint()) != primary['model_signature']['model']:
            raise ValueError('Measured primary model fingerprint is stale')
        names = model.correction_names

        def channels(t):
            terms = model.correction_rates(t)
            return np.column_stack((baseline(t), *(terms[name] for name in names)))

        major = SegmentedIntegrator.load('results/cache/phase2_gauss_legendre24_halfday.npz',
                  channels, fingerprint=primary['model_signature'], anchor=anchor,
                  method='gauss_legendre', order=24, channels=2+len(names), step_seconds=43200.)
        sb_signature = small.fingerprint()
        small_integrator = SegmentedIntegrator.load('results/cache/phase2_sb441_gauss_legendre24.npz',
                  small, fingerprint=sb_signature, anchor=anchor,
                  method='gauss_legendre', order=24, channels=2, step_seconds=43200.)
        major_values = major.integral(times)
        primary_delta = TDB0+A_SCALE*major_values.sum(axis=-1)
        small_values = A_SCALE*small_integrator.integral(times)
        mixed_delta = primary_delta+small_values.sum(axis=-1)
        spin, ecliptic = model.solar_axis, spice.pxform('ECLIPJ2000', 'J2000', 0.)[:, 2]

        def latitude_difference_rate(t):
            state = e.states(t)
            kwargs = dict(earth_position=state.earth_position, sun_position=state.body_positions['sun'],
                          sun_gm=state.gm['sun'], j2=model.solar_j2, radius=model.solar_radius)
            return (solar_j2_potential(**kwargs, axis=ecliptic)-
                    solar_j2_potential(**kwargs, axis=spin))/C_LIGHT**2

        axis_integrator = SegmentedIntegrator(latitude_difference_rate, anchor=anchor,
                method='gauss_legendre', order=24, channels=1,
                fingerprint={'purpose': 'source_predeclared_solar_latitude_difference',
                             'model': model.fingerprint(), 'ecliptic_axis': ecliptic.tolist()})
        axis_difference = A_SCALE*axis_integrator.integral(times)[:, 0]
        with JPLTimeEphemeris() as reference:
            jpl_delta, _ = reference.evaluate(times)
        rows = [{
            'label': label, 'tdb_j2000_seconds': float(t),
            'primary_raw_residual_seconds': float(p-j),
            'mixed_SB441_spin_raw_residual_seconds': float(m-j),
            'mixed_SB441_literal_ecliptic_raw_residual_seconds': float(m+d-j),
            'asteroid_anchor_integral_seconds': float(s[0]),
            'kbo_anchor_integral_seconds': float(s[1]),
            'literal_ecliptic_minus_spin_J2_seconds': float(d)}
            for label, t, p, m, d, j, s in zip(labels, times, primary_delta, mixed_delta,
                                              axis_difference, jpl_delta, small_values)]
    _common.write_json('results/tables/phase2_history_profile.json', {
        'kind': 'few_historical_fixed_original_anchor_raw_checks; no fitted offsets',
        'not_full_long_interval_acceptance': True, 'primary_signature': primary['model_signature'],
        'SB441_signature': sb_signature, 'epochs': rows})
    _common.csv_table('results/tables/phase2_history_profile.csv',
          ['tdb_j2000_seconds', 'primary_raw_seconds', 'SB441_spin_raw_seconds',
           'SB441_literal_ecliptic_raw_seconds'],
          np.column_stack((times, primary_delta-jpl_delta, mixed_delta-jpl_delta,
                           mixed_delta+axis_difference-jpl_delta)))
    figures = _common.make_plots_directory()
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(2, 1, figsize=(10, 7), sharex=True)
    x = [1977]+years
    axes[0].plot(x, (primary_delta-jpl_delta)*1e9, 'o-', label='DE440 major bodies + ring')
    axes[0].set_ylabel('Raw residual (ns)')
    axes[1].plot(x, (mixed_delta-jpl_delta)*1e12, 'o-', label='SB441 diagnostic, solar equator')
    axes[1].plot(x, (mixed_delta+axis_difference-jpl_delta)*1e12, 'o-', label='SB441 diagnostic, literal ecliptic')
    axes[1].set_ylabel('Raw residual (ps)')
    axes[1].set_xlabel('Selected January 1 TDB epochs; 1977 point is defining anchor')
    for ax in axes:
        ax.grid(alpha=.3); ax.legend()
    fig.suptitle('Historical profile from original anchor; sparse diagnostic only')
    fig.tight_layout(); fig.savefig(figures/'phase2_history_profile.png', dpi=170); plt.close(fig)
    print(json.dumps(rows, indent=2))


if __name__ == '__main__':
    main()
