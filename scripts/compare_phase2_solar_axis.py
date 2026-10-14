"""Source-defined solar J2 latitude diagnostic; never select a fitted axis.

Folkner 2014 Eq. (7) specifies latitude relative to the solar equator.
Park 2021 Eq. (5) calls it heliocentric ecliptic latitude. Both fixed definitions
are evaluated and reported without modifying the default model.
"""
import hashlib
import json
from pathlib import Path
import _common
import numpy as np
import spiceypy as spice
from tt_tdb.constants import A_SCALE, C_LIGHT, T0_TDB
from tt_tdb.ephemeris import DE440Ephemeris
from tt_tdb.integration import SegmentedIntegrator
from tt_tdb.jpl_model import read_de440_constants, solar_j2_potential
from tt_tdb.reference import difference_metrics


def main():
    header_path = Path('docs/sources/de440_header.440')
    header = read_de440_constants(header_path)
    ra, dec = np.deg2rad([286.13, 63.87])
    spin = np.array([np.cos(dec)*np.cos(ra), np.cos(dec)*np.sin(ra), np.sin(dec)])
    ecliptic = spice.pxform('ECLIPJ2000', 'J2000', 0.)[:, 2]
    times = 757339200. + np.arange(732)*86400.
    anchor = T0_TDB.seconds_since_j2000()
    with DE440Ephemeris() as eph:
        def rates(t):
            state = eph.states(t)
            return np.column_stack([
                solar_j2_potential(state.earth_position, state.body_positions['sun'],
                                   state.gm['sun'], header['J2SUN'], header['ASUN']*1000, axis)
                / C_LIGHT**2 for axis in (spin, ecliptic)])
        signature = {'purpose': 'two_predeclared_original_source_latitude_conventions',
                     'bsp_sha256': hashlib.sha256(eph.path.read_bytes()).hexdigest(),
                     'header_sha256': hashlib.sha256(header_path.read_bytes()).hexdigest(),
                     'axes': [spin.tolist(), ecliptic.tolist()]}
        values = {}
        for method, order in (('gauss_legendre', 24), ('romberg', 8)):
            integrator = SegmentedIntegrator(rates, anchor=anchor, method=method, order=order,
                                            fingerprint=signature)
            values[method] = A_SCALE*integrator.integral(times)
        rate_values = A_SCALE*rates(times)
    difference = values['gauss_legendre'][:, 1]-values['gauss_legendre'][:, 0]
    result = {
        'purpose': 'published latitude conventions, not calibration; default model unchanged',
        'axes_predeclared_from_sources': True, 'reference_used_to_choose_axis': False,
        'sources': {'spin_equator': 'https://ipnpr.jpl.nasa.gov/progress_report/42-196/196C.pdf#page=7',
                    'literal_ecliptic': 'https://ssd.jpl.nasa.gov/doc/Park.2021.AJ.DE440.pdf#page=3'},
        'source_equations': ['Folkner 2014 Eq. (7)', 'Park 2021 Eq. (5)'],
        'signature': signature, 'anchor_seconds': anchor, 'sample_count': len(times),
        'ecliptic_minus_spin_first_seconds': float(difference[0]),
        'ecliptic_minus_spin_last_seconds': float(difference[-1]),
        'modern_max_abs_rate_difference': float(abs(rate_values[:, 1]-rate_values[:, 0]).max()),
        'ecliptic_minus_spin_time_metrics': difference_metrics(times, difference),
        'independent_quadrature_difference': difference_metrics(times, (values['gauss_legendre']-values['romberg']).sum(axis=-1)),
        'independent_quadrature_max_each_channel_seconds': abs(values['gauss_legendre']-values['romberg']).max(axis=0).tolist(),
        'independent_ecliptic_minus_spin_quadrature_difference': difference_metrics(
            times, difference-(values['romberg'][:, 1]-values['romberg'][:, 0])),
        'not_used_in_primary_acceptance': True}
    # If available, show both interpretations for the separately labelled
    # mixed SB441 diagnostic. Both values are retained; no axis is selected.
    diagnostic = Path('results/tables/phase2_sb441_time_series.csv')
    if diagnostic.exists():
        table = np.genfromtxt(diagnostic, delimiter=',', names=True)
        if not np.array_equal(table['tdb_j2000_seconds'], times):
            raise ValueError('Solar-axis and SB441 diagnostic epochs disagree')
        raw = table['sb441_delta_s']-table['jpl_delta_s']
        result['SB441_diagnostic_original_spin_residual'] = difference_metrics(times, raw)
        result['SB441_diagnostic_literal_ecliptic_residual'] = difference_metrics(times, raw+difference)
    _common.write_json('results/tables/phase2_solar_axis_diagnostic.json', result)
    _common.csv_table('results/tables/phase2_solar_axis_diagnostic.csv',
                      ['tdb_j2000_seconds', 'spin_J2_delta_s', 'ecliptic_J2_delta_s',
                       'spin_J2_rate', 'ecliptic_J2_rate'],
                      np.column_stack((times, values['gauss_legendre'], rate_values)))
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
