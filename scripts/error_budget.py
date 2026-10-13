#!/usr/bin/env python3
"""Quantify separately evidenced bounds, assumptions, and rounding sensitivity.

Run from the repository root with .venv/bin/python scripts/error_budget.py.
No estimates below alter or calibrate the core Eq. (21) integration.
"""
from pathlib import Path
import argparse
from decimal import Decimal, localcontext
import hashlib
import json
import re
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'src'))
import numpy as np
from scipy.integrate import simpson
from tt_tdb.constants import A_SCALE, C_LIGHT, T0_TDB, decimal_constants
from tt_tdb.ephemeris import DE440Ephemeris, BODY_IDS
from tt_tdb.timebase import TwoPartJD

ROOT = Path(__file__).resolve().parents[1]
PARK = 'https://ssd.jpl.nasa.gov/doc/Park.2021.AJ.DE440.pdf'
HEADER_URL = 'https://ssd.jpl.nasa.gov/ftp/eph/planets/ascii/de440/header.440'


def header_constants(path):
    text = path.read_text()
    names = text.split('GROUP   1040')[1].split('GROUP   1041')[0].split()
    values = text.split('GROUP   1041')[1].split('GROUP   1050')[0].split()
    if int(names[0]) != len(names)-1 or len(names) != len(values):
        raise ValueError('Invalid JPL ASCII constant counts')
    return dict(zip(names[1:], (Decimal(v.replace('D', 'E')) for v in values[1:])))


def constant_potential_effect(gm, distance_au, au, duration):
    potential = gm/(distance_au*au)
    return {'assumed_distance_au': distance_au, 'potential_m2_s2': potential,
            'leading_rate_U_over_c2': potential/C_LIGHT**2,
            'integrated_leading_correction_seconds': A_SCALE*potential/C_LIGHT**2*duration,
            'interpretation': 'hypothetical fixed distance; not actual trajectory or rigorous bound'}


def gm_last_places(path):
    """Half one printed decimal last place; not a fit covariance/uncertainty."""
    text = path.read_text()
    result = {}
    for name, body in BODY_IDS.items():
        token = re.search(rf'BODY{body}_GM\s*=\s*\(\s*([+\-\d.DEded]+)', text)[1]
        decimal = Decimal(token.replace('D', 'E').replace('d', 'e'))
        quantum_si = float(Decimal(10)**decimal.as_tuple().exponent * Decimal('1e9'))
        result[name] = {'literal_km3_s2': token, 'half_printed_last_place_m3_s2': quantum_si/2}
    return result


def real_geometry_effects(eph, start, end, step_days, constants, last_places):
    count = int(np.ceil((end-start)/(step_days*86400.)))
    seconds = np.linspace(start, end, count+1)
    state = eph.states(seconds)
    relative = state.earth_position - state.body_positions['sun']
    radius = np.linalg.norm(relative, axis=-1)
    # Pole and reference radius have explicit provenance, Park2021 Sec3.4;
    # J2 is the actual DE440 header value, not an invented fit.
    ra, dec = np.deg2rad([286.13, 63.87])
    pole = np.array([np.cos(dec)*np.cos(ra), np.cos(dec)*np.sin(ra), np.sin(dec)])
    sine_latitude = np.einsum('ij,j->i', relative, pole)/radius
    p2 = (3*sine_latitude*sine_latitude-1)/2
    j2 = float(constants['J2SUN'])
    solar_radius = float(constants['ASUN'])*1000.
    quadrupole_potential = -eph.gm['sun']*j2*solar_radius**2/radius**3*p2
    v2 = np.sum(state.earth_velocity**2, axis=-1)
    total = np.zeros_like(seconds)
    distances = {}
    for name, body_pos in state.body_positions.items():
        distances[name] = np.linalg.norm(body_pos-state.earth_position, axis=-1)
        total += eph.gm[name]/distances[name]
    # The change of F2+F4 for a small additional potential dU is exact to
    # first order: (1/c² + (1.5v²-U)/c⁴)*dU. Quadratic dU²/c⁴ is negligible
    # here but retained for the signed J2 diagnostic.
    derivative_u = 1/C_LIGHT**2 + (1.5*v2-total)/C_LIGHT**4
    delta_rate_j2 = derivative_u*quadrupole_potential-0.5*quadrupole_potential**2/C_LIGHT**4
    rounding = {}
    for name, spec in last_places.items():
        half_decimal = spec['half_printed_last_place_m3_s2']
        half_binary = float(np.spacing(eph.gm[name]))/2.
        integral_inverse_r = float(simpson(A_SCALE*derivative_u/distances[name], x=seconds))
        rounding[name] = dict(spec, half_float64_ulp_m3_s2=half_binary,
                             fixed_orbit_d_delta_d_gm=integral_inverse_r,
                             half_decimal_place_effect_seconds=integral_inverse_r*half_decimal,
                             half_float64_ulp_effect_seconds=integral_inverse_r*half_binary)
    return {'sample_count':len(seconds), 'step_days_max':step_days,
            'earth_sun_distance_au_min':float(radius.min()/float(constants['AU']*1000)),
            'earth_sun_distance_au_max':float(radius.max()/float(constants['AU']*1000)),
            'solar_J2_signed_missing_potential_min_m2_s2':float(quadrupole_potential.min()),
            'solar_J2_signed_missing_potential_max_m2_s2':float(quadrupole_potential.max()),
            'solar_J2_fixed_orbit_integrated_correction_seconds':float(simpson(A_SCALE*delta_rate_j2, x=seconds)),
            'gm_rounding_by_body':rounding,
            'gm_decimal_places_sum_abs_effect_seconds':sum(abs(v['half_decimal_place_effect_seconds']) for v in rounding.values()),
            'gm_float64_half_ulp_sum_abs_effect_seconds':sum(abs(v['half_float64_ulp_effect_seconds']) for v in rounding.values())}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=ROOT/'results/tables/model_error_budget.json')
    args = parser.parse_args()
    header_path = ROOT/'data/header.440'
    constants = header_constants(header_path)
    small = json.loads((ROOT/'data/small_body_gm.json').read_text())
    # Verify the previously derived mass summary against the genuine header.
    with localcontext() as ctx:
        ctx.prec = 60
        au = float(constants['AU']*1000)
        factor = float((constants['AU']*1000)**3/Decimal(86400)**2)
    for group, record in small['groups'].items():
        actual = sum(float(constants[name])*factor for name in record['gm_m3_s2'])
        if not np.isclose(actual, record['total_gm_m3_s2'], rtol=2e-15):
            raise ValueError(f'Derived {group} masses do not match source header')
    exact = decimal_constants()
    end = TwoPartJD.from_decimal_jd('2461041.5').seconds_since_j2000()
    intervals = {'1977_anchor_to_2026_Jan01': (T0_TDB.seconds_since_j2000(), end),
                 '2024_Jan01_to_2026_Jan01': (TwoPartJD.from_decimal_jd('2460310.5').seconds_since_j2000(), end)}
    output = {'scope':'Separate error mechanisms, never an additive certified total accuracy',
              'sources':{'header':HEADER_URL, 'header_sha256':hashlib.sha256(header_path.read_bytes()).hexdigest(), 'Park2021':PARK},
              'constants':{'A':A_SCALE,'C_minus_LC_paper':str(exact['C_minus_L_C']),
                           'solar_J2':float(constants['J2SUN']),'solar_radius_m':float(constants['ASUN'])*1000,
                           'solar_pole_ra_deg':286.13,'solar_pole_dec_deg':63.87},
              'small_body_mass_sums':{group:{'count':record['count'],'GM_m3_s2':record['total_gm_m3_s2']} for group,record in small['groups'].items()},
              'limitations':['No asteroid/KBO trajectories are available: scenarios cannot establish a rigorous total missing-potential bound.',
                             'Ring geometry bounds require Earth heliocentric distance <=1.02 AU for the whole interval; finite samples verify plausibility, not a proof between samples.',
                             'Solar J2 uses genuine fixed DE440 orbit and a constant pole; it is an isolated additional potential diagnostic, not a complete JPL model discrepancy.',
                             'GM printed digits and binary ULP are representation sensitivity only, not physical parameter covariance.',
                             'Header GMSDOT=-9e-14 and JDEPOC=2440400.5 merit a separate mass-variation convention audit; no undocumented interpretation is applied.',
                             'DE440 observational range residuals do not directly provide an Earth BCRS state covariance or a certified TT−TDB physical error.'],
              'intervals':{}}
    with DE440Ephemeris() as eph:
        places = gm_last_places(eph.gm_path)
        for label, (start, stop) in intervals.items():
            duration = stop-start
            geometry = real_geometry_effects(eph,start,stop,0.5,constants,places)
            coarse = real_geometry_effects(eph,start,stop,1.0,constants,places)
            ring_gm = small['groups']['ring_points']['total_gm_m3_s2']
            ring_low, ring_high = ring_gm/(45.02*au), ring_gm/(42.98*au)
            asteroid_gm = small['groups']['asteroids']['total_gm_m3_s2']
            kbo_gm = small['groups']['individual_kbos']['total_gm_m3_s2']
            cpv_gm = sum(small['groups']['asteroids']['gm_m3_s2'][k] for k in ['MA0001','MA0002','MA0004'])
            # r >= 0.98 AU and |P2| <= 1 give a conditional absolute bound.
            j2_bound = eph.gm['sun']*abs(float(constants['J2SUN']))*(float(constants['ASUN'])*1000)**2/(0.98*au)**3
            output['intervals'][label] = {'start_tdb_seconds_J2000':start,'end_tdb_seconds_J2000':stop,
                'duration_seconds':duration, 'duration_days':duration/86400.,
                'paper_2_14e_19_times_elapsed_seconds':2.14e-19*duration,
                'paper_position_remainder_geocenter_seconds':0.,
                'paper_position_remainder_generic_bound_seconds':1.91e-16,
                'LC_paper_unrestored_linear_error_seconds':float(exact['A']*exact['C_minus_L_C'])*duration,
                'LC_paper_error_note':'Positive if subtracting smaller paper LC instead of C_exact and omitting restoration; not an error in the implemented Eq21.',
                'ring_conditional_geometry_bound':{'condition':'44AU heliocentric ring; Earth heliocentric radius <=1.02AU throughout',
                    'potential_m2_s2_range':[ring_low,ring_high],
                    'leading_integrated_correction_seconds_range':[A_SCALE*ring_low/C_LIGHT**2*duration,A_SCALE*ring_high/C_LIGHT**2*duration]},
                'solar_J2_conditional_absolute_bound':{'condition':'Earth-Sun r>=0.98AU throughout; |P2|<=1; fixed GM,J2,R',
                    'abs_potential_m2_s2_max':j2_bound,'leading_abs_integrated_correction_seconds_max':A_SCALE*j2_bound/C_LIGHT**2*duration},
                'asteroid_distance_scenarios':[constant_potential_effect(asteroid_gm,r,au,duration) for r in (0.5,1.,2.,3.)],
                'Ceres_Pallas_Vesta_distance_scenarios':[constant_potential_effect(cpv_gm,r,au,duration) for r in (1.,2.)],
                'KBO_distance_scenarios':[constant_potential_effect(kbo_gm,r,au,duration) for r in (30.,44.,100.)],
                'real_DE440_geometry_diagnostics':geometry,
                'solar_J2_step_half_vs_one_day_difference_seconds':geometry['solar_J2_fixed_orbit_integrated_correction_seconds']-coarse['solar_J2_fixed_orbit_integrated_correction_seconds']}
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(output,indent=2)+'\n')
    print(args.output)
    for label,record in output['intervals'].items():
        print(label, 'paper_remainder_ps=',record['paper_2_14e_19_times_elapsed_seconds']*1e12,
              'LC_error_ns=',record['LC_paper_unrestored_linear_error_seconds']*1e9,
              'ring_ns=',[v*1e9 for v in record['ring_conditional_geometry_bound']['leading_integrated_correction_seconds_range']],
              'solar_J2_ps=',record['real_DE440_geometry_diagnostics']['solar_J2_fixed_orbit_integrated_correction_seconds']*1e12)


if __name__ == '__main__':
    main()
