"""Independent moon-center TCL−TDB model with the original 1977 boundary.

No LTE/JPL time-reference module is imported. Park (2021) Eqs. (3)--(7)
are evaluated at the Moon, excluding lunar self-potential and adding Earth.
Earth's external J2 uses a source-defined long-term pole diagnostic; the
unresolved fitted frame offsets and time-varying gravity are explicit.
"""
from dataclasses import replace
from pathlib import Path
from types import SimpleNamespace
import numpy as np
from .constants import C_LIGHT, L_B, TDB0, T0_TDB
from .extended_ephemeris import ExtendedDE440Ephemeris
from .jpl_model import JPLModel, read_de440_constants
from .small_body import SmallBodyPotential, file_sha256

# Fixed computational centering constant, NOT a measured/fitted LTE rate.
LUNAR_CENTERING_RATE = 1.4825e-8
ARCSEC = np.pi/(180.*3600.)

# Vondrak, Capitaine & Wallace (2011, corrected 2012), equator/ecliptic
# pole coefficient tables; evaluated independently, ERFA is test-only.
_XYPOL = np.array([[5453.282155, .4252841, -.00037173, -.000000152],
                   [-73750.930350, -.7675452, -.00018725, .000000231]])
_XYPER = np.array([
    [256.75,-819.940624,75004.344875,81491.287984,1558.515853],
    [708.15,-8444.676815,624.033993,787.163481,7774.939698],
    [274.20,2600.009459,1251.136893,1251.296102,-2219.534038],
    [241.45,2755.175630,-1102.212834,-1257.950837,-2523.969396],
    [2309.,-167.659835,-2660.664980,-2966.799730,247.850422],
    [492.20,871.855056,699.291817,639.744522,-846.485643],
    [396.10,44.769698,153.167220,131.600209,-1393.124055],
    [288.90,-512.313065,-950.865637,-445.040117,368.526116],
    [231.10,-819.415595,499.754645,584.522874,749.045012],
    [1610.,-538.071099,-145.188210,-89.756563,444.704518],
    [620.,-189.793622,558.116553,524.429630,235.934465],
    [157.87,-402.922932,-23.923029,-13.549067,374.049623],
    [220.30,179.516345,-165.405086,-210.157124,-171.330180],
    [1200.,-9.814756,9.344131,-44.919798,-22.899655]])
_PQPOL = np.array([[5851.607687,-.1189000,-.00028913,.000000101],
                   [-1600.886300,1.1689818,-.00000020,-.000000437]])
_PQPER = np.array([
    [708.15,-5486.751211,-684.661560,667.666730,-5523.863691],
    [2309.,-17.127623,2446.283880,-2354.886252,-549.747450],
    [1620.,-617.517403,399.671049,-428.152441,-310.998056],
    [492.20,413.442940,-356.652376,376.202861,421.535876],
    [1183.,78.614193,-186.387003,184.778874,-36.776172],
    [622.,-180.732815,-316.800070,335.321713,-145.278396],
    [882.,-87.676083,198.296701,-185.138669,-34.744450],
    [547.,46.140315,101.135679,-120.972830,22.885731]])


def _pole_series(t, polynomial, periodic):
    result = np.stack([np.polynomial.polynomial.polyval(t, p) for p in polynomial], axis=-1)
    for period, xc, yc, xs, ys in periodic:
        a = 2*np.pi*t/period
        result += np.stack((xc*np.cos(a)+xs*np.sin(a), yc*np.cos(a)+ys*np.sin(a)), axis=-1)
    return result*ARCSEC


def vondrak_poles(seconds):
    """Vondrak equator/ecliptic poles in mean J2000 axes, no ERFA input."""
    t = np.asarray(seconds, float)/(36525.*86400.)
    x, y = np.moveaxis(_pole_series(t, _XYPOL, _XYPER), -1, 0)
    equator = np.stack((x, y, np.sqrt(1-x*x-y*y)), axis=-1)
    p, q = np.moveaxis(_pole_series(t, _PQPOL, _PQPER), -1, 0)
    w = np.sqrt(1-p*p-q*q)
    eps0 = 84381.406*ARCSEC
    ecliptic = np.stack((p, -q*np.cos(eps0)-w*np.sin(eps0),
                         -q*np.sin(eps0)+w*np.cos(eps0)), axis=-1)
    return equator, ecliptic


def earth_pole_park(seconds, *, x_offset_arcsec=0., y_offset_arcsec=0.):
    """Park Eqs. (20)--(26), with explicit externally declared offsets.

    The Eq. (26) row basis maps inertial to date; its transpose maps the
    Eq. (23) date pole into inertial axes. Arcsecond units for ROTEX/Y are
    explicitly traced through Folkner2014 Table10 + matching DE430 header;
    applying that field convention to DE440 is a documented version inference.
    """
    t = np.asarray(seconds, float)/(36525.*86400.)
    omega = ((125*3600+2*60+40.280)-(1934*3600+8*60+10.539)*t+7.455*t*t+.008*t**3)*ARCSEC
    dpsi, deps = -17.1996*np.sin(omega)*ARCSEC, 9.2025*np.cos(omega)*ARCSEC
    eps = (84381.448-46.815*t-.00059*t*t+.001813*t**3)*ARCSEC
    se, ce = np.sin(eps), np.cos(eps)
    st, ct = np.sin(eps+deps), np.cos(eps+deps)
    pd = np.stack((np.sin(dpsi)*st, np.cos(dpsi)*st*ce-ct*se,
                   np.cos(dpsi)*st*se+ct*ce), axis=-1)
    fx, fy = x_offset_arcsec*ARCSEC, y_offset_arcsec*ARCSEC
    # Passive Park Rx(-Fx) Ry(-Fy), in precisely the published order.
    xp = np.cos(fy)*pd[..., 0]+np.sin(fy)*pd[..., 2]
    yp = pd[..., 1]
    zp = -np.sin(fy)*pd[..., 0]+np.cos(fy)*pd[..., 2]
    pd = np.stack((xp, np.cos(fx)*yp-np.sin(fx)*zp,
                   np.sin(fx)*yp+np.cos(fx)*zp), axis=-1)
    n, k = vondrak_poles(seconds)
    x = np.cross(n, k)
    x /= np.linalg.norm(x, axis=-1)[..., None]
    y = np.cross(n, x)
    # Columns are date basis vectors in inertial axes: date -> inertial.
    return x*pd[..., 0, None]+y*pd[..., 1, None]+n*pd[..., 2, None]


def external_j2_potential(displacement, gm, radius, j2, pole):
    r = np.asarray(displacement, float)
    distance = np.linalg.norm(r, axis=-1)
    n = np.asarray(pole, float)
    if np.any(distance <= radius) or not np.all(np.isfinite(r)) or not np.allclose(np.linalg.norm(n, axis=-1), 1., rtol=0, atol=2e-14):
        raise ValueError('External J2 requires finite exterior positions and unit poles')
    sine = np.sum(r*n, axis=-1)/distance
    return -gm*j2*radius**2*(3*sine*sine-1)/(2*distance**3)


class MoonObserverEphemeris(ExtendedDE440Ephemeris):
    """Compatibility adapter: earth_* fields represent the observer Moon.

    The generic Park algebra uses these legacy field names. The external
    Earth key is 'earth_external' to avoid the original explicitly geocentric
    potential guard. GM in the PN inner self-body slot is the Moon's GM.
    """
    def states_jd(self, jd1, jd2=0.):
        original = super().states_jd(jd1, jd2)
        position = {k: v for k, v in original.body_positions.items() if k != 'moon'}
        velocity = {k: v for k, v in original.body_velocities.items() if k != 'moon'}
        acceleration = {k: v for k, v in original.body_accelerations.items() if k != 'moon'}
        # Earth acceleration from the same analytic SPK polynomial.
        from .extended_ephemeris import type2_state_derivatives
        emb = type2_state_derivatives(self.kernel[0, 3], jd1, jd2)[2]
        offset = type2_state_derivatives(self.kernel[3, 399], jd1, jd2)[2]
        position['earth_external'], velocity['earth_external'] = original.earth_position, original.earth_velocity
        acceleration['earth_external'] = (emb+offset)*1000.
        gm = {k: v for k, v in original.gm.items() if k != 'moon'}
        gm['earth_external'] = original.earth_gm
        return replace(original, earth_position=original.body_positions['moon'],
                       earth_velocity=original.body_velocities['moon'], body_positions=position,
                       body_velocities=velocity, body_accelerations=acceleration, gm=gm,
                       earth_gm=original.gm['moon'])


class LunarModel:
    channel_names = ('F2_minus_fixed_center', 'major_pointmass_F4', 'solar_J2_c2',
                     'kuiper_ring_c2', 'major_vector_c4', 'major_PN_velocity_c4',
                     'major_PN_body_potential_c4', 'major_PN_radial_c4',
                     'major_PN_acceleration_c4', 'earth_static_J2_Park_pole_c2',
                     'earth_J2_zero_offset_minus_primary_DIAGNOSTIC')
    physical_channel_count = 10

    def __init__(self, ephemeris, *, header_path='docs/sources/de440_header.440',
                 centering_rate=LUNAR_CENTERING_RATE):
        self.ephemeris = ephemeris
        self.header_path = Path(header_path)
        self.header = read_de440_constants(self.header_path)
        self.centering_rate = float(centering_rate)
        self.park = JPLModel(ephemeris, header_path=header_path)

    def __call__(self, seconds):
        d = self.park.diagnostics(seconds)
        s = self.ephemeris.states(seconds)
        displacement = s.earth_position-s.body_positions['earth_external']
        u = external_j2_potential(displacement,
                                  s.gm['earth_external'], self.header['RE']*1000.,
                                  self.header['J2E'], earth_pole_park(seconds,
                                      x_offset_arcsec=self.header['ROTEX'],
                                      y_offset_arcsec=self.header['ROTEY']))
        zero_offset = external_j2_potential(displacement,
                                  s.gm['earth_external'], self.header['RE']*1000.,
                                  self.header['J2E'], earth_pole_park(seconds))
        term = d['terms']
        return np.stack((d['leading_F2']-self.centering_rate,
                         sum(term[n] for n in ('kinetic_c4', 'mixed_c4', 'potential_squared_c4')),
                         term['solar_j2_c2'], term['kuiper_ring_c2'], term['vector_c4'],
                         term['pn_velocity_c4'], term['pn_body_potential_c4'],
                         term['pn_radial_velocity_c4'], term['pn_acceleration_c4'],
                         u/C_LIGHT**2, (zero_offset-u)/C_LIGHT**2), axis=-1)

    def rate(self, seconds, small_body_rates=None):
        centered = self(seconds)[..., :self.physical_channel_count].sum(axis=-1)
        if small_body_rates is not None:
            centered += np.asarray(small_body_rates).sum(axis=-1)
        return ((L_B-self.centering_rate)-centered)/(1-L_B)

    def restore(self, seconds, centered_integrals, small_body_integrals=None):
        """Restore exact affine term, without any reference/fitted correction."""
        t = np.asarray(seconds, dtype=np.longdouble)
        elapsed = t-np.longdouble(T0_TDB.seconds_since_j2000())
        integral = np.asarray(centered_integrals, dtype=np.longdouble)[..., :self.physical_channel_count].sum(axis=-1)
        if small_body_integrals is not None:
            integral += np.asarray(small_body_integrals, dtype=np.longdouble).sum(axis=-1)
        return -np.longdouble(TDB0)+((np.longdouble(L_B)-np.longdouble(self.centering_rate))*elapsed-integral)/(1-np.longdouble(L_B))

    def fingerprint(self):
        files = ('lunar_model.py', 'constants.py', 'timebase.py', 'integration.py',
                 'extended_ephemeris.py', 'ephemeris.py', 'pn_rate.py', 'potential.py',
                 'jpl_model.py', 'small_body.py')
        return {'model': 'moon_center_Park_major_PN_plus_ring_Earth_staticJ2_v1',
                'source_sha256': {n: file_sha256(Path(__file__).with_name(n)) for n in files},
                'earth_moon_bsp_sha256': file_sha256(self.ephemeris.path),
                'gm_kernel_sha256': file_sha256(self.ephemeris.gm_path),
                'header_sha256': file_sha256(self.header_path),
                'C_LIGHT': C_LIGHT, 'L_B': L_B, 'TDB0': TDB0,
                'anchor_jd_parts': T0_TDB.jd_parts(), 'centering_rate': self.centering_rate,
                'formula': 'Q=-TDB0+[(LB-Fstar)*(s-s0)-integral(F-Fstar)ds]/(1-LB)',
                'channels': self.channel_names,
                'observer': 'Moon center; no self potential; external Earth; observer Moon GM in PN inner sum',
                'earth_J2': {'source': 'Park2021 Sec2.5 Eq20-26 and header J2E/RE',
                             'J2': self.header['J2E'], 'radius_m': self.header['RE']*1000.,
                             'pole': 'Vondrak2011 corrected2012 + 18.6-year nutation + Park25 frame rotations',
                             'frame_offset_arcsec': [self.header['ROTEX'], self.header['ROTEY']],
                             'offset_units_evidence': 'Folkner2014 Table10 explicitly arcsec and exact match to header430 ROTEX/Y; same-field convention interpreted for DE440',
                             'zero_offset_diagnostic_channel_excluded_from_primary': 10,
                             'static_only': True},
                'excluded': ['Earth J2 time dependence; original DE440 pole integration code unavailable',
                             'higher Earth harmonics and tides', 'small-body/ring PN and vector terms',
                             'solar and other direct spin terms'],
                'reference_used_in_integrand': False}


class LunarSmallBodyPotential(SmallBodyPotential):
    """Unmodified real SB441 direct potential, with observer adapter Moon."""
    def __init__(self, ephemeris, *args, **kwargs):
        class MoonPositionView:
            path = ephemeris.path

            def states_jd(self, first, second):
                # This provider needs positions only, not the expensive
                # analytic velocity/acceleration of every major body.
                emb = ephemeris.kernel[0, 3].compute(first, second)
                moon = ephemeris.kernel[3, 301].compute(first, second)
                sun = ephemeris.kernel[0, 10].compute(first, second)
                shape = np.broadcast_arrays(first, second)[0].shape+(3,)
                return SimpleNamespace(earth_position=((emb+moon).T*1000.).reshape(shape),
                                       body_positions={'sun': (sun.T*1000.).reshape(shape)})
        super().__init__(MoonPositionView(), *args, **kwargs)

    def fingerprint(self):
        result = super().fingerprint()
        result.update(model='SB441373_direct_scalar_Moon_observer_diagnostic_v1',
                      observer='Moon center', lunar_source_sha256=file_sha256(__file__),
                      frame='ICRF heliocentric geometric; subtract DE440 Moon-Sun vector; SI TDB-compatible')
        return result
