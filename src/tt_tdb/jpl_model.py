"""Independent geocentric Park et al. (2021) Eqs. (3)–(7) candidate.

All rates use SI, geometric ICRF states and consistent TDB-compatible GMs.
No time-ephemeris reference is imported. The default major-body input is
INCOMPLETE: direct asteroid/KBO potentials require unavailable states.
The published Kuiper ring leading scalar potential is implemented separately;
its c^-4 mixed terms, vector potential and Eq. (7) terms remain omitted.
"""
from dataclasses import dataclass, asdict
from pathlib import Path
import hashlib
import json
import re
import numpy as np
import spiceypy as spice
from .constants import C_LIGHT, C_EXACT, A_SCALE, TDB0, T0_JD_DECIMAL
from .pn_rate import rate_c2, rate_c4
from .potential import external_potential


def read_de440_constants(path):
    """Parse paired GROUP 1040/1041 names/values; refuse malformed headers."""
    text = Path(path).read_text()
    groups = dict(re.findall(r'GROUP\s+(\d+)\s*\n(.*?)(?=GROUP\s+\d+|\Z)',text,re.S))
    names = groups['1040'].split(); values = groups['1041'].split()
    count = int(names.pop(0))
    if int(values.pop(0)) != count or len(names) != count or len(values) != count:
        raise ValueError('DE440 header constant counts disagree')
    result = dict(zip(names,(float(x.replace('D','E')) for x in values)))
    if result['DENUM'] != 440:
        raise ValueError('DE440 constants required')
    return result


def solar_j2_potential(earth_position, sun_position, sun_gm, j2, radius, axis):
    """Park Eq. (5): −GM J2 R² P2(axis·n)/r³, m²/s².

    ``axis`` is explicit to expose the paper's latitude-convention ambiguity.
    Default model uses the solar spin axis given in Park Sec. 3.5; Eq. (5)
    calls the latitude 'ecliptic', which is retained as an unresolved wording.
    """
    axis = np.asarray(axis,float)
    if axis.shape != (3,) or not np.all(np.isfinite(axis)) or not np.isclose(np.linalg.norm(axis),1,rtol=0,atol=1e-12):
        raise ValueError('Solar symmetry axis must be a finite unit vector')
    if not np.isfinite(j2) or not np.isfinite(radius) or radius <= 0 or sun_gm <= 0:
        raise ValueError('Invalid solar J2 parameters')
    displacement = np.asarray(earth_position)-np.asarray(sun_position)
    distance = np.linalg.norm(displacement,axis=-1)
    if np.any(distance <= 0):
        raise ValueError('Earth coincides with Sun')
    latitude_sine = np.sum(displacement*axis,axis=-1)/distance
    return -sun_gm*j2*radius**2*(3*latitude_sine**2-1)/(2*distance**3)


def kuiper_ring_potential(earth_position,sun_position,total_gm,*,radius=44*149597870700.,phase=0.):
    """Published Park Sec. 2.2 ring: 36 equal masses, heliocentric 44 AU.

    Ecliptic plane is explicitly SPICE ECLIPJ2000. The unreported orbital
    phase is set to zero as a geometry convention, never fitted. Its leading
    potential dependence begins at harmonic 36, bounded in the Earth regime
    by 4(GM/R) q**36/(1-q) for q=|r_E-r_S|/R < 1 (see model mapping).
    No unavailable ring orbital velocities or PN potentials are synthesized.
    """
    if not np.isfinite(total_gm) or total_gm <= 0 or radius <= 0 or not np.isfinite(phase):
        raise ValueError('Positive finite ring GM/radius and finite phase required')
    relative=np.asarray(earth_position)-np.asarray(sun_position)
    if relative.shape[-1:] != (3,) or not np.all(np.isfinite(relative)):
        raise ValueError('Finite Earth/Sun Cartesian positions required')
    radius_ratio=np.linalg.norm(relative,axis=-1)/radius
    if np.any(radius_ratio >= 1):
        raise ValueError('Ring calculation restricted to interior geocentric regime')
    transform=spice.pxform('ECLIPJ2000','J2000',0.)
    phases=phase+2*np.pi*np.arange(36)/36
    ecliptic=np.stack((np.cos(phases),np.sin(phases),np.zeros(36)),axis=-1)
    points=radius*(ecliptic@transform.T)
    distances=np.linalg.norm(relative[...,None,:]-points,axis=-1)
    return (total_gm/36)*np.sum(1/distances,axis=-1)


def vector_potential(earth_position, body_positions, body_velocities, gm):
    """Park Eq. (6): Σ GM_i v_i/r_iE, m³/s³."""
    _, potentials = external_potential(earth_position,body_positions,gm,breakdown=True)
    if body_velocities.keys() != gm.keys():
        raise ValueError('Body velocities and GMs must name the same bodies')
    return np.sum(np.stack([potentials[n][...,None]*body_velocities[n] for n in gm]),axis=0)


def pn_scalar_components(state):
    """Four additive parts of Park Eq. (7), each in m⁴/s⁴.

    Outer sum excludes Earth, inner scalar potential includes Earth's GM.
    Accelerations must be actual barycentric ephemeris derivatives. All
    supplied bodies contribute; completeness is the provider's responsibility.
    """
    names = state.gm.keys()
    if state.body_velocities.keys() != names or state.body_accelerations.keys() != names:
        raise ValueError('Body positions, velocities, accelerations and GMs must agree')
    if not np.isfinite(state.earth_gm) or state.earth_gm <= 0:
        raise ValueError('Earth GM required by inner Eq. (7) sum')
    _, outer = external_potential(state.earth_position,state.body_positions,state.gm,breakdown=True)
    components = {k:np.zeros_like(next(iter(outer.values()))) for k in
                  ('velocity','body_potential','radial_velocity','acceleration')}
    for name in names:
        position = state.body_positions[name]
        velocity = np.asarray(state.body_velocities[name],float)
        acceleration = np.asarray(state.body_accelerations[name],float)
        if velocity.shape != np.shape(position) or acceleration.shape != np.shape(position) or not np.all(np.isfinite(velocity)) or not np.all(np.isfinite(acceleration)):
            raise ValueError('Finite velocity and acceleration shapes must match positions')
        relative = state.earth_position-position
        radius = np.linalg.norm(relative,axis=-1)
        inner = state.earth_gm/radius
        for other in names:
            if other == name: continue
            separation = np.linalg.norm(state.body_positions[other]-position,axis=-1)
            if np.any(separation <= 0): raise ValueError('Distinct external bodies coincide')
            inner = inner+state.gm[other]/separation
        components['velocity'] += outer[name]*(-2*np.sum(velocity**2,axis=-1))
        components['body_potential'] += outer[name]*inner
        components['radial_velocity'] += outer[name]*0.5*(np.sum(velocity*relative,axis=-1)/radius)**2
        components['acceleration'] += outer[name]*0.5*np.sum(acceleration*relative,axis=-1)
    return components


@dataclass(frozen=True)
class JPLTermSwitches:
    """Independent diagnostic ablations; disabled terms never recalibrate C."""
    kinetic_c2: bool = True
    monopole_c2: bool = True
    kinetic_c4: bool = True
    mixed_c4: bool = True
    potential_squared_c4: bool = True
    solar_j2_c2: bool = True
    kuiper_ring_c2: bool = True
    vector_c4: bool = True
    pn_velocity_c4: bool = True
    pn_body_potential_c4: bool = True
    pn_radial_velocity_c4: bool = True
    pn_acceleration_c4: bool = True


class JPLModel:
    """Park candidate with explicit provenance, options and missing-body status.

    ``__call__`` returns [F2−C_exact,F4] for the unchanged integrator/converter.
    ``correction_rates`` supplies stable additive differences from Phase 1.
    ``contribution_rates`` supplies the full physical rate decomposition.
    """
    correction_names = ('solar_j2_c2','kuiper_ring_c2','vector_c4','pn_velocity_c4','pn_body_potential_c4',
                        'pn_radial_velocity_c4','pn_acceleration_c4')
    def __init__(self,ephemeris,*,header_path='docs/sources/de440_header.440',switches=None,solar_axis=None):
        self.ephemeris = ephemeris
        self.header_path = Path(header_path)
        self.header_constants = read_de440_constants(header_path)
        self.switches = switches or JPLTermSwitches()
        self.solar_radius = self.header_constants['ASUN']*1000  # header length unit km
        self.solar_j2 = self.header_constants['J2SUN']
        ra, dec = np.deg2rad([286.13,63.87])  # Park Sec. 3.5, Archinal et al. 2018
        self.solar_axis = np.asarray(solar_axis if solar_axis is not None else
                        [np.cos(dec)*np.cos(ra),np.cos(dec)*np.sin(ra),np.sin(dec)],float)
        self.latitude_convention = 'Park 3.5 solar spin axis' if solar_axis is None else 'explicit diagnostic axis'
        ring_names=['MA'+str(n) for n in range(8201,8237)]
        ring_values=[self.header_constants[n] for n in ring_names]
        if len(set(ring_values)) != 1: raise ValueError('Published equal-mass ring requires equal header GMs')
        self.ring_gm=sum(ring_values)*(149597870700.)**3/86400**2
        self.ring_radius=44*149597870700.
        self.ring_phase=0.
        self.missing_potential_sources = tuple(n for n in ephemeris.missing_potential_sources if n != 'Kuiper-belt ring')+(
            'Kuiper ring velocity/acceleration, mixed c^-4 and PN potentials',) 

    def diagnostics(self,seconds):
        state = self.ephemeris.states(seconds)
        if state.scale != 'TDB-compatible':
            raise ValueError('Consistently TDB-compatible states and GMs required')
        u,bodies = external_potential(state.earth_position,state.body_positions,state.gm,breakdown=True)
        v2 = np.sum(state.earth_velocity**2,axis=-1)
        w = vector_potential(state.earth_position,state.body_positions,state.body_velocities,state.gm)
        scalar = pn_scalar_components(state)
        c2,c4 = C_LIGHT**2,C_LIGHT**4
        terms = {'kinetic_c2':v2/(2*c2),'monopole_c2':u/c2,
                 'kinetic_c4':v2*v2/(8*c4),'mixed_c4':1.5*v2*u/c4,
                 'potential_squared_c4':-0.5*u*u/c4,
                 'kuiper_ring_c2':kuiper_ring_potential(state.earth_position,state.body_positions['sun'],self.ring_gm,radius=self.ring_radius,phase=self.ring_phase)/c2,
                 'solar_j2_c2':solar_j2_potential(state.earth_position,state.body_positions['sun'],state.gm['sun'],self.solar_j2,self.solar_radius,self.solar_axis)/c2,
                 'vector_c4':-4*np.sum(state.earth_velocity*w,axis=-1)/c4}
        terms.update({'pn_'+name+'_c4':-value/c4 for name,value in scalar.items()})
        terms = {name:(value if getattr(self.switches,name) else np.zeros_like(value)) for name,value in terms.items()}
        # Keep Phase 1 arithmetic for the default leading terms, avoiding an
        # extra rounding error from separate kinetic and monopole division.
        f2 = rate_c2(v2,u) if self.switches.kinetic_c2 and self.switches.monopole_c2 else terms['kinetic_c2']+terms['monopole_c2']
        leading_f2 = f2
        f2 = leading_f2+terms['solar_j2_c2']+terms['kuiper_ring_c2']
        baseline4 = rate_c4(v2,u) if all((self.switches.kinetic_c4,self.switches.mixed_c4,self.switches.potential_squared_c4)) else sum(terms[n] for n in ('kinetic_c4','mixed_c4','potential_squared_c4'))
        f4 = baseline4+sum(terms[n] for n in self.correction_names if n not in ('solar_j2_c2','kuiper_ring_c2'))
        return {'F2':f2,'F4':f4,'leading_F2':leading_f2,
                'residual_rate':leading_f2-C_EXACT+terms['solar_j2_c2']+terms['kuiper_ring_c2']+f4,
                'terms':terms,'potential':u,'potential_by_body':bodies,'velocity_squared':v2,
                'vector_potential':w,'pn_scalar_components':scalar,
                'baseline_F2':rate_c2(v2,u),'baseline_F4':rate_c4(v2,u),
                'missing_potential_sources':self.missing_potential_sources,
                'latitude_convention':self.latitude_convention,'observer_position_term_seconds':0.0}

    def contribution_rates(self,seconds):
        return self.diagnostics(seconds)['terms']

    def correction_rates(self,seconds):
        terms = self.contribution_rates(seconds)
        return {name:terms[name] for name in self.correction_names}

    def __call__(self,seconds):
        d = self.diagnostics(seconds)
        # Add J2 after subtracting C so its small contribution is preserved.
        f2_residual = d['leading_F2']-C_EXACT+d['terms']['solar_j2_c2']+d['terms']['kuiper_ring_c2']
        return np.stack((f2_residual,d['F4']),axis=-1)

    def fingerprint(self):
        def sha(path):
            digest = hashlib.sha256()
            with Path(path).open('rb') as stream:
                for chunk in iter(lambda:stream.read(1024*1024),b''): digest.update(chunk)
            return digest.hexdigest()
        return json.dumps({'model':'Park2021_eq3_7_major_bodies_plus_ring_scalar_candidate_v2',
          'sources':{name:sha(Path(__file__).with_name(name)) for name in
                     ('jpl_model.py','extended_ephemeris.py','pn_rate.py','potential.py','ephemeris.py','constants.py','timebase.py')},
          'bsp_sha256':sha(self.ephemeris.path),'gm_sha256':sha(self.ephemeris.gm_path),
          'header_sha256':sha(self.header_path),'gm':self.ephemeris.gm,
          'switches':asdict(self.switches),'solar_axis':self.solar_axis.tolist(),
          'latitude_convention':self.latitude_convention,'ring':{'gm_si':self.ring_gm,'radius_m':self.ring_radius,
          'frame':'SPICE ECLIPJ2000','phase_radians':self.ring_phase,'count':36,'channels':'leading scalar only'},
          'missing':self.missing_potential_sources,
          'scale':'TDB-compatible SI ICRF geometric','A':A_SCALE,'C_exact':C_EXACT,
          'T0':T0_JD_DECIMAL,'TDB0':TDB0},sort_keys=True)
