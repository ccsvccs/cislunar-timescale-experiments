"""Independent direct-solar-spin diagnostic, not part of the Park candidate.

Turyshev arXiv:2507.22145v2 p35 (A5), (A8), (A9), p40 (A33),
(A34), (A39): F_spin = (-4 v_E dot w_spin - D_spin)/c**4.
Both scalar and vector terms are required. This module never imports a time
ephemeris or alters the Earth trajectory (which already includes LT dynamics).
"""
from pathlib import Path
import hashlib
import numpy as np
from .constants import C_LIGHT, L_B


G_CODATA2022 = 6.67430e-11  # m^3 kg^-1 s^-2; NIST SP961 (May 2024)
G_STANDARD_UNCERTAINTY = 1.5e-15
SOLAR_ANGULAR_MOMENTUM_SI = 1.8838e41  # kg m^2/s; Turyshev v2 p40 A39
SOLAR_POLE_RA_DEG = 286.13  # Park2021 Sec3.4, p7, following Eq39
SOLAR_POLE_DEC_DEG = 63.87


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024*1024), b''):
            digest.update(block)
    return digest.hexdigest()


def solar_axis():
    """Fixed J2000/ICRF north spin pole from Park2021 Sec3.4, p7."""
    ra, dec = np.deg2rad([SOLAR_POLE_RA_DEG, SOLAR_POLE_DEC_DEG])
    return np.array([np.cos(dec)*np.cos(ra), np.cos(dec)*np.sin(ra), np.sin(dec)])


def solar_spin_parameter(*, tdb_compatible=True):
    """Return GS vector, m^5/s^3, with an explicitly declared scale.

    We interpret the published SI angular momentum as an unscaled physical
    parameter. Under t'=k t+const, x'=k x, v'=v, invariance of GS/r^2
    requires (GS)'=k^2 GS. This is a scale interpretation for the diagnostic,
    not proof of the scale used by the source's original solar-spin estimate.
    The CLI also integrates the unscaled-parameter sensitivity separately.
    """
    k = 1-L_B
    return G_CODATA2022*SOLAR_ANGULAR_MOMENTUM_SI*solar_axis()*(k*k if tdb_compatible else 1.)


def spin_rate_components(displacement, earth_velocity, sun_velocity, gs_vector, *, c=C_LIGHT):
    """Return unscaled delta-F scalar/vector/combined channels (dimensionless).

    displacement=x_E-x_sun, all states SI and consistently scaled; gs_vector
    has units m^5/s^3. Multiply by A_SCALE only when forming d(TDB-TT)/dTDB.
    ``relative`` is the independent algebraic recombination, useful for closure.
    """
    r, ve, vs = (np.asarray(a, float) for a in (displacement, earth_velocity, sun_velocity))
    gs = np.asarray(gs_vector, float)
    if r.shape[-1:] != (3,) or ve.shape != r.shape or vs.shape != r.shape:
        raise ValueError('Displacement and velocities must have matching Cartesian shapes')
    if gs.shape != (3,) or not all(np.all(np.isfinite(a)) for a in (r, ve, vs, gs)):
        raise ValueError('Finite Cartesian states and GS vector required')
    if not np.isfinite(c) or c <= 0:
        raise ValueError('Positive finite speed of light required')
    radius = np.linalg.norm(r, axis=-1)
    if np.any(radius <= 0):
        raise ValueError('Earth and Sun must be distinct')
    cross = np.cross(r, gs)
    # A9: w_spin = -(r cross GS)/(2r^3); A8: D_spin = 2v_s dot(r cross GS)/r^3.
    w_spin = -cross/(2*radius[..., None]**3)
    d_spin = 2*np.sum(vs*cross, axis=-1)/radius**3
    vector = -4*np.sum(ve*w_spin, axis=-1)/c**4
    scalar = -d_spin/c**4
    relative = -2*np.sum(gs*np.cross(r, ve-vs), axis=-1)/(c**4*radius**3)
    return {'vector': vector, 'scalar': scalar, 'combined': vector+scalar, 'relative': relative}


class SolarSpinDiagnostic:
    """Real DE440 geometry; additive channels, no constants or reference fit."""
    channel_names = ('vector', 'scalar', 'combined', 'relative', 'unscaled_GS')

    def __init__(self, ephemeris):
        self.ephemeris = ephemeris
        self.gs = solar_spin_parameter()

    def __call__(self, seconds):
        state = self.ephemeris.states(seconds)
        if state.scale != 'TDB-compatible':
            raise ValueError('TDB-compatible states required')
        r = state.earth_position-state.body_positions['sun']
        components = spin_rate_components(r, state.earth_velocity, state.body_velocities['sun'], self.gs)
        # This channel only measures the explicitly ambiguous scale convention.
        unscaled = components['relative']/(1-L_B)**2
        return np.stack([components[n] for n in self.channel_names[:-1]]+[unscaled], axis=-1)

    def fingerprint(self):
        return {
            'model': 'Turyshev_v2_A5_A8_A9_A33_A34_solar_direct_spin_diagnostic_v1',
            'source_equations': ['Turyshev2507.22145v2 A5,A8,A9,A33,A34,A39'],
            'paper_sha256': sha256('docs/sources/turyshev_2507.22145v2.pdf'),
            'axis_source_sha256': sha256('docs/sources/park_2021_de440_de441.pdf'),
            'axis_source': 'Park2021 Sec3.4 p7 following Eq39',
            'G_source': 'https://physics.nist.gov/cuu/pdf/wall_2022.pdf',
            'G_SI': G_CODATA2022, 'G_standard_uncertainty_SI': G_STANDARD_UNCERTAINTY,
            'S_SI_kg_m2_per_s': SOLAR_ANGULAR_MOMENTUM_SI,
            'S_source': 'Turyshev v2 p40 A39 quotes its Ref45; underlying Ref45 not independently checked',
            'S_uncertainty': 'not supplied by the cited passage; no absolute physical precision claim',
            'axis_radec_degrees': [SOLAR_POLE_RA_DEG, SOLAR_POLE_DEC_DEG],
            'axis_ICRF': solar_axis().tolist(), 'GS_TDB_vector_m5_per_s3': self.gs.tolist(),
            'spin_scaling': '(GS)_TDB=(1-L_B)^2 G_SI S_SI; published S interpreted as unscaled physical parameter',
            'scale': 'TDB-compatible SI J2000/ICRF geometric', 'channels': self.channel_names,
            'bsp_sha256': sha256(self.ephemeris.path), 'gm_sha256': sha256(self.ephemeris.gm_path),
            'sources': {n: sha256(Path(__file__).with_name(n)) for n in
                        ('solar_spin.py', 'constants.py', 'timebase.py', 'ephemeris.py', 'extended_ephemeris.py', 'integration.py')},
            'limitations': ['not established as a term in DE440 time ephemeris', 'fixed solar pole and angular momentum',
                            'other external-body spins omitted', 'main Park/Phase1 models unchanged']}
