"""Turyshev Eq. (21) geocentric c^-2 and c^-4 rates, no SOFA dependency."""
from pathlib import Path
import hashlib
import json
import numpy as np
from .constants import C_LIGHT, C_EXACT, L_C_PAPER, A_SCALE, TDB0, T0_JD_DECIMAL
from .potential import external_potential


def rate_c2(velocity_squared, potential):
    """F2=(v²/2+U)/c², dimensionless; arguments use SI."""
    return (0.5*np.asarray(velocity_squared)+np.asarray(potential))/C_LIGHT**2


def rate_c4(velocity_squared, potential):
    """F4=(v⁴/8+3v²U/2−U²/2)/c⁴, this paper's truncated model."""
    v2, u = np.asarray(velocity_squared), np.asarray(potential)
    return (v2*v2/8.0+1.5*v2*u-0.5*u*u)/C_LIGHT**4


def relativistic_rate(velocity_squared, potential, *, include_c4=True):
    return rate_c2(velocity_squared, potential)+(rate_c4(velocity_squared,potential) if include_c4 else 0.0)


class PNModel:
    """Geocentric finite-monopole Eq. (21) model on TDB-compatible DE440.

    Observer displacement from Earth is exactly zero; Earth's barycentric
    position and velocity are the nonzero real ephemeris states. Missing
    asteroid/KBO potentials and non-monopole terms are documented explicitly.
    """
    def __init__(self, ephemeris):
        self.ephemeris = ephemeris

    def diagnostics(self, seconds):
        state = self.ephemeris.states(seconds)
        if state.scale != 'TDB-compatible':
            raise ValueError('Model requires consistently TDB-compatible states and GMs')
        u, bodies = external_potential(state.earth_position, state.body_positions, state.gm, breakdown=True)
        v2 = np.sum(state.earth_velocity**2, axis=-1)
        f2, f4 = rate_c2(v2,u), rate_c4(v2,u)
        return {'velocity_squared':v2, 'potential':u, 'potential_by_body':bodies,
                'sun_potential': bodies.get('sun',np.zeros_like(u)),
                'other_potential': u-bodies.get('sun',np.zeros_like(u)),
                'F2':f2, 'F4':f4, 'F4_over_F2':f4/f2,
                'C_exact':C_EXACT, 'L_C_paper':L_C_PAPER,
                'residual_rate':f2-C_EXACT+f4, 'observer_position_term_seconds':0.0}

    def external_potential(self, seconds, *, breakdown=False):
        state = self.ephemeris.states(seconds)
        return external_potential(state.earth_position,state.body_positions,state.gm,breakdown=breakdown)

    def __call__(self, seconds):
        """Independent channels [F2-C_exact,F4], integrated with respect to TDB seconds."""
        d = self.diagnostics(seconds)
        return np.stack((d['F2']-C_EXACT,d['F4']),axis=-1)

    def fingerprint(self):
        """Actual source and data content signatures, not just path names."""
        def sha(path):
            digest=hashlib.sha256()
            with Path(path).open('rb') as f:
                for chunk in iter(lambda:f.read(1024*1024),b''):
                    digest.update(chunk)
            return digest.hexdigest()
        sources = {name:sha(Path(__file__).with_name(name)) for name in
                   ('pn_rate.py','potential.py','constants.py','ephemeris.py','timebase.py')}
        return json.dumps({'model':'Turyshev2507.22145v2_eq21_geocenter',
                           'bsp_sha256':sha(self.ephemeris.path),
                           'gm_sha256':sha(self.ephemeris.gm_path),'gm':self.ephemeris.gm,
                           'sources':sources, 'scale':'TDB-compatible SI ICRF geometric',
                           'C_exact':C_EXACT,'A':A_SCALE,'T0':T0_JD_DECIMAL,'TDB0':TDB0,
                           'missing':self.ephemeris.missing_potential_sources},sort_keys=True)
