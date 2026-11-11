"""Separate leading lunar degree-2 external-potential diagnostic at Earth.

J2M/C22M and AM come from the DE440 header, with the actual DE440 PA
orientation from NAIF's binary PCK. This is a STATIC coefficient experiment,
not the full time-dependent DE440 lunar gravity model (Park 2021 Eqs. 30--34,
54). C21M/S21M/S22M are exposed only as literal-header sensitivity channels:
their static interpretation is not established by the public header.
No time-ephemeris reference or fitted axis enters this module.
"""
from pathlib import Path
import numpy as np
from jplephem.pck import PCK
from .constants import C_LIGHT, A_SCALE
from .jpl_model import read_de440_constants
from .small_body import file_sha256


def degree2_potential(displacement_pa, gm, radius, *, j2, c21=0., s21=0.,
                      c22=0., s22=0.):
    """Return five unnormalized degree-2 potential channels, in m²/s².

    Order: J2, C21, S21, C22, S22. C20=-J2; geodesy Legendre convention
    without Condon--Shortley phase: P21=3 sin(phi) cos(phi), P22=3 cos²(phi).
    The Cartesian expression has no longitude singularity at the poles.
    """
    r = np.asarray(displacement_pa, dtype=float)
    parameters = np.array([gm, radius, j2, c21, s21, c22, s22], dtype=float)
    if r.shape[-1:] != (3,) or not np.all(np.isfinite(r)):
        raise ValueError('Finite Cartesian displacement required')
    distance = np.linalg.norm(r, axis=-1)
    if not np.all(np.isfinite(parameters)) or gm <= 0 or radius <= 0:
        raise ValueError('Positive finite GM/radius and finite coefficients required')
    if np.any(distance <= radius):
        raise ValueError('Lunar external-potential diagnostic requires r > R')
    x, y, z = np.moveaxis(r / distance[..., None], -1, 0)
    angular = np.stack((-j2 * (3*z*z-1)/2, 3*c21*x*z, 3*s21*y*z,
                        3*c22*(x*x-y*y), 6*s22*x*y), axis=-1)
    return (gm * radius**2 / distance**3)[..., None] * angular


class LunarPAOrientation:
    """Vectorized original type-2 PCK Euler angles, checked against CSPICE.

    For this DE440 kernel, R(J2000 -> PA) = R3(psi) R1(theta) R3(phi),
    where the R matrices are passive rotations as in Park Eq. (8)--(11).
    Validation independently uses the official FK and spice.pxform.
    """
    frame_name = 'MOON_PA_DE440'
    frame_id = 31008

    def __init__(self, path='data/moon_pa_de440_200625.bpc'):
        self.path = Path(path)
        if not self.path.is_file():
            raise FileNotFoundError(f'Required authentic DE440 lunar PA PCK missing: {self.path}')
        self.kernel = PCK.open(str(self.path))
        self.segments = self.kernel.segments
        if not self.segments or any((s.body, s.frame, s.data_type) != (31008, 1, 2)
                                    for s in self.segments):
            self.close()
            raise ValueError('Expected DE440 PA frame 31008, ICRF/J2000 base, type-2 PCK')
        ordered = sorted(self.segments, key=lambda s: s.initial_second)
        if any(b.initial_second > a.final_second for a, b in zip(ordered, ordered[1:])):
            self.close()
            raise ValueError('Gap in DE440 PA orientation coverage')
        self.start = min(s.initial_second for s in self.segments)
        self.end = max(s.final_second for s in self.segments)

    def matrices(self, seconds):
        times = np.asarray(seconds, dtype=float)
        flat = times.ravel()
        if np.any(~np.isfinite(flat)) or np.any(flat < self.start) or np.any(flat > self.end):
            raise ValueError('Epoch outside authentic DE440 PA orientation coverage')
        angles = np.empty((len(flat), 3))
        covered = np.zeros(len(flat), dtype=bool)
        days = np.floor(flat/86400.)
        first, second = 2451545.+days, (flat-days*86400.)/86400.
        for segment in reversed(self.segments):
            use = (~covered) & (flat >= segment.initial_second) & (flat <= segment.final_second)
            if np.any(use):
                angles[use] = segment.compute(first[use], second[use], derivative=False).T
                covered[use] = True
        if not np.all(covered):
            raise ValueError('Uncovered DE440 PA orientation epoch')
        phi, theta, psi = angles.T
        cp, sp = np.cos(phi), np.sin(phi)
        ct, st = np.cos(theta), np.sin(theta)
        cw, sw = np.cos(psi), np.sin(psi)
        matrices = np.empty((len(flat), 3, 3))
        matrices[:, 0, :] = np.stack((cw*cp-sw*ct*sp, cw*sp+sw*ct*cp, sw*st), axis=-1)
        matrices[:, 1, :] = np.stack((-sw*cp-cw*ct*sp, -sw*sp+cw*ct*cp, cw*st), axis=-1)
        matrices[:, 2, :] = np.stack((st*sp, -st*cp, ct), axis=-1)
        return matrices.reshape(times.shape+(3, 3))

    def close(self):
        self.kernel.close()


class LunarMultipolePotential:
    """Five additive U_2/c² channels; apply A_SCALE only after integration.

    The primary documented static diagnostic is J2+C22. The other three
    channels are explicit literal-header sensitivities. No default model is
    changed and no tide/rotation correction is fabricated.
    """
    channel_names = ('lunar_J2_c2', 'lunar_C21_literal_c2', 'lunar_S21_literal_c2',
                     'lunar_C22_c2', 'lunar_S22_literal_c2')
    primary_indices = (0, 3)

    def __init__(self, ephemeris, *, header_path='docs/sources/de440_header.440',
                 pck_path='data/moon_pa_de440_200625.bpc',
                 fk_path='data/moon_de440_250416.tf'):
        self.ephemeris = ephemeris
        self.header_path, self.fk_path = Path(header_path), Path(fk_path)
        if not self.fk_path.is_file():
            raise FileNotFoundError(f'Required official DE440 frame kernel missing: {self.fk_path}')
        self.orientation = LunarPAOrientation(pck_path)
        h = read_de440_constants(self.header_path)
        self.parameters = {name: h[name] for name in ('AM', 'J2M', 'C21M', 'S21M', 'C22M', 'S22M')}
        self.radius = h['AM']*1000
        self.coefficients = {k: h[v] for k, v in
                             (('j2', 'J2M'), ('c21', 'C21M'), ('s21', 'S21M'),
                              ('c22', 'C22M'), ('s22', 'S22M'))}
        self.gm = self.ephemeris.gm['moon']

    def displacement_pa(self, seconds):
        times = np.asarray(seconds, dtype=float)
        flat = times.ravel()
        days = np.floor(flat/86400.)
        first, second = 2451545.+days, (flat-days*86400.)/86400.
        # Subtract two EMB-relative segments directly: avoid cancellation
        # from subtracting two ~AU-sized barycentric position vectors.
        r = (self.ephemeris.kernel[3, 399].compute(first, second)
             - self.ephemeris.kernel[3, 301].compute(first, second)).T*1000.
        r = r.reshape(times.shape+(3,))
        return np.einsum('...ij,...j->...i', self.orientation.matrices(times), r)

    def __call__(self, seconds):
        return degree2_potential(self.displacement_pa(seconds), self.gm, self.radius,
                                 **self.coefficients)/C_LIGHT**2

    def fingerprint(self):
        return {'model': 'DE440_PA_static_lunar_degree2_diagnostic_v1',
                'source_sha256': {name: file_sha256(Path(__file__).with_name(name)) for name in
                                 ('lunar_multipole.py', 'constants.py', 'timebase.py',
                                  'integration.py', 'ephemeris.py', 'jpl_model.py', 'small_body.py')},
                'c_light_si': C_LIGHT, 'A_SCALE': A_SCALE,
                'earth_moon_bsp_sha256': file_sha256(self.ephemeris.path),
                'gm_kernel_sha256': file_sha256(self.ephemeris.gm_path),
                'header_sha256': file_sha256(self.header_path),
                'pck_sha256': file_sha256(self.orientation.path),
                'fk_sha256': file_sha256(self.fk_path),
                'parameters': self.parameters, 'moon_gm_si': self.gm,
                'radius_m': self.radius, 'channels': self.channel_names,
                'primary_channels': ['lunar_J2_c2', 'lunar_C22_c2'],
                'normalization': 'unnormalized; C20=-J2; no Condon-Shortley phase',
                'coefficient_scope': 'constant header J2M/C22M; C21M/S21M/S22M literal-header sensitivities only',
                'frame': 'ICRF/J2000 geometric -> MOON_PA_DE440 (31008); never ME or IAU_MOON',
                'scale': 'SI, TDB-compatible DE440 GM/positions; header AM in km; radius unscaled diagnostic convention',
                'radius_scale_effect': 'an O(LB) relative radius scaling ambiguity is below 1e-7 ps for this experiment',
                'orientation_coverage_tdb_seconds': [self.orientation.start, self.orientation.end],
                'excluded': ['time-dependent tidal/rotational degree-2 distortion',
                             'degree >=3', 'multipole PN scalar/vector/mixed c^-4 terms'],
                'reference_used_in_integrand': False}

    def close(self):
        self.orientation.close()

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()
