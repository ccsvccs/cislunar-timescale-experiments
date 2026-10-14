"""DE440 state derivatives for Park (2021) Eqs. (6)–(7).

The Phase 1 reader is retained unchanged. Acceleration is the analytic second
TDB derivative of the same Type-2 SPK position polynomial, not a finite
-difference or a replacement Newtonian force law. Units are SI throughout.
"""
from dataclasses import dataclass
import numpy as np
from .ephemeris import DE440Ephemeris, BODY_IDS, DAY_SECONDS, J2000_JD


@dataclass(frozen=True)
class ExtendedEphemerisState:
    earth_position: np.ndarray
    earth_velocity: np.ndarray
    body_positions: dict
    gm: dict
    body_velocities: dict
    body_accelerations: dict
    earth_gm: float
    scale: str = 'TDB-compatible'


def type2_state_derivatives(segment, jd1, jd2=0.0):
    """Return position, velocity, acceleration in km, km/s, km/s².

    Public ``load_array`` exposes increasing-order coefficients with shape
    (3, records, degree+1). Record lookup keeps the input date in two parts.
    A Clenshaw recurrence differentiated twice evaluates all derivatives.
    """
    if segment.data_type != 2:
        raise ValueError('Analytic derivative reader requires an SPK Type 2 segment')
    jd1, jd2 = np.broadcast_arrays(np.asarray(jd1, float), np.asarray(jd2, float))
    shape = jd1.shape
    first, second = jd1.ravel(), jd2.ravel()
    initial_jd, interval_days, coefficients = segment.load_array()
    whole, fractional = np.divmod(first-initial_jd, interval_days)
    carry, offset = np.divmod(fractional+second, interval_days)
    indices = (whole+carry).astype(np.int64)
    records = coefficients.shape[1]
    at_end = (indices == records) & (offset == 0)
    indices[at_end] -= 1
    offset[at_end] = interval_days
    if np.any(indices < 0) or np.any(indices >= records):
        raise ValueError('Epoch outside SPK segment coverage')
    x = 2*offset/interval_days-1
    selected = coefficients[:, indices, :]
    b1 = np.zeros(selected.shape[:2]); b2 = b1.copy()
    d1 = b1.copy(); d2 = b1.copy(); dd1 = b1.copy(); dd2 = b1.copy()
    for k in range(selected.shape[-1]-1, 0, -1):
        b = selected[..., k] + 2*x*b1-b2
        d = 2*b1 + 2*x*d1-d2
        dd = 4*d1 + 2*x*dd1-dd2
        b2, b1 = b1, b
        d2, d1 = d1, d
        dd2, dd1 = dd1, dd
    position = selected[..., 0] + x*b1-b2
    velocity = b1 + x*d1-d2
    acceleration = 2*d1 + x*dd1-dd2
    factor = 2/(interval_days*DAY_SECONDS)
    return tuple(value.T.reshape(shape+(3,))*scale for value, scale in
                 ((position, 1.), (velocity, factor), (acceleration, factor**2)))


class ExtendedDE440Ephemeris(DE440Ephemeris):
    """Known major-body Park candidate input; missing small bodies explicit."""
    def states_jd(self, jd1, jd2=0.0):
        jd1, jd2 = np.broadcast_arrays(np.asarray(jd1, float), np.asarray(jd2, float))
        elapsed = (jd1-self.start_jd)+jd2
        if np.any(~np.isfinite(elapsed)) or np.any(elapsed < 0) or np.any(elapsed > self.end_jd-self.start_jd):
            raise ValueError('TDB epoch outside DE440 coverage')
        def state(center, target):
            return tuple(x*1000 for x in type2_state_derivatives(self.kernel[center,target], jd1, jd2))
        emb = state(0,3)
        earth = tuple(a+b for a,b in zip(emb,state(3,399)))
        bodies = {name: state(0,body) for name,body in BODY_IDS.items() if name != 'moon'}
        bodies['moon'] = tuple(a+b for a,b in zip(emb,state(3,301)))
        return ExtendedEphemerisState(earth[0],earth[1],{n:s[0] for n,s in bodies.items()},
              self.gm,{n:s[1] for n,s in bodies.items()},{n:s[2] for n,s in bodies.items()},self.all_gm[399])
