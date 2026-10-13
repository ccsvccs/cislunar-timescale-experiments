"""Read JPL's genuine numerical TT−TDB reference from de440t.bsp.

The special SPK X coordinate is seconds, not kilometres. See the official
https://ssd.jpl.nasa.gov/ftp/eph/planets/bsp/README.txt . This module is strictly
a validation reference and never supplies the core relativistic integrand.
"""
from pathlib import Path
import os
import numpy as np
from jplephem.spk import SPK


class JPLTimeEphemeris:
    """Special SPK segment 1000000000→1000000001; output sign is TDB−TT."""
    def __init__(self, path=None):
        self.path = Path(path or os.environ.get('JPL_TIME_BSP', 'data/de440t.bsp'))
        if not self.path.is_file():
            raise FileNotFoundError(f'Genuine numerical JPL time reference missing: {self.path}')
        self.kernel = SPK.open(str(self.path))
        try:
            self.segment = self.kernel[1000000000, 1000000001]
        except KeyError:
            self.kernel.close()
            raise ValueError('This SPK does not contain the JPL TT−TDB time segment') from None
        self.start_jd = self.segment.start_jd
        self.end_jd = self.segment.end_jd

    def evaluate_jd(self, jd1, jd2=0.):
        """Return (TDB−TT seconds, d(TDB−TT)/dTDB dimensionless)."""
        first, second = np.broadcast_arrays(np.asarray(jd1, dtype=float), np.asarray(jd2, dtype=float))
        elapsed = first - self.start_jd + second
        if np.any(~np.isfinite(elapsed)) or np.any(elapsed < 0.) or np.any(elapsed > self.end_jd-self.start_jd):
            raise ValueError('Epoch outside numerical JPL time reference coverage')
        shape = first.shape
        value, derivative = self.segment.compute_and_differentiate(first.reshape(-1), second.reshape(-1))
        # jplephem differentiates with respect to JD (days); SPICE VX is /s.
        delta = -value[0].reshape(shape)
        rate = -derivative[0].reshape(shape) / 86400.
        if shape == ():
            return float(delta), float(rate)
        return delta, rate

    def evaluate(self, seconds_since_j2000):
        seconds = np.asarray(seconds_since_j2000, dtype=float)
        days = np.floor(seconds/86400.)
        return self.evaluate_jd(2451545.+days, (seconds-days*86400.)/86400.)

    def tdb_minus_tt(self, epoch):
        """Accept project TwoPartJD or an explicit (jd1, jd2) pair."""
        parts = epoch.jd_parts() if hasattr(epoch, 'jd_parts') else epoch
        return self.evaluate_jd(*parts)[0]

    def close(self):
        self.kernel.close()

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()
