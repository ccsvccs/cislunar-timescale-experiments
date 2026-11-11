"""Read the fixed 20250620 LTE440 reference product, for validation only.

The special SPK X is seconds, not km. Its companion TPC linear term is
mandatory. This module is never imported by the physical lunar integrand.
Author product: Lu Xu, https://github.com/xlucn/LTE440/releases/tag/20250620 .
"""
from pathlib import Path
import re
import numpy as np
from jplephem.spk import SPK
from .constants import T0_TDB
from .small_body import file_sha256


class LTE440Reference:
    bsp_sha256 = 'f48c36b0c1b872f51616790db753538b4895bc2ceae73386a72eb0150e44f1e8'
    tpc_sha256 = 'aa897dd4aeba15dfeb6f80883ad6b3e2d7c4f098690b02cdeb021167f5734c92'

    def __init__(self, path='data/lte440.bsp', tpc_path='data/lte440.tpc', *, verify_hash=True):
        self.path, self.tpc_path = Path(path), Path(tpc_path)
        for p in (self.path, self.tpc_path):
            if not p.is_file():
                raise FileNotFoundError(f'Fixed original LTE440 reference missing: {p}')
        self.actual_bsp_sha256 = file_sha256(self.path)
        self.actual_tpc_sha256 = file_sha256(self.tpc_path)
        if verify_hash and (self.actual_bsp_sha256 != self.bsp_sha256 or self.actual_tpc_sha256 != self.tpc_sha256):
            raise ValueError('LTE440 reference does not match fixed 20250620 release SHA256')
        active, data = False, []
        for line in self.tpc_path.read_text().splitlines():
            if '\\begindata' in line:
                active = True
            elif '\\begintext' in line:
                active = False
            elif active:
                data.append(line)
        rates = re.findall(r'BODY1000000005_RATE\s*=\s*\(\s*([+\-\d.DEde]+)\s*\)', '\n'.join(data))
        if len(rates) != 1:
            raise ValueError('Exactly one original LTE440 companion TPC rate is required')
        self.rate = float(rates[0].replace('D', 'E').replace('d', 'e'))
        self.kernel = SPK.open(str(self.path))
        if len(self.kernel.segments) != 1:
            self.close()
            raise ValueError('Fixed LTE440 release must contain one time segment')
        self.segment = self.kernel.segments[0]
        if (self.segment.center, self.segment.target, self.segment.frame, self.segment.data_type) != (1000000000, 1000000005, 1, 2):
            self.close()
            raise ValueError('Invalid LTE440 TCL-TDB segment identity/frame/type')
        self.start_jd, self.end_jd = self.segment.start_jd, self.segment.end_jd

    def evaluate_jd(self, jd1, jd2=0.):
        """Return full (TCL−TDB seconds, d(TCL−TDB)/dTDB)."""
        first, second = np.broadcast_arrays(np.asarray(jd1, float), np.asarray(jd2, float))
        offset = (first-self.start_jd)+second
        if np.any(~np.isfinite(offset)) or np.any(offset < 0) or np.any(offset > self.end_jd-self.start_jd):
            raise ValueError('Epoch outside fixed LTE440 reference coverage')
        values, derivative = self.segment.compute_and_differentiate(first.ravel(), second.ravel())
        # Exact project anchor in two parts; the BSP's polynomial is seconds.
        a, b = T0_TDB.jd_parts()
        elapsed = ((first.astype(np.longdouble)-np.longdouble(a))
                   +(second.astype(np.longdouble)-np.longdouble(b)))*np.longdouble(86400)
        value = values[0].reshape(first.shape).astype(np.longdouble)+np.longdouble(self.rate)*elapsed
        rate = derivative[0].reshape(first.shape)/86400.+self.rate
        return value, rate

    def evaluate(self, seconds):
        times = np.asarray(seconds, float)
        days = np.floor(times/86400.)
        return self.evaluate_jd(2451545.+days, (times-days*86400.)/86400.)

    def fingerprint(self):
        return {'classification': 'independent published reference; never integrand input',
                'release_url': 'https://github.com/xlucn/LTE440/releases/tag/20250620',
                'bsp_sha256': self.actual_bsp_sha256, 'tpc_sha256': self.actual_tpc_sha256,
                'source_sha256': file_sha256(__file__), 'rate': self.rate,
                'anchor_jd_parts': T0_TDB.jd_parts(),
                'dependency_source_sha256': {n: file_sha256(Path(__file__).with_name(n)) for n in
                                            ('constants.py', 'timebase.py', 'small_body.py')},
                'segment': {'center': self.segment.center, 'target': self.segment.target,
                            'frame': self.segment.frame, 'type': self.segment.data_type,
                            'start_jd': self.start_jd, 'end_jd': self.end_jd},
                'unit': 'SPK X seconds; derivative /day converted /second; TPC linear term restored'}

    def close(self):
        self.kernel.close()

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()
