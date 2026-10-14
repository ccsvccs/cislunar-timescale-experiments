"""Independent SB441 small-body direct-potential diagnostic (NOT DE440 states).

Farnocchia IOM 392R-21-005 supplies the 373 target identities and masses.
These orbits were separately fitted and integrated with DE441 planets; the
module never labels them as the original DE440 integration trajectories.
"""
from pathlib import Path
import hashlib
import json
import re
import numpy as np
from jplephem.spk import SPK
from .constants import C_LIGHT
from .jpl_model import read_de440_constants
from .ephemeris import DAY_SECONDS, J2000_JD


def file_sha256(path):
    digest=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda:stream.read(1024*1024),b''):digest.update(chunk)
    return digest.hexdigest()


def build_mapping(iom_text_path='docs/sources/sb441_perturbers.txt',header_path='docs/sources/de440_header.440'):
    """Map all actual IOM identities to independently matching DE440 GMs.

    Asteroid MA#### identifiers match the numbered NAIF target minus 2000000.
    KBO MA8001–8030 lack names in the ASCII header, so the IOM's 30 distinct
    full-precision masses establish a one-to-one mapping, validated for
    uniqueness and relative agreement <2e-15. This is input matching, not a
    reference-time fit or a GM adjustment. Every returned GM uses the header.
    """
    h=read_de440_constants(header_path)
    entries=[]
    for line in Path(iom_text_path).read_text().splitlines():
        match=re.match(r'^\s*(.*?)\s+(\d{7})\s+(\d\.\d+[eE][+-]\d+)\s*$',line)
        if match:
            label,target,gm=match.groups()
            entries.append({'iom_label':label.strip(),'target':int(target),'iom_gm_au3_day2':float(gm)})
    if len(entries) != 373 or len({r['target'] for r in entries}) != 373:
        raise ValueError('IOM Table 1 must provide exactly 373 unique targets')
    ast_names={n for n in h if re.fullmatch(r'MA\d{4}',n) and 0<int(n[2:])<8000}
    kbo_names={f'MA{n}' for n in range(8001,8031)}
    if len(ast_names) != 343: raise ValueError('DE440 header must contain 343 asteroid masses')
    factor=149597870700.**3/DAY_SECONDS**2
    used=set()
    for row in entries:
        candidate=f"MA{row['target']-2000000:04d}"
        if candidate in ast_names:
            name=candidate;group='asteroids'
        else:
            matches=[n for n in kbo_names if abs(h[n]-row['iom_gm_au3_day2']) <= 2e-15*abs(h[n])]
            if len(matches) != 1: raise ValueError('KBO mass must uniquely match one header constant')
            name=matches[0];group='kbos'
        if name in used:raise ValueError('Header constant matched to multiple targets')
        if abs(h[name]-row['iom_gm_au3_day2']) > 2e-15*abs(h[name]):
            raise ValueError('IOM and DE440 header mass parameters disagree')
        used.add(name)
        row.update(header_name=name,group=group,gm_si=h[name]*factor)
    if used != ast_names|kbo_names:raise ValueError('Incomplete asteroid/KBO mapping')
    return {'source':'JPL Farnocchia IOM 392R-21-005 Table 1; DE440 GROUP1040/1041',
            'orbit_provenance':'SB441 separate fit/integration relative to DE441 planets; not original DE440 states',
            'header_sha256':file_sha256(header_path),'iom_text_sha256':file_sha256(iom_text_path),
            'units':'TDB-compatible AU^3/day^2 converted to SI with exact AU149597870700m',
            'entries':entries}


class SmallBodyPotential:
    """Real heliocentric SB441 positions plus DE440 geocenter; two direct rates.

    ``__call__`` returns [asteroid_U/c²,KBO_U/c²], no C or IAU scaling. Apply
    A_SCALE to its independent anchor integral when adding to TDB−TT. This
    is an explicitly separate diagnostic input model, not an automatic
    fallback in the Park candidate. Small-body vector/PN terms are omitted.
    """
    channel_names=('asteroids_c2','kbos_c2')
    def __init__(self,earth_ephemeris,path='data/sb441-n373s.bsp',*,mapping=None):
        self.ephemeris=earth_ephemeris;self.path=Path(path)
        if not self.path.is_file(): raise FileNotFoundError(f'Required actual SB441 kernel missing: {self.path}')
        self.mapping=mapping or build_mapping()
        self.kernel=SPK.open(str(path))
        self.segments={r['target']:[] for r in self.mapping['entries']}
        for segment in self.kernel.segments:
            if segment.center != 10 or segment.frame != 1:
                self.kernel.close();raise ValueError('SB441 input must be Sun-relative ICRF geometric segments')
            if segment.target in self.segments:self.segments[segment.target].append(segment)
        if any(not pieces for pieces in self.segments.values()):
            self.kernel.close();raise ValueError('SB441 kernel must contain all 373 IOM targets')
        bounds=[]
        for pieces in self.segments.values():
            ordered=sorted(pieces,key=lambda s:s.start_jd)
            if any(b.start_jd > a.end_jd for a,b in zip(ordered,ordered[1:])):
                self.kernel.close();raise ValueError('Gap in a small-body split-segment coverage')
            bounds.append((ordered[0].start_jd,max(s.end_jd for s in ordered)))
        self.start_jd=max(a for a,b in bounds)
        self.end_jd=min(b for a,b in bounds)
        self.gm=np.array([r['gm_si'] for r in self.mapping['entries']])
        self.groups=np.array([r['group']=='kbos' for r in self.mapping['entries']])

    def diagnostics(self,seconds,*,by_body=False):
        seconds=np.asarray(seconds,float)
        days=np.floor(seconds/DAY_SECONDS)
        jd1=J2000_JD+days;jd2=(seconds-days*DAY_SECONDS)/DAY_SECONDS
        state=self.ephemeris.states_jd(jd1,jd2)
        earth_heliocentric=state.earth_position-state.body_positions['sun']
        elapsed=(jd1-self.start_jd)+jd2
        if np.any(elapsed < 0) or np.any(elapsed > self.end_jd-self.start_jd):
            raise ValueError('Epoch outside common SB441 coverage')
        totals=[np.zeros(seconds.shape),np.zeros(seconds.shape)]
        bodies={}
        for row in self.mapping['entries']:
            xyz=self._heliocentric_position(row['target'],jd1,jd2)
            distance=np.linalg.norm(xyz-earth_heliocentric,axis=-1)
            if np.any(distance <= 0) or np.any(~np.isfinite(distance)):
                raise ValueError('Invalid small-body distance')
            value=row['gm_si']/distance
            totals[int(row['group']=='kbos')]+=value
            if by_body:bodies[row['header_name']]=value
        result={'potential_asteroids':totals[0],'potential_kbos':totals[1],
                'rates':np.stack(totals,axis=-1)/C_LIGHT**2,
                'provenance':self.mapping['orbit_provenance']}
        if by_body:result['potential_by_body']=bodies
        return result

    def _heliocentric_position(self,target,jd1,jd2):
        first,second=jd1.ravel(),jd2.ravel()
        output=np.empty((len(first),3));covered=np.zeros(len(first),bool)
        # Later SPK segments have priority, including a shared boundary.
        for segment in reversed(self.segments[target]):
            offset=(first-segment.start_jd)+second
            use=(~covered)&(offset>=0)&(offset<=segment.end_jd-segment.start_jd)
            if np.any(use):
                output[use]=segment.compute(first[use],second[use])[:3].T*1000
                covered[use]=True
        if not np.all(covered):raise ValueError('Epoch outside small-body segment coverage')
        return output.reshape(jd1.shape+(3,))

    def __call__(self,seconds):return self.diagnostics(seconds)['rates']

    def fingerprint(self):
        return {'model':'SB441373_direct_scalar_diagnostic_v1',
                'bsp_sha256':file_sha256(self.path),'mapping':self.mapping,
                'sources':{name:file_sha256(Path(__file__).with_name(name)) for name in
                    ('small_body.py','ephemeris.py','constants.py','timebase.py')},
                'c_light_si':C_LIGHT,'earth_bsp_sha256':file_sha256(self.ephemeris.path),
                'missing':'small-body vector, PN scalar, mixed c^-4 terms; DE440 original small-body trajectories',
                'channels':self.channel_names,'frame':'ICRF heliocentric geometric; subtract DE440 Earth-Sun vector; SI TDB-compatible'}

    def close(self):self.kernel.close()
    def __enter__(self):return self
    def __exit__(self,*args):self.close()
