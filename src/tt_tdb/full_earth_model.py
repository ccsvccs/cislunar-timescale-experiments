"""Explicit full-coverage Earth candidate composed from unchanged audited cores.

Park major bodies + published ring leading potential + actual SB441 leading
potentials + static DE440 lunar J2/C22. Direct solar spin is excluded.
The solar-equator J2 axis is primary; a separate ecliptic-axis difference is
diagnostic only. SB441 is a DE441-refitted input, not original DE440 orbits.
"""
import json
from pathlib import Path
import numpy as np
import spiceypy as spice
from .constants import C_LIGHT
from .ephemeris import DE440Ephemeris
from .extended_ephemeris import ExtendedDE440Ephemeris
from .jpl_model import JPLModel, solar_j2_potential
from .lunar_multipole import LunarMultipolePotential
from .pn_rate import PNModel
from .small_body import SmallBodyPotential, file_sha256


def anchored_prefix(increments, origin_index, anchor_to_origin):
    """Neumaier accumulation left/right from a fixed-anchor split point.

    increments[i] is the actual model integral from grid[i] to grid[i+1].
    anchor_to_origin is the independently integrated physical anchor ->
    grid[origin_index]. No reference value or inferred offset is accepted.
    """
    increments=np.asarray(increments,float)
    offset=np.asarray(anchor_to_origin,float)
    if increments.ndim!=2 or offset.shape!=(increments.shape[1],):
        raise ValueError('Matched two-dimensional increments and channel offset required')
    if not 0<=origin_index<=len(increments) or not np.all(np.isfinite(increments)) or not np.all(np.isfinite(offset)):
        raise ValueError('Finite integrals and an in-grid origin required')
    values=np.empty((len(increments)+1,increments.shape[1]))
    values[origin_index]=offset
    for direction in (1,-1):
        total=offset.copy();compensation=np.zeros_like(total)
        indices=range(origin_index,len(increments)) if direction==1 else range(origin_index-1,-1,-1)
        for index in indices:
            increment=direction*increments[index]
            updated=total+increment
            compensation+=np.where(abs(total)>=abs(increment),(total-updated)+increment,(increment-updated)+total)
            total=updated
            values[index+(1 if direction==1 else 0)]=total+compensation
    return values


class FullEarthModel:
    """Independent additive integrand channels; no reference-time imports.

    ``fast_channels`` and ``small_channels`` can use independently verified
    quadratures. Their combined integral uses exactly the original IAU anchor.
    Rates are unscaled: multiply by A_SCALE after integration/evaluation.
    """
    fast_names=('baseline_F2_minus_C', 'baseline_F4', *JPLModel.correction_names,
                'lunar_J2_c2', 'lunar_C22_c2', 'ecliptic_minus_equator_J2_c2')
    small_names=SmallBodyPotential.channel_names
    channel_names=fast_names+small_names
    axis_index=fast_names.index('ecliptic_minus_equator_J2_c2')
    primary_indices=tuple(i for i,n in enumerate(channel_names) if n!='ecliptic_minus_equator_J2_c2')

    def __init__(self):
        self.ephemeris=ExtendedDE440Ephemeris()
        self.baseline_ephemeris=DE440Ephemeris()
        self.baseline=PNModel(self.baseline_ephemeris)
        self.park=JPLModel(self.ephemeris)
        self.small=SmallBodyPotential(self.baseline_ephemeris)
        self.lunar=LunarMultipolePotential(self.baseline_ephemeris)
        self.ecliptic_axis=spice.pxform('ECLIPJ2000','J2000',0.)[:,2]
        self.start_jd=max(self.ephemeris.start_jd,self.small.start_jd,
                          2451545.+self.lunar.orientation.start/86400.)
        self.end_jd=min(self.ephemeris.end_jd,self.small.end_jd,
                        2451545.+self.lunar.orientation.end/86400.)

    def fast_channels(self, seconds):
        times=np.asarray(seconds,float)
        terms=self.park.correction_rates(times)
        lunar=self.lunar(times)
        state=self.baseline_ephemeris.states(times)
        ecliptic=solar_j2_potential(state.earth_position,state.body_positions['sun'],state.gm['sun'],
                                  self.park.solar_j2,self.park.solar_radius,self.ecliptic_axis)/C_LIGHT**2
        axis_difference=ecliptic-terms['solar_j2_c2']
        return np.column_stack((self.baseline(times),*(terms[n] for n in self.park.correction_names),
                                lunar[...,0],lunar[...,3],axis_difference))

    def small_channels(self, seconds):
        return self.small(seconds)

    def __call__(self, seconds):
        return np.column_stack((self.fast_channels(seconds),self.small_channels(seconds)))

    def fingerprint(self):
        return {'model':'Earth_Park_ring_SB441_lunar_static_J2_C22_no_solar_spin_v1',
                'park':json.loads(self.park.fingerprint()), 'baseline':json.loads(self.baseline.fingerprint()),
                'small_body':self.small.fingerprint(), 'lunar_static':self.lunar.fingerprint(),
                'source_sha256':file_sha256(Path(__file__)),
                'channels':self.channel_names,'primary_indices':self.primary_indices,'axis_index':self.axis_index,
                'primary_solar_J2_axis':'solar equator, Park2021 Sec3.4 pole; no selection using residuals',
                'secondary_axis':'SPICE ECLIPJ2000; separately labelled diagnostic only',
                'solar_direct_spin_included':False,
                'lunar_static_primary_coefficients':['J2M','C22M'],
                'limitations':['SB441 separate DE441-referenced small-body fit, not original DE440 trajectories',
                              'published ring geometry only; no original ring dynamical states',
                              'small-body/ring vector, PN scalar and mixed c^-4 terms omitted',
                              'lunar time-variable distortion, other degree-2 coefficients and degree>=3 omitted',
                              'GMSDOT time-potential convention unresolved'],
                'reference_time_used_in_integrand':False}

    def close(self):
        self.lunar.close();self.small.close();self.ephemeris.close();self.baseline_ephemeris.close()

    def __enter__(self):return self
    def __exit__(self,*args):self.close()
