"""Geocentric TDB ↔ TT, Eq. (21), using an independent residual integral."""
from dataclasses import dataclass, asdict
import numpy as np
from .constants import A_SCALE,TDB0,T0_TDB
from .timebase import TwoPartJD


@dataclass(frozen=True)
class IterationInfo:
    initial_delta_seconds: float
    iterations: int
    tolerance_seconds: float
    final_delta_seconds: float
    implicit_residual_seconds: float
    final_step_seconds: float


class TTTDBConverter:
    """Two-part external epochs, with small delta evaluated on TDB arguments.

    ``integrator`` must start at T0_TDB and integrate [F2-C_exact,F4].
    No reference-derived offset or linear trend is introduced.
    """
    def __init__(self, integrator, *, include_c4=True):
        if integrator.anchor != T0_TDB.seconds_since_j2000():
            raise ValueError('Eq. (21) requires the exact defined 1977 integral anchor')
        if integrator.channels != 2:
            raise ValueError('Expected independent channels [F2-C_exact,F4]')
        self.integrator=integrator
        self.include_c4=include_c4
        self.last_iteration=None

    def deltas_at_seconds(self, seconds, *, include_c4=None):
        channels=self.integrator.integral(seconds)
        enabled=self.include_c4 if include_c4 is None else include_c4
        return TDB0+A_SCALE*(channels[...,0]+(channels[...,1] if enabled else 0.))

    def tdb_minus_tt(self,tdb_epoch:TwoPartJD, *, include_c4=None)->float:
        if not isinstance(tdb_epoch,TwoPartJD):
            raise TypeError('Use TwoPartJD, never a summed complete float JD')
        return float(self.deltas_at_seconds(tdb_epoch.seconds_since_j2000(),include_c4=include_c4))

    def tdb_to_tt(self,tdb_epoch:TwoPartJD)->TwoPartJD:
        """TT=TDB−Δ(TDB), apply the seconds correction to the small fraction."""
        return tdb_epoch.add_seconds(-self.tdb_minus_tt(tdb_epoch))

    def tt_to_tdb(self,tt_epoch:TwoPartJD, *, tolerance_seconds=1e-14,max_iterations=12)->TwoPartJD:
        """Solve TDB=TT+Δ(TDB) with measured fixed-point convergence."""
        if tolerance_seconds<=0 or max_iterations<1:
            raise ValueError('Positive tolerance and iteration limit required')
        initial=self.tdb_minus_tt(tt_epoch)
        candidate=tt_epoch.add_seconds(initial)
        for iteration in range(1,max_iterations+1):
            delta=self.tdb_minus_tt(candidate)
            updated=tt_epoch.add_seconds(delta)
            step=updated.difference_seconds(candidate)
            candidate=updated
            if abs(step)<=tolerance_seconds:
                final=self.tdb_minus_tt(candidate)
                residual=candidate.difference_seconds(tt_epoch)-final
                self.last_iteration=IterationInfo(initial,iteration,tolerance_seconds,final,residual,step)
                if abs(residual)>tolerance_seconds:
                    raise ArithmeticError('Iteration step small but implicit equation residual too large')
                return candidate
        raise ArithmeticError('TT→TDB iteration failed to converge')

    def iteration_diagnostics(self):
        return asdict(self.last_iteration) if self.last_iteration else None
