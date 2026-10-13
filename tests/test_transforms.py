import numpy as np
import pytest
from tt_tdb.constants import T0_TDB,TDB0,A_SCALE
from tt_tdb.timebase import TwoPartJD
from tt_tdb.transforms import TTTDBConverter


class KnownIntegral:
    anchor=T0_TDB.seconds_since_j2000()
    channels=2
    def integral(self,t):
        # Analytic synthetic residual with the correct anchor. High enough
        # derivative to ensure a single plug-in TT argument is insufficient.
        x=np.asarray(t)-self.anchor
        return np.stack((.0015*np.sin(x/3e6)/A_SCALE,np.zeros_like(x)),axis=-1)


def test_anchor_sign_and_implicit_roundtrips():
    converter=TTTDBConverter(KnownIntegral())
    assert converter.tdb_minus_tt(T0_TDB)==TDB0
    for date in ('2443144.5003725','2451545.0','2460310.500123456789'):
        t=TwoPartJD.from_decimal_jd(date)
        tt=converter.tdb_to_tt(t)
        recovered=converter.tt_to_tdb(tt)
        assert abs(recovered.difference_seconds(t))<1e-15
        tdb=converter.tt_to_tdb(t)
        assert abs(converter.tdb_to_tt(tdb).difference_seconds(t))<1e-15
        assert abs(converter.last_iteration.implicit_residual_seconds)<1e-15
        assert converter.last_iteration.iterations<=4
    with pytest.raises(TypeError):
        converter.tdb_minus_tt(2451545.)


def test_nonconverging_iteration_limit():
    c=TTTDBConverter(KnownIntegral())
    with pytest.raises(ArithmeticError):
        c.tt_to_tdb(TwoPartJD.from_decimal_jd('2460310.5'),max_iterations=1,tolerance_seconds=1e-20)
