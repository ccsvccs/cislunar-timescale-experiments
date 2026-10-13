"""Real model, real 1977 anchor, no calibration or long cached prerequisites."""
import numpy as np
import pytest
from tt_tdb.ephemeris import DE440Ephemeris
from tt_tdb.pn_rate import PNModel
from tt_tdb.integration import SegmentedIntegrator
from tt_tdb.constants import T0_TDB,TDB0,A_SCALE
from tt_tdb.transforms import TTTDBConverter


@pytest.mark.de440
def test_real_anchor_bidirectional_integration_and_roundtrips(de440_paths):
    with DE440Ephemeris(*de440_paths) as eph:
        model=PNModel(eph)
        integral=SegmentedIntegrator(model,anchor=T0_TDB.seconds_since_j2000(),
                                     fingerprint='nonpersistent-real-DE440-test',order=8)
        converter=TTTDBConverter(integral)
        assert converter.tdb_minus_tt(T0_TDB)==TDB0
        for seconds in (-86400,-.001,0,.001,86400):
            epoch=T0_TDB.add_seconds(seconds)
            tt=converter.tdb_to_tt(epoch)
            recovered=converter.tt_to_tdb(tt)
            assert abs(recovered.difference_seconds(epoch))<1e-15
            assert abs(converter.tdb_to_tt(converter.tt_to_tdb(epoch)).difference_seconds(epoch))<1e-15
        # An infinitesimal forward derivative checks the sign of the residual
        # at the physically fixed anchor, not just an annual amplitude.
        minus=converter.tdb_minus_tt(T0_TDB.add_seconds(-10))
        plus=converter.tdb_minus_tt(T0_TDB.add_seconds(10))
        rate=A_SCALE*model.diagnostics(T0_TDB.seconds_since_j2000())['residual_rate']
        assert abs((plus-minus)/20-rate)<1e-20
