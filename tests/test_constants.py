from decimal import Decimal, localcontext
from tt_tdb.constants import A_SCALE,C_EXACT,L_C_PAPER,T0,T0_TDB,TDB0,decimal_constants


def test_exact_decimal_definitions_and_geocentric_identity():
    d=decimal_constants()
    with localcontext() as ctx:
        ctx.prec=60
        assert abs(d['A']*(1-d['C_exact'])-1) < Decimal('1e-58')
        assert abs(d['omitted_drift_seconds_per_day']-Decimal('11.3547151318e-12')) < Decimal('1e-22')
    assert C_EXACT != L_C_PAPER
    assert abs(T0_TDB.difference_seconds(T0)-TDB0) < 3e-17
    assert T0.seconds==43232
    assert abs(T0.fraction-.184)<1e-16
    assert A_SCALE>1
