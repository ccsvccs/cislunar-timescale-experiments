"""SI and IAU defining constants, kept separate from an empirical paper value.

Sources and decimal derivation: docs/model_audit.md, IAU 2000 B1.9,
IAU 2006 B3, Turyshev arXiv:2507.22145v2 Eqs. (15), (19)-(21).
"""
from decimal import Decimal, localcontext
from .timebase import TwoPartJD

C_LIGHT = 299792458.0
L_G = 6.969290134e-10
L_B = 1.550519768e-8
L_C_PAPER = 1.48082685455e-8
TDB0 = -65.5e-6
T0_JD_DECIMAL = '2443144.5003725'
T0 = TwoPartJD.from_decimal_jd(T0_JD_DECIMAL)
T0_TDB = T0.add_seconds(TDB0)


def decimal_constants():
    """Recompute defining linear quantities from exact decimal literals."""
    with localcontext() as context:
        context.prec = 60
        lg, lb, lc = map(Decimal, ('6.969290134e-10','1.550519768e-8','1.48082685455e-8'))
        a = (1-lg)/(1-lb)
        c = (lb-lg)/(1-lg)
        return {'A': a, 'C_exact': c, 'L_C_paper': lc, 'C_minus_L_C': c-lc,
                'omitted_drift_seconds_per_day': a*(c-lc)*86400,
                'omitted_drift_seconds_per_julian_year': a*(c-lc)*Decimal('31557600')}


_DECIMAL = decimal_constants()
A_SCALE = float(_DECIMAL['A'])
C_EXACT = float(_DECIMAL['C_exact'])
K_LINEAR = float(_DECIMAL['A']*_DECIMAL['C_exact'])
