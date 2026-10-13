"""Time arithmetic preserving small corrections separately from large epochs.

An integral Julian day starts at noon. ``day`` is the integer part of JD;
``seconds`` and ``fraction`` represent the elapsed SI seconds within that day.
The third component is necessary: a float containing 43200 seconds already has
a 7.3 ps spacing. Ephemeris evaluation may use a float offset, but conversion
corrections and round-trip differences must use the split representation.
"""
from __future__ import annotations

from dataclasses import dataclass
from decimal import Decimal, localcontext
import math

DAY_SECONDS = 86400
J2000_JD = 2451545


@dataclass(frozen=True)
class TwoPartJD:
    """Split Julian date (integer day, integer seconds, fractional seconds).

    The name describes the external two-double JD interface; internally the
    extra integer seconds ensure sub-picosecond corrections survive at noon.
    Pass decimal-defined epochs through ``from_decimal_jd``, not a float JD.
    """

    day: int
    seconds: int | float = 0
    fraction: float = 0.0

    def __post_init__(self) -> None:
        if isinstance(self.day, bool) or int(self.day) != self.day:
            raise ValueError("day must be an integer Julian day")
        if not math.isfinite(self.seconds) or not math.isfinite(self.fraction):
            raise ValueError("time components must be finite")
        integer_seconds = math.floor(self.seconds)
        fractional = math.fsum((self.seconds - integer_seconds, self.fraction))
        carry = math.floor(fractional)
        fractional -= carry
        day_carry, seconds = divmod(integer_seconds + carry, DAY_SECONDS)
        object.__setattr__(self, "day", int(self.day) + day_carry)
        object.__setattr__(self, "seconds", seconds)
        object.__setattr__(self, "fraction", fractional)

    @classmethod
    def from_decimal_jd(cls, jd: str | Decimal) -> TwoPartJD:
        """Construct from a decimal-defined JD without first rounding to float."""
        with localcontext() as ctx:
            ctx.prec = 60
            value = Decimal(jd)
            if not value.is_finite():
                raise ValueError("Julian date must be finite")
            day = int(value // 1)
            elapsed = (value - day) * DAY_SECONDS
            seconds = int(elapsed // 1)
            return cls(day, seconds, float(elapsed - seconds))

    @classmethod
    def from_jd(cls, jd1: float, jd2: float = 0.0) -> TwoPartJD:
        """Construct from two binary doubles, preserving their exact sum."""
        if not math.isfinite(jd1) or not math.isfinite(jd2):
            raise ValueError("Julian date parts must be finite")
        with localcontext() as ctx:
            ctx.prec = 60
            return cls.from_decimal_jd(Decimal.from_float(float(jd1))
                                       + Decimal.from_float(float(jd2)))

    @classmethod
    def from_seconds_since_j2000(cls, seconds: float) -> TwoPartJD:
        return cls(J2000_JD, seconds)

    def add_seconds(self, delta_seconds: float) -> TwoPartJD:
        """Apply a correction without rounding it into an absolute epoch."""
        if not math.isfinite(delta_seconds):
            raise ValueError("time correction must be finite")
        # For small negative corrections, floor would introduce 1 + delta;
        # truncation keeps the small quantity intact instead.
        integer = math.trunc(delta_seconds)
        remainder = delta_seconds - integer
        return TwoPartJD(self.day, self.seconds + integer,
                         math.fsum((self.fraction, remainder)))

    def difference_seconds(self, other: TwoPartJD) -> float:
        """Return self - other in seconds using compensated small differences."""
        whole = (self.day - other.day) * DAY_SECONDS + self.seconds - other.seconds
        return math.fsum((float(whole), self.fraction, -other.fraction))

    def seconds_since_j2000(self) -> float:
        """Float64 ephemeris argument, not a representation of tiny corrections."""
        return self.difference_seconds(TwoPartJD(J2000_JD))

    def decimal_jd(self) -> Decimal:
        with localcontext() as ctx:
            ctx.prec = 60
            return (Decimal(self.day)
                    + (Decimal(self.seconds) + Decimal.from_float(self.fraction))
                    / Decimal(DAY_SECONDS))

    def jd_parts(self) -> tuple[float, float]:
        """Return a coarse complete JD and its accurately computed tiny residual.

        Do not add the returned doubles before passing them to a two-part API.
        A conventional integer-day + fractional-day split can itself lose
        several picoseconds; the small residual used here avoids that loss.
        """
        with localcontext() as ctx:
            ctx.prec = 60
            exact = self.decimal_jd()
            coarse = float(exact)
            residual = float(exact - Decimal.from_float(coarse))
            return coarse, residual

    @property
    def tdb1(self) -> float:
        return self.jd_parts()[0]

    @property
    def tdb2(self) -> float:
        return self.jd_parts()[1]


def time_precision_experiment() -> list[dict[str, float]]:
    """Measure actual correction retention; no ephemeris or NumPy required."""
    epoch = TwoPartJD.from_decimal_jd("2460676.123456789")
    coarse = float(epoch.decimal_jd())
    elapsed = float(epoch.seconds + epoch.fraction)
    offset = epoch.seconds_since_j2000()
    rows = []
    for correction in (1e-6, 1e-9, 1e-12, 1e-13):
        updated = epoch.add_seconds(correction)
        jd1, jd2 = updated.jd_parts()
        roundtrip = TwoPartJD.from_jd(jd1, jd2)
        rows.append({
            "requested_seconds": correction,
            "single_jd_recovered_seconds": (coarse + correction / DAY_SECONDS - coarse) * DAY_SECONDS,
            "day_float_seconds_recovered_seconds": (elapsed + correction) - elapsed,
            "j2000_float_seconds_recovered_seconds": (offset + correction) - offset,
            "split_recovered_seconds": updated.difference_seconds(epoch),
            "split_error_seconds": updated.difference_seconds(epoch) - correction,
            "two_double_jd_roundtrip_error_seconds": roundtrip.difference_seconds(updated),
        })
    return rows
