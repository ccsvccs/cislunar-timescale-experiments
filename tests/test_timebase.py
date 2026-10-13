from decimal import Decimal
import math

import pytest

from tt_tdb.timebase import TwoPartJD, time_precision_experiment


def test_decimal_defined_reference_epoch():
    epoch = TwoPartJD.from_decimal_jd("2443144.5003725")
    assert epoch.day == 2443144
    assert epoch.seconds == 43232
    assert epoch.fraction == pytest.approx(0.184, abs=1e-17)
    assert abs(epoch.decimal_jd() - Decimal("2443144.5003725")) < Decimal("1e-21")


@pytest.mark.parametrize("correction", [1e-6, 1e-9, 1e-12, 1e-13, -1e-13])
def test_subpicosecond_correction(correction):
    original = TwoPartJD(2460000, 43232, 0.184)
    corrected = original.add_seconds(correction)
    assert abs(corrected.difference_seconds(original) - correction) < 3e-17
    assert abs(corrected.add_seconds(-correction).difference_seconds(original)) < 3e-17


def test_carry_and_negative_correction():
    epoch = TwoPartJD(2451545, 86399, 0.75)
    assert epoch.add_seconds(0.5) == TwoPartJD(2451546, 0, 0.25)
    assert TwoPartJD(2451545).add_seconds(-0.25) == TwoPartJD(2451544, 86399, 0.75)
    assert TwoPartJD(2451545, -86400.25) == TwoPartJD(2451543, 86399, 0.75)


def test_j2000_epoch_and_long_difference():
    epoch = TwoPartJD.from_seconds_since_j2000(-1234567.125)
    assert epoch.seconds_since_j2000() == -1234567.125
    assert TwoPartJD(2451546).difference_seconds(TwoPartJD(2451545)) == 86400


def test_two_double_jd_retains_small_correction():
    epoch = TwoPartJD.from_decimal_jd("2460676.123456789").add_seconds(1e-13)
    reconstructed = TwoPartJD.from_jd(*epoch.jd_parts())
    assert abs(reconstructed.difference_seconds(epoch)) < 3e-17


def test_single_jd_actually_loses_nanoseconds():
    rows = time_precision_experiment()
    for row in rows:
        assert row["single_jd_recovered_seconds"] == 0
        assert abs(row["split_error_seconds"]) < 6e-17
        assert abs(row["two_double_jd_roundtrip_error_seconds"]) < 3e-17


@pytest.mark.parametrize("args", [(1.5, 0), (2451545, math.nan), (2451545, 0, math.inf)])
def test_invalid_time(args):
    with pytest.raises(ValueError):
        TwoPartJD(*args)
