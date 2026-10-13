"""Independent reference comparisons; never imported by the physical model."""
from __future__ import annotations
import numpy as np


def sofa_geocentric(tdb1, tdb2):
    """ERFA dtdb: TDB-TT in seconds, geocenter u=v=0, no station terms.

    Source: https://github.com/liberfa/erfa/blob/master/src/dtdb.c
    This is a comparison model, not a DE440 numerical time ephemeris.
    """
    import erfa
    return erfa.dtdb(tdb1, tdb2, 0.0, 0.0, 0.0, 0.0)


def difference_metrics(seconds, difference, *, value_unit='seconds'):
    """Raw differences; time abscissae are seconds, slopes are per day.

    Rates use ``value_unit='dimensionless'``. Fits are diagnostic only.
    """
    x, d = np.asarray(seconds, dtype=float), np.asarray(difference, dtype=float)
    if x.shape != d.shape or x.size < 2 or not np.all(np.isfinite(d)):
        raise ValueError("Need at least two matched finite observations")
    centered_days = (x - np.mean(x)) / 86400.0
    centered_difference = d - np.mean(d)
    slope = float(np.dot(centered_days, centered_difference) / np.dot(centered_days, centered_days))
    fit = np.mean(d) + slope * centered_days
    result = {
        "count": int(d.size), "max_abs_seconds": float(np.max(np.abs(d))),
        "rms_seconds": float(np.sqrt(np.mean(d*d))), "mean_seconds": float(np.mean(d)),
        "min_seconds": float(np.min(d)), "max_seconds": float(np.max(d)),
        "p50_abs_seconds": float(np.quantile(np.abs(d), .5)),
        "p95_abs_seconds": float(np.quantile(np.abs(d), .95)),
        "linear_slope_seconds_per_day": slope,
        "rms_after_diagnostic_linear_fit_seconds": float(np.sqrt(np.mean((d-fit)**2))),
    }
    return {key.replace('_seconds', '_'+value_unit): value for key, value in result.items()}
