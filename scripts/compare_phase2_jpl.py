"""Independent Park Eq. (3)–(7) candidate: raw time, rate and term attribution.

The reference is only read after integration. No fitted offset, rate or reference
checkpoint enters the model. Missing small-body data are reported by the model.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import time

import _common
import numpy as np

from tt_tdb.constants import A_SCALE, TDB0, T0_TDB
from tt_tdb.ephemeris import DE440Ephemeris
from tt_tdb.extended_ephemeris import ExtendedDE440Ephemeris
from tt_tdb.integration import SegmentedIntegrator
from tt_tdb.jpl_model import JPLModel
from tt_tdb.jpl_reference import JPLTimeEphemeris
from tt_tdb.pn_rate import PNModel
from tt_tdb.reference import difference_metrics
from tt_tdb.timebase import TwoPartJD


def metrics(times, difference, unit="seconds"):
    """Acceptance uses only raw max/RMS; any fitted fields are diagnostics."""
    return difference_metrics(times, difference, value_unit=unit)


def integrate(rate, times, signature, method, order, channels, no_cache):
    path = Path(f"results/cache/phase2_{method}{order}_halfday.npz")
    parameters = dict(fingerprint=signature,
                      anchor=T0_TDB.seconds_since_j2000(), step_seconds=43200.,
                      method=method, order=order, channels=channels, batch_segments=32)
    if path.exists() and not no_cache:
        try:
            integrator = SegmentedIntegrator.load(path, rate, **parameters)
            print(f"Validated phase2 cache: {path}", flush=True)
        except ValueError as error:
            print(f"Rebuilding phase2 cache: {error}", flush=True)
            integrator = SegmentedIntegrator(rate, **parameters)
    else:
        integrator = SegmentedIntegrator(rate, **parameters)
    started = time.perf_counter()
    values = integrator.integral(times)
    integrator.save(path)
    info = {"method": method, "order": order, "step_days": .5,
            "elapsed_wall_seconds": time.perf_counter() - started,
            "stats": integrator.stats.copy(),
            "sum_segment_error_indicators_seconds": integrator.cumulative_error_diagnostic.tolist(),
            "max_segment_error_indicators_seconds": integrator.segment_errors.max(axis=0).tolist(),
            "cache_fingerprint": integrator.fingerprint}
    print(f"{method} completed: {info['elapsed_wall_seconds']:.2f} s", flush=True)
    return values, info


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--start-jd", default="2460310.5")
    parser.add_argument("--end-jd", default="2461041.5")
    parser.add_argument("--no-cache", action="store_true")
    parser.add_argument("--reference", default=None)
    args = parser.parse_args()
    started = time.perf_counter()
    start = TwoPartJD.from_decimal_jd(args.start_jd).seconds_since_j2000()
    end = TwoPartJD.from_decimal_jd(args.end_jd).seconds_since_j2000()
    if end <= start:
        parser.error("end must follow start")
    times = np.arange(start, end + 1., 86400.)
    if times[-1] != end:
        times = np.append(times, end)
    anchor = T0_TDB.seconds_since_j2000()
    with ExtendedDE440Ephemeris() as ephemeris, DE440Ephemeris() as baseline_ephemeris:
        baseline = PNModel(baseline_ephemeris)
        model = JPLModel(ephemeris)
        term_names = tuple(model.correction_rates(np.array([start])).keys())

        def channels(t):
            terms = model.correction_rates(t)
            return np.column_stack((baseline(t), *(terms[name] for name in term_names)))

        signature = {"model": json.loads(model.fingerprint()),
                     "baseline": json.loads(baseline.fingerprint()),
                     "channels": ["baseline_F2_minus_C", "baseline_F4", *term_names],
                     "integration_source_sha256": hashlib.sha256(Path("src/tt_tdb/integration.py").read_bytes()).hexdigest(),
                     "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                     "versions": _common.versions()}
        # Short modern-window convergence before integrating from 1977.
        short_times = start + np.arange(33) * 86400.
        short = {}
        reference_short = None
        for method, order, step in (("gauss_legendre", 24, 21600.),
                                    ("gauss_legendre", 16, 43200.),
                                    ("romberg", 10, 43200.),
                                    ("romberg", 10, 86400.)):
            check = SegmentedIntegrator(model, anchor=start, method=method, order=order,
                                        step_seconds=step, fingerprint=signature)
            values = check.integral(short_times)
            if reference_short is None:
                reference_short = values
            short[f"{method}{order}_{step/86400:g}day"] = {
                "delta_difference": metrics(short_times, A_SCALE * (values - reference_short).sum(axis=-1)),
                "stats": check.stats}
        probe = start + np.arange(33) * 1234.567
        closure = channels(probe).sum(axis=-1) - model(probe).sum(axis=-1)
        print(f"Term reconstruction max rate discrepancy: {abs(closure).max():.3e}", flush=True)
        query = np.concatenate(([anchor, 0.], times))
        channel_count = 2 + len(term_names)
        raw, romberg_info = integrate(channels, query, signature, "romberg", 10,
                                      channel_count, args.no_cache)
        independent, gauss_info = integrate(channels, query, signature, "gauss_legendre", 24,
                                            channel_count, args.no_cache)
        phase1_delta = TDB0 + A_SCALE * raw[2:, :2].sum(axis=-1)
        contributions = A_SCALE * raw[2:, 2:]
        delta = phase1_delta + contributions.sum(axis=-1)
        delta_gauss = TDB0 + A_SCALE * independent[2:].sum(axis=-1)
        instantaneous = channels(times)
        phase1_rate = A_SCALE * instantaneous[:, :2].sum(axis=-1)
        added_rate = A_SCALE * instantaneous[:, 2:]
        model_rate = phase1_rate + added_rate.sum(axis=-1)
        with JPLTimeEphemeris(args.reference) as reference:
            jpl_delta, jpl_rate = reference.evaluate(times)
            anchor_jpl, _ = reference.evaluate(anchor)
            coverage = [max(ephemeris.start_jd, reference.start_jd),
                        min(ephemeris.end_jd, reference.end_jd)]
            # This is ONLY a sparse derivative diagnostic; no 1100-year
            # integrated residual or longer-interval acceptance is claimed.
            wide_jds = np.linspace(coverage[0]+.5, coverage[1]-.5, 221)
            wide_times = (wide_jds - 2451545.) * 86400.
            _, wide_ref_rate = reference.evaluate(wide_times)
            wide_model_rate = A_SCALE * model(wide_times).sum(axis=-1)
            wide_baseline_rate = A_SCALE * baseline(wide_times).sum(axis=-1)
            ref_sha = hashlib.sha256(reference.path.read_bytes()).hexdigest()
        residual, rate_residual = delta-jpl_delta, model_rate-jpl_rate
        attribution = {}
        for column, name in enumerate(term_names):
            term, term_rate = contributions[:, column], added_rate[:, column]
            without = residual-term
            attribution[name] = {
                "max_abs_instantaneous_delta_rate": float(np.max(np.abs(term_rate))),
                "anchor_to_first_epoch_seconds": float(term[0]),
                "anchor_to_last_epoch_seconds": float(term[-1]),
                "modern_peak_to_peak_seconds": float(np.ptp(term)),
                "net_modern_accumulation_seconds": float(term[-1]-term[0]),
                "raw_residual_without_term": metrics(times, without),
                "max_abs_residual_change_when_added_seconds": float(abs(without).max()-abs(residual).max()),
                "rate_residual_without_term": metrics(times, rate_residual-term_rate, "dimensionless")}
        document = {
            "command": " ".join(__import__("sys").argv), "versions": _common.versions(),
            "model_signature": signature, "sample_count": len(times),
            "start_jd": args.start_jd, "end_jd": args.end_jd,
            "anchor_jd_parts": T0_TDB.jd_parts(), "anchor_j2000_seconds": anchor,
            "anchor_delta_seconds": TDB0,
            "anchor_model_minus_jpl_seconds": float(TDB0-anchor_jpl),
            "sign": "TDB-TT; SI seconds; d/dTDB dimensionless",
            "reference": {"path": str(reference.path), "sha256": ref_sha,
                          "segment": [1000000000, 1000000001], "common_coverage_jd": coverage,
                          "derivative": "analytic Chebyshev derivative / 86400, sign reversed; independently checked against SPICE"},
            "acceptance": {"threshold_seconds": 1e-12,
                           "max_abs_raw_seconds": float(np.max(np.abs(residual))),
                           "passed": bool(np.max(np.abs(residual)) < 1e-12),
                           "zero_point_fitted": False, "linear_rate_fitted": False},
            "baseline_time_residual": metrics(times, phase1_delta-jpl_delta),
            "baseline_rate_residual": metrics(times, phase1_rate-jpl_rate, "dimensionless"),
            "time_residual": metrics(times, residual),
            "rate_residual": metrics(times, rate_residual, "dimensionless"),
            "residual_first_seconds": float(residual[0]), "residual_last_seconds": float(residual[-1]),
            "residual_peak_to_peak_seconds": float(np.ptp(residual)),
            "independent_quadratures": metrics(times, delta-delta_gauss),
            "max_channel_quadrature_difference_seconds": abs(raw-independent).max(axis=0).tolist(),
            "term_reconstruction_max_abs_rate": float(abs(closure).max()),
            "short_convergence": short, "romberg": romberg_info, "gauss": gauss_info,
            "attribution": attribution,
            "sparse_full_coverage_rate_only": {"count": len(wide_times),
                "start_jd": float(wide_jds[0]), "end_jd": float(wide_jds[-1]),
                "extended_rate_residual": metrics(wide_times, wide_model_rate-wide_ref_rate, "dimensionless"),
                "baseline_rate_residual": metrics(wide_times, wide_baseline_rate-wide_ref_rate, "dimensionless"),
                "integrated_long_interval_test_performed": False},
            "elapsed_wall_seconds": time.perf_counter()-started}
    _common.write_json("results/tables/phase2_jpl_comparison.json", document)
    headers = ["tdb_j2000_seconds", "phase1_delta_s", "extended_delta_s", "gauss_delta_s",
               "jpl_delta_s", "phase1_rate", "extended_rate", "jpl_rate"]
    values = [times, phase1_delta, delta, delta_gauss, jpl_delta,
              phase1_rate, model_rate, jpl_rate]
    for index, name in enumerate(term_names):
        headers.extend((name+"_delta_s", name+"_rate"))
        values.extend((contributions[:, index], added_rate[:, index]))
    _common.csv_table("results/tables/phase2_time_series.csv", headers, np.column_stack(values))
    _common.csv_table("results/tables/phase2_full_coverage_rates.csv",
                      ["tdb_j2000_seconds", "extended_rate_residual", "baseline_rate_residual"],
                      np.column_stack((wide_times, wide_model_rate-wide_ref_rate, wide_baseline_rate-wide_ref_rate)))
    figures = _common.make_plots_directory()
    import matplotlib.pyplot as plt
    days = (times-times[0])/86400.
    fig, axes = plt.subplots(2, 1, figsize=(10, 7), sharex=True)
    axes[0].plot(days, (phase1_delta-jpl_delta)*1e9, label="Phase 1")
    axes[0].plot(days, residual*1e9, label="Park candidate (available bodies)")
    axes[0].set_ylabel("Raw TDB-TT residual (ns)")
    axes[1].plot(days, (phase1_rate-jpl_rate)*1e17, label="Phase 1")
    axes[1].plot(days, rate_residual*1e17, label="Park candidate")
    axes[1].set_ylabel("Raw rate residual (10⁻¹⁷)")
    axes[1].set_xlabel("Days since first TDB epoch")
    for ax in axes:
        ax.grid(alpha=.3); ax.legend()
    fig.tight_layout(); fig.savefig(figures/"phase2_raw_residuals.png", dpi=170); plt.close(fig)
    fig, axes = plt.subplots(2, 1, figsize=(11, 8), sharex=True)
    for index, name in enumerate(term_names):
        axes[0].plot(days, contributions[:, index]*1e12, label=name)
        axes[1].plot(days, added_rate[:, index]*1e19, label=name)
    axes[0].set_ylabel("Term integral from 1977 (ps)")
    axes[1].set_ylabel("Term rate (10⁻¹⁹)")
    axes[1].set_xlabel("Days since first TDB epoch")
    for ax in axes:
        ax.grid(alpha=.3); ax.legend(fontsize=8, ncol=2)
    fig.tight_layout(); fig.savefig(figures/"phase2_term_contributions.png", dpi=170); plt.close(fig)
    fig, ax = plt.subplots(figsize=(10, 4))
    ax.plot(wide_jds, (wide_model_rate-wide_ref_rate)*1e17, label="Park candidate")
    ax.plot(wide_jds, (wide_baseline_rate-wide_ref_rate)*1e17, label="Phase 1")
    ax.set_xlabel("TDB Julian day; sparse full-coverage rate samples")
    ax.set_ylabel("Raw rate residual (10⁻¹⁷)")
    ax.set_title("Derivative diagnostic only; no long-interval time acceptance")
    ax.grid(alpha=.3); ax.legend(); fig.tight_layout()
    fig.savefig(figures/"phase2_full_coverage_rates.png", dpi=170); plt.close(fig)
    print(json.dumps({key: document[key] for key in
                      ("acceptance", "time_residual", "rate_residual", "independent_quadratures", "elapsed_wall_seconds")}, indent=2), flush=True)


if __name__ == "__main__":
    main()
