#!/usr/bin/env python3
"""Real DE440 convergence, fixed-input precision, and explicit Kepler checks.

Run from the repository root with .venv/bin/python scripts/convergence_study.py.
Missing real data raises an error; synthetic results are labelled separately.
"""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
os.environ.setdefault("MPLCONFIGDIR", "/tmp/tt_tdb_matplotlib")

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter
import mpmath as mp
import numpy as np
from scipy.integrate import quad

from tt_tdb.constants import A_SCALE, C_EXACT, C_LIGHT
from tt_tdb.ephemeris import DE440Ephemeris
from tt_tdb.integration import SegmentedIntegrator, integrate_rate
from tt_tdb.pn_rate import PNModel
from tt_tdb.timebase import time_precision_experiment
from tt_tdb.validation import KeplerReference


def metric(values):
    values = np.asarray(values, dtype=float)
    return {"max_abs_seconds": float(np.max(np.abs(values))),
            "rms_seconds": float(np.sqrt(np.mean(values ** 2))),
            "p95_abs_seconds": float(np.percentile(np.abs(values), 95)),
            "mean_seconds": float(np.mean(values))}


def comparison(values, baseline):
    difference = values - baseline
    return {"F2_minus_C_integral": metric(difference[:, 0]),
            "F4_integral": metric(difference[:, 1]),
            "sum_integral": metric(difference.sum(axis=1)),
            "delta_equivalent": metric(A_SCALE * difference.sum(axis=1))}


def run_case(model, times, fingerprint, step, method, order):
    cache = SegmentedIntegrator(model, anchor=float(times[0]), step_seconds=step,
                                 method=method, order=order, fingerprint=fingerprint)
    started = time.perf_counter()
    result = cache.integral(times)
    info = {"method": method, "order": order, "step_days": step / 86400,
            "runtime_seconds": time.perf_counter() - started,
            "stats": cache.stats,
            "sum_segment_error_indicators_seconds": cache.cumulative_error_diagnostic.tolist(),
            "max_segment_error_indicator_seconds": cache.segment_errors.max(axis=0).tolist(),
            "max_final_partial_error_indicator_seconds": cache.query_errors.max(axis=0).tolist()}
    return result, info


def fixed_input_precision(model, start):
    """Separate arithmetic checks from the accuracy of double ephemeris inputs."""
    c = mp.mpf(int(C_LIGHT))
    state = model.ephemeris.states(start)
    diagnostics = model.diagnostics(start)
    earth_position = [mp.mpf(float(x)) for x in state.earth_position]
    velocity = [mp.mpf(float(x)) for x in state.earth_velocity]
    v2_exact = mp.fsum(x * x for x in velocity)
    potential_exact = mp.fsum(
        mp.mpf(float(state.gm[body])) / mp.sqrt(mp.fsum(
            (mp.mpf(float(position[k])) - earth_position[k]) ** 2 for k in range(3)))
        for body, position in state.body_positions.items())
    v2_double = mp.mpf(float(diagnostics["velocity_squared"]))
    u_double = mp.mpf(float(diagnostics["potential"]))

    def rates(v2, u):
        return ((v2 / 2 + u) / c ** 2,
                (v2 ** 2 / 8 + mp.mpf("1.5") * v2 * u - u ** 2 / 2) / c ** 4)

    formula_f2, formula_f4 = rates(v2_double, u_double)
    state_f2, state_f4 = rates(v2_exact, potential_exact)
    actual_f2, actual_f4 = float(diagnostics["F2"]), float(diagnostics["F4"])

    # Fixed binary64 nodes, weights, times and real DE440 samples. The 60-digit
    # calculation isolates weighted-sum arithmetic, not ephemeris evaluation.
    degree = 24
    width = 43200.0
    nodes, weights = np.polynomial.legendre.leggauss(degree)
    sampled_times = start + width * (nodes + 1) / 2
    samples = model(sampled_times)
    double_integral = width / 2 * np.einsum("n,nk->k", weights, samples)
    high_integral = [mp.mpf(width) / 2 * mp.fsum(
        mp.mpf(float(w)) * mp.mpf(float(y)) for w, y in zip(weights, samples[:, k]))
        for k in range(2)]
    high_nodes, high_weights = mp.gauss_quadrature(degree, "legendre")
    exact_weight_sum = [mp.mpf(width) / 2 * mp.fsum(
        high_weights[j] * mp.mpf(float(samples[j, k])) for j in range(degree))
        for k in range(2)]
    return {
        "kind": "real_DE440_states_and_samples_fixed_at_binary64_inputs",
        "mpmath_digits": mp.mp.dps,
        "epoch_seconds_since_j2000": start,
        "limitations": [
            "DE440 Chebyshev evaluation, coordinates, velocities and GMs remain binary64 inputs.",
            "High precision does not recover information lost when computing those inputs.",
            "High precision GL weights reuse the binary64-node samples; their nodes are not reevaluated.",
            "These tests measure arithmetic discrepancies, not physical or absolute ephemeris error."],
        "potential_binary64_minus_fixed_state_mp_SI": float(u_double - potential_exact),
        "v2_binary64_minus_fixed_state_mp_SI": float(v2_double - v2_exact),
        "F2": actual_f2, "F4": actual_f4,
        "formula_F2_binary64_minus_mp": float(mp.mpf(actual_f2) - formula_f2),
        "formula_F4_binary64_minus_mp": float(mp.mpf(actual_f4) - formula_f4),
        "fixed_state_F2_binary64_minus_mp": float(mp.mpf(actual_f2) - state_f2),
        "fixed_state_F4_binary64_minus_mp": float(mp.mpf(actual_f4) - state_f4),
        "residual_binary64_minus_fixed_state_mp": float(
            mp.mpf(float(actual_f2 - C_EXACT)) - (state_f2 - mp.mpf(C_EXACT))),
        "gauss_fixed_inputs": {
            "degree": degree, "interval_seconds": [start, start + width],
            "double_integral_seconds": double_integral.tolist(),
            "mp_fixed_weights_integral_seconds": [str(x) for x in high_integral],
            "weighted_sum_binary64_minus_mp_seconds": [float(mp.mpf(float(x)) - y)
                                                        for x, y in zip(double_integral, high_integral)],
            "double_weights_minus_mp_weights_seconds": [float(x - y)
                                                         for x, y in zip(high_integral, exact_weight_sum)],
            "max_double_node_minus_mp_node": float(max(abs(mp.mpf(float(nodes[i])) - high_nodes[i])
                                                       for i in range(degree)))}}


def independent_quad(model, start):
    end = start + 43200
    quad_values, quad_errors, details = [], [], []
    for channel in range(2):
        result = quad(lambda t: float(model(np.array([t]))[0, channel]), start, end,
                      epsabs=1e-21, epsrel=1e-12, limit=200, full_output=1)
        value, error, info = result[:3]
        quad_values.append(value)
        quad_errors.append(error)
        details.append({"channel": channel, "evaluations": info["neval"],
                        "subintervals": info["last"],
                        "message": result[3] if len(result) > 3 else "converged"})
    romberg = integrate_rate(model, start, end, order=10)
    gauss = integrate_rate(model, start, end, method="gauss_legendre", order=24)
    return {"kind": "real_DE440_independent_scipy_QUADPACK_short_interval",
            "interval_seconds": [start, end],
            "quad_integral_seconds": quad_values,
            "quad_reported_errors_seconds": quad_errors, "quad_details": details,
            "romberg_minus_quad_seconds": (romberg.value - quad_values).tolist(),
            "gauss24_minus_quad_seconds": (gauss.value - quad_values).tolist(),
            "romberg_evaluations": romberg.evaluations,
            "gauss_evaluations": gauss.evaluations}


def kepler_experiment(figures):
    model = KeplerReference()
    times = np.linspace(0, model.period, 129)
    expected = model.integral(0, times)
    results = {}
    computed = {}
    for method, order in (("romberg", 10), ("gauss_legendre", 24)):
        cache = SegmentedIntegrator(model.rate, anchor=0, channels=1,
                                     method=method, order=order,
                                     fingerprint="explicit-synthetic-Newtonian-Kepler")
        started = time.perf_counter()
        values = cache.integral(times)[:, 0]
        results[method] = {"difference_from_exact": metric(values - expected),
                           "runtime_seconds": time.perf_counter() - started,
                           "stats": cache.stats}
        computed[method] = values
    fig, axes = plt.subplots(2, 1, figsize=(8, 6), sharex=True)
    axes[0].plot(times / 86400, expected * 1e3, label="Exact Kepler integral", color="black")
    axes[0].plot(times / 86400, computed["romberg"] * 1e3, linestyle="--", label="Romberg")
    axes[0].set_ylabel("Residual integral (ms)")
    axes[0].legend()
    for method, values in computed.items():
        axes[1].plot(times / 86400, (values - expected) * 1e15, label=method)
    axes[1].set_xlabel("Days since perihelion (synthetic)")
    axes[1].set_ylabel("Numerical difference (fs)")
    axes[1].legend()
    for axis in axes:
        axis.grid(alpha=0.3)
    fig.suptitle("Synthetic Kepler validation; not a DE440 accuracy test")
    fig.tight_layout()
    fig.savefig(figures / "kepler_validation.png", dpi=160)
    plt.close(fig)
    return {"kind": "explicit_synthetic_exact_Newtonian_ellipse_c2_only",
            "period_days": model.period / 86400,
            "analytic_amplitude_seconds": model.amplitude, "eccentricity": model.e,
            "derivation": "docs/kepler_reference.md", "methods": results}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bsp", default="data/de440.bsp")
    parser.add_argument("--gm", default="data/gm_de440.tpc")
    parser.add_argument("--start-seconds", type=float, default=0.0)
    parser.add_argument("--days", type=int, default=32)
    parser.add_argument("--output", default="results")
    args = parser.parse_args()
    if args.days < 1:
        parser.error("--days must be positive")
    output = Path(args.output)
    tables, figures = output / "tables", output / "figures"
    tables.mkdir(parents=True, exist_ok=True)
    figures.mkdir(parents=True, exist_ok=True)
    times = args.start_seconds + np.arange(args.days + 1) * 86400.0
    with DE440Ephemeris(args.bsp, args.gm) as ephemeris:
        model = PNModel(ephemeris)
        fingerprint = model.fingerprint()
        reference, reference_info = run_case(model, times, fingerprint, 21600,
                                             "gauss_legendre", 24)
        cases, all_values = [], {}
        for step in (86400, 43200, 21600):
            for order in (6, 8, 10):
                values, info = run_case(model, times, fingerprint, step, "romberg", order)
                info["difference_from_gauss24_quarterday"] = comparison(values, reference)
                cases.append(info)
                all_values[(step, order)] = values
                print(f"Romberg {order}, step {step / 86400:g} day: "
                      f"max total difference {info['difference_from_gauss24_quarterday']['sum_integral']['max_abs_seconds']:.6g} s",
                      flush=True)
        independent, independent_info = run_case(model, times, fingerprint, 43200,
                                                  "gauss_legendre", 16)
        independent_info["difference_from_gauss24_quarterday"] = comparison(independent, reference)
        document = {
            "kind": "real_DE440_Turyshev_eq21_finite_monopole_geocenter_short_window",
            "anchor_seconds_since_j2000": args.start_seconds,
            "duration_days": args.days, "sample_count": len(times),
            "channels": ["integral(F2-C_exact) dt", "integral F4 dt"],
            "time_units": "seconds", "rate_input_scale": "TDB-compatible",
            "romberg_order_definition": "number of table rows; order 10 uses 513 points per segment",
            "reference": reference_info, "reference_channels_seconds": reference.tolist(),
            "model_fingerprint": json.loads(fingerprint), "cases": cases,
            "independent_gauss16": independent_info,
            "step_comparisons_at_order10": {
                "one_day_minus_half_day": comparison(all_values[(86400, 10)], all_values[(43200, 10)]),
                "half_day_minus_quarter_day": comparison(all_values[(43200, 10)], all_values[(21600, 10)])},
            "limitations": "Convergence checks numerical quadrature, not absolute physical-model accuracy. "
                           "No zero-point adjustment or detrending was applied."}
        (tables / "convergence.json").write_text(json.dumps(document, indent=2, allow_nan=False) + "\n")
        with mp.workdps(60):
            precision = fixed_input_precision(model, args.start_seconds)
        precision["independent_scipy_quad"] = independent_quad(model, args.start_seconds)
    precision["time_representation"] = time_precision_experiment()
    precision["kepler_validation"] = kepler_experiment(figures)
    (tables / "precision_validation.json").write_text(json.dumps(precision, indent=2, allow_nan=False) + "\n")
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    for order in (6, 8, 10):
        selected = [row for row in cases if row["order"] == order]
        x = [row["step_days"] for row in selected]
        for axis, channel in zip(axes, ("sum_integral", "F4_integral")):
            y = [row["difference_from_gauss24_quarterday"][channel]["max_abs_seconds"]
                 for row in selected]
            axis.loglog(x, np.maximum(y, 1e-30), marker="o", label=f"Romberg {order} rows")
            axis.set_xlabel("Segment width (day)")
            axis.set_xticks([0.25, 0.5, 1.0], ["0.25", "0.5", "1"])
            axis.xaxis.set_minor_formatter(NullFormatter())
            axis.set_ylabel("Maximum absolute difference (s)")
            axis.set_title("Total residual integral" if channel == "sum_integral" else "c^-4 integral")
            axis.grid(alpha=0.3)
            axis.legend()
    fig.suptitle(f"Real DE440: {args.days} days, vs independent Gauss24 / 0.25 day")
    fig.tight_layout()
    fig.savefig(figures / "convergence.png", dpi=160)
    plt.close(fig)
    print(f"Saved {tables / 'convergence.json'} and {tables / 'precision_validation.json'}", flush=True)


if __name__ == "__main__":
    main()
