"""Independent, vectorized quadrature and anchored cumulative checkpoints.

For Eq. (21), callers supply columns [F2-C_exact, F4] in dimensionless
rates. This module does not know physical constants and never calls ERFA.
The LTE440-inspired default is ten Romberg rows per half-day segment.
Last-diagonal/high-low-order differences are empirical diagnostics, not
rigorous error bounds or estimates of the physical model error.
"""
from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
import math
from pathlib import Path
import time
from typing import Callable

import numpy as np

Rate = Callable[[np.ndarray], np.ndarray]


@dataclass
class IntegrationResult:
    value: np.ndarray
    error: np.ndarray
    evaluations: int
    levels: int
    method: str


def _evaluate(rate: Rate, times: np.ndarray, channels: int) -> np.ndarray:
    values = np.asarray(rate(times.ravel()), dtype=np.float64)
    if values.ndim == 1 and channels == 1:
        values = values[:, None]
    if values.shape != (times.size, channels):
        raise ValueError(f"rate must return ({times.size}, {channels}), got {values.shape}")
    if not np.all(np.isfinite(values)):
        raise ValueError("rate returned nonfinite values")
    return values.reshape(*times.shape, channels)


def _quadrature_intervals(rate: Rate, starts: np.ndarray, ends: np.ndarray,
                          method: str, order: int, channels: int) -> IntegrationResult:
    """Integrate independent intervals simultaneously; no per-node Python loop."""
    starts = np.atleast_1d(np.asarray(starts, dtype=np.float64))
    ends = np.atleast_1d(np.asarray(ends, dtype=np.float64))
    if starts.shape != ends.shape or starts.ndim != 1:
        raise ValueError("interval endpoint arrays must have the same one-dimensional shape")
    if not np.all(np.isfinite(starts)) or not np.all(np.isfinite(ends)):
        raise ValueError("integration endpoints must be finite")
    if channels < 1:
        raise ValueError("channels must be positive")
    if method not in ("romberg", "gauss_legendre"):
        raise ValueError("method must be romberg or gauss_legendre")
    if not isinstance(order, int) or order < 2:
        raise ValueError("order must be an integer >= 2")
    if method == "romberg" and order > 18:
        raise ValueError("Romberg order exceeds the supported memory limit (18 rows)")
    widths = ends - starts
    if method == "romberg":
        endpoints = np.column_stack((starts, ends))
        endpoint_values = _evaluate(rate, endpoints, channels)
        trapezoid = widths[:, None] * endpoint_values.sum(axis=1) / 2
        previous = [trapezoid]
        evaluations = 2 * len(starts)
        error = np.zeros_like(trapezoid)
        for level in range(1, order):
            count = 2 ** (level - 1)
            fractions = (2 * np.arange(count) + 1) / 2 ** level
            nodes = starts[:, None] + widths[:, None] * fractions
            values = _evaluate(rate, nodes, channels)
            trapezoid = (previous[0] / 2
                         + widths[:, None] / 2 ** level * values.sum(axis=1))
            row = [trapezoid]
            for column in range(1, level + 1):
                # Algebraically equivalent to Richardson's expression, with
                # less subtraction of large nearly equal scaled values.
                row.append(row[-1] + (row[-1] - previous[column - 1])
                           / (4 ** column - 1))
            error = np.abs(row[-1] - previous[-1])
            previous = row
            evaluations += count * len(starts)
        return IntegrationResult(previous[-1], error, evaluations, order, method)

    def gauss(degree: int) -> np.ndarray:
        nodes, weights = np.polynomial.legendre.leggauss(degree)
        times = starts[:, None] + widths[:, None] * (nodes + 1) / 2
        values = _evaluate(rate, times, channels)
        return widths[:, None] / 2 * np.einsum("n,mnk->mk", weights, values)

    high = gauss(order)
    low_order = max(1, order // 2)
    low = gauss(low_order)
    return IntegrationResult(high, np.abs(high - low),
                             (order + low_order) * len(starts), 2, method)


def integrate_rate(rate: Rate, start: float, end: float, *, method: str = "romberg",
                   order: int = 10, channels: int = 2) -> IntegrationResult:
    """Single-interval integral with seconds as the integration variable.

    ``order=10`` means ten Romberg table rows, 512 trapezoid subintervals
    and 513 total function evaluation points. For Gauss-Legendre it means
    ten high-order nodes (plus five independent diagnostic nodes).
    """
    result = _quadrature_intervals(rate, np.array([start]), np.array([end]),
                                   method, order, channels)
    return IntegrationResult(result.value[0], result.error[0], result.evaluations,
                             result.levels, result.method)


class SegmentedIntegrator:
    """Integrate from a fixed anchor and reuse compensated grid checkpoints.

    ``fingerprint`` identifies the complete physical model: ephemeris SHA256,
    GM, constants, included bodies and PN options. The caller is responsible
    for those contents; algorithm settings and the anchor are added here.
    Supplying an empty identity is prohibited for persistent-cache safety.
    Grid construction never changes the anchor or substitutes a calibration.
    """

    def __init__(self, rate: Rate, *, anchor: float, step_seconds: float = 43200,
                 method: str = "romberg", order: int = 10, channels: int = 2,
                 fingerprint: str | dict, batch_segments: int = 32) -> None:
        if not math.isfinite(anchor) or not math.isfinite(step_seconds) or step_seconds <= 0:
            raise ValueError("anchor must be finite and step_seconds positive")
        if not fingerprint:
            raise ValueError("a nonempty model fingerprint is required")
        if method not in ("romberg", "gauss_legendre"):
            raise ValueError("unknown integration method")
        if not isinstance(order, int) or order < 2 or (method == "romberg" and order > 18):
            raise ValueError("invalid integration order")
        if channels < 1 or batch_segments < 1:
            raise ValueError("channels and batch_segments must be positive")
        self.rate = rate
        self.anchor = float(anchor)
        self.step_seconds = float(step_seconds)
        self.method = method
        self.order = order
        self.channels = channels
        self.batch_segments = batch_segments
        self.model_fingerprint = fingerprint
        self.configuration = {
            "schema": 1, "anchor": self.anchor, "step_seconds": self.step_seconds,
            "method": method, "order": order, "channels": channels,
            "model": fingerprint,
        }
        self.fingerprint = hashlib.sha256(
            json.dumps(self.configuration, sort_keys=True, separators=(",", ":"),
                       allow_nan=False).encode()).hexdigest()
        zero = np.zeros(channels)
        self._sum = {0: zero.copy()}
        self._compensation = {0: zero.copy()}
        self._errors = {0: zero.copy()}
        self._queries: dict[float, np.ndarray] = {}
        self._query_errors: dict[float, np.ndarray] = {}
        self.stats = {"segments": 0, "evaluations": 0, "batches": 0,
                      "elapsed_seconds": 0.0, "query_intervals": 0, "rate_calls": 0}

    def _grid_time(self, index: int | np.ndarray) -> float | np.ndarray:
        return self.anchor + index * self.step_seconds

    def _integrate(self, starts: np.ndarray, ends: np.ndarray) -> IntegrationResult:
        began = time.perf_counter()
        result = _quadrature_intervals(self.rate, starts, ends, self.method,
                                       self.order, self.channels)
        self.stats["elapsed_seconds"] += time.perf_counter() - began
        self.stats["evaluations"] += result.evaluations
        self.stats["batches"] += 1
        self.stats["rate_calls"] += self.order if self.method == "romberg" else 2
        return result

    def _extend(self, target: int, direction: int) -> None:
        current = max(self._sum) if direction == 1 else min(self._sum)
        while (target - current) * direction > 0:
            count = min(abs(target - current), self.batch_segments)
            indices = current + direction * np.arange(count)
            result = self._integrate(self._grid_time(indices),
                                     self._grid_time(indices + direction))
            total = self._sum[current].copy()
            compensation = self._compensation[current].copy()
            for j, increment in enumerate(result.value):
                updated = total + increment
                correction = np.where(np.abs(total) >= np.abs(increment),
                                      (total - updated) + increment,
                                      (increment - updated) + total)
                compensation += correction
                total = updated
                index = int(indices[j] + direction)
                self._sum[index] = total.copy()
                self._compensation[index] = compensation.copy()
                self._errors[index] = result.error[j].copy()
            current += direction * count
            self.stats["segments"] += count

    def build(self, start: float, end: float | None = None) -> SegmentedIntegrator:
        """Extend anchor→start, or cover both start/end, preserving existing work.

        Only complete anchor-aligned segments are constructed. A query's last
        partial interval is integrated directly, avoiding extrapolation or a
        segment beyond ephemeris coverage near its endpoint.
        """
        endpoints = [float(start)] if end is None else [float(start), float(end)]
        if not all(math.isfinite(t) for t in endpoints):
            raise ValueError("build endpoints must be finite")
        offsets = [(t - self.anchor) / self.step_seconds for t in endpoints]
        right = math.floor(max(0.0, *offsets))
        left = math.ceil(min(0.0, *offsets))
        self._extend(right, 1)
        self._extend(left, -1)
        return self

    def integral(self, times: float | np.ndarray) -> np.ndarray:
        """Return anchor-based channel integrals, shape ``times.shape + (K,)``."""
        values = np.asarray(times, dtype=np.float64)
        flat = values.ravel()
        if not np.all(np.isfinite(flat)):
            raise ValueError("query times must be finite")
        if not len(flat):
            return np.empty((*values.shape, self.channels))
        self.build(float(flat.min()), float(flat.max()))
        missing = np.array(sorted(set(float(t) for t in flat) - self._queries.keys()))
        for start in range(0, len(missing), self.batch_segments):
            batch = missing[start:start + self.batch_segments]
            # Truncation picks the checkpoint towards the anchor in both directions.
            indices = np.trunc((batch - self.anchor) / self.step_seconds).astype(np.int64)
            checkpoints = self._grid_time(indices)
            for j in np.flatnonzero(batch == checkpoints):
                index = int(indices[j])
                self._queries[float(batch[j])] = self._sum[index] + self._compensation[index]
                self._query_errors[float(batch[j])] = np.zeros(self.channels)
            pending = np.flatnonzero(batch != checkpoints)
            if len(pending):
                result = self._integrate(checkpoints[pending], batch[pending])
                self.stats["query_intervals"] += len(pending)
                for j, row, error in zip(pending, result.value, result.error):
                    index = int(indices[j])
                    self._queries[float(batch[j])] = np.array([
                        math.fsum((self._sum[index][k], self._compensation[index][k], row[k]))
                        for k in range(self.channels)])
                    self._query_errors[float(batch[j])] = error.copy()
        return np.stack([self._queries[float(t)] for t in flat]).reshape(*values.shape, self.channels)

    @property
    def checkpoint_times(self) -> np.ndarray:
        return np.asarray(self._grid_time(np.array(sorted(self._sum))))

    @property
    def segment_errors(self) -> np.ndarray:
        """Absolute last-table differences per checkpoint, sorted in time."""
        return np.stack([self._errors[i] for i in sorted(self._errors)])

    @property
    def cumulative_error_diagnostic(self) -> np.ndarray:
        """Sum of per-segment indicators, not a rigorous error bound."""
        return self.segment_errors.sum(axis=0)

    @property
    def query_times(self) -> np.ndarray:
        return np.array(sorted(self._queries))

    @property
    def query_errors(self) -> np.ndarray:
        """Error indicators of final partial intervals, aligned with query_times."""
        if not self._queries:
            return np.empty((0, self.channels))
        return np.stack([self._query_errors[t] for t in sorted(self._queries)])

    @staticmethod
    def _payload_digest(arrays: tuple[np.ndarray, ...]) -> str:
        digest = hashlib.sha256()
        for array in arrays:
            digest.update(str(array.shape).encode())
            digest.update(str(array.dtype).encode())
            digest.update(array.tobytes())
        return digest.hexdigest()

    def save(self, path: str | Path) -> None:
        """Persist compensated checkpoints and a content/hash-checked identity."""
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        indices = np.array(sorted(self._sum), dtype=np.int64)
        totals = np.stack([self._sum[i] for i in indices])
        compensations = np.stack([self._compensation[i] for i in indices])
        errors = np.stack([self._errors[i] for i in indices])
        metadata = {"configuration": self.configuration, "fingerprint": self.fingerprint,
                    "stats": self.stats,
                    "payload_sha256": self._payload_digest((indices, totals, compensations, errors))}
        # Writing via a file handle respects the supplied filename, without
        # numpy silently appending another .npz suffix.
        temporary = path.with_name(path.name + ".tmp")
        with temporary.open("wb") as stream:
            np.savez_compressed(stream, indices=indices, totals=totals,
                                compensations=compensations, errors=errors,
                                metadata=json.dumps(metadata, sort_keys=True, allow_nan=False))
        temporary.replace(path)

    @classmethod
    def load(cls, path: str | Path, rate: Rate, *, fingerprint: str | dict,
             anchor: float | None = None, step_seconds: float | None = None,
             method: str | None = None, order: int | None = None,
             channels: int | None = None, batch_segments: int = 32) -> SegmentedIntegrator:
        """Validate model/algorithm identity, then resume stored checkpoints.

        Explicit algorithm overrides are checked against the stored signature.
        Omitting an override intentionally accepts that cached configuration.
        """
        with np.load(path, allow_pickle=False) as data:
            metadata = json.loads(str(data["metadata"]))
            indices = data["indices"].copy()
            totals = data["totals"].copy()
            compensations = data["compensations"].copy()
            errors = data["errors"].copy()
        config = metadata["configuration"]
        instance = cls(rate, anchor=config["anchor"] if anchor is None else anchor,
                       step_seconds=config["step_seconds"] if step_seconds is None else step_seconds,
                       method=config["method"] if method is None else method,
                       order=config["order"] if order is None else order,
                       channels=config["channels"] if channels is None else channels,
                       fingerprint=fingerprint, batch_segments=batch_segments)
        if config.get("schema") != 1 or instance.fingerprint != metadata["fingerprint"]:
            raise ValueError("cache fingerprint mismatch: model, anchor or algorithm changed")
        if cls._payload_digest((indices, totals, compensations, errors)) != metadata["payload_sha256"]:
            raise ValueError("cache payload checksum mismatch")
        expected_shape = (len(indices), instance.channels)
        if (indices.ndim != 1 or not len(indices) or 0 not in indices
                or not np.array_equal(indices, np.arange(indices[0], indices[-1] + 1))
                or any(a.shape != expected_shape or not np.all(np.isfinite(a))
                       for a in (totals, compensations, errors)) or np.any(errors < 0)):
            raise ValueError("cache checkpoint structure is invalid")
        zero_index = int(np.flatnonzero(indices == 0)[0])
        if np.any(totals[zero_index]) or np.any(compensations[zero_index]):
            raise ValueError("cache anchor is not zero")
        instance._sum = {int(i): row for i, row in zip(indices, totals)}
        instance._compensation = {int(i): row for i, row in zip(indices, compensations)}
        instance._errors = {int(i): row for i, row in zip(indices, errors)}
        instance.stats = metadata["stats"]
        return instance
