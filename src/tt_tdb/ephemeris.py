"""Geometric DE440 states, in SI, on the TDB-compatible J2000/ICRF system.

No aberration corrections, Earth self-potential, or time-scale rescaling are
applied here. The physical model must explicitly choose its scale convention.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import os
import re

import numpy as np
from jplephem.spk import SPK

DAY_SECONDS = 86400.0
J2000_JD = 2451545.0
BODY_IDS = {"sun": 10, "moon": 301, "mercury": 1, "venus": 2,
            "mars": 4, "jupiter": 5, "saturn": 6, "uranus": 7,
            "neptune": 8, "pluto": 9}


def read_gm_kernel(path: str | Path) -> dict[int, float]:
    """Read BODY*_GM in a NAIF text PCK; convert km³/s² to m³/s².

    Only \\begindata blocks are parsed. Other satellite and asteroid GMs in
    the supplied PCK are exposed but never silently included in the model.
    """
    values = {}
    active = False
    for line in Path(path).read_text().splitlines():
        if "\\begindata" in line:
            active = True
        elif "\\begintext" in line:
            active = False
        elif active:
            match = re.search(r"BODY(\d+)_GM\s*=\s*\(\s*([+\-\d.DEded]+)\s*\)", line)
            if match:
                values[int(match[1])] = float(match[2].replace("D", "E").replace("d", "e")) * 1e9
    required = set(BODY_IDS.values()) | {3, 399}
    if not required.issubset(values):
        raise ValueError(f"GM kernel missing required body IDs {sorted(required - values.keys())}")
    return values


@dataclass(frozen=True)
class EphemerisState:
    """Earth state and external-body positions; last axis is Cartesian xyz."""
    earth_position: np.ndarray
    earth_velocity: np.ndarray
    body_positions: dict[str, np.ndarray]
    gm: dict[str, float]
    scale: str = "TDB-compatible"


class DE440Ephemeris:
    """Memory-mapped, vectorized reader of an authentic DE440 SPK.

    ``states(t)`` accepts TDB seconds since JD 2451545.0. Scalar results are
    shape (3,), and array inputs of shape S produce states of shape S+(3,).
    ``states_jd(jd1, jd2)`` retains a two-part date without first summing it.
    Missing files raise FileNotFoundError; there is no synthetic fallback.
    """

    def __init__(self, path: str | Path | None = None, gm_path: str | Path | None = None):
        self.path = Path(path or os.environ.get("DE440_BSP", "data/de440.bsp"))
        self.gm_path = Path(gm_path or os.environ.get("DE440_GM", "data/gm_de440.tpc"))
        for file in (self.path, self.gm_path):
            if not file.is_file():
                raise FileNotFoundError(f"Required DE440 data missing: {file}; run scripts/download_data.py")
        self.all_gm = read_gm_kernel(self.gm_path)
        self.gm = {name: self.all_gm[body] for name, body in BODY_IDS.items()}
        self.kernel = SPK.open(str(self.path))
        required = {(0, n) for n in range(1, 11)} | {(3, 301), (3, 399)}
        actual = {(s.center, s.target) for s in self.kernel.segments}
        if not required.issubset(actual):
            self.kernel.close()
            raise ValueError(f"SPK missing DE440 segments {sorted(required - actual)}")
        self.start_jd = max(self.kernel[pair].start_jd for pair in required)
        self.end_jd = min(self.kernel[pair].end_jd for pair in required)
        self.frame = "J2000/ICRF, geometric"
        self.missing_potential_sources = ("343 individual asteroids", "30 KBOs", "Kuiper-belt ring", "resolved planetary satellite multipoles")

    def states(self, seconds_since_j2000) -> EphemerisState:
        seconds = np.asarray(seconds_since_j2000, dtype=float)
        days = np.floor(seconds / DAY_SECONDS)
        # Avoid rounding decades of seconds divided by one day to a large
        # fractional JD (that can move a state by millimetres).
        remainder = seconds - days * DAY_SECONDS
        return self.states_jd(J2000_JD + days, remainder / DAY_SECONDS)

    def states_jd(self, jd1, jd2=0.0) -> EphemerisState:
        jd1, jd2 = np.broadcast_arrays(np.asarray(jd1, dtype=float), np.asarray(jd2, dtype=float))
        shape = jd1.shape
        first, second = jd1.reshape(-1), jd2.reshape(-1)
        elapsed = (first - self.start_jd) + second
        duration = self.end_jd - self.start_jd
        if np.any(~np.isfinite(elapsed)) or np.any(elapsed < 0) or np.any(elapsed > duration):
            raise ValueError(f"TDB epoch outside DE440 coverage JD {self.start_jd}..{self.end_jd}")

        def position(center, target):
            value = self.kernel[center, target].compute(first, second)
            return (value.T * 1000.0).reshape(shape + (3,))

        emb, emb_v = self.kernel[0, 3].compute_and_differentiate(first, second)
        offset, offset_v = self.kernel[3, 399].compute_and_differentiate(first, second)
        earth = ((emb + offset).T * 1000.0).reshape(shape + (3,))
        velocity = ((emb_v + offset_v).T * (1000.0 / DAY_SECONDS)).reshape(shape + (3,))
        positions = {name: position(0, body) for name, body in BODY_IDS.items() if name != "moon"}
        # Moon and Earth are separate; EMB is never added as an external mass.
        positions["moon"] = (emb.T * 1000.0).reshape(shape + (3,)) + position(3, 301)
        return EphemerisState(earth, velocity, positions, self.gm)

    def close(self):
        self.kernel.close()

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()
