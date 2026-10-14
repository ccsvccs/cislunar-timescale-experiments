# Phase 2 — Independent DE440 geocentric TT↔TDB time ephemeris reproduction

**Status:** Implemented and measured, target not yet achieved (2026-10-14). The separate Park candidate with major-body PN terms, solar J2 and the published ring's leading potential has a maximum raw residual of **11.282037595 ns** over the 732 modern epochs. See the [term mapping](phase2_model_mapping.md) and [measured report](phase2_jpl_reproduction_report.md) for additional small-body diagnostics, provenance and unresolved model inputs. **Scope:** geocenter only; retain Phase 1 Turyshev Eq. (21) exactly as the unchanged baseline. **Goal:** reproduce JPL `de440t.bsp` TT−TDB independently with a model matching Park et al. (2021), rather than improve the already stable quadrature in isolation.

Measured follow-up: the independent SB441 direct-potential diagnostic reaches
**2.818949569 ps** maximum raw residual; the separately specified Park literal
ecliptic-latitude diagnostic reaches **1.953519690 ps**. Both fail the 1 ps gate.
Neither introduces a fitted offset or linear correction. Original DE440
small-body trajectories and complete all-body PN inputs remain unverified.

## Reference sources

- [Park 2021 Sec. 2.3 Eq. (3)–(7)](https://ssd.jpl.nasa.gov/doc/Park.2021.AJ.DE440.pdf): higher-completeness JPL time ephemeris formulation; acquire and commit PDF per [Park provenance note](sources/park_2021_de440_de441.reference.md).
- [LTE440 Sec. 3 / Fig. 1](sources/lte440_2509.18511v1.pdf): independent numerical integration and DE440 geocentric validation (<1 ps reported over DE440 span).
- [Turyshev Eq. (21)](sources/turyshev_2507.22145v2.pdf), [Phase 1 model audit](model_audit.md), [Phase 1 reproducibility report](tt_tdb_eq21_report.md).
- [Reference index](references.md).

## Baseline (existing repository results; not newly rerun)

- 2024-01-01 to 2026-01-01, 732 daily samples: current truncated Eq. (21) vs true `de440t.bsp`: maximum absolute time difference **34.585 ns**; derivative residual about **2.30e-17** peak (dimensionless).
- Independent Romberg–Gauss-Legendre integration agreement ~**6.07e-18 s** maximum over those samples; this only verifies internal numerical convergence of the same physical model.
- Original epoch defined by IAU 2006 B3: JD 2443144.5003725 plus TDB0 = −65.5 μs; no reference-based offset fit or artificial detrending.
- `de440.bsp` and `de440t.bsp` are different artifacts; original `de440t.bsp` includes true TT−TDB time pseudo-target 1000000001 centered at 1000000000.
- See `results/tables/jpl_comparison.json`, `experiment_summary_fresh.json`, `model_error_budget.json`.

## Work sequence

1. Synchronize current `main` safely with GitHub CLI, inspect `git status`; never overwrite uncommitted work. Verify LTE440 PDF, obtain Park PDF through verified official source, record SHA256.
2. Write `docs/phase2_model_mapping.md`: line-by-line *read-the-actual-PDF* mapping of Park Eq. (3)–(7) to Turyshev Eq. (21) and code. Explicitly define potential conventions, geocenter vs ground, signs, time scales, TDB-compatible GM/positions, vector-potential and extra PN terms, boundary time, and expected magnitudes. No guessed coefficients.
3. Implement a **separate** `JPL-complete candidate` integrand with explicit optional contributions and per-channel diagnostics. Reuse `TwoPartJD`, ephemeris interface, integrators and independent references; never import JPL time results into the integrand. Include explicit direct gravitational potentials from modeled asteroids, Kuiper belt objects, and ring, with source/geometry/GMs verified; do not assume DE440 state perturbations also contribute to the direct potential. If full body state sources or model constants unavailable, explicitly classify as blocked and implement only verified terms.
4. Validate **rate** residual `d(TDB−TT)/dTDB` and **raw time** residual `(TDB−TT)model − (TDB−TT)de440t`, without zero-point alignment, fit removal, or trend removal in acceptance. Diagnostic detrended plots may be separate and labeled.
5. Build targeted diagnostic experiments for each omitted contribution, attribution table, error-budget breakdown, ablation tests, independent quadrature checks, short vs long epochs, and full regression of Phase 1.
6. Short-range integration and precision checks first; then 1977 fixed anchor to modern epoch; optionally extend only after valid interval and original anchor handling are clear. Store deterministic scripts, data checksums, environment, reference segment/units, tables, figures and explanation.

## Acceptance gates

- **Primary target:** over the fixed 2024-01-01 through 2026-01-01 732 daily TDB samples, **max absolute raw delta < 1e-12 s (1 ps)** versus official JPL `de440t.bsp`, preserving the original defined epoch and sign conventions; don't claim success if not actually measured.
- Rate residual, numerical quadrature error, and physical mismatch reported separately; report extended date ranges only when actually tested. No automatic calibration to reference.
- Existing Phase 1 reproducibility, baseline results, files and tests unchanged except clearly separated documented phase-2 additions.
- If 1 ps cannot be reached, deliver highest quality partial implementation plus verified residual plots and ranked missing physics, **not fabricated achievement**.

## Initial focus

Investigate the role of consistent dynamical direct potentials of asteroids/KBO/ring and all terms from Park Eq. (3)–(7), including solar asphericity and vector/pN potentials, alongside mass and coordinate-time conventions. Existing error budget's 44 AU ring estimate is a conditional scenario only, not proof of a 23 ns JPL correction. Reuse local cached data, avoid committing huge ephemeris kernels.
