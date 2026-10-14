# Phase 2 small-body input provenance

Acquired on 2026-10-14 from the actual JPL directories. The recorded sources are snapshots, not proof that all inputs use the identical DE440 dynamical solution.

| Artifact | Official source | Local verification |
|---|---|---|
| Park 2021 | https://ssd.jpl.nasa.gov/doc/Park.2021.AJ.DE440.pdf | 15 pages, 5,501,569 bytes; SHA256 `96c3945deb12e0f0843c2bf0703d3542d950cf85ff903b1fe11a501634f68ee0` |
| Farnocchia, IOM392R-21-005, March 30 2021 | https://ssd.jpl.nasa.gov/ftp/eph/small_bodies/asteroids_de441/SB441_IOM392R-21-005_perturbers.pdf | `docs/sources/sb441_perturbers.pdf`, 12 pages, 993,300 bytes; SHA256 `124fe0fb48f8e71f3f1dba9607cfe5b47a016b95aa5429b61bf7cba301ed4b22`; author/Caltech copyright preserved |
| Actual short 373-body kernel | https://ssd.jpl.nasa.gov/ftp/eph/small_bodies/asteroids_de441/sb441-n373s.bsp | `data/sb441-n373s.bsp`, 982,106,112 bytes; SHA256 `2143113282bfc2b2a0b0b4626125d4f84362339b5a8ae7eea40f4120ca8da10b`; genuine DAF/SPK; not committed |
| JPL directory evidence | https://ssd.jpl.nasa.gov/ftp/eph/small_bodies/ and https://ssd.jpl.nasa.gov/ftp/eph/small_bodies/asteroids_de441/ | `jpl_small_bodies_index.html`, `jpl_asteroids_de441_index.html`; DE430/441 available, no DE440 subdirectory advertised |

`sb441-n373s.bsp` has 374 SPK Type 2 segments for 373 distinct targets. All centers are Sun (NAIF 10), all frames are J2000; common **per-target union** coverage is TDB JD 2287184.5–2688976.5, the DE440 interval. Target 2001036 has two adjacent segments, JD 2287184.5–2399474.5 and 2399474.5–2688976.5. A reader must select by epoch; a dictionary that keeps only the last segment gives incorrect coverage. Other targets have one segment. Add the same DE440 Sun barycentric position before constructing small-body geocentric distances; do not interpret heliocentric SPK positions as barycentric.

The actual embedded comments identify the short file as a 2021-09-02 SPKMERGE extraction from `sb441-n373.bsp`. The original full kernel is advertised as 14.13 GB, the short kernel as 936.6 MB (MiB display). The short file fully includes the original 1977 anchor and modern test interval. Download used a `.part` file; only after successful HTTP completion, binary read/segment verification and coverage checks was it renamed.

The IOM explicitly says the **same set** of 343 asteroids + 30 KBOs is used for DE440/441. However its trajectories were re-integrated relative to DE441 reference planetary ephemerides and fitted iteratively to optical/Gaia/radar/Dawn data through March 2021. That does not establish identity with the original DE440 time-integrator small-body trajectories. Its Table 1 GM values use the DE441 solution and can be checked against each DE440 header value. Asteroid `MA####` labels correspond to asteroid numbers; KBO `MA8001..8030` labels must be mapped using Table 1, not assigned fabricated NAIF IDs. The independent implementation writes `sb441_mapping.json` with that traceable correspondence.

Consequently SB441 + DE440 is an explicit **diagnostic combination**, not certified complete DE440 input. Rate and accumulated leading-potential effects from real distances are useful evidence; a remaining picosecond difference cannot automatically be assigned to that version difference. SB441 orbit convergence differences are not formal uncertainty bounds. The ring's 36 bodies are absent from this kernel; its published geometry remains a separate leading-scalar candidate.

The archived source text is for audit, while the kernel is a large ignored input. Do not treat the memo's original full-kernel -8001..9000-year coverage as the measured short-file coverage. Public NAIF 2010 300-asteroid kernels and current Horizons orbit solutions were examined as alternatives and rejected as matching DE440 inputs.
