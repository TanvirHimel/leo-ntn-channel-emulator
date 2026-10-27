# Measured results

Every performance or accuracy number cited in the README must appear here first, together with
how it was measured. Each entry records the commit, date and command. Tolerances follow SPEC.md §8;
if a tolerance is ever loosened, this file records why.

## Machine
Captured with the M1a measurements (`out/machine_m1a.txt`, 2026-10-21T02:58:51Z).

| Item | Value |
|---|---|
| CPU | AMD Ryzen 5 5500U with Radeon Graphics, 6 cores / 12 threads, max clock 2100 MHz |
| RAM | 5.9 GiB |
| OS | Microsoft Windows 11 Home 10.0.26300 build 26300 |
| Power plan | Balanced |
| Libraries | Python 3.11.0; numpy 2.4.6; scipy 1.17.1; skyfield 1.55; sgp4 2.27 (C++ accelerated: True); matplotlib 3.11.2 |

MSVC: not used in M1a (Python only).

## M1a — Python golden slice
Code commit: `6feb7fc` (clean tree; full hash `6feb7fcca2048a07387e81f59966195d5d438909`).
Header of the measurement run: `M1a golden acceptance measurements  (commit 6feb7fc, 2026-10-21 04:47 W. Europe Daylight Time)`.

Commands (repo root, PowerShell):
- `.\.venv\Scripts\python.exe scripts\measure_m1a.py --slow` -> `out/measure_m1a_slow.txt` (313.2 s)
- `.\.venv\Scripts\python.exe scripts\demo_m1a.py` -> `out/demo_m1a.txt`
- `.\.venv\Scripts\python.exe -m pytest` and `... -m slow` -> `out/pytest_fast_m1a.txt`, `out/pytest_slow_m1a.txt`

Test suites at this commit: fast `74 passed, 7 deselected in 62.62s (0:01:02)`; slow `7 passed, 74 deselected in 234.07s (0:03:54)`.

### Acceptance tests 1-4: golden geometric channel and fractional delay
| Test | Measured | Limit |
|---|---|---|
| 1 phase, rho0=1e+06 m, v=-7000 m/s, N=1e+06 | 6.42e-09 rad (max \|\|y\|-1\| 2.2e-16) | 1e-7 rad |
| 1 FFT peak, rho0=1e+06 m, v=-7000 m/s, N=1e+06 | 1.52e-04 Hz | 0.1 Hz |
| 1 phase, rho0=1.23457e+06 m, v=6853.37 m/s, N=1e+06 | 1.30e-08 rad (max \|\|y\|-1\| 2.2e-16) | 1e-7 rad |
| 1 FFT peak, rho0=1.23457e+06 m, v=6853.37 m/s, N=1e+06 | 7.64e-07 Hz | 0.1 Hz |
| 1 phase, rho0=1e+06 m, v=-7000 m/s, N=1e+07 | 6.42e-09 rad (max \|\|y\|-1\| 2.2e-16) | 1e-7 rad |
| 1 FFT peak, rho0=1e+06 m, v=-7000 m/s, N=1e+07 | 1.37e-05 Hz | 0.1 Hz |
| 1 phase, rho0=1.23457e+06 m, v=6853.37 m/s, N=1e+07 | 1.30e-08 rad (max \|\|y\|-1\| 2.2e-16) | 1e-7 rad |
| 1 FFT peak, rho0=1.23457e+06 m, v=6853.37 m/s, N=1e+07 | 2.13e-04 Hz | 0.1 Hz |
| 2 integer shift, k=0 (snap path) | exact (0) | 1e-12 rel |
| 2 integer shift, k=1 (snap path) | exact (0) | 1e-12 rel |
| 2 integer shift, k=37 (snap path) | exact (0) | 1e-12 rel |
| 2 integer shift, k=70001 (snap path) | exact (0) | 1e-12 rel |
| 3 fractional delay, d=0.1 | -156.8 dB | -80 dB |
| 3 fractional delay, d=0.25 | -150.7 dB | -80 dB |
| 3 fractional delay, d=0.5 | -148.4 dB | -80 dB |
| 3 fractional delay, d=0.73 | -150.3 dB | -80 dB |
| 3 fractional delay, d=0.95 | -162.4 dB | -80 dB |
| 3 fractional delay, d=5.3 | -149.7 dB | -80 dB |
| 3 fractional delay, d=-2.4 | -148.7 dB | -80 dB |
| 3 fractional delay, d=0.999 | -195.9 dB | -80 dB |
| 3b integer delay, snap off (kernel path), k=0 | exact (0) | -120 dB |
| 3b integer delay, snap off (kernel path), k=3 | exact (0) | -120 dB |
| 3b integer delay, snap off (kernel path), k=70001 | exact (0) | -120 dB |
| 3b integer delay, snap off (kernel path), k=0 + 1e-10 | -195.1 dB | -120 dB |
| 3b integer delay, snap off (kernel path), k=3 + 1e-10 | -195.1 dB | -120 dB |
| 3b integer delay, snap off (kernel path), k=70001 + 1e-10 | -194.9 dB | -120 dB |
| 4 quadratic phase, inst. frequency | 1.46e-04 Hz | 0.01 Hz |

Notes:
- Test 2 runs the exact-integer ("snap") path; 3b runs the same integer delays through the
  windowed-sinc kernel with snapping disabled, plus a 1e-10-sample offset that makes all 128
  taps nonzero.
- Test 3 delays are >= 0.05 samples from an integer, except the separate near-integer case 0.999.

### Acceptance test 7: SigMF
Write -> read round trip is byte-identical (including +-0, subnormals, +-inf, NaN);
`core:datatype` is `cf32_le`; `core:frequency` is in `captures[0]` and absent from `global`.
Pass/fail check with no numeric metric (`python/tests/test_sigmf_io.py`, part of the fast suite above).

### Pinned scenario (inputs, not expected values)
Main TLE `tle/starlink_53934.txt` (STARLINK-4775, epoch 2026-10-20T04:06:17.748840Z); UE
52.5 N, 13.4 E, 50.0 m (rounded); mask 10 deg.
Selection and fetch details: `tests/reference/SOURCES.md`.

| Pass | Max elevation (deg) | AOS (UTC) | TCA (UTC) | LOS (UTC) |
|---|---|---|---|---|
| low | 34.763254 | 2026-10-20T01:09:22.282911Z | 2026-10-20T01:13:05.546313Z | 2026-10-20T01:16:50.083969Z |
| mid | 46.424768 | 2026-10-20T06:06:56.014650Z | 2026-10-20T06:10:52.194184Z | 2026-10-20T06:14:47.283690Z |
| high | 84.544908 | 2026-10-20T04:27:31.609378Z | 2026-10-20T04:31:38.206034Z | 2026-10-20T04:35:44.482198Z |

### Acceptance tests 5-6: segments and light-time range (per pass)
| Test | Measured | Limit |
|---|---|---|
| 5 segment fit, high pass, TCA +-60 s: rho | 1.14e-04 m | 1e-3 m |
| 5 segment fit, high pass, TCA +-60 s: rho_dot | 3.51e-04 m/s | 1e-3 m/s |
| 5 segment fit, low pass, full pass: rho | 2.74e-05 m | 1e-3 m |
| 5 segment fit, low pass, full pass: rho_dot | 8.46e-05 m/s | 1e-3 m/s |
| 5 segment fit, mid pass, full pass: rho | 4.96e-05 m | 1e-3 m |
| 5 segment fit, mid pass, full pass: rho_dot | 1.53e-04 m/s | 1e-3 m/s |
| 5 segment fit, high pass, full pass: rho | 1.14e-04 m | 1e-3 m |
| 5 segment fit, high pass, full pass: rho_dot | 3.52e-04 m/s | 1e-3 m/s |
| 6 light-time first order, low pass | 9.73e-04 m | 1e-2 m |
| 6 light-time first order, mid pass | 9.71e-04 m | 1e-2 m |
| 6 light-time first order, high pass | 1.02e-03 m | 1e-2 m |

Test 5 in the fast suite uses only the high pass, TCA +-60 s (SPEC.md M1a); the slow suite runs
all three full passes. Both are listed above.

### Range-rate derivation (rho_dot)
rho_dot is the 4th-order central difference (h = 0.2 s) of the position-based range track,
never SGP4 velocity (CLAUDE.md).

| Test | Measured | Limit |
|---|---|---|
| rho_dot routine, circular orbit (closed form) | 1.29e-07 m/s | 1e-6 m/s |
| rho_dot step consistency h=0.1 vs 0.2 s, high pass, geometric | max 2.27e-05 m/s, median 5.07e-07 m/s | max 1e-04 m/s |
| rho_dot step consistency h=0.1 vs 0.2 s, high pass, lighttime | max 2.33e-05 m/s, median 5.07e-07 m/s | max 1e-04 m/s |

The step-consistency limit is budget-derived: test 5 allows 1 mm/s, i.e. 6.7e-3 Hz at 2 GHz, and
this check may use 10% of it (1e-4 m/s). It is a new check, not a loosened SPEC tolerance. A
1 um/s target was tried first and is not reachable on the real track because of the two
effects in the next table; the derivative routine itself is held to 1e-6 m/s in closed form.

### SGP4 track properties (informational)
| Quantity | Measured | Limit |
|---|---|---|
| SGP4 velocity vs d(position)/dt, high pass TCA +-200 s (informational) | 1.08e-02 m/s | - |
| SGP4 solver step, mid pass, AOS+18.96 s (informational) | +5.95 um range; carrier phase 2.5e-04 rad; rho_dot artefact <= 2.0e-05 m/s = 1.3e-04 Hz | - |
| SGP4 solver step, high pass, AOS+23.09 s (informational) | +6.29 um range; carrier phase 2.6e-04 rad; rho_dot artefact <= 2.1e-05 m/s = 1.4e-04 Hz | - |
| SGP4 solver step, high pass, AOS+259.81 s (informational) | -1.12 um range; carrier phase 4.7e-05 rad; rho_dot artefact <= 3.7e-06 m/s = 2.5e-05 Hz | - |
| Time-quantization range noise, high pass AOS+60 s (informational) | rms 7.9e-08 m, max 2.3e-07 m | - |

- SGP4 velocity disagrees with the derivative of SGP4 position by up to the value above; this is
  why rho_dot is position-derived.
- Solver steps: the SGP4 Kepler solver stops at a fixed 1e-12 tolerance, so range jumps by a few
  um when its iteration count changes. None were found on the low pass above the 1.5e-6 m
  detection threshold.

### End-to-end demo (`scripts/demo_m1a.py`)
- `pass 'high': max el 84.54 deg, AOS 2026-10-20T04:27:31.609378Z, TCA 2026-10-20T04:31:38.206034Z`
- `window: 2026-10-20T04:29:34.407706Z + 1.0 s (centred at (AOS+TCA)/2)`
- `segments: 1; rho(0) = 1024.452 km, tau0 = 26244 samples (3.4172 ms)`
- `golden channel: 7680000 samples in 105.8 s (73 ksps)`

Input: NR-numerology CP-OFDM (QPSK, 5 MHz) plus a test tone at +2.8 MHz (-10 dB). Each 0.1 s
sub-window's tone frequency is compared with the prediction at that sub-window's centre,
-(fc + f_tone)*rho_dot/c (the tone's RF frequency; see the README design note). The last column
is the gap to the carrier-only value -fc*rho_dot/c, shown for reference.

That gap is about 55 Hz (+54.9865 to +55.1002 Hz over the ten sub-windows): the extra
Doppler of the tone's baseband frequency, f_tone*rho_dot/c.

| Sub | Centre (s) | rho_dot (m/s) | Measured tone offset (Hz) | Predicted -(fc+f_tone)rho_dot/c (Hz) | Error (Hz) | Measured minus -fc*rho_dot/c (Hz) |
|---|---|---|---|---|---|---|
| 0 | 0.050 | -5899.5026 | 39412.3450 | 39412.3449 | +0.0001 | +55.1002 |
| 1 | 0.150 | -5898.1597 | 39403.3743 | 39403.3738 | +0.0004 | +55.0880 |
| 2 | 0.250 | -5896.8145 | 39394.3868 | 39394.3866 | +0.0002 | +55.0752 |
| 3 | 0.350 | -5895.4668 | 39385.3846 | 39385.3834 | +0.0012 | +55.0637 |
| 4 | 0.450 | -5894.1167 | 39376.3626 | 39376.3641 | -0.0015 | +55.0484 |
| 5 | 0.550 | -5892.7642 | 39367.3283 | 39367.3287 | -0.0003 | +55.0369 |
| 6 | 0.650 | -5891.4094 | 39358.2766 | 39358.2772 | -0.0006 | +55.0240 |
| 7 | 0.750 | -5890.0521 | 39349.2095 | 39349.2096 | -0.0001 | +55.0118 |
| 8 | 0.850 | -5888.6924 | 39340.1250 | 39340.1260 | -0.0009 | +54.9982 |
| 9 | 0.950 | -5887.3303 | 39331.0263 | 39331.0262 | +0.0000 | +54.9865 |

Result: `max |measured - predicted| = 0.0015 Hz over 10 sub-windows (limit 1.0 Hz)` -> PASS. Figure: `docs/figures/m1a_doppler.png`.

## M1b — C++ scalar engine
Header of the measurement run: `M1b measurements (base commit a3ef802, uncommitted M1b working tree, 2026-10-23 12:13 W. Europe Daylight Time)`.

Commands (repo root, PowerShell; build commands inside the VsDevCmd environment, see CLAUDE.md):
- `cmake --build build -j 2`, then `ctest --test-dir build --output-on-failure` -> `out/ctest_fast_m1b.txt`;
  `set "NTN_SLOW=1" && ctest --test-dir build -L slow --output-on-failure` -> `out/ctest_slow_m1b.txt`
- `set "CMAKE_BUILD_PARALLEL_LEVEL=2" && .\.venv\Scripts\python.exe -m pip install -e . --no-build-isolation`
- `.\.venv\Scripts\python.exe -m pytest python/tests` and `... -m slow` -> `out/pytest_fast_m1b.txt`, `out/pytest_slow_m1b.txt`
- `.\.venv\Scripts\python.exe scripts\demo_m1b.py` -> `out/demo_m1b.txt`, `out/m1b_results.json`
- `.\.venv\Scripts\python.exe scripts\measure_m1b.py --slow` -> `out/measure_m1b_slow.txt` (everything below)

Notes:
- Scalar kernels (`frac_delay`, `rotator`, `cpu_check`) are a separate CMake target (`ntn_scalar`)
  compiled without `/arch:AVX2`. No AVX2 code exists yet (M2).
- Fractional-delay table: provisional (32 taps x 256 phases, Kaiser beta 9.7, linear
  interpolation between phases); M2 finalises it.
- Latency: `L = max(0, 15 - floor(min D - 1e-6))`, bounded over the whole range track, so it
  differs between windows cut from different tracks. Minimality is tested (impulse lines below).
- Acceptance test 3 stimulus: random QPSK OFDM band-limited to |f| < 0.4 fs (zero-phase Kaiser
  FIR). The unshaped generator leaks power into 0.45-0.5 fs, outside the C++ design band
  (|f| < 0.4 fs, SPEC §2) and near Nyquist, where the golden kernel is not a reference either;
  the diagnostic table shows the band split. The -55 dB limit is unchanged (not loosened).
- The demo keeps the M1a input unchanged (unshaped OFDM + tone), so it records the full-band
  error (not asserted) and asserts the error below 0.38 fs against -55 dB.
- gtest FFT-peak errors are below the peak search's 1e-5 Hz resolution.
- Wall times are informational, not benchmarks (M2 measures throughput).
- Figure: `docs/figures/m1b_golden_vs_cpp.png` (output spectra overlaid with their difference;
  error power per 10 ms, full band and below 0.38 fs).

### Toolchain
| Item | Value |
|---|---|
| MSVC (cl) | 19.50.35724.0 (from build/CMakeFiles/*/CMakeCXXCompiler.cmake) |
| CMake | cmake version 4.3.0-rc1 |
| Ninja | 1.13.2 |
| Python | 3.11.0 |
| pybind11 | 3.1.0 |
| scikit_build_core | 1.1.1 |
| numpy | 2.4.6 |
| scipy | 1.17.1 |
| skyfield | 1.55 |
| ntn_emu._core | _core.cp311-win_amd64.pyd; cpu_has_avx2() = True |

### Test suites (from out/ logs)
| Suite | Result |
|---|---|
| pytest fast | `104 passed, 10 deselected in 82.61s (0:01:22)` |
| pytest slow | `10 passed, 104 deselected in 387.88s (0:06:27)` |
| ctest fast | `100% tests passed, 0 tests failed out of 55`, `Total Test time (real) =   8.72 sec` |
| ctest slow (NTN_SLOW=1, -L slow) | `100% tests passed, 0 tests failed out of 1`, `Total Test time (real) =   2.22 sec` |

### C++ closed-form tests (gtest, scalar path; `[ measure ]` lines, NTN_SLOW=1)
```
frac delay d = 0.1       error -103.8 dB
frac delay d = 0.25      error -102.9 dB
frac delay d = 0.5       error -100.8 dB
frac delay d = 0.73      error -102.2 dB
frac delay d = 0.95      error -107.7 dB
frac delay d = 5.3       error -101.3 dB
frac delay d = -2.4      error -99.9 dB
frac delay d = 0.999     error -107.5 dB
frac delay d = 0.00195312 error -105.2 dB
frac delay d = 0.0013    error -106.2 dB
frac delay worst -99.9 dB (limit -90 dB)
tone rho0=1000000 v=-7000 N=1000000: max phase err 1.637e-07 rad, last 1e5 1.627e-07 rad, max ||y|-1| 4.29e-08
tone rho0=1234567.89 v=6853.37 N=1000000: max phase err 1.651e-07 rad, last 1e5 1.602e-07 rad, max ||y|-1| 4.30e-08
tone rho0=1000000 v=-7000 N=10000000: max phase err 1.648e-07 rad, last 1e5 1.619e-07 rad, max ||y|-1| 4.31e-08
tone rho0=1234567.89 v=6853.37 N=10000000: max phase err 1.698e-07 rad, last 1e5 1.632e-07 rad, max ||y|-1| 4.31e-08
tone v=-7000: peak 46698.973328 Hz, f_d 46698.973328 Hz, error 0.00e+00 Hz
tone v=6853.37: peak -45720.763262 Hz, f_d -45720.763262 Hz, error 7.28e-12 Hz
integer shift k=0: latency L = 16
integer shift k=1: latency L = 15
integer shift k=37: latency L = 0
integer shift k=70001: latency L = 0
quadratic phase: max |f_inst - (-fc rho_dot / c)| = 1.725e-03 Hz
approaching track (v=-7000 m/s, 1 s): min D = -176.3240, latency L = 192
D=2.3, L=13, block 1: |out[n+L]| with impulse at n+L = 3.403e-05, at n+L+1 = 0.000e+00
D=2.3, L=13, block 4096: |out[n+L]| with impulse at n+L = 3.403e-05, at n+L+1 = 0.000e+00
approaching, D(N-1)=-22.9495, L=38: |last output| with impulse at last input = 9.784e-06
time-varying delay: min D = -22.9795, latency L = 38
time-varying delay vs closed form: -103.6 dB
block invariance track: D 0.1302 -> -7.7403 (8 integer slips), knot at 200000, L = 23
cpu_has_avx2() = 1, IsProcessorFeaturePresent(AVX2) = 1
[  PASSED  ] 55 tests.
```

### Acceptance test 1: M1a tests 1, 2, 4 through `_core` (pytest, test_core.py)
| Test | Measured | Limit |
|---|---|---|
| 1 phase, rho0=1e+06 m, v=-7000 m/s, N=1e+06 | 1.61e-07 rad (last 1e5: 1.60e-07; max \|\|y\|-1\| 4.3e-08) | 2e-6 rad |
| 1 phase, rho0=1.23457e+06 m, v=6853.37 m/s, N=1e+06 | 1.64e-07 rad (last 1e5: 1.62e-07; max \|\|y\|-1\| 4.3e-08) | 2e-6 rad |
| 1 phase, rho0=1e+06 m, v=-7000 m/s, N=1e+07 | 1.63e-07 rad (last 1e5: 1.59e-07; max \|\|y\|-1\| 4.3e-08) | 2e-6 rad |
| 1 phase, rho0=1.23457e+06 m, v=6853.37 m/s, N=1e+07 | 1.68e-07 rad (last 1e5: 1.65e-07; max \|\|y\|-1\| 4.3e-08) | 2e-6 rad |
| 1 FFT peak, rho0=1e+06 m, v=-7000 m/s, N=1e+06 | 1.52e-04 Hz | 0.1 Hz |
| 1 FFT peak, rho0=1.23457e+06 m, v=6853.37 m/s, N=1e+06 | 7.64e-07 Hz | 0.1 Hz |
| 2 integer shift, k=0 | L = 16; 0 differing samples (bit-exact) | bit-exact |
| 2 integer shift, k=1 | L = 15; 0 differing samples (bit-exact) | bit-exact |
| 2 integer shift, k=37 | L = 0; 0 differing samples (bit-exact) | bit-exact |
| 2 integer shift, k=70001 | L = 0; 0 differing samples (bit-exact) | bit-exact |
| 4 quadratic phase, inst. frequency | 1.72e-03 Hz | 0.01 Hz |

### Acceptance test 3: C++ vs golden, real pass (pytest, test_core_vs_golden.py)
Stimulus: band-limited QPSK OFDM; power share at |f| >= 0.4 fs: -139.7 dB (0.1 s, Welch Kaiser(14)).

| Window (crosses the knot) | Error power | L | tau0 | Limit |
|---|---|---|---|---|
| 0.1 s | -105.30 dB | 174 | 26176 | -55 dB |
| 1 s | -106.04 dB | 242 | 26244 | -55 dB |
| negative control: 0.01 s, comparison off by 1 sample | -0.21 dB | | | must be > -20 dB |

Diagnostic, unshaped OFDM (0.1 s, as generated): full-band error -50.75 dB.
| Band | Input power share | Error power, rel. golden output |
|---|---|---|
| 0.00-0.30 fs | -0.0 dB | -106.1 dB |
| 0.30-0.40 fs | -30.4 dB | -128.8 dB |
| 0.40-0.45 fs | -37.9 dB | -69.5 dB |
| 0.45-0.50 fs | -38.6 dB | -50.7 dB |

### End-to-end demo (out/m1b_results.json, out/demo_m1b.txt)
| Quantity | Value |
|---|---|
| Error power, full band (recorded, not asserted) | -51.24 dB |
| Error power, \|f\| < 0.38 fs (zero-phase Kaiser low-pass, 401 taps, beta 12, cutoff 0.39 fs), 2000 edge samples excluded | -104.20 dB (limit -55 dB) |
| Per 10 ms, full band (min / median / max) | -53.2 / -51.3 / -49.9 dB |
| Per 10 ms, in band (min / median / max) | -104.7 / -104.2 / -103.7 dB |
| Tone Doppler on the C++ output, max \|measured - predicted\| | 0.0015 Hz (limit 1 Hz) |
| Latency L / tau0 | 166 / 26244 samples |
| Wall time, C++ scalar / golden (informational, not a benchmark) | 2.14 s / 125.7 s |

_measure_m1b.py ran in 194.0 s._

## M2 — Accuracy, AVX2, benchmarks
Code commit: `c215b25` + uncommitted changes (full hash `c215b2543897e2e4453f64fbfe44fcae52968d36`), measurements started 2026-10-27T22:38:53+01:00; benchmarks at commit `c215b25` + uncommitted changes.

Commands (repo root, PowerShell; build and C++ commands inside the VsDevCmd environment, see CLAUDE.md):
- `cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release -DNTN_BUILD_BENCH=ON`, `cmake --build build -j 2`
- `.\.venv\Scripts\python.exe scripts\bench_m2.py` -> `out/bench_m2.txt`, `out/bench_m2.json`, `out/bench_m2_summary.json`
- `.\.venv\Scripts\python.exe scripts\measure_m2.py` -> suites, gtest measure lines, ISA/FMA checks, design sweep, demo (`out/*_m2*`; 406 s)
- `.\.venv\Scripts\python.exe scripts\write_results_m2.py` -> this section

Notes:
- Scalar kernels (`frac_delay`, `rotator`) are a separate CMake target (`ntn_scalar`) compiled without `/arch:AVX2`; the AVX2 kernels (xsimd) are in `ntn_avx2`, compiled with `/arch:AVX2`. Proof: compile flags and disassembly below (`out/isa_check_m2.txt`).
- Departure from SPEC §6 (recorded in SPEC.md): the bindings (`_core`) and `ntn_scalar` (which holds the CPU check) are built without `/arch:AVX2`, and the `.pyd` without LTCG, so the import-time check runs before any AVX2 instruction.
- CPU requirement: AVX2, FMA3, BMI1, BMI2 and OS support for AVX state (XGETBV). `ntn_core` is compiled with `/arch:AVX2` (MSVC emits BMI1/BMI2 there; xsimd uses FMA3), so there is no runtime fallback: `import ntn_emu._core` raises ImportError and `GeometricChannel`/`Pipeline` throw for either kernel when a feature is missing. The scalar kernels are a test oracle and benchmark baseline.
- Fractional delay: 32 taps x 256 phases kept (decision from the design sweep below); Kaiser beta from the window design formula (not tuned).
- Derived limit (new check, not a loosened tolerance): the vectorised float32 sincos alone must stay below 1e-6 per sample (|(c + js) - e^{ja}| and phase), half of the 2e-6 rad tone-phase budget of M1a test 1 on the C++ path; the float32 angle rounding and the complex multiply take the rest.
- "Full chain" in the benchmarks is the geometric stage only (TDL arrives in M3, AWGN in M4).
- Throughput is not pass/fail (SPEC §3).

### Toolchain and dependencies
| Item | Value |
|---|---|
| MSVC (cl) | 19.50.35724.0 |
| CMake | cmake version 4.3.0-rc1 |
| Ninja | 1.13.2 |
| Python | 3.11.0 |
| pybind11 | 3.1.0 |
| scikit_build_core | 1.1.1 |
| numpy | 2.4.6 |
| scipy | 1.17.1 |
| skyfield | 1.55 |
| googletest (FetchContent) | `v1.17.0` (https://github.com/google/googletest.git) |
| xsimd (FetchContent) | `13.2.0` (https://github.com/xtensor-stack/xsimd.git) |
| benchmark (FetchContent) | `v1.9.4` (https://github.com/google/benchmark.git) |
| ntn_emu._core | _core.cp311-win_amd64.pyd; cpu_has_avx2() = True |

### Test suites (from out/ logs)
| Suite | Result |
|---|---|
| pytest fast | `127 passed, 13 deselected in 41.37s` |
| pytest slow | `13 passed, 127 deselected in 188.94s (0:03:08)` |
| ctest fast | `100% tests passed, 0 tests failed out of 95`, `Total Test time (real) =  10.01 sec` |
| ctest slow (NTN_SLOW=1, -L slow) | `100% tests passed, 0 tests failed out of 4`, `Total Test time (real) =  58.24 sec` |
| gtest, all tests with NTN_SLOW=1 (source of the measure lines below) | `[  PASSED  ] 95 tests.` |

### Fractional-delay design sweep (`out/fd_design_sweep_m2.txt`)
200 delays (100 grid fractions, 100 random in [-8, 8)), band-limited periodic noise |f| <= 0.4 fs. Worst/median: error against the FFT-domain reference. Peak: max |H(f) - exp(-j2 pi f d)| over |f| <= 0.4 fs and all delays, measured through the kernel.

| Taps | Phases | Kaiser beta | Design A (dB) | Worst (dB) | Median (dB) | Peak \|H - ideal\| | Peak (dB) |
|---|---|---|---|---|---|---|---|
| 16 | 128 | 4.6693 | 51.07 | -54.3 | -56.3 | 5.282e-03 | -45.5 |
| 16 | 256 | 4.6693 | 51.07 | -54.3 | -56.3 | 5.291e-03 | -45.5 |
| 24 | 128 | 7.2008 | 74.04 | -78.2 | -80.0 | 4.215e-04 | -67.5 |
| 24 | 256 | 7.2008 | 74.04 | -78.3 | -80.0 | 4.230e-04 | -67.5 |
| 32 | 128 | 9.7322 | 97.01 | -92.7 | -95.2 | 7.170e-05 | -82.9 |
| 32 (production) | 256 | 9.7322 | 97.01 | -99.8 | -102.1 | 3.860e-05 | -88.3 |

Figure: `docs/figures/m2_fd_design_sweep.png` (worst-case error vs taps, 128 and 256 phases, -60 dB target).

![Design sweep](figures/m2_fd_design_sweep.png)

### Acceptance tests 1-2: fractional-delay accuracy and slip continuity (gtest, scalar kernel)
Limits: -60 dB for each of 50 random delays; -60 dB overall and within +-64 samples of every slip. The negative control must exceed -60 dB.
```
FracDelayM2.RandomDelaysBandLimitedNoise: M2 test 1: 50 random delays, worst -99.8 dB (limit -60 dB)
SlipContinuity.FastRampKernelIncreasingDelay: kernel ramp +2.5e-4: 25 slips, overall -101.8 dB, worst +-64 window -105.9 dB (limit -60 dB)
SlipContinuity.FastRampKernelDecreasingDelay: kernel ramp -2.5e-4: 25 slips, overall -101.8 dB, worst +-64 window -105.5 dB (limit -60 dB)
SlipContinuity.NegativeControlShortDesignFails: negative control, 16 x 128 design, same ramp: overall -55.7 dB (must be > -60 dB)
SlipContinuity.RealisticRampThroughGeometricChannel: GeometricChannel ramp v=-7000 m/s (-23.35 ppm), L = 21: 23 slips, overall -101.9 dB, worst +-64 window -108.7 dB (limit -60 dB)
```

### Acceptance test 3: AVX2 vs scalar (limit -100 dB)
Lengths 0, 1, 7, 8, 9, 31, 33, 4095, 4097; buffers aligned and misaligned by 8 bytes; integer and fractional delays; NaN sentinels next to the allowed read range; output guard cells.
```
FracDelayAvx2.IsNotTheScalarKernel: AVX2 vs scalar, 4096 fractional delays: 3842 outputs differ in bits, error -138.7 dB
FracDelayAvx2.MatchesScalarLengthsAlignmentDelayMix: AVX2 vs scalar frac delay: worst -136.4 dB over lengths 0..4097, offsets 0/8 B, integer+fractional delay mix (limit -100 dB)
FracDelayAvx2.SingleSampleExactTapSpan: AVX2 vs scalar, single samples on the exact tap span: worst -123.2 dB
RotatorAvx2.IsNotTheScalarRotator: AVX2 vs scalar rotator, 4096 samples: 897 outputs differ in bits, error -149.4 dB
RotatorAvx2.MatchesScalarLengthsAlignment: AVX2 vs scalar rotator: worst -146.0 dB over lengths 0..4097, offsets 0/8 B (limit -100 dB)
```

### Vectorised sincos accuracy (derived limit 1e-6, see notes)
Reference: float64 sin/cos of the same float32 argument; range [-pi_f, pi_f] (the rotator's wrapped argument). Fast: 2^22 + 1 grid points plus special values; slow: every float32 in the range.
```
SincosAvx2.AccuracyOverWrappedRange: sincos_avx2 (xsimd), 4194317 args in [-pi_f, pi_f]: max |sin err| 6.892e-08, |cos err| 7.714e-08, |e^ja err| 9.172e-08 (at a = -2.41851687), phase 8.922e-08 rad, ||.|-1| 6.689e-08
SincosAvx2.AccuracyOverWrappedRange: CRT sinf/cosf (scalar rotator, for reference), 4194317 args in [-pi_f, pi_f]: max |sin err| 3.149e-08, |cos err| 3.148e-08, |e^ja err| 4.291e-08 (at a = -2.35024142), phase 4.291e-08 rad, ||.|-1| 4.279e-08
SincosAvx2.AccuracyExhaustiveSlow: sincos_avx2 (xsimd), exhaustive, 2157060024 args in [-pi_f, pi_f]: max |sin err| 6.957e-08, |cos err| 7.773e-08, |e^ja err| 9.337e-08 (at a = 2.41084123), phase 9.178e-08 rad, ||.|-1| 6.911e-08
```

### M1a closed-form tests on both kernels (gtest, channel level; same limits as M1b)
Phase < 2e-6 rad with no drift (last 1e5 samples), FFT peak within 0.1 Hz, integer shift bit-exact, quadratic phase within 0.01 Hz, time-varying delay < -90 dB. Kernel in the test name (`/Scalar`, `/Avx2`); `RotatorAvx2.TonePhase*` is the same tone test on the AVX2 rotator alone.
```
RotatorAvx2.TonePhase: AVX2 rotator tone rho0=1000000 v=-7000 N=1000000: max phase err 2.017e-07 rad, first 1e5 1.968e-07, last 1e5 1.953e-07, max ||y|-1| 6.91e-08
RotatorAvx2.TonePhase: AVX2 rotator tone rho0=1234567.89 v=6853.37 N=1000000: max phase err 2.075e-07 rad, first 1e5 2.048e-07, last 1e5 2.019e-07, max ||y|-1| 6.79e-08
RotatorAvx2.TonePhaseSlow: AVX2 rotator tone rho0=1000000 v=-7000 N=10000000: max phase err 2.094e-07 rad, first 1e5 1.968e-07, last 1e5 1.918e-07, max ||y|-1| 6.91e-08
RotatorAvx2.TonePhaseSlow: AVX2 rotator tone rho0=1234567.89 v=6853.37 N=10000000: max phase err 2.117e-07 rad, first 1e5 2.048e-07, last 1e5 1.918e-07, max ||y|-1| 6.91e-08
Kernels/GeoTone.Phase/Scalar: tone rho0=1000000 v=-7000 N=1000000: max phase err 1.637e-07 rad, last 1e5 1.627e-07 rad, max ||y|-1| 4.29e-08
Kernels/GeoTone.Phase/Scalar: tone rho0=1234567.89 v=6853.37 N=1000000: max phase err 1.651e-07 rad, last 1e5 1.602e-07 rad, max ||y|-1| 4.30e-08
Kernels/GeoTone.Phase/Avx2: tone rho0=1000000 v=-7000 N=1000000: max phase err 2.017e-07 rad, last 1e5 1.953e-07 rad, max ||y|-1| 6.91e-08
Kernels/GeoTone.Phase/Avx2: tone rho0=1234567.89 v=6853.37 N=1000000: max phase err 2.075e-07 rad, last 1e5 2.019e-07 rad, max ||y|-1| 6.79e-08
Kernels/GeoTone.PhaseSlow/Scalar: tone rho0=1000000 v=-7000 N=10000000: max phase err 1.648e-07 rad, last 1e5 1.619e-07 rad, max ||y|-1| 4.31e-08
Kernels/GeoTone.PhaseSlow/Scalar: tone rho0=1234567.89 v=6853.37 N=10000000: max phase err 1.698e-07 rad, last 1e5 1.632e-07 rad, max ||y|-1| 4.31e-08
Kernels/GeoTone.PhaseSlow/Avx2: tone rho0=1000000 v=-7000 N=10000000: max phase err 2.094e-07 rad, last 1e5 1.918e-07 rad, max ||y|-1| 6.91e-08
Kernels/GeoTone.PhaseSlow/Avx2: tone rho0=1234567.89 v=6853.37 N=10000000: max phase err 2.117e-07 rad, last 1e5 1.918e-07 rad, max ||y|-1| 6.91e-08
Kernels/GeoTone.FftPeakAtDoppler/Scalar: tone v=-7000: peak 46698.973328 Hz, f_d 46698.973328 Hz, error 0.00e+00 Hz
Kernels/GeoTone.FftPeakAtDoppler/Scalar: tone v=6853.37: peak -45720.763262 Hz, f_d -45720.763262 Hz, error 7.28e-12 Hz
Kernels/GeoTone.FftPeakAtDoppler/Avx2: tone v=-7000: peak 46698.973328 Hz, f_d 46698.973328 Hz, error 7.28e-12 Hz
Kernels/GeoTone.FftPeakAtDoppler/Avx2: tone v=6853.37: peak -45720.763262 Hz, f_d -45720.763262 Hz, error 7.28e-12 Hz
Kernels/GeoShift.IntegerShiftBitExact/Scalar: integer shift k=0: latency L = 16
Kernels/GeoShift.IntegerShiftBitExact/Scalar: integer shift k=1: latency L = 15
Kernels/GeoShift.IntegerShiftBitExact/Scalar: integer shift k=37: latency L = 0
Kernels/GeoShift.IntegerShiftBitExact/Scalar: integer shift k=70001: latency L = 0
Kernels/GeoShift.IntegerShiftBitExact/Avx2: integer shift k=0: latency L = 16
Kernels/GeoShift.IntegerShiftBitExact/Avx2: integer shift k=1: latency L = 15
Kernels/GeoShift.IntegerShiftBitExact/Avx2: integer shift k=37: latency L = 0
Kernels/GeoShift.IntegerShiftBitExact/Avx2: integer shift k=70001: latency L = 0
Kernels/GeoQuadratic.InstantaneousFrequency/Scalar: quadratic phase: max |f_inst - (-fc rho_dot / c)| = 1.725e-03 Hz
Kernels/GeoQuadratic.InstantaneousFrequency/Avx2: quadratic phase: max |f_inst - (-fc rho_dot / c)| = 2.032e-03 Hz
Kernels/GeoLatency.ApproachingTrack/Scalar: approaching track (v=-7000 m/s, 1 s): min D = -176.3240, latency L = 192
Kernels/GeoLatency.ApproachingTrack/Avx2: approaching track (v=-7000 m/s, 1 s): min D = -176.3240, latency L = 192
Kernels/GeoLatency.MinimalForConstantFractionalDelay/Scalar: D=2.3, L=13, block 1: |out[n+L]| with impulse at n+L = 3.326e-05, at n+L+1 = 0.000e+00
Kernels/GeoLatency.MinimalForConstantFractionalDelay/Scalar: D=2.3, L=13, block 4096: |out[n+L]| with impulse at n+L = 3.326e-05, at n+L+1 = 0.000e+00
Kernels/GeoLatency.MinimalForConstantFractionalDelay/Avx2: D=2.3, L=13, block 1: |out[n+L]| with impulse at n+L = 3.326e-05, at n+L+1 = 0.000e+00
Kernels/GeoLatency.MinimalForConstantFractionalDelay/Avx2: D=2.3, L=13, block 4096: |out[n+L]| with impulse at n+L = 3.326e-05, at n+L+1 = 0.000e+00
Kernels/GeoLatency.MinimalForApproachingTrack/Scalar: approaching, D(N-1)=-22.9495, L=38: |last output| with impulse at last input = 9.576e-06
Kernels/GeoLatency.MinimalForApproachingTrack/Avx2: approaching, D(N-1)=-22.9495, L=38: |last output| with impulse at last input = 9.576e-06
Kernels/GeoDelay.TimeVaryingDelayMatchesClosedForm/Scalar: time-varying delay: min D = -22.9795, latency L = 38
Kernels/GeoDelay.TimeVaryingDelayMatchesClosedForm/Scalar: time-varying delay vs closed form: -103.8 dB
Kernels/GeoDelay.TimeVaryingDelayMatchesClosedForm/Avx2: time-varying delay: min D = -22.9795, latency L = 38
Kernels/GeoDelay.TimeVaryingDelayMatchesClosedForm/Avx2: time-varying delay vs closed form: -103.8 dB
```

### Acceptance test 4: AVX2 vs golden, block invariance
Golden: float64 model on the M1b band-limited OFDM stimulus, real pass, window crossing a segment knot.

| Kernel | Window | Error power | Limit | Log |
|---|---|---|---|---|
| Scalar | 0.1 s | -105.48 dB | -55.0 dB | `pytest_fast_m2.txt` |
| Avx2 | 0.1 s | -105.48 dB | -55.0 dB | `pytest_fast_m2.txt` |
| Scalar | 1.0 s | -106.25 dB | -55.0 dB | `pytest_slow_m2.txt` |
| Avx2 | 1.0 s | -106.25 dB | -55.0 dB | `pytest_slow_m2.txt` |

Block sizes 1, 7 and 4096 (and a mixed pattern); limits: scalar bit-identical, AVX2 < -120 dB.
```
PipelineBlocks.ScalarBlockSizeInvarianceBitIdentical: block invariance track: D 0.1302 -> -7.7403 (8 integer slips), knot at 200000, L = 23
PipelineBlocks.Avx2BlockSizeInvariance: AVX2 block invariance, blocks {1} vs {4096}: -inf dB, 0 of 400023 samples differ in bits (limit -120 dB)
PipelineBlocks.Avx2BlockSizeInvariance: AVX2 block invariance, blocks {7} vs {4096}: -inf dB, 0 of 400023 samples differ in bits (limit -120 dB)
PipelineBlocks.Avx2BlockSizeInvariance: AVX2 block invariance, blocks {1,7,4096,3,1000,13} vs {4096}: -inf dB, 0 of 400023 samples differ in bits (limit -120 dB)
```

### CPU feature check
```
CpuCheck.MatchesWindowsFeatureQuery: cpu_has_avx2() = 1, IsProcessorFeaturePresent(AVX2) = 1
CpuCheck.Avx2AndOsBitsMatchWindowsQuery: cpu features: AVX2 1, FMA3 1, BMI1 1, BMI2 1, OS AVX state 1; IsProcessorFeaturePresent(AVX2) 1
CpuCheck.HideSwitchEachFeatureIsReportedByName: NTN_TEST_HIDE_AVX2=1 -> "ntn: this CPU/OS lacks instruction-set features the engine is built for (missing: AVX2; required: AVX2, FMA3, BMI1, BMI2, OS AVX state (XGETBV))"
CpuCheck.HideSwitchEachFeatureIsReportedByName: NTN_TEST_HIDE_AVX2=AVX2 -> "ntn: this CPU/OS lacks instruction-set features the engine is built for (missing: AVX2; required: AVX2, FMA3, BMI1, BMI2, OS AVX state (XGETBV))"
CpuCheck.HideSwitchEachFeatureIsReportedByName: NTN_TEST_HIDE_AVX2=FMA3 -> "ntn: this CPU/OS lacks instruction-set features the engine is built for (missing: FMA3; required: AVX2, FMA3, BMI1, BMI2, OS AVX state (XGETBV))"
CpuCheck.HideSwitchEachFeatureIsReportedByName: NTN_TEST_HIDE_AVX2=BMI1 -> "ntn: this CPU/OS lacks instruction-set features the engine is built for (missing: BMI1; required: AVX2, FMA3, BMI1, BMI2, OS AVX state (XGETBV))"
CpuCheck.HideSwitchEachFeatureIsReportedByName: NTN_TEST_HIDE_AVX2=BMI2 -> "ntn: this CPU/OS lacks instruction-set features the engine is built for (missing: BMI2; required: AVX2, FMA3, BMI1, BMI2, OS AVX state (XGETBV))"
CpuCheck.HideSwitchEachFeatureIsReportedByName: NTN_TEST_HIDE_AVX2=OSAVX -> "ntn: this CPU/OS lacks instruction-set features the engine is built for (missing: OS AVX state (XGETBV); required: AVX2, FMA3, BMI1, BMI2, OS AVX state (XGETBV))"
```
Import-time check: `python/tests/test_cpu_guard.py` (ImportError naming each hidden feature, in a fresh interpreter; part of the fast pytest suite above).

### Instruction-set checks (`out/isa_check_m2.txt`, `out/fma_contraction_m2.txt`)
ISA check (object files of our sources; the full report also covers the `.pyd` build tree):
```
ISA check, M2 (2026-10-27 22:45 W. Europe Standard Time)
git: c215b25 + uncommitted changes
scope: inspects the object files (.obj) of our sources, NOT the final linked _core.pyd or
       ntn_tests.exe; the link step is covered by check 4 (functions shared across the boundary).
[build] hello          src/hello.cpp                      /arch:AVX2=no
[build] ntn_avx2       src/ntn/frac_delay_avx2.cpp        /arch:AVX2=yes ok
[build] ntn_avx2       src/ntn/rotator_avx2.cpp           /arch:AVX2=yes ok
[build] ntn_bench      bench/bench_m2.cpp                 /arch:AVX2=no
[build] ntn_core       src/ntn/geometric_channel.cpp      /arch:AVX2=yes ok
[build] ntn_core       src/ntn/pipeline.cpp               /arch:AVX2=yes ok
[build] ntn_core       src/ntn/range_track.cpp            /arch:AVX2=yes ok
[build] ntn_isa_probe  tests/cpp/isa_probe_avx2.cpp       /arch:AVX2=yes ok
[build] ntn_scalar     src/ntn/cpu_check.cpp              /arch:AVX2=no  ok
[build] ntn_scalar     src/ntn/frac_delay.cpp             /arch:AVX2=no  ok
[build] ntn_scalar     src/ntn/rotator.cpp                /arch:AVX2=no  ok
[build] ntn_tests      tests/cpp/test_avx2_kernels.cpp    /arch:AVX2=no
[build] ntn_tests      tests/cpp/test_cpu_check.cpp       /arch:AVX2=no
[build] ntn_tests      tests/cpp/test_frac_delay.cpp      /arch:AVX2=no
[build] ntn_tests      tests/cpp/test_geometric_channel.cpp /arch:AVX2=no
[build] ntn_tests      tests/cpp/test_pipeline.cpp        /arch:AVX2=no
[build] ntn_tests      tests/cpp/test_range_track.cpp     /arch:AVX2=no
[build] ntn_tests      tests/cpp/test_rotator.cpp         /arch:AVX2=no
[build] ntn_tests      tests/cpp/test_slip_continuity.cpp /arch:AVX2=no
[build] ntn_tests      tests/cpp/test_smoke.cpp           /arch:AVX2=no
== 2-3. Disassembly per object (dumpbin /disasm) ==
build dir                  target         object                               insns    VEX    ymm  zmm
build                      hello          hello.cpp.obj                           53      0      0    0
build                      ntn_avx2       frac_delay_avx2.cpp.obj                275    145     61    0
build                      ntn_avx2       rotator_avx2.cpp.obj                  8215   3587    991    0
build                      ntn_bench      bench_m2.cpp.obj                     10940      0      0    0
build                      ntn_core       geometric_channel.cpp.obj            27732    761     53    0
build                      ntn_core       pipeline.cpp.obj                       475     11      2    0
build                      ntn_core       range_track.cpp.obj                  27314    957     30    0
build                      ntn_isa_probe  isa_probe_avx2.cpp.obj                6089   2676    679    0
build                      ntn_scalar     cpu_check.cpp.obj                     3001      0      0    0
build                      ntn_scalar     frac_delay.cpp.obj                    5616      0      0    0
build                      ntn_scalar     rotator.cpp.obj                         77      0      0    0
build                      ntn_tests      test_avx2_kernels.cpp.obj            23194      0      0    0
build                      ntn_tests      test_cpu_check.cpp.obj               12227      0      0    0
build                      ntn_tests      test_frac_delay.cpp.obj              23533      0      0    0
build                      ntn_tests      test_geometric_channel.cpp.obj       49984      0      0    0
build                      ntn_tests      test_pipeline.cpp.obj                19256      0      0    0
build                      ntn_tests      test_range_track.cpp.obj             22737      0      0    0
build                      ntn_tests      test_rotator.cpp.obj                  8024      0      0    0
build                      ntn_tests      test_slip_continuity.cpp.obj         13831      0      0    0
build                      ntn_tests      test_smoke.cpp.obj                    3864      0      0    0
== 4. Functions defined in both a non-AVX2 and an AVX2 object of one link unit ==
(the linker keeps one copy, chosen arbitrarily. FAIL if the AVX2 copy has VEX code, the
 non-AVX2 side is a scalar kernel or the bindings, and some object references the function;
 an unreferenced copy can never execute. Test/hello objects are reported as INFO.)
[build] relocation targets found: 10471
[build] INFO fabsf: VEX code in ntn_core:geometric_channel.cpp, ntn_core:range_track.cpp; also defined in ntn_bench:bench_m2.cpp, ntn_tests:test_avx2_kernels.cpp, ntn_tests:test_frac_delay.cpp, ntn_tests:test_geometric_channel.cpp, ntn_tests:test_pipeline.cpp; referenced by no object in this link unit (dead copy)
[build] shared functions: 282 (with ntn_scalar/_core: 100); with VEX code in the AVX2 copy: 1 (referenced and shared with ntn_scalar/_core: 0)
== 5. Dynamic initialisers (.CRT$XC*) in AVX2 objects ==
build                      ntn_avx2       frac_delay_avx2.cpp.obj            0
build                      ntn_avx2       rotator_avx2.cpp.obj               0
build                      ntn_core       range_track.cpp.obj                0
build                      ntn_core       geometric_channel.cpp.obj          0
build                      ntn_core       pipeline.cpp.obj                   0
== 6. Instruction families beyond AVX/AVX2 in /arch:AVX2 objects ==
(cpu_check tests: AVX2, OS AVX state, BMI1, BMI2, FMA3; any other family FAILs)
build                      ntn_avx2       frac_delay_avx2.cpp.obj            FMA3=10
build                      ntn_avx2       rotator_avx2.cpp.obj               BMI2=24 FMA3=57
build                      ntn_core       geometric_channel.cpp.obj          BMI1=3 BMI2=39
build                      ntn_core       pipeline.cpp.obj                   -
build                      ntn_core       range_track.cpp.obj                BMI1=3 BMI2=39
build                      ntn_isa_probe  isa_probe_avx2.cpp.obj             BMI2=18 FMA3=42
INFO list: shared functions whose AVX2 copy has VEX code, not counted as hazards:
  - [build] fabsf (VEX copy in ntn_core:geometric_channel.cpp, ntn_core:range_track.cpp; unreferenced)
  - [py-cp311-cp311-win_amd64] fabsf (VEX copy in ntn_core:geometric_channel.cpp, ntn_core:range_track.cpp; unreferenced)
  - [py-cp311-cp311-win_amd64] fabsl (VEX copy in ntn_core:geometric_channel.cpp, ntn_core:range_track.cpp; unreferenced)
RESULT: PASS
```

FMA contraction under `/arch:AVX2 /fp:precise` (each FMA attributed via a `/FAs` listing to an xsimd header or to our own code):
```
FMA contraction check, M2 (2026-10-27 22:45 W. Europe Standard Time)
git: c215b25 + uncommitted changes

== FMAs per /arch:AVX2 object ==
target         source                              total  xsimd  own (contracted)  build obj
ntn_avx2       src/ntn/frac_delay_avx2.cpp            10     10                 0         10
ntn_core       src/ntn/geometric_channel.cpp           0      0                 0          0
ntn_core       src/ntn/pipeline.cpp                    0      0                 0          0
ntn_core       src/ntn/range_track.cpp                 0      0                 0          0
ntn_avx2       src/ntn/rotator_avx2.cpp               57     57                 0         57
ntn_isa_probe  tests/cpp/isa_probe_avx2.cpp           42     42                 0         42

== Probe functions (tests/cpp/isa_probe_avx2.cpp) ==
ntn_isa_probe_contract_f32     plain C++ float a*b + c                  FMAs:   0 (xsimd 0, own 0)
ntn_isa_probe_contract_f64     plain C++ double a*b + c                 FMAs:   0 (xsimd 0, own 0)
ntn_isa_probe_contract_loop    plain C++ loop acc = a*b + acc           FMAs:   0 (xsimd 0, own 0)
ntn_isa_probe_xsimd_mul_add    xsimd operators a*b + c                  FMAs:   0 (xsimd 0, own 0)
ntn_isa_probe_xsimd_fma        explicit xsimd::fma (positive control)   FMAs:   1 (xsimd 1, own 0)
ntn_isa_probe_xsimd_sincos     xsimd::sincos                            FMAs:   8 (xsimd 8, own 0)

== Conclusion ==
MSVC contracts a*b + c to FMA under /arch:AVX2 /fp:precise: NO (all four probe forms stay mul + add)
FMAs attributed to our own code in ntn_avx2/ntn_core (contractions): 0
All FMAs in ntn_avx2/ntn_core come from xsimd (explicit xsimd::fma, xsimd math routines).

RESULT: PASS (consistency checks)
```

### Benchmarks (`out/bench_m2.txt`, `out/bench_m2_summary.json`)
| Item | Value |
|---|---|
| Commit | `c215b25` + uncommitted changes |
| Time | 2026-10-27T22:36:40+01:00 |
| CPU | AMD Ryzen 5 5500U with Radeon Graphics, 6 cores / 12 threads, max clock 2100 MHz |
| Caches | L2 3072 KiB total (512 KiB per core), L3 8192 KiB |
| Clock before / after | 2659 MHz effective (2100 MHz nominal x 126.6% processor performance, 2026-10-27T22:36:47+01:00) / 3679 MHz effective (2100 MHz nominal x 175.2% processor performance, 2026-10-27T22:38:38+01:00) |
| Power plan (powercfg) | Balanced (`381b4222-f694-41f0-9685-ff5bb260df2e`) |
| Power mode (Windows overlay) | Best performance (`ded574b5-45a0-4f42-8737-46345c09c238`) |
| AC power | True (battery 35%) |
| Build | Release; Google Benchmark library build: release |
| cpu_has_avx2() / default kernel | True / Avx2 |
| Flags | `--benchmark_repetitions=5 --benchmark_enable_random_interleaving=true --benchmark_min_time=1s` |

Msps = 1e6 complex samples per second; median of the repetitions, spread = (max - min) / median in parentheses (⚠ marks spread > 5%). Single-threaded. Throughput is reported, not pass/fail.

#### Per stage, block 4096 samples

| Stage | Scalar Msps | AVX2 Msps | AVX2 / scalar | Working set (est.), scalar / AVX2 |
|---|---|---|---|---|
| Fractional delay (kernel) | 12.30 (1.2%) | 80.17 (2.1%) | 6.52x | 129 KiB (L2) / 225 KiB (L2) |
| Rotator (kernel) | 69.53 (0.9%) | 265.32 (0.3%) | 3.82x | 64 KiB (L2) / 64 KiB (L2) |
| Geometric channel (fs 7.68 Msps) | 10.26 (0.2%) | 55.18 (0.6%) | 5.38x | 227 KiB (L2) / 323 KiB (L2) |

#### Full chain (M2: geometric stage only; TDL arrives in M3, AWGN in M4)

| fs (Msps) | Block | Scalar Msps | Scalar RTF | AVX2 Msps | AVX2 RTF | AVX2 / scalar | Working set (est.), AVX2 |
|---|---|---|---|---|---|---|---|
| 7.68 | 4096 | 10.24 (1.1%) | 1.33 | 55.02 (0.4%) | 7.16 | 5.37x | 323 KiB (L2) |
| 7.68 | 65536 | 10.20 (1.2%) | 1.33 | 54.84 (2.6%) | 7.14 | 5.38x | 2243 KiB (L3) |
| 30.72 | 4096 | 10.27 (0.5%) | 0.33 | 55.02 (0.5%) | 1.79 | 5.36x | 330 KiB (L2) |
| 30.72 | 65536 | 10.23 (1.6%) | 0.33 | 54.81 (1.3%) | 1.78 | 5.36x | 2250 KiB (L3) |

Notes:
- Real-time factor (RTF) = Msps / fs; below 1 means slower than real time.
- "Full chain" is `Pipeline::process`, which in M2 contains only the geometric stage (range evaluation, fractional delay, carrier rotation, input history). TDL fading (M3) and AWGN (M4) will lower these numbers.
- Each benchmark processes the same block repeatedly, so the data stay in cache as far as the working set allows. The working set is an estimate (input, output and history buffers, per-4096-sample scratch, filter table), compared with the per-core L2 and the shared L3; the 65536-sample block shows the effect of a working set beyond L2.
- Scalar and AVX2 repetitions were interleaved in random order (`--benchmark_enable_random_interleaving`).
- Rates use wall-clock time (`UseRealTime`). A first run used Google Benchmark's default basis, process CPU time, which on Windows advances in 15.625 ms ticks; its medians and spreads were snapped to that grid (kept as `out/bench_m2_cputime_quantized.*`, not used here).
- Spread > 5%: none.
- The scalar kernels are a separate target (`ntn_scalar`) compiled without `/arch:AVX2` (verified by `scripts/check_isa.py`), so the comparison is not inflated by compiler auto-vectorisation of the scalar path.

### End-to-end demo (`out/m2_results.json`, `out/demo_m2.txt`)
Window: pass 'high', 2026-10-20T04:29:34.407706Z + 1.0 s (same as M1a/M1b); kernel Avx2; cpu_has_avx2() = True; latency L = 166, tau0 = 26244 samples.

| Sub | Centre (s) | rho_dot (m/s) | Measured tone offset (Hz) | Predicted -(fc+f_tone)rho_dot/c (Hz) | Error (Hz) |
|---|---|---|---|---|---|
| 0 | 0.050 | -5899.5026 | 39412.3450 | 39412.3449 | +0.0001 |
| 1 | 0.150 | -5898.1597 | 39403.3743 | 39403.3738 | +0.0004 |
| 2 | 0.250 | -5896.8145 | 39394.3868 | 39394.3866 | +0.0002 |
| 3 | 0.350 | -5895.4668 | 39385.3846 | 39385.3834 | +0.0012 |
| 4 | 0.450 | -5894.1167 | 39376.3626 | 39376.3641 | -0.0015 |
| 5 | 0.550 | -5892.7642 | 39367.3283 | 39367.3287 | -0.0003 |
| 6 | 0.650 | -5891.4094 | 39358.2766 | 39358.2772 | -0.0006 |
| 7 | 0.750 | -5890.0521 | 39349.2095 | 39349.2096 | -0.0001 |
| 8 | 0.850 | -5888.6924 | 39340.1250 | 39340.1260 | -0.0009 |
| 9 | 0.950 | -5887.3303 | 39331.0263 | 39331.0262 | +0.0000 |

Result: max \|measured - predicted\| = 0.0015 Hz over 10 sub-windows (limit 1 Hz) -> PASS.
AVX2 vs scalar on the same 1 s input: -138.0 dB, 6882148 of 7680000 samples differ in bits. Wall time AVX2 / scalar (informational, not a benchmark): 0.18 s / 0.81 s.

## M3 — NTN-TDL-C/D fading
_Pending._

## M4 — Receiver harness
_Pending._

## M5 — Cross-checks
_Pending._
