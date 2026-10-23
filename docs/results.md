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
_Pending._

## M3 — NTN-TDL-C/D fading
_Pending._

## M4 — Receiver harness
_Pending._

## M5 — Cross-checks
_Pending._
