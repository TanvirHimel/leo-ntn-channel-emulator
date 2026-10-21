# NTN Channel Emulator — Specification

Portfolio-scale LEO non-terrestrial-network channel emulator. S-band (~2 GHz, `f_c` runtime),
3GPP Rel-17 (TR 38.811) channel models, service-link **downlink** only. C++20 engine (MSVC x64,
Ninja), Python 3.11 control plane. Conservative scope: when in doubt, it is out.

## 1. Key numbers (2 GHz, 600 km LEO, first-order estimates)
| Quantity | Value | Consequence |
|---|---|---|
| One-way delay | 2–7 ms | relative-delay mode (tau − tau0) keeps buffers small |
| Max Doppler | ~±46 kHz ≈ ±3 SCS @15 kHz | Stage-0 integer-CFO search required |
| Max Doppler rate | ~630 Hz/s (overhead TCA) | receiver needs a 2nd-order FLL |
| Delay drift | ~23 ppm → 1 sample slip / ~43k samples | integer slip handling in frac-delay |
| f_c·tau | ~1.4e7 cycles | phase evaluated in float64, wrapped before float32 |
| Full pass @30.72 Msps | ~150 GB cf32 | emulate short windows only |

## 2. Decisions (fixed)
- **Captures**: short windows (default 0.5–1 s) at chosen epochs (AOS / mid / TCA / LOS).
  The engine streams with continuous state; full passes are possible but never stored by default.
- **UE location**: a generic, rounded city-centre coordinate (e.g. 52.5° N, 13.4° E, 50 m;
  0.1° rounding). No personal address anywhere in the repo.
- **Geometry**: `rho(t_rx) = |r_sat(t_rx − tau) − r_ue(t_rx)|`, light-time solved iteratively
  in GCRS via Skyfield (converged < 1 µm). Satellite from SGP4 using a pinned TLE. A geometric
  (no light-time) mode exists for cross-checks.
- **Bulk delay**: relative. Output delay = `tau(t) − tau0`, where `tau0` is the integer-sample
  delay at window start, recorded in SigMF. **Phase always uses the full `tau(t)`.**
- **Segments**: cubic Hermite on `rho` with knots every `T_k` = 1 s (configurable), built
  from `rho`, `rho_dot` at the knots, so C1 holds by construction. Coeffs use local time `u = (n − n_k)/fs`.
- **Time base**: exact int64 global sample index `n`, segment lookup by `n`, local `u` in double.
- **Phase**: `phi(n) = −2π f_c rho(u)/c`. Carrier cycles `f_c·rho(u)/c` are computed in float64
  (Horner), reduced to fractional cycles in [−0.5, 0.5), and *then* cast to float32 for SIMD
  sincos. Phase is never accumulated. `f_d = −f_c rho_dot/c` is derived, never applied separately.
- **Fractional delay (C++)**: polyphase windowed-sinc (32 taps × 256 phases, Kaiser window,
  linear interpolation between phases). Target error < −60 dB over |f| < 0.4 fs against ideal sinc.
  The integer part is an int64 read pointer, so slips are exact. The phase-0 filter is an exact unit impulse.
- **Fractional delay (golden)**: float64 windowed sinc evaluated directly at the exact delay of each
  output sample (long kernel, ≥ 128 taps, no table). It is more accurate than the C++ path, so
  "C++ vs golden" measures the C++ error. It processes input in chunks of ~64k samples to bound memory.
- **Pre-compensation**: none in the engine; full physical Doppler. An optional constant
  residual-CFO mode exists only for the EVM sweep.
- **TDL**: NTN-TDL-C and NTN-TDL-D (TR 38.811), both LOS. Sum-of-sinusoids (SoS) fading with
  angles and phases drawn in Python (seeded NumPy RNG) and passed to C++, so the engine is
  deterministic. Tap delays are band-limited windowed-sinc FIRs summed into one channel FIR,
  refreshed every 64 samples. DS, K (table default) and UE max Doppler (default 5 Hz) are runtime
  params, not elevation-dependent within a window. Satellite Doppler comes **only** from the range track.
- **Pipeline order**: TDL → geometric delay/Doppler → AWGN (SNR at the receiver).
- **Waveform**: "NR-numerology CP-OFDM (not NR-compliant)". Default: SCS 15 kHz, 5 MHz
  (25 PRB, 300 SC), N_FFT 512, fs 7.68 Msps, normal CP. The frequency-selective demo uses
  fs 30.72 Msps (see M3). NR PSS m-sequence in symbol 0 of each 10 ms frame; comb pilots (every
  6th SC, every symbol) carrying known QPSK; random QPSK/16QAM data. No SSS/PBCH/LDPC.
- **Receiver**:
  - Stage 0: PSS correlation over ±50 kHz in 2 kHz steps, giving coarse timing and CFO.
  - Stage 1: CP-correlation fine CFO and timing, tracked by a 2nd-order FLL and a timing loop.
  - Per symbol: pilot-based CPE correction, LS channel estimate with linear interpolation,
    ZF equalization, then EVM.
- **Noise**: AWGN at a given SNR, seeded in C++. No path loss, antennas or link budget.
- **SIMD**: xsimd with an AVX2 baseline. Scalar kernels are kept as test oracle and bench baseline.
- **Sources rule**: TR 38.811 tables and OAI equations are transcribed **only from fetched source
  documents** (3GPP archive, OAI GitLab at a fixed commit). Citations go in the file: document,
  version, table/section, URL, commit hash, fetch date. If a source cannot be fetched, the item is
  marked `PENDING-SOURCE`, its tests are skipped with that reason, and the user is told. Nothing
  is transcribed from memory.
- **Protection timeline**: `python/ntn_emu/golden/` and `tests/reference/` are under construction
  during M1a. They become protected (the CLAUDE.md "never edit" rule) **only once the golden model
  passes its closed-form tests**, at the M1a gate (tag `golden-v1`). CLAUDE.md already states the
  rule from golden-v1 onward. After that, any change needs the user's sign-off.

## 3. Performance budget (modest laptop)
- **Fast suite** (default `pytest`, default `ctest`): under ~2 min total on the dev laptop.
- **Slow suite**: pytest marker `slow` (`-m slow`) and ctest label `slow` (`-L slow`), run on
  demand. Fast-mode variants use reduced counts with tolerances widened to a 4σ statistical bound
  for the reduced sample size. Slow mode uses the full counts and tight tolerances.
  | Test | Fast | Slow |
  |---|---|---|
  | Phase-drift tone test | 1e6 samples | 1e7 samples |
  | SoS realizations (PDP/K/ACF) | 200 | 2000 |
  | Stage-0 trials | 30 at SNR 5 dB (0 misses) | 200 at SNR 0 dB (≥ 99%) |
  | Tracking run | 0.5 s | 2 s |
- **Golden model cost**: the golden fractional delay evaluates a ≥ 128-tap kernel per output
  sample, so a 1 s window at 7.68 Msps costs ~1e9 kernel evaluations. Fast-suite tests that run
  the golden fractional delay use windows ≤ 0.25 s.
- **Demo windows**: 0.5–1 s at 7.68 Msps (~30–60 MB), ≤ 0.25 s at 30.72 Msps.
  Generated samples go in `out/` (git-ignored), and the scripts overwrite them.
- **Build**: Release, `cmake --build build -j 2`, FetchContent with `GIT_SHALLOW`, dep
  tests/examples disabled. No sanitizers and no Debug builds by default.
- **Throughput is not pass/fail.** Benchmarks report what is measured (Msps, real-time factor,
  including < 1 at 30.72 Msps if that is the result) with CPU model, clock, power plan and commit.

## 4. Repo layout
```
CMakeLists.txt                 core lib, gtest tests, benches, _core.pyd
pyproject.toml                 scikit-build-core, wraps the same CMakeLists
cmake/                         Deps.cmake (FetchContent, pinned tags), Warnings.cmake
src/ntn/                       C++ engine (namespace ntn)
  range_track.hpp/.cpp         RangeSegment, RangeTrack (validation, lookup, rho/rho_dot)
  frac_delay.hpp/.cpp          polyphase table + scalar/AVX2 kernels
  rotator.hpp/.cpp             float64 phase eval + wrap + float32 SIMD sincos
  geometric_channel.hpp/.cpp   delay + Doppler from one RangeTrack
  tdl_channel.hpp/.cpp         SoS gains + per-block FIR            (M3)
  awgn.hpp/.cpp                seeded complex Gaussian noise         (M4)
  pipeline.hpp/.cpp            TDL → geo → AWGN, block streaming
  cpu_check.hpp/.cpp           AVX2 CPUID guard
src/py/bindings.cpp            ntn_emu._core (pybind11)
src/hello.cpp                  existing sanity target (kept)
tests/cpp/                     GoogleTest closed-form tests (label `slow` for heavy ones)
tests/reference/               pinned TLEs (tle/*.txt) + SOURCES.md citations; protected from golden-v1
bench/                         Google Benchmark
python/ntn_emu/
  orbit/{tle,geometry,passes}.py   Skyfield/SGP4, light-time range, pass finder
  segments.py                      Hermite fit, C1 check, → _core
  waveform/{ofdm,pss}.py           CP-OFDM generator
  channel/{tdl_tables,tdl_params}.py  38.811 tables (cited), SoS parameter draw
  golden/{fracdelay,geometric,tdl,pipeline}.py   float64 reference; protected from golden-v1
  rx/{pss_search,cp_sync,loops,equalize,evm}.py
  sigmf_io.py, plotting.py
python/tests/                  pytest (closed-form, golden-vs-C++, cross-checks; marker `slow`)
python/tests/oracles/oai_circular.py   OAI equations transcribed from fetched source (M5)
scripts/                       demo_m1a.py … demo_m5.py, fetch_tle.py (demo only, not tests)
out/                           generated SigMF files (git-ignored)
docs/results.md                every measured number
docs/figures/                  generated plots (small PNGs, committed)
README.md                      skeleton from M1a; claims only with numbers in results.md
```

## 5. Interfaces

### 5.1 C++ (namespace `ntn`)
```cpp
using cf32 = std::complex<float>;
struct RangeSegment {          // rho(u) = c[0] + c[1]u + c[2]u^2 + c[3]u^3; meters, u in seconds
  int64_t n0, n1;              // [n0, n1) in global samples
  std::array<double,4> c;
};
class RangeTrack {             // contiguous; C1 checked at construction (|Δrho|<1e-6 m, |Δrho_dot|<1e-6 m/s)
 public:
  RangeTrack(double fs, std::vector<RangeSegment>);
  double rho(int64_t n) const;  double rho_dot(int64_t n) const;  double tau(int64_t n) const;
};
struct GeoConfig { double fc_hz, fs_hz; int64_t tau0_samples; bool delay = true, doppler = true;
                   double residual_ref_doppler_hz = 0.0; };
struct TdlConfig {             // all arrays come from Python; the engine draws nothing random
  std::vector<double> tap_delay_s, tap_power_lin;  std::vector<uint8_t> tap_is_los;
  double k_factor_lin, fd_max_hz, los_doppler_hz;
  std::vector<double> sos_alpha, sos_phi;   // [n_taps × n_sinusoids], row-major
  int update_block = 64; int fir_half_len = 16; };
struct AwgnConfig { double snr_db; uint64_t seed; bool enabled = false; };
enum class Kernel { Scalar, Avx2 };
class Pipeline {
 public:
  Pipeline(RangeTrack, GeoConfig, std::optional<TdlConfig>, AwgnConfig, Kernel = Kernel::Avx2);
  void process(std::span<const cf32> in, std::span<cf32> out);  // any block size, same output
  int64_t sample_index() const;  void reset(int64_t n_start);
};
```
Block-size invariance is a hard requirement: on the scalar path, block sizes 1, 7 and 4096 give
bit-identical output. On the AVX2 path they agree to < −120 dB (SIMD tails and FMA can round differently).

### 5.2 Python ↔ C++ (pybind11 module `ntn_emu._core`)
- `RangeTrack(fs, n0: int64[N], n1: int64[N], coeffs: float64[N,4])`. It takes segment
  coefficients only, **never sampled rho/tau/phase vectors**.
- `GeoConfig`, `TdlConfig`, `AwgnConfig`, `Kernel`, `Pipeline(...)`,
  `Pipeline.process(x: complex64[M]) -> complex64[M]` (GIL released, zero-copy input), `cpu_has_avx2()`.
- SigMF I/O stays in Python (`sigmf_io.py`); C++ never touches files in M1–M5.

### 5.3 SigMF (Python writer/reader)
- `core:datatype = "cf32_le"`, `core:sample_rate`, `core:version`, `core:extensions` declaring `ntn`.
- `captures[0]`: `core:sample_start = 0`, **`core:frequency = f_c`**, `core:datetime` (UTC window start).
- `global.ntn:*`: `tau0_s`, `tau0_samples`, `tle_name`, `tle_epoch`, `ue_lla` (rounded),
  `tdl_model`, `delay_spread_s`, `snr_db`, waveform params, `git_sha`. Annotations mark generator
  PSS positions (ground truth for the clean waveform only).
- Files are written and read with explicit little-endian `complex64` (`<c8`).

## 6. Windows / MSVC build rules
- FetchContent with pinned tags, `GIT_SHALLOW`, and deps declared `SYSTEM` (needs CMake ≥ 3.25): xsimd, GoogleTest (`gtest_force_shared_crt ON`),
  Google Benchmark (`BENCHMARK_ENABLE_TESTING OFF`). pybind11 comes from the venv (`pybind11.get_cmake_dir()`).
- `/MD` CRT everywhere, matching python311.dll. `.pyd` is built only for Release/RelWithDebInfo
  (Debug would need python311_d). Deps are statically linked into `_core.pyd` with no extra DLLs,
  because Python ≥ 3.8 ignores PATH for DLL search.
- Our targets use `/W4 /permissive- /utf-8 /EHsc`; `/WX` is a CMake option (`NTN_WARNINGS_AS_ERRORS`),
  default OFF; deps are exempt. `/fp:precise` everywhere, never `/fp:fast`.
- `/arch:AVX2` on core/pyd targets. Scalar kernels are a separate CMake target compiled
  **without** `/arch:AVX2`, so scalar-vs-AVX2 benchmarks are honest (the compiler cannot
  auto-vectorize the scalar path with AVX2). They are selected at runtime via `Kernel`.
  `cpu_check` raises a clean error at import or startup if AVX2 is missing.
- Use `std::numbers::pi` (no `M_PI`). Define `NOMINMAX` if `<windows.h>` is ever included.
- `pip install -e .` runs inside the same `cmd /c VsDevCmd.bat … &&` environment, with
  `.\.venv\Scripts\python.exe` (3.11) for ABI match. Parallelism is capped with `CMAKE_BUILD_PARALLEL_LEVEL=2`.
- Benchmarks run only from Release.

## 7. Out of scope
- Ka-band.
- Uplink, timing advance, feeder link, regenerative payload.
- Engine-side Doppler/delay pre-compensation (only the constant residual-CFO mode exists).
- CDL, spatial models, MIMO, antenna patterns, beams, handover.
- Path loss, link budget, atmospheric, ionospheric and scintillation effects.
- Elevation-dependent DS/K within a window.
- A full NR PHY (SSS, PBCH, PDSCH, LDPC, HARQ).
- Real-time SDR streaming, multi-threading, runtime SIMD dispatch (Highway).
- Linux/macOS builds, running OAI itself, live TLE fetching in tests, sanitizer and Debug CI builds.

**Stretch (not part of done)**:
- A C++ CLI with scenario JSON, C++ SigMF I/O, and a CLI ≡ pybind11 equivalence test.
- NTN-TDL-A/B.
- Elevation-dependent DS/K.

## 8. Milestones
Each milestone ends with tests green (fast suite, plus the slow suite once at the milestone gate),
numbers in `docs/results.md`, and a commit. Tolerances are starting values. Tightening is fine;
loosening needs a note in results.md.

### M1a — Python golden slice (no C++)
Scope: `orbit/` (light-time and geometric range), `segments.py`, `waveform/`, `sigmf_io.py`,
`golden/geometric.py` + `golden/fracdelay.py`, pinned TLEs (main case: a Starlink-class satellite
at ~550 km; secondary: the ISS), fetched once, with fetch date and URL recorded in
`tests/reference/SOURCES.md`. At least 3 passes of the main satellite with different maximum
elevations are pinned in `tests/reference/` (M5 needs them, and the folder is protected after golden-v1),
`.gitignore` for `out/`, README skeleton (purpose, honest labels, build placeholder, "results pending").
Acceptance tests (pytest, on the golden model):
1. **Tone/Doppler**: DC input, linear rho = rho0 + v·t, delay off. Output phase equals
   −2π f_c (rho0 + v t)/c, |error| < 1e-7 rad (float64 resolution at ~1e7 cycles is ~1e-8 rad) over 1e6 samples (slow: 1e7). The FFT peak
   sits at −f_c v/c within 0.1 Hz (interpolated).
2. **Integer shift**: f_c = 0, tau − tau0 = k/fs exactly. Output equals input shifted by k with relative error < 1e-12.
3. **Fractional delay**: band-limited periodic noise, constant fractional delay d. Error against the
   FFT-domain reference exp(−j2πfd) is < −80 dB over |f| < 0.4 fs.
4. **Quadratic phase**: cubic segment with known coefficients. Instantaneous frequency equals
   −f_c rho_dot(t)/c within 0.01 Hz.
5. **Segment fit**: Hermite track against dense (10 ms) Skyfield light-time range.
   max |Δrho| < 1 mm, |Δrho_dot| < 1 mm/s, and C1 at all knots. Slow: all three pinned passes,
   full pass (AOS → LOS). Fast (reduced): the highest pass only, ±60 s around TCA.
6. **Light-time**: rho_lt − rho_geo matches −(v_sat,GCRS · û)·tau to < 1 cm.
7. **SigMF**: write → read round-trip is bit-identical. `core:datatype == "cf32_le"`.
   `core:frequency` is present in captures and absent from global.
Gate: on all-green, tag `golden-v1`; `golden/` and `tests/reference/` become protected.
E2E: `scripts/demo_m1a.py`. Pinned TLE + rounded UE → pass → 1 s mid-pass window → OFDM (+ test
tone) → golden channel → SigMF in `out/` → figure in `docs/figures/` (spectrogram + measured vs
predicted f_d). Asserts |measured − predicted| < 1 Hz.

### M1b — C++ scalar engine + pybind11
Scope: CMake deps, pyproject/scikit-build-core, `RangeTrack`, scalar frac-delay + rotator,
`GeometricChannel`, `Pipeline` (geo only), `cpu_check`, bindings. No C++ file I/O.
Acceptance tests:
1. gtest closed-form: M1a tests 1, 2 and 4 re-run on the C++ scalar path. Phase error < 2e-6 rad
   (float32 output), no drift over 1e6 samples (slow: 1e7). The integer shift is bit-exact (C++ phase-0 impulse path).
2. `RangeTrack` rejects non-contiguous segments and C1 violations.
3. **C++ vs golden** on random OFDM with a real-pass segment set: error power < −55 dB rel.
4. **Block invariance** (scalar): block sizes 1, 7 and 4096 give bit-identical output.
5. `import ntn_emu._core` works from a clean venv install (no DLL errors).
E2E: `scripts/demo_m1b.py`. Same window as M1a through `_core` → SigMF → figure overlaying
golden vs C++ with the error power printed and recorded in results.md.

### M2 — Accuracy, AVX2, benchmarks
Scope: final polyphase design, AVX2 kernels (frac-delay, rotator), Google Benchmark.
Acceptance tests:
1. **Frac-delay accuracy**: 50 random fractional delays, band-limited noise (|f| < 0.4 fs).
   Error against the FFT-domain reference is < −60 dB.
2. **Slip continuity**: linear delay ramp over ≥ 20 integer slips. Error against the ideal reference
   is < −60 dB everywhere, including ±64 samples around each slip.
3. **AVX2 vs scalar**: relative error < −100 dB.
4. **AVX2 vs golden** < −55 dB. AVX2 block invariance: block sizes 1, 7 and 4096 agree to < −120 dB.
Bench (reported, not pass/fail): Msps per stage, scalar vs AVX2, and real-time factor at
7.68 and 30.72 Msps, with CPU, clock and power plan. Accuracy numbers from tests 1–4 are recorded.
The results.md notes state that the scalar kernels are a separate target built without `/arch:AVX2`.
E2E: `scripts/demo_m2.py`. A 1 s pass window through the AVX2 pipeline, re-checking the M1a
tone criterion, then the benchmarks, with the table pasted into results.md.

### ⛔ Publish gate (after M2, before M3 starts)
The repo must be publishable as-is:
- README with build steps, architecture sketch, honest scope labels, and M1–M2 figures.
- `docs/results.md` with the first accuracy and benchmark numbers plus machine details.
- Fast suite green from a clean clone.
- Every README claim traced to results.md.
- A LICENSE chosen by the user.

M3 does not start until the user confirms the gate.

### M3 — NTN-TDL-C/D fading
Scope: `tdl_tables.py` transcribed from fetched TR 38.811 (cited; `PENDING-SOURCE` if not
fetched), SoS draw, `TdlChannel` (scalar + AVX2), `golden/tdl.py`, pipeline TDL → geo.
Note: golden/tdl.py is a new file added under the protection rule with user sign-off. Its
closed-form tests must pass before it is protected.
**Delay spread and bandwidth**: at 7.68 Msps the table delay spreads give a nearly flat channel
over 5 MHz, which is expected and stated. The frequency-selective demo uses fs 30.72 Msps with
delay-spread scaling (normalized table delays × DS, with DS taken from the fetched 38.811
delay-spread values for a stated scenario/elevation, or an explicitly labelled larger DS chosen
so the coherence bandwidth < occupied bandwidth). The README says which and why.
Acceptance tests:
1. **Static response**: fading frozen. Measured H(f) equals Σ a_k exp(−j2πf τ_k) within −50 dB over |f| < 0.4 fs.
2. **PDP**: mean tap powers match the table × DS scaling. Fast: 200 realizations, 4σ bound. Slow:
   2000 realizations, ±0.3 dB. Total power is 1 within 0.05 dB.
3. **Rician K** of the LOS tap matches the table within the 4σ bound (fast) or ±0.5 dB (slow).
4. **Doppler**: NLOS autocorrelation matches J0(2π f_D τ) for f_D τ ∈ [0, 2]
   (fast 4σ, slow |err| < 0.05). The LOS component frequency equals the configured LOS Doppler within 0.01 Hz.
5. **Block update**: 64-sample updates against per-sample (golden) evaluation at f_D = 5 Hz give < −60 dB.
6. **Determinism**: the same seed gives C++ ≡ golden within −55 dB.
7. **Ordering (closed form)**: static frozen TDL with response H(f), constant Doppler f_d (linear
   rho, delay off), input tone at baseband f0. TDL → geo produces a tone at f0 + f_d with complex
   amplitude H(f0) (within −50 dB), **not** H(f0 + f_d). The reversed order gives H(f0 + f_d).
   This uses a test DS where |H(f0) − H(f0 + f_d)| is ≥ 6 dB above the tolerance.
E2E: `scripts/demo_m3.py`. TDL-C and TDL-D at a mid-pass window: 7.68 Msps (flat) and 30.72 Msps
≤ 0.25 s (frequency-selective). Figures: PDP, |H(f,t)|, Doppler spectrum. Numbers and TDL FIR Msps go in results.md.

### M4 — Receiver harness
Scope: `rx/` (PSS search, CP sync, FLL and timing loop, equalization, EVM), AWGN in C++, residual-CFO mode.
Acceptance tests:
1. **AWGN EVM**: no other impairments, SNR 10–30 dB. EVM² = 1/SNR within 0.5 dB.
2. **ICI theory**: static residual CFO ε ∈ [0.01, 0.1] SCS at high SNR. EVM matches ICI power ≈ (π ε)²/3 within 1 dB.
3. **Stage 0**: CFO uniform in ±50 kHz. Fast: SNR 5 dB (miss rate essentially zero), 30 trials,
   0 misses. Slow: SNR 0 dB, 200 trials, ≥ 99%.
   |coarse CFO error| ≤ 1 kHz and timing error ≤ 1 sample.
4. **Fine CFO** (CP, SNR 20 dB): bias < 0.5% SCS, std within 2× of the closed-form CP-ML variance.
5. **Tracking**: CFO ramp 630 Hz/s + delay ramp 23 ppm from linear-rho segments. Residual CFO
   < 1% SCS RMS, timing error < 0.5 sample RMS, no cycle slips (fast 0.5 s, slow 2 s).
E2E: `scripts/demo_m4.py`. AOS/mid/TCA/LOS windows (0.5 s) through TDL-C + geo + AWGN 20 dB;
the receiver must acquire all four. EVM-vs-residual-CFO plot shows the static sweep, theory curve,
and the operating points from the pass windows.

### M5 — Cross-checks, final docs
Scope: OAI oracle, synthetic circular-orbit geometry provider, Skyfield TCA cross-check.
Acceptance tests:
1. **OAI oracle**: `oracles/oai_circular.py` is transcribed from OAI source fetched at a fixed
   commit (file path, commit hash and URL cited; `PENDING-SOURCE` if not fetched).
   - A synthetic circular orbit through *our* segment pipeline, with light-time off and a
     non-rotating Earth (OAI's assumptions), must agree within |Δtau| < 1 ns and |Δf_d| < 0.1 Hz.
   - A light-time-on run documents the expected difference.
2. **TCA vs dense Skyfield minimum range**: for 3 pinned passes, an independent reference
   t_min comes from dense (10 ms) Skyfield **geometric** topocentric range plus a parabolic fit
   around the minimum. Our geometric-mode TCA is the root of rho_dot from the Hermite segments.
   - **Tolerance |Δt| < 1 ms.** Justification: the Hermite rho_dot error (< 1 mm/s) divided by
     rho_ddot at TCA (≥ ~10 m/s² for any pass with min range < 3000 km) is < 0.1 ms. Parabolic
     interpolation of a 10 ms grid on a quartic-dominated minimum is ≪ 1 ms. A 10× margin goes on top.
   - Light-time mode: the TCA shift relative to geometric equals rho_min/c × (1 ± 0.1). The
     ±10% covers the UE's inertial velocity, ≤ ~0.47 km/s against the satellite's ~7.5 km/s.
3. **Full regression**: fast and slow suites are green from a clean clone.
E2E: fresh clone → build per README → `scripts/demo_m5.py` regenerates every figure and number
in results.md. README claims are cross-checked line by line against results.md.

## 9. Definition of done
- M1a–M5 complete and the publish gate passed.
- Accuracy numbers and benchmarks (scalar vs AVX2, real-time factor as measured) in `docs/results.md`.
- README with build steps, architecture, figures and honest labels, making only claims backed by results.md.
