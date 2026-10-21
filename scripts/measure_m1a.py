"""Re-run the M1a acceptance measurements and print a results table.

Golden model: tests 1-4. Orbit/segments: tests 5-6, derivative checks, and the measured SGP4
velocity-vs-position inconsistency (the reason rho_dot is position-derived, see CLAUDE.md).

The measurement functions are the ones the tests assert on (python/tests), so the printed
numbers are exactly what the test suite checks. Output is Markdown, for docs/results.md.

    .\\.venv\\Scripts\\python.exe scripts\\measure_m1a.py          # test 1 at 1e6 samples
    .\\.venv\\Scripts\\python.exe scripts\\measure_m1a.py --slow   # test 1 also at 1e7
"""

import argparse
import platform
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python"))
sys.path.insert(0, str(ROOT / "python" / "tests"))

import numpy as np  # noqa: E402
import scipy  # noqa: E402

import test_golden_fracdelay as tf  # noqa: E402
import test_golden_geometric as tg  # noqa: E402
import test_orbit as to  # noqa: E402
import test_segments as ts  # noqa: E402
from ntn_emu.orbit import geometry, passes, tle  # noqa: E402


def sgp4_velocity_inconsistency(sat, p):
    """max |v_SGP4 - d r/dt| in GCRS, 41 epochs over TCA +-200 s of the given pass.

    d r/dt is the 2nd-order central difference of SGP4 position with h = 10 ms (truncation
    ~h^2/6 * |r'''| ~ 1e-7 m/s, negligible against the cm/s effect being measured).
    """
    off = np.linspace(-200.0, 200.0, 41)
    h = 0.01

    def pos(o):
        return sat.at(geometry.times_from_offsets(p["tca"], o)).position.m

    v = sat.at(geometry.times_from_offsets(p["tca"], off)).velocity.m_per_s
    return np.max(np.linalg.norm(v - (pos(off + h) - pos(off - h)) / (2 * h), axis=0))


def geometric_range(sat, ue, t0, off):
    return geometry.range_geometric(sat, ue, geometry.times_from_offsets(t0, off)).rho


def sgp4_solver_steps(sat, ue, p, dt=0.005, threshold_m=1.5e-6, fit_s=2.0):
    """Range steps from SGP4 Kepler-solver iteration switches over a full pass.

    Steps are located where |3rd difference| of the 5 ms geometric range exceeds threshold_m
    (time-quantization noise gives ~2e-7 m). Each step size is the jump between quintic fits
    over fit_s on either side, evaluated at the step. Pass edges (< fit_s) are skipped.
    """
    duration = np.floor((p["los"] - p["aos"]) * 86_400.0)
    off = np.arange(0.0, duration, dt)
    rho = geometry._range(sat, ue, geometry.times_from_offsets(p["aos"], off), False)[0]
    d3 = np.abs(np.diff(rho, 3))
    w = int(fit_s / dt)
    steps = []
    for i in np.argsort(-d3):
        if d3[i] < threshold_m:
            break
        j = i + 2  # first sample after the jump (3rd difference spans i .. i+3)
        if j < w or j > len(off) - w or any(abs(off[j] - t) < 1.0 for t, _ in steps):
            continue
        tc = off[j] - dt / 2
        left, right = slice(j - w, j - 2), slice(j + 2, j + w)
        fl = np.polyval(np.polyfit(off[left] - tc, rho[left], 5), 0.0)
        fr = np.polyval(np.polyfit(off[right] - tc, rho[right], 5), 0.0)
        steps.append((off[j], fr - fl))
    return sorted(steps)


def time_quantization_noise(sat, ue, p, at_s=60.0):
    """RMS / max residual of geometric range vs a degree-7 fit over 2 s at 10 ms.

    Taken at AOS + at_s, where u_hat . v_sat is large (at TCA the time error is nearly
    perpendicular to u_hat and barely shows in range). No solver step falls in the window
    for the pinned high pass (steps at 23.1 s and 259.8 s).
    """
    x = np.arange(-1.0, 1.0 + 1e-9, 0.01)
    rho = geometric_range(sat, ue, geometry.times_from_offsets(p["aos"], at_s), x)
    res = rho - np.polyval(np.polyfit(x, rho - rho.mean(), 7), x) - rho.mean()
    return np.sqrt(np.mean(res**2)), np.max(np.abs(res))


def git_rev():
    try:
        sha = subprocess.run(["git", "rev-parse", "--short", "HEAD"], cwd=ROOT,
                             capture_output=True, text=True, check=True).stdout.strip()
        dirty = subprocess.run(["git", "status", "--porcelain"], cwd=ROOT,
                               capture_output=True, text=True, check=True).stdout.strip()
        return sha + ("-dirty" if dirty else "")
    except (OSError, subprocess.CalledProcessError):
        return "unknown"


def to_labels():
    from conftest import PASS_LABELS
    return PASS_LABELS


def db(x):
    return "exact (0)" if x == -np.inf else f"{x:.1f} dB"


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--slow", action="store_true", help="also run test 1 at 1e7 samples")
    args = ap.parse_args()
    counts = [tg.TONE_COUNT_FAST] + ([tg.TONE_COUNT_SLOW] if args.slow else [])

    rows = []
    t_start = time.perf_counter()
    for count in counts:
        for rho0, v in tg.TONE_CASES:
            ph, mag = tg.tone_phase_error_rad(rho0, v, count)
            rows.append((f"1 phase, rho0={rho0:g} m, v={v:g} m/s, N={count:.0e}",
                         f"{ph:.2e} rad (| \\|y\\|-1 | {mag:.1e})", "1e-7 rad"))
            pk = tg.tone_peak_error_hz(rho0, v, count)
            rows.append((f"1 FFT peak, rho0={rho0:g} m, v={v:g} m/s, N={count:.0e}",
                         f"{pk:.2e} Hz", "0.1 Hz"))
    for k in tg.SHIFT_KS:
        e = tg.integer_shift_rel_error(k)
        rows.append((f"2 integer shift, k={k} (snap path)",
                     "exact (0)" if e == 0 else f"{e:.2e}", "1e-12 rel"))
    for d in tf.DELAYS + tf.NEAR_INT_DELAYS:
        rows.append((f"3 fractional delay, d={d}", db(tf.fracdelay_error_db(d)), "-80 dB"))
    for eps in tf.NOSNAP_EPS:
        for k in tf.NOSNAP_KS:
            label = f"k={k}" if eps == 0 else f"k={k} + {eps:g}"
            rows.append((f"3b integer delay, snap off (kernel path), {label}",
                         db(tf.nosnap_shift_error_db(k, eps)), "-120 dB"))
    rows.append(("4 quadratic phase, inst. frequency", f"{tg.quadratic_freq_error_hz():.2e} Hz",
                 "0.01 Hz"))

    sat, ue = tle.load_pinned("starlink").sat, geometry.ue_position()
    pinned = passes.load_pinned_passes()
    fits = [("high", ts.FAST_HALF_WINDOW_S)] + ([(lb, None) for lb in to_labels()] if args.slow else [])
    for label, half in fits:
        t0, dur = ts.fit_window(pinned[label], half)
        d_rho, d_rate = ts.segment_fit_errors(sat, ue, t0, dur)
        span = f"TCA +-{half:g} s" if half else "full pass"
        rows.append((f"5 segment fit, {label} pass, {span}: rho", f"{d_rho:.2e} m", "1e-3 m"))
        rows.append((f"5 segment fit, {label} pass, {span}: rho_dot", f"{d_rate:.2e} m/s", "1e-3 m/s"))
    for label in to_labels():
        rows.append((f"6 light-time first order, {label} pass",
                     f"{to.lighttime_first_order_error(sat, ue, pinned[label]):.2e} m", "1e-2 m"))
    rows.append(("rho_dot routine, circular orbit (closed form)",
                 f"{to.circular_diff_error():.2e} m/s", "1e-6 m/s"))
    for mode in ("geometric", "lighttime"):
        d = to.rate_step_diffs(sat, ue, pinned["high"], mode)
        rows.append((f"rho_dot step consistency h=0.1 vs 0.2 s, high pass, {mode}",
                     f"max {d.max():.2e} m/s, median {np.median(d):.2e} m/s",
                     f"max {to.RATE_STEP_LIMIT_M_S:.0e} m/s"))
    rows.append(("SGP4 velocity vs d(position)/dt, high pass TCA +-200 s (informational)",
                 f"{sgp4_velocity_inconsistency(sat, pinned['high']):.2e} m/s", "-"))
    k_phase = 2.0 * np.pi * tg.FC / tg.C          # rad of carrier phase per m of range
    k_dopp = tg.FC / tg.C                          # Hz per m/s
    h = geometry.RATE_STEP_S
    for label in to_labels():
        for t_s, s in sgp4_solver_steps(sat, ue, pinned[label]):
            rows.append((f"SGP4 solver step, {label} pass, AOS+{t_s:.2f} s (informational)",
                         f"{s * 1e6:+.2f} um range; carrier phase {abs(s) * k_phase:.1e} rad; "
                         f"rho_dot artefact <= {8 * abs(s) / (12 * h):.1e} m/s "
                         f"= {8 * abs(s) / (12 * h) * k_dopp:.1e} Hz", "-"))
    rms, mx = time_quantization_noise(sat, ue, pinned["high"])
    rows.append(("Time-quantization range noise, high pass AOS+60 s (informational)",
                 f"rms {rms:.1e} m, max {mx:.1e} m", "-"))
    elapsed = time.perf_counter() - t_start

    print(f"M1a golden acceptance measurements  (commit {git_rev()}, "
          f"{time.strftime('%Y-%m-%d %H:%M %Z')})")
    print(f"Python {platform.python_version()}, NumPy {np.__version__}, SciPy {scipy.__version__}, "
          f"{platform.processor() or platform.machine()}, {platform.platform()}")
    print()
    print("| Test | Measured | Limit |")
    print("|---|---|---|")
    for name, meas, tol in rows:
        print(f"| {name} | {meas} | {tol} |")
    print()
    print(f"Measurement time: {elapsed:.1f} s")


if __name__ == "__main__":
    main()
