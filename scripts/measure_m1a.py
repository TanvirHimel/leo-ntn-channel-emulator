"""Re-run the M1a golden-model acceptance measurements (tests 1-4) and print a results table.

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


def git_rev():
    try:
        sha = subprocess.run(["git", "rev-parse", "--short", "HEAD"], cwd=ROOT,
                             capture_output=True, text=True, check=True).stdout.strip()
        dirty = subprocess.run(["git", "status", "--porcelain"], cwd=ROOT,
                               capture_output=True, text=True, check=True).stdout.strip()
        return sha + ("-dirty" if dirty else "")
    except (OSError, subprocess.CalledProcessError):
        return "unknown"


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
    for d in tf.DELAYS:
        rows.append((f"3 fractional delay, d={d}", db(tf.fracdelay_error_db(d)), "-80 dB"))
    for k in tf.NOSNAP_KS:
        rows.append((f"3b integer delay, snap off (kernel path), k={k}",
                     db(tf.nosnap_shift_error_db(k)), "-120 dB"))
    rows.append(("4 quadratic phase, inst. frequency", f"{tg.quadratic_freq_error_hz():.2e} Hz",
                 "0.01 Hz"))
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
