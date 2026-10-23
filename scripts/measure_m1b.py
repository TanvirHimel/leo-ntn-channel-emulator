"""Re-run the M1b acceptance measurements and print the docs/results.md M1b section (Markdown).

The measurement functions are the ones the tests assert on (python/tests/test_core*.py), so
the printed numbers are exactly what the test suite checks. The C++ closed-form values come
from the gtest binary's "[ measure  ]" lines. Suite results, the demo numbers and the
toolchain are read from files: out/pytest_*_m1b.txt, out/ctest_*_m1b.txt,
out/m1b_results.json, build/CMakeFiles/*/CMakeCXXCompiler.cmake. Nothing is typed in.

    .\\.venv\\Scripts\\python.exe scripts\\measure_m1b.py [--slow] > out\\measure_m1b.txt
"""

import argparse
import json
import os
import platform
import re
import subprocess
import sys
import time
from importlib import metadata
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python"))
sys.path.insert(0, str(ROOT / "python" / "tests"))

import numpy as np  # noqa: E402

import test_core as tc  # noqa: E402
import test_core_vs_golden as tv  # noqa: E402
from ntn_emu import _core  # noqa: E402
from ntn_emu.golden import geometric  # noqa: E402
from ntn_emu.orbit import geometry, passes, tle  # noqa: E402
from ntn_emu.waveform import ofdm  # noqa: E402
from test_golden_geometric import SHIFT_KS, TONE_CASES  # noqa: E402

OUT = ROOT / "out"


def git(*args):
    return subprocess.run(["git", *args], cwd=ROOT, capture_output=True, text=True).stdout.strip()


def msvc_version():
    for f in sorted((ROOT / "build" / "CMakeFiles").glob("*/CMakeCXXCompiler.cmake")):
        m = re.search(r'set\(CMAKE_CXX_COMPILER_VERSION "([^"]+)"\)', f.read_text())
        if m:
            return m.group(1)
    return "unknown"


def tool_version(cmd):
    try:
        return subprocess.run(cmd, capture_output=True, text=True).stdout.splitlines()[0].strip()
    except (OSError, IndexError):
        return "unknown"


def suite_line(path, pattern):
    """Last line of an out/ suite log matching pattern, or a note that the file is missing."""
    p = OUT / path
    if not p.exists():
        return f"(missing {p.name})"
    hits = [ln.strip() for ln in p.read_text(encoding="utf-8", errors="replace").splitlines()
            if re.search(pattern, ln)]
    return hits[-1] if hits else f"(no match in {p.name})"


def gtest_measures(slow):
    env = dict(os.environ, NTN_SLOW="1" if slow else "0")
    res = subprocess.run([str(ROOT / "build" / "ntn_tests.exe")], capture_output=True, text=True,
                         env=env)
    lines = [ln[len("[ measure  ] "):] for ln in res.stdout.splitlines() if ln.startswith("[ measure")]
    summary = [ln for ln in res.stdout.splitlines() if ln.startswith(("[  PASSED", "[  FAILED", "[  SKIPPED ] 1"))]
    return lines, summary


def unshaped_band_split(seg):
    """Diagnostic behind the stimulus choice: 0.1 s window, unshaped OFDM (as generated)."""
    from scipy import signal
    n_win = int(round(0.1 * tv.FS))
    n_start = tv.FS - n_win // 2
    x = ofdm.generate(100, ofdm.OfdmConfig(modulation="qpsk", seed=3)).samples.astype(np.complex64)
    tau0 = int(np.floor(seg.rho(np.array([n_start]))[0] / tv.C * tv.FS))
    n0, n1, co = seg.to_arrays()
    yg = geometric.apply(x, fs=tv.FS, fc=tv.FC, n0=n0, n1=n1, coeffs=co, tau0_samples=tau0,
                         n_start=n_start)
    p = _core.Pipeline(_core.RangeTrack(tv.FS, n0, n1, co),
                       _core.GeoConfig(fc_hz=tv.FC, fs_hz=tv.FS, tau0_samples=tau0))
    p.reset(n_start)
    lat = p.latency()
    a = p.process(x).astype(complex)[lat:]
    b = yg[:len(yg) - lat]
    e = a - b
    kw = {"fs": 1.0, "window": ("kaiser", 14.0), "nperseg": 4096, "return_onesided": False,
          "scaling": "spectrum"}
    f, px = signal.welch(x.astype(complex), **kw)
    _, pe = signal.welch(e, **kw)
    _, pb = signal.welch(b, **kw)
    rows = []
    for lo, hi in [(0.0, 0.3), (0.3, 0.4), (0.4, 0.45), (0.45, 0.5)]:
        m = (np.abs(f) >= lo) & (np.abs(f) < hi)
        rows.append((lo, hi, 10 * np.log10(px[m].sum() / px.sum()),
                     10 * np.log10(pe[m].sum() / pb.sum())))
    full = 10 * np.log10(np.sum(np.abs(e) ** 2) / np.sum(np.abs(b) ** 2))
    return full, rows


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--slow", action="store_true", help="also 1e7-sample tone and 1 s window")
    args = ap.parse_args()
    t_start = time.perf_counter()
    head = git("rev-parse", "--short", "HEAD")
    dirty = git("status", "--porcelain")
    print(f"M1b measurements (base commit {head}{', uncommitted M1b working tree' if dirty else ''}, "
          f"{time.strftime('%Y-%m-%d %H:%M %Z')})\n")

    print("### Toolchain")
    print("| Item | Value |\n|---|---|")
    print(f"| MSVC (cl) | {msvc_version()} (from build/CMakeFiles/*/CMakeCXXCompiler.cmake) |")
    print(f"| CMake | {tool_version(['cmake', '--version'])} |")
    print(f"| Ninja | {tool_version(['ninja', '--version'])} |")
    print(f"| Python | {platform.python_version()} |")
    for pkg in ("pybind11", "scikit_build_core", "numpy", "scipy", "skyfield"):
        print(f"| {pkg} | {metadata.version(pkg)} |")
    print(f"| ntn_emu._core | {Path(_core.__file__).name}; cpu_has_avx2() = {_core.cpu_has_avx2()} |")

    print("\n### Test suites (from out/ logs)")
    print("| Suite | Result |\n|---|---|")
    print(f"| pytest fast | `{suite_line('pytest_fast_m1b.txt', r'passed|failed')}` |")
    print(f"| pytest slow | `{suite_line('pytest_slow_m1b.txt', r'passed|failed')}` |")
    print(f"| ctest fast | `{suite_line('ctest_fast_m1b.txt', r'tests passed|tests failed')}`, "
          f"`{suite_line('ctest_fast_m1b.txt', r'Total Test time')}` |")
    print(f"| ctest slow (NTN_SLOW=1, -L slow) | `{suite_line('ctest_slow_m1b.txt', r'tests passed|tests failed')}`, "
          f"`{suite_line('ctest_slow_m1b.txt', r'Total Test time')}` |")

    lines, summary = gtest_measures(args.slow)
    print("\n### C++ closed-form tests (gtest, scalar path; `[ measure ]` lines"
          f"{', NTN_SLOW=1' if args.slow else ''})")
    print("```")
    for ln in lines + summary:
        print(ln)
    print("```")

    print("\n### Acceptance test 1: M1a tests 1, 2, 4 through `_core` (pytest, test_core.py)")
    print("| Test | Measured | Limit |\n|---|---|---|")
    counts = [1_000_000] + ([10_000_000] if args.slow else [])
    for count in counts:
        for rho0, v in TONE_CASES:
            e, et, me = tc.tone_phase_stats(rho0, v, count)
            print(f"| 1 phase, rho0={rho0:g} m, v={v:g} m/s, N={count:.0e} | {e:.2e} rad "
                  f"(last 1e5: {et:.2e}; max \\|\\|y\\|-1\\| {me:.1e}) | 2e-6 rad |")
    for rho0, v in TONE_CASES:
        print(f"| 1 FFT peak, rho0={rho0:g} m, v={v:g} m/s, N=1e+06 | "
              f"{tc.tone_peak_error_hz(rho0, v):.2e} Hz | 0.1 Hz |")
    for k in SHIFT_KS:
        lat, mism = tc.integer_shift_result(k)
        print(f"| 2 integer shift, k={k} | L = {lat}; {mism} differing samples "
              f"({'bit-exact' if mism == 0 else 'NOT exact'}) | bit-exact |")
    print(f"| 4 quadratic phase, inst. frequency | {tc.quadratic_freq_error_hz():.2e} Hz | 0.01 Hz |")

    print("\n### Acceptance test 3: C++ vs golden, real pass (pytest, test_core_vs_golden.py)")
    sat, ue = tle.load_pinned("starlink").sat, geometry.ue_position()
    seg = tv.knot_segments(sat, ue, passes.load_pinned_passes()[tv.PASS_LABEL])
    print(f"Stimulus: band-limited QPSK OFDM; power share at |f| >= 0.4 fs: "
          f"{tv.power_share_above_db(tv.stimulus(0.1)):.1f} dB (0.1 s, Welch Kaiser(14)).\n")
    print("| Window (crosses the knot) | Error power | L | tau0 | Limit |\n|---|---|---|---|---|")
    for w in [0.1] + ([1.0] if args.slow else []):
        err, lat, tau0, _ = tv.core_vs_golden(seg, w)
        print(f"| {w:g} s | {err:.2f} dB | {lat} | {tau0} | {tv.LIMIT_DB:g} dB |")
    err_shift, *_ = tv.core_vs_golden(seg, 0.01, shift=1)
    print(f"| negative control: 0.01 s, comparison off by 1 sample | {err_shift:.2f} dB | | | must be > -20 dB |")

    full, rows = unshaped_band_split(seg)
    print(f"\nDiagnostic, unshaped OFDM (0.1 s, as generated): full-band error {full:.2f} dB.")
    print("| Band | Input power share | Error power, rel. golden output |\n|---|---|---|")
    for lo, hi, ps, es in rows:
        print(f"| {lo:.2f}-{hi:.2f} fs | {ps:.1f} dB | {es:.1f} dB |")

    print("\n### End-to-end demo (out/m1b_results.json, out/demo_m1b.txt)")
    res = json.loads((OUT / "m1b_results.json").read_text(encoding="utf-8"))
    run = res["run"]
    print("| Quantity | Value |\n|---|---|")
    print(f"| Error power, full band (recorded, not asserted) | {res['error_full_db']:.2f} dB |")
    inband_def = res["inband_def"].replace("|", "\\|")  # Markdown table cell
    print(f"| Error power, {inband_def} | {res['error_inband_db']:.2f} dB (limit {res['limit_db']:g} dB) |")
    pf, pi = res["per_10ms_full_db"], res["per_10ms_inband_db"]
    print(f"| Per 10 ms, full band (min / median / max) | {pf['min']:.1f} / {pf['median']:.1f} / {pf['max']:.1f} dB |")
    print(f"| Per 10 ms, in band (min / median / max) | {pi['min']:.1f} / {pi['median']:.1f} / {pi['max']:.1f} dB |")
    print(f"| Tone Doppler on the C++ output, max \\|measured - predicted\\| | {res['doppler_max_err_hz']:.4f} Hz (limit {res['doppler_limit_hz']:g} Hz) |")
    print(f"| Latency L / tau0 | {run['latency_samples']} / {run['tau0_samples']} samples |")
    print(f"| Wall time, C++ scalar / golden (informational, not a benchmark) | "
          f"{run['cpp_wall_s']:.2f} s / {run['golden_wall_s']:.1f} s |")
    print(f"\n_measure_m1b.py ran in {time.perf_counter() - t_start:.1f} s._")


if __name__ == "__main__":
    main()
