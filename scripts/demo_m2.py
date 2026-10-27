"""M2 end-to-end demo: a 1 s pass window through the AVX2 pipeline, re-checking the M1a tone criterion.

Window and input are those of demo_m1a.py / demo_m1b.py: highest pinned pass, 1 s centred at
(AOS+TCA)/2, light-time Hermite segments (1 s knots), tau0 = floor(tau(0) fs), NR-numerology
CP-OFDM (QPSK, seed 1) plus a test tone at +2.8 MHz (-10 dB), cast to complex64.
The pipeline uses the default kernel, which must be Kernel.Avx2. Its output is flushed with
latency() zeros and the first L samples dropped, so index i is channel time i.
Checks and records:
  - M1a tone criterion on the AVX2 output: per 0.1 s sub-window, measured tone frequency vs
    -(fc + f_tone) rho_dot / c at the sub-window centre; max |error| < 1 Hz (asserted);
  - AVX2 vs scalar on the same input: error power relative to the scalar output (recorded).
Wall times are informational, not benchmarks (scripts/bench_m2.py measures throughput).

    .\\.venv\\Scripts\\python.exe scripts\\demo_m2.py   ->  out/m2_results.json (print: out/demo_m2.txt)
"""

import json
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python"))
sys.path.insert(0, str(ROOT / "scripts"))

import numpy as np  # noqa: E402

import demo_m1a as m1a  # noqa: E402  (tone-Doppler measurement and limit)
import demo_m1b as m1b  # noqa: E402  (same window, segments, tau0 and input)
from ntn_emu import _core  # noqa: E402
from ntn_emu.orbit import geometry  # noqa: E402

FC, FS = m1a.FC, m1a.FS
RESULTS = ROOT / "out" / "m2_results.json"


def run_pipeline(seg, tau0, x, kernel=None):
    n0, n1, coeffs = seg.to_arrays()
    args = (_core.RangeTrack(FS, n0, n1, coeffs), _core.GeoConfig(fc_hz=FC, fs_hz=FS, tau0_samples=tau0))
    pipe = _core.Pipeline(*args) if kernel is None else _core.Pipeline(*args, kernel)
    lat = pipe.latency()
    t1 = time.perf_counter()
    y = pipe.process(np.concatenate([x, np.zeros(lat, np.complex64)]))[lat:]
    return y, pipe.kernel, lat, time.perf_counter() - t1


def main():
    main_tle, t0, seg, tau0, frame, x = m1b.build_window()
    print(f"pass '{m1a.PASS_LABEL}', window {geometry.time_to_iso(t0)} + {m1a.WINDOW_S} s; tau0 = {tau0} samples; "
          f"{len(x)} samples (complex64)")
    print(f"cpu_has_avx2() = {_core.cpu_has_avx2()}")

    y_avx2, kernel, lat, t_avx2 = run_pipeline(seg, tau0, x)
    print(f"default kernel: {kernel.name}; latency L = {lat} samples; {t_avx2:.2f} s wall (informational)")
    assert kernel == _core.Kernel.Avx2, "the default kernel must be Avx2"
    y_scalar, _, _, t_scalar = run_pipeline(seg, tau0, x, _core.Kernel.Scalar)
    print(f"scalar kernel on the same input: {t_scalar:.2f} s wall (informational)")

    n0, n1, coeffs = seg.to_arrays()
    truth = {"n0": n0, "n1": n1, "coeffs": coeffs, "fs": FS, "fc": FC, "tone_hz": m1a.TONE_HZ}
    rows, _, _ = m1a.measure(y_avx2, truth)
    err = rows[:, 3] - rows[:, 4]
    print("\nsub | centre (s) | rho_dot (m/s) | measured (Hz) | predicted -(fc+f_tone)rho_dot/c (Hz) | error (Hz)")
    for r, e in zip(rows, err):
        print(f"{int(r[0]):3d} | {r[1]:.3f} | {r[2]:.4f} | {r[3]:.4f} | {r[4]:.4f} | {e:+.4f}")
    dmax = float(np.max(np.abs(err)))
    print(f"\ntone Doppler on the AVX2 output: max |measured - predicted| = {dmax:.4f} Hz over {len(rows)} "
          f"sub-windows (limit {m1a.LIMIT_HZ} Hz)")

    a, s = y_avx2.astype(np.complex128), y_scalar.astype(np.complex128)
    avx2_vs_scalar_db = float(10 * np.log10(np.sum(np.abs(a - s) ** 2) / np.sum(np.abs(s) ** 2)))
    differ = int(np.count_nonzero(y_avx2.view(np.uint64) != y_scalar.view(np.uint64)))
    print(f"AVX2 vs scalar, same input: {avx2_vs_scalar_db:.1f} dB; {differ} of {len(y_avx2)} samples differ in bits")

    RESULTS.parent.mkdir(exist_ok=True)
    RESULTS.write_text(json.dumps({
        "window_start_utc": geometry.time_to_iso(t0), "window_s": m1a.WINDOW_S, "pass_label": m1a.PASS_LABEL,
        "tau0_samples": tau0, "latency_samples": lat, "kernel": kernel.name,
        "cpu_has_avx2": bool(_core.cpu_has_avx2()),
        "tone_rows": [{"sub": int(r[0]), "centre_s": r[1], "rho_dot_mps": r[2], "measured_hz": r[3],
                       "predicted_hz": r[4], "error_hz": float(e)} for r, e in zip(rows, err)],
        "doppler_max_err_hz": dmax, "doppler_limit_hz": m1a.LIMIT_HZ,
        "avx2_vs_scalar_db": avx2_vs_scalar_db, "avx2_vs_scalar_differing_samples": differ,
        "n_samples": len(y_avx2), "wall_s": {"avx2": t_avx2, "scalar": t_scalar},
    }, indent=2), encoding="utf-8")
    print(f"saved {RESULTS}")
    assert dmax < m1a.LIMIT_HZ, f"tone criterion failed: {dmax:.4f} Hz"
    print("PASS")


if __name__ == "__main__":
    main()
