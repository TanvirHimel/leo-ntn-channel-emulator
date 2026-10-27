"""M1b end-to-end demo: the M1a window through the C++ engine (ntn_emu._core) and the golden
model, compared sample by sample -> SigMF in out/ -> figure in docs/figures/.

Window and input are those of demo_m1a.py: highest pinned pass, 1 s centred at (AOS+TCA)/2,
light-time Hermite segments (1 s knots), tau0 = floor(tau(0) fs), NR-numerology CP-OFDM (QPSK,
seed 1) plus a test tone at +2.8 MHz (-10 dB). The input is cast to complex64 and the same
array goes to both paths. The C++ output is flushed with latency() zeros and its first L
samples dropped, so index i of both outputs is channel time i.

Error power, C++ vs golden, relative to the golden output power:
- full band: includes the unshaped OFDM leakage at 0.45-0.5 fs (~ -38 dB of the input power),
  outside the C++ fractional-delay design band (|f| < 0.4 fs, SPEC §2); recorded, not asserted;
- in band: both outputs through the same zero-phase low-pass (passband 0.38 fs, stopband
  0.40 fs), the first and last 2000 samples excluded; asserted < -55 dB (SPEC M1b test 3 limit).
The M1a tone-Doppler check (max |measured - predicted| < 1 Hz) is re-run on the C++ output.

    .\\.venv\\Scripts\\python.exe scripts\\demo_m1b.py [--plot-only]
"""

import argparse
import json
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python"))
sys.path.insert(0, str(ROOT / "scripts"))

import numpy as np  # noqa: E402
from scipy import signal  # noqa: E402

import demo_m1a as m1a  # noqa: E402  (window constants and the tone-Doppler measurement)
from ntn_emu import _core, segments, sigmf_io  # noqa: E402
from ntn_emu.constants import C_LIGHT, UE_ALT_M, UE_LAT_DEG, UE_LON_DEG  # noqa: E402
from ntn_emu.golden import geometric  # noqa: E402
from ntn_emu.orbit import geometry, passes, tle  # noqa: E402
from ntn_emu.waveform import ofdm  # noqa: E402

FC, FS, WINDOW_S = m1a.FC, m1a.FS, m1a.WINDOW_S
TONE_HZ, TONE_REL_DB, PASS_LABEL = m1a.TONE_HZ, m1a.TONE_REL_DB, m1a.PASS_LABEL
LIMIT_DB = -55.0
INBAND_LP = signal.firwin(401, 0.39 * FS, window=("kaiser", 12.0), fs=FS)  # pass 0.38, stop 0.40 fs
EDGE = 2000  # samples excluded at each end for the in-band metric (filter transients)
BLOCK = 76_800  # 10 ms, for error vs time

OUT = ROOT / "out"
BASE_CPP = OUT / "m1b_window"
BASE_GOLD = OUT / "m1b_golden"
RESULTS = OUT / "m1b_results.json"
FIG = ROOT / "docs" / "figures" / "m1b_golden_vs_cpp.png"


def build_window():
    """Same window, segments, tau0 and input as demo_m1a.run_channel()."""
    main_tle = tle.load_pinned("starlink")
    sat, ue = main_tle.sat, geometry.ue_position()
    p = passes.load_pinned_passes()[PASS_LABEL]
    mid = geometry.times_from_offsets(p["aos"], (p["tca"] - p["aos"]) * 86_400.0 / 2.0)
    t0 = geometry.times_from_offsets(mid, -WINDOW_S / 2.0)

    def rho_fn(off):
        st = geometry.range_lighttime(sat, ue, geometry.times_from_offsets(t0, off))
        return st.rho, st.rho_dot

    seg = segments.hermite_segments(rho_fn, WINDOW_S, FS, knot_s=1.0)
    segments.check_c1(seg)
    tau0 = int(np.floor(seg.rho(np.array([0]))[0] / C_LIGHT * FS))
    frame = ofdm.generate(int(round(WINDOW_S * 1e3)), ofdm.OfdmConfig(modulation="qpsk", seed=1))
    n = np.arange(len(frame.samples))
    tone = np.sqrt(10 ** (TONE_REL_DB / 10)) * np.exp(2j * np.pi * (TONE_HZ / FS) * n)
    x = (frame.samples + tone).astype(np.complex64)
    return main_tle, t0, seg, tau0, frame, x


def run_channels():
    main_tle, t0, seg, tau0, frame, x = build_window()
    n0, n1, coeffs = seg.to_arrays()
    print(f"pass '{PASS_LABEL}', window {geometry.time_to_iso(t0)} + {WINDOW_S} s; "
          f"tau0 = {tau0} samples; {len(x)} samples (complex64 input to both paths)")

    t1 = time.perf_counter()
    pipe = _core.Pipeline(_core.RangeTrack(FS, n0, n1, coeffs),
                          _core.GeoConfig(fc_hz=FC, fs_hz=FS, tau0_samples=tau0),
                          _core.Kernel.Scalar)  # M1b record: scalar engine (default is Avx2 since M2)
    lat = pipe.latency()
    y_cpp = pipe.process(np.concatenate([x, np.zeros(lat, np.complex64)]))[lat:]
    t_cpp = time.perf_counter() - t1
    print(f"C++ scalar (_core): latency L = {lat} samples; {t_cpp:.2f} s wall (informational)")

    t1 = time.perf_counter()
    y_gold = geometric.apply(x, fs=FS, fc=FC, n0=n0, n1=n1, coeffs=coeffs, tau0_samples=tau0)
    t_gold = time.perf_counter() - t1
    print(f"golden (float64): {t_gold:.1f} s wall (informational)")

    ntn = {"tau0_s": tau0 / FS, "tau0_samples": tau0, "tle_name": main_tle.name,
           "tle_norad": main_tle.norad, "tle_epoch": geometry.time_to_iso(main_tle.sat.epoch),
           "ue_lla": [UE_LAT_DEG, UE_LON_DEG, UE_ALT_M], "pass_label": PASS_LABEL,
           "tdl_model": "none", "delay_spread_s": None, "snr_db": None,
           "range_mode": "light-time", "test_tone_hz": TONE_HZ, "test_tone_rel_db": TONE_REL_DB,
           "waveform": {"type": "NR-numerology CP-OFDM (not NR-compliant)", "scs_hz": ofdm.SCS,
                        "n_fft": ofdm.N_FFT, "n_sc": ofdm.N_SC,
                        "modulation": frame.cfg.modulation, "n_id2": frame.cfg.n_id2,
                        "seed": frame.cfg.seed, "pilot_seed": frame.cfg.pilot_seed}}
    sigmf_io.write(BASE_CPP, y_cpp, fs=FS, fc=FC, datetime_utc=t0.utc_datetime(),
                   ntn={**ntn, "engine": "C++ scalar (ntn_emu._core)", "latency_samples": lat,
                        "alignment": "flushed with latency zeros, first L samples dropped"})
    sigmf_io.write(BASE_GOLD, y_gold.astype(np.complex64), fs=FS, fc=FC,
                   datetime_utc=t0.utc_datetime(), ntn={**ntn, "engine": "golden (float64)"})
    truth = {"n0": n0, "n1": n1, "coeffs": coeffs, "fs": FS, "fc": FC, "tone_hz": TONE_HZ}
    return y_cpp, y_gold, truth, {"latency_samples": lat, "tau0_samples": tau0,
                                  "cpp_wall_s": t_cpp, "golden_wall_s": t_gold}


def db(v):
    return 10.0 * np.log10(v)


def compare(y_cpp, y_gold):
    a, b = y_cpp.astype(np.complex128), y_gold.astype(np.complex128)
    e = a - b
    full = db(np.sum(np.abs(e) ** 2) / np.sum(np.abs(b) ** 2))
    af = np.convolve(a, INBAND_LP, mode="same")[EDGE:-EDGE]
    bf = np.convolve(b, INBAND_LP, mode="same")[EDGE:-EDGE]
    inband = db(np.sum(np.abs(af - bf) ** 2) / np.sum(np.abs(bf) ** 2))
    nb = len(b) // BLOCK
    per_ms_full = [db(np.sum(np.abs(e[i * BLOCK:(i + 1) * BLOCK]) ** 2)
                      / np.sum(np.abs(b[i * BLOCK:(i + 1) * BLOCK]) ** 2)) for i in range(nb)]
    ef = np.convolve(e, INBAND_LP, mode="same")
    bf_full = np.convolve(b, INBAND_LP, mode="same")
    per_ms_in = [db(np.sum(np.abs(ef[i * BLOCK:(i + 1) * BLOCK]) ** 2)
                    / np.sum(np.abs(bf_full[i * BLOCK:(i + 1) * BLOCK]) ** 2))
                 for i in range(1, nb - 1)]  # first/last block hold the filter transients
    return full, inband, np.array(per_ms_full), np.array(per_ms_in)


def plot(y_cpp, y_gold, per_ms_full, per_ms_in, full, inband):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    # Reference palette slots 1-2 (blue, orange); text in ink tones, never series colour.
    ink, ink2, grid_c, surf = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"
    c1, c2 = "#2a78d6", "#eb6834"
    plt.rcParams.update({"font.size": 9, "axes.edgecolor": ink2, "axes.labelcolor": ink,
                         "xtick.color": ink2, "ytick.color": ink2, "text.color": ink,
                         "axes.facecolor": surf, "figure.facecolor": surf})
    fig, ax = plt.subplots(2, 1, figsize=(7.5, 7.0), constrained_layout=True)

    kw = {"fs": FS, "window": ("kaiser", 14.0), "nperseg": 8192, "return_onesided": False,
          "scaling": "density"}
    f, pg = signal.welch(y_gold.astype(np.complex128), **kw)
    _, pc = signal.welch(y_cpp.astype(np.complex128), **kw)
    _, pe = signal.welch((y_cpp.astype(np.complex128) - y_gold), **kw)
    order = np.argsort(f)
    f, pg, pc, pe = f[order] / 1e6, pg[order], pc[order], pe[order]
    ref = pg.max()
    ax[0].axvspan(0.4 * FS / 1e6, 0.5 * FS / 1e6, color=grid_c, lw=0)
    ax[0].axvspan(-0.5 * FS / 1e6, -0.4 * FS / 1e6, color=grid_c, lw=0)
    ax[0].plot(f, db(pg / ref), color=c1, lw=2.0, label="golden (float64)")
    ax[0].plot(f, db(pc / ref), color=c2, lw=1.2, ls=(0, (4, 3)), label="C++ scalar (_core)")
    ax[0].plot(f, db(pe / ref), color=ink2, lw=1.0, label="difference C++ - golden")
    ax[0].text(0.5 * FS / 1e6 - 0.06, 2, "shaded: outside C++ design band, |f| > 0.4 fs",
               ha="right", va="top", color=ink2, fontsize=8)
    ax[0].set(title="Output spectrum, M1a window (OFDM + 2.8 MHz tone), Welch, Kaiser(14)",
              xlabel="Baseband frequency (MHz)", ylabel="PSD (dB rel. golden peak)",
              xlim=(-FS / 2e6, FS / 2e6), ylim=(-180, 5))
    ax[0].legend(loc="lower left", frameon=False)

    t_full = (np.arange(len(per_ms_full)) + 0.5) * BLOCK / FS
    t_in = (np.arange(1, len(per_ms_in) + 1) + 0.5) * BLOCK / FS
    ax[1].axhline(LIMIT_DB, color=ink2, lw=1, ls="--",
                  label=f"limit {LIMIT_DB:g} dB (asserted for |f| < 0.38 fs)")
    ax[1].plot(t_full, per_ms_full, color=c2, lw=2.0, label=f"full band (total {full:.1f} dB)")
    ax[1].plot(t_in, per_ms_in, color=c1, lw=2.0, label=f"|f| < 0.38 fs (total {inband:.1f} dB)")
    ax[1].set(title="C++ vs golden error power per 10 ms, relative to the golden output",
              xlabel="Time in window (s)", ylabel="Error power (dB)")
    ax[1].legend(loc="center right", frameon=False)
    for a in ax:
        a.grid(True, color=grid_c, lw=0.6)
        a.set_axisbelow(True)
        for side in ("top", "right"):
            a.spines[side].set_visible(False)
    FIG.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG, dpi=110)
    print(f"saved {FIG}")


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--plot-only", action="store_true", help="reuse out/ (skip both channels)")
    args = ap.parse_args()
    if args.plot_only:
        res = json.loads(RESULTS.read_text(encoding="utf-8"))
        y_cpp, _ = sigmf_io.read(BASE_CPP)
        y_gold, _ = sigmf_io.read(BASE_GOLD)
        seg_arrays = np.load(OUT / "m1b_truth.npz")
        truth = {k: seg_arrays[k] for k in seg_arrays.files}
        info = res["run"]
    else:
        y_cpp, y_gold, truth, info = run_channels()
        np.savez(OUT / "m1b_truth.npz", **truth)
        y_gold = y_gold.astype(np.complex64)  # compare what is stored, so --plot-only matches

    full, inband, per_ms_full, per_ms_in = compare(y_cpp, y_gold)
    print(f"\nerror power C++ vs golden, full band: {full:.2f} dB (recorded; includes the "
          f"unshaped OFDM leakage above 0.4 fs)")
    print(f"error power C++ vs golden, |f| < 0.38 fs: {inband:.2f} dB (limit {LIMIT_DB} dB)")
    print(f"per-10-ms error, full band: min {per_ms_full.min():.1f}, median "
          f"{np.median(per_ms_full):.1f}, max {per_ms_full.max():.1f} dB")
    print(f"per-10-ms error, |f| < 0.38 fs: min {per_ms_in.min():.1f}, median "
          f"{np.median(per_ms_in):.1f}, max {per_ms_in.max():.1f} dB")

    rows, _, _ = m1a.measure(y_cpp, truth)
    dop = np.abs(rows[:, 3] - rows[:, 4])
    print(f"tone Doppler on the C++ output: max |measured - predicted| = {dop.max():.4f} Hz over "
          f"{len(rows)} sub-windows (limit {m1a.LIMIT_HZ} Hz)")

    RESULTS.write_text(json.dumps({
        "error_full_db": full, "error_inband_db": inband, "inband_def": "|f| < 0.38 fs "
        "(zero-phase Kaiser low-pass, 401 taps, beta 12, cutoff 0.39 fs), 2000 edge samples "
        "excluded", "limit_db": LIMIT_DB,
        "per_10ms_full_db": {"min": per_ms_full.min(), "median": float(np.median(per_ms_full)),
                           "max": per_ms_full.max()},
        "per_10ms_inband_db": {"min": per_ms_in.min(), "median": float(np.median(per_ms_in)),
                             "max": per_ms_in.max()},
        "doppler_max_err_hz": float(dop.max()), "doppler_limit_hz": m1a.LIMIT_HZ,
        "run": info}, indent=2), encoding="utf-8")
    print(f"saved {RESULTS}")
    plot(y_cpp, y_gold, per_ms_full, per_ms_in, full, inband)
    assert inband < LIMIT_DB, f"in-band error {inband:.2f} dB"
    assert dop.max() < m1a.LIMIT_HZ, f"Doppler check failed: {dop.max():.4f} Hz"
    print("PASS")


if __name__ == "__main__":
    main()
