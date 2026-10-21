"""M1a end-to-end demo: pinned TLE + rounded UE -> pass -> 1 s window -> OFDM + tone ->
golden geometric channel -> SigMF in out/ -> figure in docs/figures/.

Window: highest pinned pass, 1 s centred at mid = (AOS + TCA)/2 (large Doppler).
Range track: light-time Hermite segments (1 s knots); tau0 = floor(tau(0) * fs).
Input: NR-numerology CP-OFDM (QPSK) plus a test tone at +2.8 MHz, -10 dB relative to the
OFDM power. The tone stays inside 0.4 fs and clear of the OFDM band edge (+-2.25 MHz) for
Doppler up to +-46 kHz.

Doppler check: the tone is mixed to DC, low-pass filtered and decimated to 240 kHz, then its
frequency is estimated in ten 0.1 s sub-windows (FFT peak refined by |DTFT| maximization).
Each estimate is compared with the prediction at THAT sub-window's centre sample n_i:
    tone offset = -(fc + f_tone) * rho_dot(n_i) / c
The tone travels at RF frequency fc + f_tone. Its baseband frequency is Doppler-shifted by
the relative delay too (f_tone * rho_dot / c, ~50 Hz here), so -fc*rho_dot/c alone is the
carrier Doppler, not the tone's. Both are printed; the assertion uses the exact form.
Assert: max |measured - predicted| < 1 Hz over all ten.

The channel output is saved before any plotting. --plot-only reloads out/ and skips the channel.

    .\\.venv\\Scripts\\python.exe scripts\\demo_m1a.py [--plot-only]
"""

import argparse
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python"))

import numpy as np  # noqa: E402
from scipy import signal  # noqa: E402
from scipy.optimize import minimize_scalar  # noqa: E402

from ntn_emu import segments, sigmf_io  # noqa: E402
from ntn_emu.constants import C_LIGHT, UE_ALT_M, UE_LAT_DEG, UE_LON_DEG  # noqa: E402
from ntn_emu.golden import geometric  # noqa: E402
from ntn_emu.orbit import geometry, passes, tle  # noqa: E402
from ntn_emu.waveform import ofdm  # noqa: E402

FC = 2.0e9
FS = ofdm.FS
WINDOW_S = 1.0
TONE_HZ = 2.8e6
TONE_REL_DB = -10.0
PASS_LABEL = "high"
DECIM = 32               # 7.68 Msps -> 240 kHz
N_SUB = 10               # ten 0.1 s sub-windows
LIMIT_HZ = 1.0

OUT = ROOT / "out"
BASE = OUT / "m1a_window"
TRUTH = OUT / "m1a_truth.npz"
FIG = ROOT / "docs" / "figures" / "m1a_doppler.png"


def run_channel():
    t_start = time.perf_counter()
    main_tle = tle.load_pinned("starlink")
    sat, ue = main_tle.sat, geometry.ue_position()
    p = passes.load_pinned_passes()[PASS_LABEL]
    mid = geometry.times_from_offsets(p["aos"], (p["tca"] - p["aos"]) * 86_400.0 / 2.0)
    t0 = geometry.times_from_offsets(mid, -WINDOW_S / 2.0)
    print(f"pass '{PASS_LABEL}': max el {p['max_el_deg']:.2f} deg, AOS {p['aos_utc']}, "
          f"TCA {p['tca_utc']}")
    print(f"window: {geometry.time_to_iso(t0)} + {WINDOW_S} s (centred at (AOS+TCA)/2)")

    def rho_fn(off):
        st = geometry.range_lighttime(sat, ue, geometry.times_from_offsets(t0, off))
        return st.rho, st.rho_dot

    seg = segments.hermite_segments(rho_fn, WINDOW_S, FS, knot_s=1.0)
    segments.check_c1(seg)
    n0, n1, coeffs = seg.to_arrays()
    tau0 = int(np.floor(seg.rho(np.array([0]))[0] / C_LIGHT * FS))
    print(f"segments: {len(n0)}; rho(0) = {seg.rho(np.array([0]))[0] / 1e3:.3f} km, "
          f"tau0 = {tau0} samples ({tau0 / FS * 1e3:.4f} ms)")

    n_slots = int(round(WINDOW_S * 1e3))
    frame = ofdm.generate(n_slots, ofdm.OfdmConfig(modulation="qpsk", seed=1))
    n = np.arange(len(frame.samples))
    tone = np.sqrt(10 ** (TONE_REL_DB / 10)) * np.exp(2j * np.pi * (TONE_HZ / FS) * n)
    x = frame.samples + tone
    t1 = time.perf_counter()
    y = geometric.apply(x, fs=FS, fc=FC, n0=n0, n1=n1, coeffs=coeffs, tau0_samples=tau0)
    t_chan = time.perf_counter() - t1
    print(f"golden channel: {len(x)} samples in {t_chan:.1f} s ({len(x) / t_chan / 1e3:.0f} ksps)")

    annotations = [{"core:sample_start": int(frame.symbol_starts[s]),
                    "core:sample_count": int(frame.cp_lengths[s] + ofdm.N_FFT),
                    "core:label": f"PSS N_ID2={frame.cfg.n_id2} (generator position, pre-channel)"}
                   for s in frame.pss_symbols]
    sigmf_io.write(BASE, y, fs=FS, fc=FC, datetime_utc=t0.utc_datetime(), annotations=annotations,
                   ntn={"tau0_s": tau0 / FS, "tau0_samples": tau0,
                        "tle_name": main_tle.name, "tle_norad": main_tle.norad,
                        "tle_epoch": geometry.time_to_iso(sat.epoch),
                        "ue_lla": [UE_LAT_DEG, UE_LON_DEG, UE_ALT_M],
                        "pass_label": PASS_LABEL, "tdl_model": "none", "delay_spread_s": None,
                        "snr_db": None, "range_mode": "light-time",
                        "test_tone_hz": TONE_HZ, "test_tone_rel_db": TONE_REL_DB,
                        "waveform": {"type": "NR-numerology CP-OFDM (not NR-compliant)",
                                     "scs_hz": ofdm.SCS, "n_fft": ofdm.N_FFT,
                                     "n_sc": ofdm.N_SC, "modulation": frame.cfg.modulation,
                                     "n_id2": frame.cfg.n_id2, "seed": frame.cfg.seed,
                                     "pilot_seed": frame.cfg.pilot_seed}})
    np.savez(TRUTH, n0=n0, n1=n1, coeffs=coeffs, fs=FS, fc=FC, tone_hz=TONE_HZ, tau0=tau0,
             window_start_utc=geometry.time_to_iso(t0), pass_label=PASS_LABEL)
    print(f"saved {BASE}.sigmf-data/-meta and {TRUTH} ({time.perf_counter() - t_start:.1f} s total)")


def dtft_peak(z, fs_z, f_lo, f_hi):
    n = np.arange(len(z))

    def neg(f):
        return -np.abs(np.dot(z, np.exp(-2j * np.pi * (f / fs_z) * n)))

    return minimize_scalar(neg, bounds=(f_lo, f_hi), method="bounded",
                           options={"xatol": 1e-5}).x


def measure(y, truth):
    fs, fc, f_tone = float(truth["fs"]), float(truth["fc"]), float(truth["tone_hz"])
    seg = segments.RangeSegments(fs, truth["n0"], truth["n1"], truth["coeffs"])
    n = np.arange(len(y))
    mixed = y.astype(np.complex128) * np.exp(-2j * np.pi * (f_tone / fs) * n)
    h = signal.firwin(1601, 100e3, fs=fs, window=("kaiser", 10.0))
    z = signal.oaconvolve(mixed, h, mode="same")[::DECIM]
    fs_z = fs / DECIM
    sub = len(z) // N_SUB
    rows = []
    for i in range(N_SUB):
        zi = z[i * sub:(i + 1) * sub]
        nfft = 1 << 20
        spec = np.abs(np.fft.fft(zi, nfft))
        f0 = np.fft.fftfreq(nfft, 1.0 / fs_z)[np.argmax(spec)]
        f_meas = dtft_peak(zi, fs_z, f0 - 2.0, f0 + 2.0)
        n_c = int(round((i + 0.5) * sub * DECIM))  # sub-window centre, in input samples
        rate = seg.rho_dot(np.array([n_c]))[0]
        pred_tone = -(fc + f_tone) * rate / C_LIGHT
        pred_carrier = -fc * rate / C_LIGHT
        rows.append((i, n_c / fs, rate, f_meas, pred_tone, pred_carrier))
    return np.array(rows), z, fs_z


def plot(y, truth, rows, z, fs_z):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    ink, ink2, grid_c, surf = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"
    c_pred, c_meas = "#2a78d6", "#eb6834"
    plt.rcParams.update({"font.size": 9, "axes.edgecolor": ink2, "axes.labelcolor": ink,
                         "xtick.color": ink2, "ytick.color": ink2, "text.color": ink,
                         "axes.facecolor": surf, "figure.facecolor": surf})
    fs = float(truth["fs"])
    fig, ax = plt.subplots(3, 1, figsize=(7.5, 8.2), constrained_layout=True)

    f, t, s = signal.spectrogram(y, fs=fs, nperseg=2048, noverlap=0, return_onesided=False,
                                 scaling="spectrum")
    f, s = np.fft.fftshift(f), np.fft.fftshift(s, axes=0)
    db = 10 * np.log10(s + 1e-20)
    im = ax[0].pcolormesh(t, f / 1e6, db, shading="auto", cmap="Blues",
                          vmin=db.max() - 70, vmax=db.max(), rasterized=True)
    ax[0].set(title="Channel output, full band (golden model, light-time range)",
              xlabel="Time in window (s)", ylabel="Baseband frequency (MHz)")
    ax[0].annotate("test tone (+2.8 MHz)", xy=(0.5, 2.83), xytext=(0.62, 3.35), color=ink,
                   arrowprops={"arrowstyle": "-", "color": ink2, "lw": 0.8})
    fig.colorbar(im, ax=ax[0], label="Power (dB, rel.)", pad=0.01)

    fz, tz, sz = signal.spectrogram(z, fs=fs_z, nperseg=16384, noverlap=15360,
                                    return_onesided=False, scaling="spectrum")
    fz, sz = np.fft.fftshift(fz), np.fft.fftshift(sz, axes=0)
    dbz = 10 * np.log10(sz + 1e-20)
    centre = np.mean(rows[:, 4])
    sel = np.abs(fz - centre) < 300.0
    ax[1].pcolormesh(tz, (fz[sel] - centre), dbz[sel], shading="auto", cmap="Blues",
                     vmin=dbz.max() - 60, vmax=dbz.max(), rasterized=True)
    ax[1].plot(rows[:, 1], rows[:, 4] - centre, color=c_pred, lw=2, label="predicted")
    ax[1].plot(rows[:, 1], rows[:, 3] - centre, "o", ms=7, mfc=c_meas, mec=surf, mew=1.5,
               label="measured (0.1 s sub-windows)")
    ax[1].set(title=f"Test tone Doppler, offset from {centre / 1e3:.3f} kHz",
              xlabel="Time in window (s)", ylabel="Frequency offset (Hz)")
    ax[1].legend(loc="upper right", frameon=False)

    res = rows[:, 3] - rows[:, 4]
    ax[2].axhline(0.0, color=grid_c, lw=1)
    for lim in (-LIMIT_HZ, LIMIT_HZ):
        ax[2].axhline(lim * 1e3, color=ink2, lw=1, ls="--")
    ax[2].plot(rows[:, 1], res * 1e3, "o", ms=7, mfc=c_meas, mec=surf, mew=1.5)
    ax[2].set(title=f"Measured minus predicted (max |error| {np.max(np.abs(res)) * 1e3:.2f} mHz; "
                    f"limit +-{LIMIT_HZ:g} Hz, dashed)",
              xlabel="Sub-window centre time (s)", ylabel="Error (mHz)",
              ylim=(-LIMIT_HZ * 1e3 * 1.15, LIMIT_HZ * 1e3 * 1.15))
    for a in ax:
        for side in ("top", "right"):
            a.spines[side].set_visible(False)
    FIG.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG, dpi=110)
    print(f"saved {FIG}")


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--plot-only", action="store_true", help="reuse out/ (skip the channel)")
    args = ap.parse_args()
    if not args.plot_only:
        run_channel()
    y, meta = sigmf_io.read(BASE)
    truth = np.load(TRUTH)
    rows, z, fs_z = measure(y, truth)

    print("\nsub | t_centre (s) | rho_dot (m/s) | measured tone offset (Hz) | predicted "
          "-(fc+f_tone)rho_dot/c (Hz) | error (Hz) | -fc*rho_dot/c (Hz) | measured - that (Hz)")
    for i, tc, rate, fm, pt, pc in rows:
        print(f"{int(i):3d} | {tc:12.3f} | {rate:13.4f} | {fm:25.4f} | {pt:37.4f} | "
              f"{fm - pt:+10.4f} | {pc:18.4f} | {fm - pc:+19.4f}")
    err = np.abs(rows[:, 3] - rows[:, 4])
    print(f"\nmax |measured - predicted| = {err.max():.4f} Hz over {N_SUB} sub-windows "
          f"(limit {LIMIT_HZ} Hz)")
    plot(y, truth, rows, z, fs_z)
    assert err.max() < LIMIT_HZ, f"Doppler check failed: {err.max():.4f} Hz"
    print("PASS")


if __name__ == "__main__":
    main()
