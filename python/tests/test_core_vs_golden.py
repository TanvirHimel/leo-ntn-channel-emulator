"""M1b acceptance test 3 (C++ scalar engine) and M2 acceptance test 4 (AVX2 engine) vs the
float64 golden model on a real pass. Both kernels run the same comparison and limit.

Range track: highest pinned Starlink pass, light-time Hermite segments with 1 s knots over 2 s,
the knot placed at (AOS + TCA)/2. The window is centred on the knot, so it crosses it.
Input: random QPSK CP-OFDM (waveform/ofdm.py), band-limited to |f| < 0.4 fs, cast to complex64
and given to both paths. The unshaped generator leaks ~ -38 dB of its power into 0.4-0.5 fs
(rectangular symbol edges). That band is outside the C++ fractional-delay design band
(|f| < 0.4 fs, SPEC §2), and near Nyquist even the golden kernel is not a reference. Measured on
the unshaped input: -50.8 dB full band, from 0.45-0.5 fs only; -106 dB below 0.38 fs. So the
stimulus is shaped like a real transmitter's would be: a fixed zero-phase Kaiser low-pass
(301 taps, cutoff 0.35 fs, beta 10: stopband from ~0.36 fs). The -55 dB limit is unchanged.
Both paths use the same tau0 = floor(tau(n_start) * fs). Golden output index i is channel time
n_start + i; C++ output index i is channel time n_start + i - L. So C++ out[L:] is compared with
golden y[:len - L]. Error power relative to the golden output must be < -55 dB.

Fast: 0.1 s window (golden model ~10 s). Slow: 1 s window.
"""

import numpy as np
import pytest
from scipy import signal

from ntn_emu import _core, segments
from ntn_emu.golden import geometric
from ntn_emu.orbit import geometry
from ntn_emu.waveform import ofdm

C = 299_792_458.0
FS = 7_680_000
FC = 2.0e9
KNOT_S = 1.0
PASS_LABEL = "high"
LIMIT_DB = -55.0
WINDOWS = [0.1, pytest.param(1.0, marks=pytest.mark.slow)]
LOWPASS = signal.firwin(301, 0.35 * FS, window=("kaiser", 10.0), fs=FS)  # stimulus shaping


def stimulus(window_s, seed=3):
    """Random QPSK OFDM, band-limited to |f| < 0.4 fs (zero-phase FIR), as complex64."""
    n_slots = int(round(window_s * 1e3))
    x = ofdm.generate(n_slots, ofdm.OfdmConfig(modulation="qpsk", seed=seed)).samples
    return np.convolve(x, LOWPASS, mode="same").astype(np.complex64)


def power_share_above_db(x, f_rel=0.4):
    """Share of x's power at |f| >= f_rel * fs, in dB (Welch, Kaiser(14) window: a plain FFT of
    a non-periodic block leaks ~ -70 dB by itself and would measure its own edges, not x)."""
    f, p = signal.welch(x.astype(np.complex128), fs=1.0, window=("kaiser", 14.0), nperseg=4096,
                        return_onesided=False, scaling="spectrum")
    return 10.0 * np.log10(p[np.abs(f) >= f_rel].sum() / p.sum())


def knot_segments(sat, ue, p):
    """Two 1 s light-time Hermite segments with the knot (global sample FS) at (AOS+TCA)/2."""
    mid = geometry.times_from_offsets(p["aos"], (p["tca"] - p["aos"]) * 86_400.0 / 2.0)
    t_start = geometry.times_from_offsets(mid, -KNOT_S)

    def rho_fn(off):
        st = geometry.range_lighttime(sat, ue, geometry.times_from_offsets(t_start, off))
        return st.rho, st.rho_dot

    seg = segments.hermite_segments(rho_fn, 2.0 * KNOT_S, FS, knot_s=KNOT_S)
    segments.check_c1(seg)
    return seg


# Golden output per (segment set, window, seed): the float64 model is the slow part (~10 s for
# 0.1 s), and both kernels are compared against the same reference.
_GOLDEN_CACHE = {}


def golden_output(seg, x, n_start, tau0, window_s, seed):
    key = (id(seg), window_s, seed, n_start, tau0)
    if key not in _GOLDEN_CACHE:
        n0, n1, coeffs = seg.to_arrays()
        _GOLDEN_CACHE[key] = (seg, geometric.apply(x, fs=FS, fc=FC, n0=n0, n1=n1, coeffs=coeffs,
                                                   tau0_samples=tau0, n_start=n_start))
    return _GOLDEN_CACHE[key][1]


def core_vs_golden(seg, window_s, seed=3, shift=0, kernel=_core.Kernel.Scalar):
    """Returns (error dB, latency L, tau0, n_start). `shift` misaligns the comparison by that
    many samples (negative control only; 0 for the acceptance measurement). The kernel defaults
    to Scalar (scripts/measure_m1b.py imports this for the M1b record)."""
    n_win = int(round(window_s * FS))
    n_start = FS - n_win // 2  # window centred on the knot
    x = stimulus(window_s, seed)
    assert len(x) == n_win
    tau0 = int(np.floor(seg.rho(np.array([n_start]))[0] / C * FS))
    n0, n1, coeffs = seg.to_arrays()

    y_gold = golden_output(seg, x, n_start, tau0, window_s, seed)
    p = _core.Pipeline(_core.RangeTrack(FS, n0, n1, coeffs),
                       _core.GeoConfig(fc_hz=FC, fs_hz=FS, tau0_samples=tau0), kernel)
    p.reset(n_start)
    lat = p.latency()
    y_cpp = p.process(x).astype(np.complex128)

    a = y_cpp[lat + shift:]
    b = y_gold[:len(y_gold) - lat - shift]
    err_db = 10.0 * np.log10(np.sum(np.abs(a - b) ** 2) / np.sum(np.abs(b) ** 2))
    return err_db, lat, tau0, n_start


@pytest.fixture(scope="module")
def seg(starlink, ue, pinned_passes):
    return knot_segments(starlink.sat, ue, pinned_passes[PASS_LABEL])


@pytest.mark.parametrize("kernel", [_core.Kernel.Scalar, _core.Kernel.Avx2], ids=["Scalar", "Avx2"])
@pytest.mark.parametrize("window_s", WINDOWS)
def test_core_vs_golden_error_power(seg, window_s, kernel):
    err_db, lat, _, n_start = core_vs_golden(seg, window_s, kernel=kernel)
    print(f"[ measure ] core vs golden, {kernel.name}, {window_s} s: {err_db:.2f} dB (limit {LIMIT_DB} dB)")
    n_win = int(round(window_s * FS))
    assert n_start < FS < n_start + n_win - lat  # the compared span crosses the knot
    assert err_db < LIMIT_DB


def test_stimulus_is_band_limited():
    assert power_share_above_db(stimulus(0.02)) < -80.0  # 0.02 s: 37 Welch segments


def test_misaligned_comparison_fails(seg):
    """Negative control: the same comparison off by one sample is far above the limit, so the
    -55 dB criterion discriminates (it is not met by any two OFDM-like signals)."""
    err_db, *_ = core_vs_golden(seg, 0.01, shift=1)
    assert err_db > -20.0
