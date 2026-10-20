"""M1a acceptance test 3: golden fractional delay against the FFT-domain reference.

Measurement functions return the metric; tests assert on it, and scripts/measure_m1a.py prints it.
"""

import numpy as np
import pytest

from ntn_emu.golden.fracdelay import fracdelay

N = 8192  # period of the test signal
# Every delay is >= 0.05 samples from an integer, so the snap path is never used; two are > 1 in magnitude.
DELAYS = [0.1, 0.25, 0.5, 0.73, 0.95, 5.3, -2.4]
NOSNAP_KS = [0, 3, 70_001]


def bandlimited_spectrum(seed):
    """Random complex spectrum, nonzero only for |f| < 0.4 fs."""
    rng = np.random.default_rng(seed)
    f = np.fft.fftfreq(N)
    return (rng.standard_normal(N) + 1j * rng.standard_normal(N)) * (np.abs(f) < 0.4)


def fracdelay_error_db(d):
    """Constant delay d on periodic band-limited noise; error power in |f| < 0.4 fs, dB rel."""
    X = bandlimited_spectrum(3)
    x = np.fft.ifft(X)
    y = fracdelay(np.tile(x, 3), d)[N:2 * N]  # middle period: kernel never sees the edges
    f = np.fft.fftfreq(N)
    ref = X * np.exp(-2j * np.pi * f * d)
    band = np.abs(f) < 0.4
    err = np.fft.fft(y)[band] - ref[band]
    return 10.0 * np.log10(np.sum(np.abs(err) ** 2) / np.sum(np.abs(ref[band]) ** 2))


def nosnap_shift_error_db(k):
    """Integer delay k through the kernel path (snap off) on white noise; error dB rel."""
    rng = np.random.default_rng(5)
    n = 200_000
    x = rng.standard_normal(n) + 1j * rng.standard_normal(n)
    y = fracdelay(x, float(k), snap=False)
    ref = np.zeros_like(x)
    ref[k:] = x[:n - k]
    err_pow = np.sum(np.abs(y - ref) ** 2)
    return -np.inf if err_pow == 0 else 10.0 * np.log10(err_pow / np.sum(np.abs(ref) ** 2))


def test_delays_avoid_snap():
    assert all(abs(d - round(d)) >= 0.05 for d in DELAYS)
    assert any(abs(d) > 1 for d in DELAYS)


@pytest.mark.parametrize("d", DELAYS)
def test_fractional_delay(d):
    """Test 3: error against exp(-j 2 pi f d) < -80 dB over |f| < 0.4 fs."""
    assert fracdelay_error_db(d) < -80.0


@pytest.mark.parametrize("k", NOSNAP_KS)
def test_integer_delay_kernel_path(k):
    """Snap disabled: the kernel itself reproduces an integer shift to < -120 dB."""
    assert nosnap_shift_error_db(k) < -120.0


def test_chunk_and_thread_invariance():
    """Each output sample is computed independently: chunking and threads cannot change any bit."""
    rng = np.random.default_rng(4)
    n = 20_000
    x = rng.standard_normal(n) + 1j * rng.standard_normal(n)
    d = 3.0 + np.linspace(0.0, 2.5, n)  # ramp crossing integer delays
    ref = fracdelay(x, d, workers=1)
    np.testing.assert_array_equal(ref, fracdelay(x, d, chunk=777, workers=1))
    np.testing.assert_array_equal(ref, fracdelay(x, d, chunk=777, workers=2))
    np.testing.assert_array_equal(ref, fracdelay(x, d, chunk=777, workers=4))
