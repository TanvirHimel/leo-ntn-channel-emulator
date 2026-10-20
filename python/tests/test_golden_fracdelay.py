"""M1a acceptance test 3: golden fractional delay against the FFT-domain reference."""

import numpy as np
import pytest

from ntn_emu.golden.fracdelay import fracdelay

N = 8192  # period of the test signal


def bandlimited_spectrum(seed):
    """Random complex spectrum, nonzero only for |f| < 0.4 fs."""
    rng = np.random.default_rng(seed)
    f = np.fft.fftfreq(N)
    return (rng.standard_normal(N) + 1j * rng.standard_normal(N)) * (np.abs(f) < 0.4)


@pytest.mark.parametrize("d", [0.1, 0.25, 0.5, 0.73, 0.999, 5.3, -2.4])
def test_fractional_delay(d):
    """Test 3: constant delay d on periodic band-limited noise; error < -80 dB in |f| < 0.4 fs."""
    X = bandlimited_spectrum(3)
    x = np.fft.ifft(X)
    y = fracdelay(np.tile(x, 3), d)[N:2 * N]  # middle period: kernel never sees the edges
    f = np.fft.fftfreq(N)
    ref = X * np.exp(-2j * np.pi * f * d)
    band = np.abs(f) < 0.4
    err = np.fft.fft(y)[band] - ref[band]
    err_db = 10.0 * np.log10(np.sum(np.abs(err) ** 2) / np.sum(np.abs(ref[band]) ** 2))
    assert err_db < -80.0


def test_chunk_size_invariance():
    """Each output sample is computed independently, so chunking cannot change any bit."""
    rng = np.random.default_rng(4)
    n = 20_000
    x = rng.standard_normal(n) + 1j * rng.standard_normal(n)
    d = 3.0 + np.linspace(0.0, 2.5, n)  # ramp crossing integer delays
    np.testing.assert_array_equal(fracdelay(x, d), fracdelay(x, d, chunk=777))
