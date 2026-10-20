"""M1a acceptance tests 1, 2 and 4: golden geometric channel against closed-form references."""

import numpy as np
import pytest
from scipy.optimize import minimize_scalar

from ntn_emu.golden import geometric

C = 299_792_458  # m/s, exact (SI), kept as int for the exact reference
FS = 7_680_000
FC = 2_000_000_000


def exact_cycles_frac(rho0, v, count, fs=FS, fc=FC, chunk=1_000_000):
    """fc*(rho0 + v*n/fs)/c reduced to [-0.5, 0.5), computed exactly for n = 0..count-1.

    Every float64 is a dyadic rational: rho0 = a/A, v = b/B. The cycle count is then
    (fc*a*B*fs + fc*b*A*n) / (A*B*fs*c), and its remainder is computed in Python ints.
    Only the final remainder/denominator division rounds (correctly, to float64).
    """
    a, A = float(rho0).as_integer_ratio()
    b, B = float(v).as_integer_ratio()
    den = A * B * fs * C
    k0 = (fc * a * B * fs) % den
    k1 = (fc * b * A) % den
    out = np.empty(count)
    for s in range(0, count, chunk):
        e = min(s + chunk, count)
        n = np.arange(s, e, dtype=np.int64).astype(object)
        r = (k0 + k1 * n) % den
        out[s:e] = (r / den).astype(np.float64)
    out[out >= 0.5] -= 1.0
    return out


def dtft_peak_hz(y, fs):
    """Coarse FFT peak, refined by maximizing |DTFT| within +-1 bin (exact for one complex tone)."""
    n_samp = len(y)
    k = int(np.argmax(np.abs(np.fft.fft(y))))
    f0 = np.fft.fftfreq(n_samp, 1.0 / fs)[k]
    df = fs / n_samp
    n = np.arange(n_samp)

    def neg_mag(f):
        return -np.abs(np.dot(y, np.exp(-2j * np.pi * (f / fs) * n)))

    res = minimize_scalar(neg_mag, bounds=(f0 - df, f0 + df), method="bounded",
                          options={"xatol": 1e-5})
    return res.x


# (rho0 [m], v [m/s]): integer case, and non-integer case exercising Horner rounding.
TONE_CASES = [(1_000_000.0, -7000.0), (1_234_567.891, 6853.37)]
COUNTS = [1_000_000, pytest.param(10_000_000, marks=pytest.mark.slow)]


def linear_track(rho0, v, count):
    return dict(n0=[0], n1=[count], coeffs=[[rho0, v, 0.0, 0.0]])


@pytest.mark.parametrize("count", COUNTS)
@pytest.mark.parametrize("rho0,v", TONE_CASES)
def test_tone_phase(rho0, v, count):
    """Test 1: DC in, linear rho, delay off -> phase = -2*pi*fc*(rho0 + v t)/c within 1e-7 rad."""
    y = geometric.apply(np.ones(count, complex), fs=FS, fc=FC, tau0_samples=0,
                        delay=False, **linear_track(rho0, v, count))
    ref = -2.0 * np.pi * exact_cycles_frac(rho0, v, count)
    err = np.angle(y * np.exp(-1j * ref))
    assert np.max(np.abs(err)) < 1e-7
    assert np.max(np.abs(np.abs(y) - 1.0)) < 1e-12


@pytest.mark.parametrize("count", COUNTS)
@pytest.mark.parametrize("rho0,v", TONE_CASES)
def test_tone_fft_peak(rho0, v, count):
    """Test 1: the output tone sits at f_d = -fc*v/c within 0.1 Hz (interpolated peak)."""
    y = geometric.apply(np.ones(count, complex), fs=FS, fc=FC, tau0_samples=0,
                        delay=False, **linear_track(rho0, v, count))
    f_expected = -FC * v / C
    assert abs(dtft_peak_hz(y, FS) - f_expected) < 0.1


@pytest.mark.parametrize("k", [0, 1, 37, 70_001])
def test_integer_shift(k):
    """Test 2: fc = 0, tau - tau0 = k/fs exactly -> output is input shifted by k (rel err < 1e-12).

    200k samples, so the 64k-sample chunks of the fractional delay are crossed.
    """
    rng = np.random.default_rng(2)
    n_samp = 200_000
    x = rng.standard_normal(n_samp) + 1j * rng.standard_normal(n_samp)
    tau0 = 15_360  # 2 ms at 7.68 Msps
    rho = (tau0 + k) * C / FS
    y = geometric.apply(x, fs=FS, fc=0.0, tau0_samples=tau0,
                        n0=[0], n1=[n_samp], coeffs=[[rho, 0.0, 0.0, 0.0]])
    ref = np.zeros_like(x)
    ref[k:] = x[:n_samp - k]
    assert np.linalg.norm(y - ref) / np.linalg.norm(ref) < 1e-12


def test_quadratic_phase():
    """Test 4: cubic segment -> instantaneous frequency = -fc*rho_dot(t)/c within 0.01 Hz.

    Window starts 5 s into the segment so all coefficients matter. Instantaneous frequency is
    measured by a +-M central difference of the unwrapped output phase: the per-sample
    difference would amplify float64 phase quantization (~1e-9 cycles * fs ~ 0.01 Hz),
    while the central-difference bias fc/c * rho''' * (M/fs)^2 / 6 is ~1e-9 Hz.
    """
    c0, c1, c2, c3 = 1_200_000.0, -5000.0, 40.0, -0.5
    n_start = 5 * FS
    n_samp = 200_000
    y = geometric.apply(np.ones(n_samp, complex), fs=FS, fc=FC, tau0_samples=0,
                        n_start=n_start, delay=False,
                        n0=[0], n1=[10 * FS], coeffs=[[c0, c1, c2, c3]])
    phi = np.unwrap(np.angle(y))
    m = 100
    f_meas = (phi[2 * m:] - phi[:-2 * m]) * FS / (2.0 * np.pi * 2 * m)
    u = (n_start + np.arange(m, n_samp - m)) / FS
    f_true = -FC * (c1 + 2.0 * c2 * u + 3.0 * c3 * u**2) / C
    assert np.max(np.abs(f_meas - f_true)) < 0.01
