"""ntn_emu._core (pybind11): import, input validation, zero-copy input contract, and M1a
acceptance tests 1, 2 and 4 re-run through the C++ scalar path (SPEC §8 M1b tests 1 and 5).

Expected values are the same closed forms as the golden tests (test_golden_geometric.py),
with the C++ float32-output limits: phase < 2e-6 rad, integer shift bit-exact.
"""

import numpy as np
import pytest

from ntn_emu import _core
from test_golden_geometric import (  # closed-form helpers, shared with the golden tests
    C,
    COUNTS,
    FC,
    FS,
    QUAD_COEFFS,
    SHIFT_KS,
    TONE_CASES,
    dtft_peak_hz,
    exact_cycles_frac,
)

HALF = 16  # frac-delay half length (taps x[m0-15 .. m0+16]); latency L = max(0, 15 - floor(min D))


def track(n0, n1, coeffs, fs=FS):
    return _core.RangeTrack(fs, np.asarray(n0, np.int64), np.asarray(n1, np.int64),
                            np.asarray(coeffs, np.float64).reshape(-1, 4))


def pipeline(n0, n1, coeffs, *, fc, tau0=0, delay=True, doppler=True):
    geo = _core.GeoConfig(fc_hz=fc, fs_hz=FS, tau0_samples=tau0, delay=delay, doppler=doppler)
    return _core.Pipeline(track(n0, n1, coeffs), geo)


def run(p, x, block=1 << 20):
    y = np.empty_like(x)
    for s in range(0, len(x), block):
        y[s:s + block] = p.process(x[s:s + block])
    return y


# ---- Import and module surface --------------------------------------------------------------

def test_import_and_surface():
    assert _core.__file__.endswith(".pyd")
    for name in ("RangeTrack", "GeoConfig", "Kernel", "Pipeline", "cpu_has_avx2"):
        assert hasattr(_core, name), name
    assert isinstance(_core.cpu_has_avx2(), bool)
    assert _core.Kernel.Scalar != _core.Kernel.Avx2


def test_range_track_evaluates_closed_form():
    # fs = 1024 and dyadic coefficients: Horner is exact (same case as the gtest).
    t = track([5000], [9096], [[1048576.0, -7000.5, 40.25, -0.5]], fs=1024.0)
    k = 1000
    u = k / 1024.0
    assert t.rho(5000 + k) == 1048576.0 - 7000.5 * u + 40.25 * u**2 - 0.5 * u**3
    assert t.rho_dot(5000 + k) == -7000.5 + 80.5 * u - 1.5 * u**2
    assert t.tau(5000 + k) == t.rho(5000 + k) / C
    assert (t.n_begin, t.n_end, t.fs) == (5000, 9096, 1024.0)


# ---- Invalid segments raise ValueError --------------------------------------------------------

SEG_A = [1048576.0, -7000.5, 40.25, -0.5]  # end at u = 1 s: rho 1041615.25, rho_dot -6921.5


@pytest.mark.parametrize("n0,n1,coeffs", [
    ([0, 1025], [1024, 2048], [SEG_A, [1041615.25, -6921.5, 0, 0]]),         # gap
    ([0, 1000], [1024, 2048], [SEG_A, [1041615.25, -6921.5, 0, 0]]),         # overlap
    ([0, 1024], [1024, 2048], [SEG_A, [1041615.25 + 2e-6, -6921.5, 0, 0]]),  # rho jump
    ([0, 1024], [1024, 2048], [SEG_A, [1041615.25, -6921.5 + 2e-6, 0, 0]]),  # rho_dot jump
    ([0], [0], [SEG_A]),                                                      # empty segment
    ([], [], np.zeros((0, 4))),                                               # no segments
    ([0, 1024], [1024], [SEG_A]),                                             # length mismatch
])
def test_invalid_segments_raise_value_error(n0, n1, coeffs):
    with pytest.raises(ValueError):
        track(n0, n1, coeffs, fs=1024.0)


def test_coeffs_must_have_four_columns():
    with pytest.raises(ValueError):
        _core.RangeTrack(1024.0, np.array([0], np.int64), np.array([1024], np.int64), np.zeros((1, 3)))


def test_c1_within_tolerance_accepted():
    track([0, 1024], [1024, 2048], [SEG_A, [1041615.25 + 0.5e-6, -6921.5, 0, 0]], fs=1024.0)


def test_invalid_config_raises_value_error():
    t = track([0], [1000], [[1e6, 0, 0, 0]])
    with pytest.raises(ValueError):
        _core.Pipeline(t, _core.GeoConfig(fc_hz=FC, fs_hz=FS + 1.0))
    with pytest.raises(ValueError):
        _core.Pipeline(t, _core.GeoConfig(fc_hz=FC, fs_hz=FS), _core.Kernel.Avx2)


# ---- process(): zero-copy complex64 input contract --------------------------------------------

def test_complex128_input_raises_type_error():
    p = pipeline([0], [1000], [[1e6, 0, 0, 0]], fc=FC, delay=False)
    with pytest.raises(TypeError):
        p.process(np.ones(10, np.complex128))
    assert p.sample_index() == 0  # nothing consumed


def test_non_contiguous_input_raises_type_error():
    # A strided view would need a copy; noconvert refuses it instead of copying silently.
    p = pipeline([0], [1000], [[1e6, 0, 0, 0]], fc=FC, delay=False)
    with pytest.raises(TypeError):
        p.process(np.ones(20, np.complex64)[::2])


def test_input_buffer_unchanged_and_read_only_accepted():
    rng = np.random.default_rng(1)
    x = (rng.standard_normal(5000) + 1j * rng.standard_normal(5000)).astype(np.complex64)
    before = x.copy()
    p = pipeline([0], [10_000], [[1_024_452.123, -5899.5, 13.4, -0.02]], fc=FC, tau0=26244)
    y = p.process(x)
    assert x.tobytes() == before.tobytes()
    assert y.dtype == np.complex64 and y.shape == x.shape and not np.shares_memory(x, y)
    x.setflags(write=False)
    p.process(x)  # read-only input is fine: the engine only reads it


def test_two_dimensional_input_raises_value_error():
    p = pipeline([0], [1000], [[1e6, 0, 0, 0]], fc=FC, delay=False)
    with pytest.raises(ValueError):
        p.process(np.ones((2, 5), np.complex64))


def test_past_track_end_raises_index_error():
    p = pipeline([0], [100], [[1e6, 0, 0, 0]], fc=FC, delay=False)
    p.process(np.ones(100, np.complex64))
    with pytest.raises(IndexError):
        p.process(np.ones(1, np.complex64))
    assert p.sample_index() == 100


def test_state_methods():
    p = pipeline([0], [100_000], [[1e6, -7000.0, 0, 0]], fc=FC, tau0=int(1e6 * FS / C) - 3)
    assert p.sample_index() == 0
    p.process(np.zeros(7, np.complex64))
    assert p.sample_index() == 7
    p.reset(500)
    assert p.sample_index() == 500
    assert p.latency() >= 0


# ---- M1a test 1: tone / Doppler phase ---------------------------------------------------------

def tone_output(rho0, v, count):
    p = pipeline([0], [count], [[rho0, v, 0.0, 0.0]], fc=FC, delay=False)
    assert p.latency() == 0
    return run(p, np.ones(count, np.complex64))


@pytest.mark.parametrize("count", COUNTS)
@pytest.mark.parametrize("rho0,v", TONE_CASES)
def test_tone_phase(rho0, v, count):
    """Phase error < 2e-6 rad against the exact (big-int) reference; no drift to the end."""
    y = tone_output(rho0, v, count).astype(np.complex128)
    ref = -2.0 * np.pi * exact_cycles_frac(rho0, v, count)
    err = np.abs(np.angle(y * np.exp(-1j * ref)))
    assert err.max() < 2e-6
    assert err[-100_000:].max() < 2e-6
    assert np.max(np.abs(np.abs(y) - 1.0)) < 1e-6


@pytest.mark.parametrize("rho0,v", TONE_CASES)
def test_tone_fft_peak(rho0, v):
    """The output tone sits at f_d = -fc*v/c within 0.1 Hz (interpolated peak)."""
    y = tone_output(rho0, v, 1_000_000).astype(np.complex128)
    assert abs(dtft_peak_hz(y, FS) - (-FC * v / C)) < 0.1


# ---- M1a test 2: integer shift, bit-exact -----------------------------------------------------

@pytest.mark.parametrize("k", SHIFT_KS)
def test_integer_shift_bit_exact(k):
    """fc = 0, tau - tau0 = k/fs: output equals the input shifted by k, bit for bit."""
    rng = np.random.default_rng(2)
    n_samp = 200_000
    x = (rng.standard_normal(n_samp) + 1j * rng.standard_normal(n_samp)).astype(np.complex64)
    tau0 = 15_360
    p = pipeline([0], [n_samp], [[(tau0 + k) * C / FS, 0.0, 0.0, 0.0]], fc=0.0, tau0=tau0)
    lat = p.latency()
    assert lat == max(0, HALF - k)  # D = k: L = 15 - floor(k - 1e-6) = 16 - k
    y = run(p, np.concatenate([x, np.zeros(lat, np.complex64)]), block=65_536)
    ref = np.zeros(n_samp, np.complex64)
    ref[k:] = x[:n_samp - k]
    assert np.array_equal(y[:lat].view(np.uint64), np.zeros(lat, np.uint64))
    assert np.array_equal(y[lat:].view(np.uint64), ref.view(np.uint64))


# ---- M1a test 4: quadratic phase --------------------------------------------------------------

def test_quadratic_phase():
    """Instantaneous frequency = -fc*rho_dot(t)/c within 0.01 Hz (window 5 s into the segment)."""
    c0, c1, c2, c3 = QUAD_COEFFS
    n_start = 5 * FS
    n_samp = 200_000
    p = pipeline([0], [10 * FS], [QUAD_COEFFS], fc=FC, delay=False)
    p.reset(n_start)
    y = p.process(np.ones(n_samp, np.complex64)).astype(np.complex128)
    phi = np.unwrap(np.angle(y))
    m = 100
    f_meas = (phi[2 * m:] - phi[:-2 * m]) * FS / (2.0 * np.pi * 2 * m)
    u = (n_start + np.arange(m, n_samp - m)) / FS
    f_true = -FC * (c1 + 2.0 * c2 * u + 3.0 * c3 * u**2) / C
    assert np.max(np.abs(f_meas - f_true)) < 0.01
