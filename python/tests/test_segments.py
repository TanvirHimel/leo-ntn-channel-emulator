"""Range segments: acceptance test 5 (Hermite fit) plus closed-form checks.

Test 5 fast variant: highest pinned pass, +-60 s around TCA. Slow: all three passes, AOS -> LOS.
"""

import numpy as np
import pytest

from conftest import PASS_LABELS
from ntn_emu import segments
from ntn_emu.orbit import geometry

FS = 7_680_000.0
DENSE_STEP_SAMPLES = 76_800  # 10 ms at FS, exact
FAST_HALF_WINDOW_S = 60.0


def lighttime_fn(sat, ue, t0):
    def fn(offsets_s):
        st = geometry.range_lighttime(sat, ue, geometry.times_from_offsets(t0, offsets_s))
        return st.rho, st.rho_dot
    return fn


def fit_window(p, half_window_s):
    """(t0, duration_s): +-half_window_s around TCA, or AOS -> LOS (whole seconds) if None."""
    if half_window_s is None:
        return p["aos"], float(np.floor((p["los"] - p["aos"]) * 86_400.0))
    return geometry.times_from_offsets(p["tca"], -half_window_s), 2.0 * half_window_s


def segment_fit_errors(sat, ue, t0, duration_s):
    """Test 5: Hermite (1 s knots) vs dense 10 ms light-time range on [0, duration).

    Returns (max |d rho| [m], max |d rho_dot| [m/s]); check_c1 must pass first.
    """
    fn = lighttime_fn(sat, ue, t0)
    seg = segments.hermite_segments(fn, duration_s, FS, knot_s=1.0)
    segments.check_c1(seg)
    n = np.arange(0, int(duration_s * FS), DENSE_STEP_SAMPLES, dtype=np.int64)
    ref_rho, ref_rate = fn(n / FS)
    return np.max(np.abs(seg.rho(n) - ref_rho)), np.max(np.abs(seg.rho_dot(n) - ref_rate))


FIT_CASES = [pytest.param("high", FAST_HALF_WINDOW_S, id="high-tca+-60s")] + [
    pytest.param(label, None, marks=pytest.mark.slow, id=f"{label}-full") for label in PASS_LABELS]


@pytest.mark.parametrize("label,half_window_s", FIT_CASES)
def test_segment_fit(label, half_window_s, starlink, ue, pinned_passes):
    """Test 5: max |d rho| < 1 mm, max |d rho_dot| < 1 mm/s, C1 at every knot."""
    t0, duration = fit_window(pinned_passes[label], half_window_s)
    d_rho, d_rate = segment_fit_errors(starlink.sat, ue, t0, duration)
    assert d_rho < 1e-3
    assert d_rate < 1e-3


def test_hermite_reproduces_cubic():
    """Cubic Hermite is exact for a cubic: rho = a + b t + c t^2 + d t^3."""
    a, b, c, d = 1.1e6, -6500.0, 35.0, -0.2

    def fn(t):
        return a + b * t + c * t**2 + d * t**3, b + 2 * c * t + 3 * d * t**2

    seg = segments.hermite_segments(fn, 20.0, FS, knot_s=1.0)
    segments.check_c1(seg)
    n = np.arange(0, 20 * 7_680_000, 99_991)
    rho, rate = fn(n / FS)
    assert np.max(np.abs(seg.rho(n) - rho)) < 1e-6
    assert np.max(np.abs(seg.rho_dot(n) - rate)) < 1e-9


def test_segment_layout():
    seg = segments.hermite_segments(lambda t: (1e6 + 0 * t, 0 * t), 3.5, FS, knot_s=1.0)
    n0, n1, coeffs = seg.to_arrays()
    assert n0.dtype == np.int64 and n1.dtype == np.int64 and coeffs.shape == (4, 4)
    assert n0[0] == 0 and n1[-1] == int(3.5 * FS)
    assert np.array_equal(n0[1:], n1[:-1])


def two_segments(d_rho=0.0, d_rate=0.0, gap=0):
    h = 1.0
    c_a = [1e6, -5000.0, 10.0, 0.5]
    end_rho = c_a[0] + c_a[1] * h + c_a[2] * h**2 + c_a[3] * h**3
    end_rate = c_a[1] + 2 * c_a[2] * h + 3 * c_a[3] * h**2
    c_b = [end_rho + d_rho, end_rate + d_rate, 0.0, 0.0]
    k = int(FS)
    return segments.RangeSegments(FS, np.array([0, k + gap]), np.array([k, 2 * k]),
                                  np.array([c_a, c_b]))


def test_check_c1_accepts_continuous():
    segments.check_c1(two_segments())


@pytest.mark.parametrize("kw", [{"d_rho": 1e-5}, {"d_rate": 1e-5}, {"gap": 1}])
def test_check_c1_rejects(kw):
    with pytest.raises(ValueError):
        segments.check_c1(two_segments(**kw))
