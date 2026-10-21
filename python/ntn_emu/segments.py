"""Cubic Hermite range segments (SPEC §2 "Segments"), the control-plane -> engine format.

Each segment k covers global samples [n0, n1) with local time u = (n - n0)/fs and
    rho(u) = c0 + c1 u + c2 u^2 + c3 u^3.
Coefficients come from (rho, rho_dot) at the knots, so value and slope match at every knot
(C1 by construction). to_arrays() gives the (n0, n1, coeffs) layout of _core.RangeTrack.
"""

from dataclasses import dataclass

import numpy as np

C1_TOL_RHO_M = 1e-6
C1_TOL_RATE_M_S = 1e-6


@dataclass
class RangeSegments:
    fs: float
    n0: np.ndarray      # int64 [N]
    n1: np.ndarray      # int64 [N]
    coeffs: np.ndarray  # float64 [N, 4]

    def _lookup(self, n):
        n = np.asarray(n, dtype=np.int64)
        k = np.searchsorted(self.n0, n, side="right") - 1
        if np.any(k < 0) or np.any(n >= self.n1[np.clip(k, 0, None)]):
            raise ValueError("sample index outside the segments")
        return k, (n - self.n0[k]) / self.fs

    def rho(self, n):
        k, u = self._lookup(n)
        c = self.coeffs[k]
        return ((c[:, 3] * u + c[:, 2]) * u + c[:, 1]) * u + c[:, 0]

    def rho_dot(self, n):
        k, u = self._lookup(n)
        c = self.coeffs[k]
        return (3.0 * c[:, 3] * u + 2.0 * c[:, 2]) * u + c[:, 1]

    def to_arrays(self):
        return self.n0.copy(), self.n1.copy(), self.coeffs.copy()


def hermite_segments(rho_fn, duration_s, fs, knot_s=1.0):
    """Hermite segments over [0, duration_s) from rho_fn(offsets_s) -> (rho, rho_dot).

    Knots every knot_s; knot_s*fs and duration_s*fs must be integers. The last segment may be
    shorter than knot_s.
    """
    step = knot_s * fs
    total = duration_s * fs
    if step != round(step) or total != round(total) or total <= 0:
        raise ValueError("knot_s*fs and duration_s*fs must be positive integers")
    step, total = int(round(step)), int(round(total))
    knots = np.append(np.arange(0, total, step, dtype=np.int64), np.int64(total))
    rho, rate = rho_fn(knots / fs)
    rho, rate = np.asarray(rho, dtype=np.float64), np.asarray(rate, dtype=np.float64)

    h = np.diff(knots) / fs
    p0, p1, m0, m1 = rho[:-1], rho[1:], rate[:-1], rate[1:]
    slope = (p1 - p0) / h
    c2 = (3.0 * slope - 2.0 * m0 - m1) / h
    c3 = (m0 + m1 - 2.0 * slope) / h**2
    coeffs = np.column_stack([p0, m0, c2, c3])
    return RangeSegments(float(fs), knots[:-1].copy(), knots[1:].copy(), coeffs)


def check_c1(seg, tol_rho=C1_TOL_RHO_M, tol_rate=C1_TOL_RATE_M_S):
    """Raise ValueError unless segments are contiguous and C1 at every interior knot."""
    if not np.array_equal(seg.n0[1:], seg.n1[:-1]) or np.any(seg.n1 <= seg.n0):
        raise ValueError("segments are not contiguous")
    h = (seg.n1 - seg.n0)[:-1] / seg.fs
    c = seg.coeffs
    end_rho = ((c[:-1, 3] * h + c[:-1, 2]) * h + c[:-1, 1]) * h + c[:-1, 0]
    end_rate = (3.0 * c[:-1, 3] * h + 2.0 * c[:-1, 2]) * h + c[:-1, 1]
    d_rho = np.abs(end_rho - c[1:, 0])
    d_rate = np.abs(end_rate - c[1:, 1])
    if np.any(d_rho >= tol_rho) or np.any(d_rate >= tol_rate):
        raise ValueError(f"C1 violated: max |d rho| {d_rho.max(initial=0):.3e} m, "
                         f"max |d rho_dot| {d_rate.max(initial=0):.3e} m/s")
