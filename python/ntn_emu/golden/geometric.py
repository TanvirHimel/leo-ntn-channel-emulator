"""Golden geometric channel: delay and carrier phase from one range track rho(n).

The range track is a set of cubic segments:
    rho(u) = c0 + c1*u + c2*u^2 + c3*u^3,   u = (n - n0) / fs  [s],   n in [n0, n1).
It is evaluated by Horner in float64 at the exact int64 global sample index n.

Channel (baseband):
    y(n) = x(n - (tau(n)*fs - tau0)) * exp(-j*2*pi*fc*tau(n)),   tau = rho/c.
- The delay is relative to tau0_samples (integer-sample delay at window start).
- The phase uses the full tau(n). Carrier cycles fc*rho/c are computed in float64 and
  reduced to [-0.5, 0.5] before the exponential. Phase is never accumulated.
- Doppler f_d = -fc*rho_dot/c is never applied separately; it follows from the phase.
"""

import numpy as np

from .fracdelay import fracdelay

C = 299_792_458.0  # m/s, exact (SI); golden/ does not import ntn_emu.constants


def eval_rho(n, fs, n0, n1, coeffs):
    """rho at global sample indices n (int64) from contiguous segments [n0, n1)."""
    n = np.asarray(n, dtype=np.int64)
    n0 = np.asarray(n0, dtype=np.int64)
    n1 = np.asarray(n1, dtype=np.int64)
    coeffs = np.asarray(coeffs, dtype=np.float64).reshape(-1, 4)
    k = np.searchsorted(n0, n, side="right") - 1
    if np.any(k < 0) or np.any(n >= n1[np.clip(k, 0, None)]):
        raise ValueError("sample index outside the range track")
    u = (n - n0[k]) / fs
    c = coeffs[k]
    return ((c[:, 3] * u + c[:, 2]) * u + c[:, 1]) * u + c[:, 0]


def apply(x, *, fs, fc, n0, n1, coeffs, tau0_samples, n_start=0, delay=True, doppler=True):
    """Apply geometric delay and carrier phase to x, whose first sample is global index n_start."""
    x = np.asarray(x, dtype=np.complex128)
    n = n_start + np.arange(len(x), dtype=np.int64)
    rho = eval_rho(n, fs, n0, n1, coeffs)
    y = fracdelay(x, rho * fs / C - tau0_samples) if delay else x.copy()
    if doppler:
        cycles = fc * rho / C
        frac = cycles - np.round(cycles)
        y *= np.exp(-2j * np.pi * frac)
    return y
