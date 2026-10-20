"""Golden fractional delay: float64 Kaiser-windowed sinc evaluated at each sample's exact delay.

y[n] = sum_m x[m] * h(n - d[n] - m), with h(t) = sinc(t) * kaiser(t) and support |t| <= L/2.
There is no polyphase table: the kernel is evaluated directly at every output sample's exact
fractional offset. That makes this path more accurate than the C++ engine (SPEC §2), so
"C++ vs golden" measures the C++ error.

Design: L = 128 taps, Kaiser beta = 14. This is about 136 dB stopband attenuation with a
transition of about +-0.035 fs around 0.5 fs. The passband |f| < 0.4 fs sits well inside it.

With snap=True (default), delays within SNAP of an integer are snapped to that integer and
applied as a pure gather. This is the golden counterpart of the C++ rule "phase-0 filter =
exact unit impulse". With snap=False every sample goes through the kernel.
Input samples outside [0, len(x)) are zero.

Output samples are independent of each other, so chunks (~64k samples, bounding memory) run
on a thread pool (NumPy releases the GIL). Chunk size and worker count do not change any bit.
"""

from concurrent.futures import ThreadPoolExecutor

import numpy as np
from scipy.special import i0

TAPS = 128
HALF = TAPS // 2
BETA = 14.0
SNAP = 1e-9  # samples
CHUNK = 65_536  # output samples per chunk
WORKERS = 2

_I0_BETA = i0(BETA)


def kaiser(t):
    """Kaiser window of half-width HALF at offsets t (samples); zero for |t| > HALF."""
    r = np.clip(1.0 - (t / HALF) ** 2, 0.0, None)
    return i0(BETA * np.sqrt(r)) / _I0_BETA


def fracdelay(x, d, chunk=CHUNK, workers=WORKERS, snap=True):
    """Delay x by d samples (scalar or one value per output sample); returns complex128."""
    x = np.asarray(x, dtype=np.complex128)
    n_samp = len(x)
    d = np.broadcast_to(np.asarray(d, dtype=np.float64), (n_samp,))
    # Zero padding: any read index outside [0, n_samp) is clipped onto a zero.
    xp = np.concatenate([np.zeros(HALF + 1), x, np.zeros(HALF + 1)])
    y = np.empty(n_samp, dtype=np.complex128)

    def run(s):
        e = min(s + chunk, n_samp)
        y[s:e] = _chunk(xp, d[s:e], np.arange(s, e, dtype=np.int64), snap)

    starts = range(0, n_samp, chunk)
    if workers > 1 and len(starts) > 1:
        with ThreadPoolExecutor(workers) as pool:
            list(pool.map(run, starts))
    else:
        for s in starts:
            run(s)
    return y


def _gather(xp, idx):
    """x[idx] with zeros outside [0, len(x)); xp is x padded by HALF + 1 zeros per side."""
    return xp[np.clip(idx + (HALF + 1), 0, len(xp) - 1)]


def _chunk(xp, d, n, snap):
    # Split d into integer and fractional parts before forming n - d, so precision is set
    # by |d| rather than by the sample index n.
    d_int = np.floor(d)
    d_frac = d - d_int
    if snap:
        up = d_frac > 1.0 - SNAP
        d_int[up] += 1.0
        d_frac[up] = 0.0
        exact = d_frac < SNAP
        d_frac[exact] = 0.0
    else:
        exact = np.zeros(len(d), dtype=bool)

    # Read position p = n - d = m0 + f with integer m0 and f in [0, 1).
    m0 = n - d_int.astype(np.int64)
    has = d_frac > 0.0
    m0 = np.where(has, m0 - 1, m0)
    f = np.where(has, 1.0 - d_frac, 0.0)
    wrap = f >= 1.0  # 1 - d_frac rounds to 1 for d_frac < 2^-53
    f[wrap] -= 1.0
    m0[wrap] += 1

    out = np.empty(len(n), dtype=np.complex128)
    out[exact] = _gather(xp, m0[exact])

    k = ~exact
    if np.any(k):
        m0k = m0[k]
        fk = f[k]
        sin_pf = np.sin(np.pi * fk)
        acc = np.zeros(len(fk), dtype=np.complex128)
        for j in range(TAPS):
            off = j - (HALF - 1)  # taps m0-63 .. m0+64
            t = fk - off
            # sinc(t) = sin(pi t)/(pi t), with sin(pi (f - off)) = (-1)^off sin(pi f).
            # t = 0 only for off = 0 and f = 0, where sinc = 1.
            num = sin_pf if off % 2 == 0 else -sin_pf
            if off == 0:
                zero = t == 0.0
                sinc = num / (np.pi * np.where(zero, 1.0, t))
                sinc[zero] = 1.0
            else:
                sinc = num / (np.pi * t)
            acc += _gather(xp, m0k + off) * (sinc * kaiser(t))
        out[k] = acc
    return out
