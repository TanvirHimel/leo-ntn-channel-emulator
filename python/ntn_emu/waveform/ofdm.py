"""NR-numerology CP-OFDM test waveform (not NR-compliant), SPEC §2 "Waveform".

Numerology: SCS 15 kHz, N_FFT 512, fs 7.68 Msps, 300 subcarriers (25 PRB, 5 MHz), normal CP
(40 samples for symbols 0 and 7 of each 1 ms slot, 36 otherwise; 7680 samples per slot).

Subcarrier indexing: signed k = -150..149 maps to FFT bin k mod 512; grid index g = k + 150
runs 0..299. There is no DC null. Subcarrier k sits at k * 15 kHz.

Content:
- PSS (TS 38.211 §7.4.2.2, see pss.py) in symbol 0 of every 10 ms frame (every 140 symbols),
  on k = -63..+63 (grid 87..213): 127 subcarriers symmetric about the centre (DC) bin, with
  d_PSS(n) on k = n - 63, so DC carries d_PSS(63). Not the NR SSB mapping.
- Comb pilots on every symbol, every 6th subcarrier (g % 6 == 0, 50 pilots), known seeded QPSK.
  In PSS symbols the pilots inside the PSS band are replaced by the PSS; other non-PSS REs of
  a PSS symbol are empty.
- Data: seeded QPSK or 16QAM on all other REs, unit average RE power.

Time-domain scaling N_FFT/sqrt(N_SC) gives mean sample power = mean RE power (1 for QPSK).
"""

from dataclasses import dataclass

import numpy as np

from ntn_emu.waveform.pss import pss_sequence

N_FFT = 512
N_SC = 300
SCS = 15e3
FS = N_FFT * SCS  # 7.68e6
SYMBOLS_PER_SLOT = 14
CP_PATTERN = np.array([40] + [36] * 6 + [40] + [36] * 6, dtype=np.int64)
SLOT_SAMPLES = int(CP_PATTERN.sum() + SYMBOLS_PER_SLOT * N_FFT)  # 7680
SYMBOLS_PER_FRAME = 140  # 10 ms
K_RANGE = (-150, 149)
PSS_K_RANGE = (-63, 63)
PSS_GRID_RANGE = (PSS_K_RANGE[0] + 150, PSS_K_RANGE[1] + 150)  # (87, 213)
PILOT_GRID = np.arange(0, N_SC, 6)

_BINS = np.arange(K_RANGE[0], K_RANGE[1] + 1) % N_FFT


@dataclass(frozen=True)
class OfdmConfig:
    modulation: str = "qpsk"  # "qpsk" or "16qam"
    n_id2: int = 0
    seed: int = 0
    pilot_seed: int = 0x5EED


@dataclass
class OfdmFrame:
    samples: np.ndarray        # complex128, len = n_slots * 7680
    grid: np.ndarray           # complex128 [n_symbols, 300], grid index g = k + 150
    symbol_starts: np.ndarray  # int64, sample index of each symbol's CP start
    cp_lengths: np.ndarray     # int64 per symbol
    pss_symbols: np.ndarray    # symbol indices carrying the PSS
    cfg: OfdmConfig


def _qpsk(rng, shape):
    b = rng.integers(0, 2, size=shape + (2,))
    return ((1 - 2 * b[..., 0]) + 1j * (1 - 2 * b[..., 1])) / np.sqrt(2.0)


def _qam16(rng, shape):
    b = rng.integers(0, 2, size=shape + (4,))
    i = (1 - 2 * b[..., 0]) * (2 - (1 - 2 * b[..., 2]))
    q = (1 - 2 * b[..., 1]) * (2 - (1 - 2 * b[..., 3]))
    return (i + 1j * q) / np.sqrt(10.0)


def pilot_values(cfg, n_symbols):
    """Known pilot QPSK, [n_symbols, 50], for grid indices PILOT_GRID."""
    return _qpsk(np.random.default_rng(cfg.pilot_seed), (n_symbols, len(PILOT_GRID)))


def generate(n_slots, cfg=OfdmConfig()):
    n_sym = SYMBOLS_PER_SLOT * n_slots
    cp = np.tile(CP_PATTERN, n_slots)
    starts = np.concatenate([[0], np.cumsum(cp + N_FFT)[:-1]]).astype(np.int64)
    pss_syms = np.arange(0, n_sym, SYMBOLS_PER_FRAME)

    rng = np.random.default_rng(cfg.seed)
    if cfg.modulation == "qpsk":
        grid = _qpsk(rng, (n_sym, N_SC))
    elif cfg.modulation == "16qam":
        grid = _qam16(rng, (n_sym, N_SC))
    else:
        raise ValueError(f"unknown modulation {cfg.modulation!r}")
    grid[:, PILOT_GRID] = pilot_values(cfg, n_sym)

    g_lo, g_hi = PSS_GRID_RANGE
    in_band = np.zeros(N_SC, dtype=bool)
    in_band[g_lo:g_hi + 1] = True
    keep_pilot = np.zeros(N_SC, dtype=bool)
    keep_pilot[PILOT_GRID] = ~in_band[PILOT_GRID]
    d = pss_sequence(cfg.n_id2)
    for s in pss_syms:
        row = np.where(keep_pilot, grid[s], 0.0)
        row[g_lo:g_hi + 1] = d
        grid[s] = row

    bins = np.zeros((n_sym, N_FFT), dtype=np.complex128)
    bins[:, _BINS] = grid
    body = np.fft.ifft(bins, axis=1) * (N_FFT / np.sqrt(N_SC))
    samples = np.empty(n_slots * SLOT_SAMPLES, dtype=np.complex128)
    for s in range(n_sym):
        c, p = cp[s], starts[s]
        samples[p:p + c] = body[s, -c:]
        samples[p + c:p + c + N_FFT] = body[s]
    return OfdmFrame(samples, grid, starts, cp, pss_syms, cfg)
