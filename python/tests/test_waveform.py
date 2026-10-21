"""Waveform: NR PSS (TS 38.211 §7.4.2.2) and NR-numerology CP-OFDM, closed-form checks."""

import numpy as np
import pytest

from ntn_emu.waveform import ofdm, pss

# --- PSS -------------------------------------------------------------------------------


@pytest.mark.parametrize("n_id2", [0, 1, 2])
def test_pss_is_bpsk_length_127(n_id2):
    d = pss.pss_sequence(n_id2)
    assert d.shape == (127,)
    assert set(np.unique(d)) == {-1.0, 1.0}


def test_pss_first_elements_from_initial_state():
    """N_ID2 = 0: m = n, so d(0..6) = 1 - 2*[x(0)..x(6)] with [x(6)..x(0)] = [1 1 1 0 1 1 0]."""
    x0_to_x6 = np.array([0, 1, 1, 0, 1, 1, 1])
    np.testing.assert_array_equal(pss.pss_sequence(0)[:7], 1 - 2 * x0_to_x6)


def test_pss_recursion():
    """x(i+7) = (x(i+4) + x(i)) mod 2 holds for the whole N_ID2 = 0 sequence (m = n)."""
    x = (1 - pss.pss_sequence(0)) / 2
    i = np.arange(127 - 7)
    np.testing.assert_array_equal(x[i + 7], (x[i + 4] + x[i]) % 2)


@pytest.mark.parametrize("n_id2", [0, 1, 2])
def test_pss_periodic_autocorrelation(n_id2):
    """m-sequence: periodic autocorrelation 127 at lag 0, -1 at every other lag."""
    d = pss.pss_sequence(n_id2)
    r = np.array([np.dot(d, np.roll(d, -k)) for k in range(127)])
    assert r[0] == 127
    assert np.all(r[1:] == -1)


@pytest.mark.parametrize("n_id2", [1, 2])
def test_pss_cyclic_shift_43(n_id2):
    """d_N(n) = d_0((n + 43 N) mod 127)."""
    np.testing.assert_array_equal(pss.pss_sequence(n_id2),
                                  pss.pss_sequence(0)[(np.arange(127) + 43 * n_id2) % 127])


# --- OFDM ------------------------------------------------------------------------------


@pytest.fixture(scope="module")
def frame():
    return ofdm.generate(n_slots=12, cfg=ofdm.OfdmConfig(modulation="qpsk", seed=1))


def test_numerology_constants():
    assert ofdm.N_FFT == 512 and ofdm.N_SC == 300 and ofdm.FS == 7.68e6 and ofdm.SCS == 15e3
    assert ofdm.FS == ofdm.N_FFT * ofdm.SCS
    assert ofdm.SYMBOLS_PER_SLOT == 14 and ofdm.SLOT_SAMPLES == 7680


def test_cp_lengths_and_slot_length(frame):
    cp = frame.cp_lengths[:14]
    np.testing.assert_array_equal(cp, [40] + [36] * 6 + [40] + [36] * 6)
    assert np.sum(cp) + 14 * 512 == 7680
    assert len(frame.samples) == 12 * 7680
    np.testing.assert_array_equal(np.diff(frame.symbol_starts), (frame.cp_lengths + 512)[:-1])


def test_cp_equals_symbol_tail(frame):
    for s, cp in zip(frame.symbol_starts, frame.cp_lengths):
        sym = frame.samples[s:s + cp + 512]
        np.testing.assert_array_equal(sym[:cp], sym[-cp:])


def demodulate(frame):
    """CP removal + FFT, undoing the generator's N_FFT/sqrt(N_SC) scaling; bins k = -150..149."""
    k = np.arange(-150, 150) % ofdm.N_FFT
    out = []
    for s, cp in zip(frame.symbol_starts, frame.cp_lengths):
        spec = np.fft.fft(frame.samples[s + cp:s + cp + 512]) * np.sqrt(ofdm.N_SC) / ofdm.N_FFT
        out.append(spec[k])
    return np.array(out)


def test_demodulation_recovers_grid(frame):
    assert np.max(np.abs(demodulate(frame) - frame.grid)) < 1e-12


def test_unused_bins_are_empty(frame):
    s, cp = frame.symbol_starts[1], frame.cp_lengths[1]
    spec = np.fft.fft(frame.samples[s + cp:s + cp + 512])
    used = np.zeros(512, bool)
    used[np.arange(-150, 150) % 512] = True
    assert np.max(np.abs(spec[~used])) < 1e-12 * np.max(np.abs(spec))


def test_pss_placement_symmetric_about_centre(frame):
    """PSS on signed SC k = -63..63 (grid 87..213), symmetric about DC; DC carries d(63)."""
    lo, hi = ofdm.PSS_K_RANGE
    assert (lo, hi) == (-63, 63) and lo == -hi and hi - lo + 1 == 127
    g_lo, g_hi = ofdm.PSS_GRID_RANGE
    assert (g_lo, g_hi) == (87, 213) and (g_lo + g_hi) / 2 == ofdm.N_SC // 2
    row = frame.grid[0]
    d = pss.pss_sequence(frame.cfg.n_id2)
    np.testing.assert_array_equal(row[87:214], d)
    assert row[150] == d[63]  # grid index 150 is k = 0 (DC)
    non_pss = np.ones(300, bool)
    non_pss[87:214] = False
    non_pss[ofdm.PILOT_GRID] = False
    assert np.all(row[non_pss] == 0)


def test_pss_once_per_frame():
    """Symbol 0 of each 10 ms frame (140 symbols = 10 slots = 76 800 samples)."""
    short = ofdm.generate(n_slots=9, cfg=ofdm.OfdmConfig(seed=1))     # 126 symbols
    assert list(short.pss_symbols) == [0]
    long = ofdm.generate(n_slots=21, cfg=ofdm.OfdmConfig(seed=1))     # 294 symbols
    assert list(long.pss_symbols) == [0, 140, 280]
    assert list(long.symbol_starts[long.pss_symbols]) == [0, 76_800, 153_600]


def test_pilots_on_every_symbol(frame):
    """Comb pilots on every 6th SC (grid index g % 6 == 0) in every symbol, known QPSK."""
    assert np.array_equal(ofdm.PILOT_GRID, np.arange(0, 300, 6))
    pilots = ofdm.pilot_values(frame.cfg, frame.grid.shape[0])
    for sym in range(frame.grid.shape[0]):
        keep = ofdm.PILOT_GRID
        if sym in frame.pss_symbols:
            keep = keep[(keep < 87) | (keep > 213)]
        idx = np.searchsorted(ofdm.PILOT_GRID, keep)
        np.testing.assert_array_equal(frame.grid[sym, keep], pilots[sym, idx])
    assert np.allclose(np.abs(pilots), 1.0, atol=1e-15)


def test_unit_power_qpsk_symbol(frame):
    """All 300 RE unit-modulus -> useful-part mean power exactly 1 (Parseval)."""
    s, cp = frame.symbol_starts[1], frame.cp_lengths[1]
    assert abs(np.mean(np.abs(frame.samples[s + cp:s + cp + 512]) ** 2) - 1.0) < 1e-12


def test_16qam_levels():
    f = ofdm.generate(n_slots=1, cfg=ofdm.OfdmConfig(modulation="16qam", seed=2))
    data = f.grid[1][~np.isin(np.arange(300), ofdm.PILOT_GRID)]
    levels = np.unique(np.round(np.concatenate([data.real, data.imag]) * np.sqrt(10), 12))
    np.testing.assert_array_equal(levels, [-3, -1, 1, 3])


def test_seeded_determinism():
    a = ofdm.generate(n_slots=2, cfg=ofdm.OfdmConfig(seed=5)).samples
    b = ofdm.generate(n_slots=2, cfg=ofdm.OfdmConfig(seed=5)).samples
    c = ofdm.generate(n_slots=2, cfg=ofdm.OfdmConfig(seed=6)).samples
    assert np.array_equal(a, b) and not np.array_equal(a, c)
