"""NR primary synchronization signal (PSS) sequence.

Source (transcribed from the fetched document, not from memory):
    ETSI TS 138 211 V17.11.0 (2026-09) = 3GPP TS 38.211 version 17.11.0 Release 17,
    "NR; Physical channels and modulation", clause 7.4.2.2.1 "Sequence generation" (p. 118).
    URL: https://www.etsi.org/deliver/etsi_ts/138200_138299/138211/17.11.00_60/ts_138211v171100p.pdf
    Fetched 2026-10-21T01:50:36Z, sha256 a4612390e7af933379baaaaf96b771a4e2ace109afe2c05233d494c6a207c90d.

    d_PSS(n) = 1 - 2 x(m),   m = (n + 43 N_ID^(2)) mod 127,   0 <= n < 127
    x(i + 7) = (x(i + 4) + x(i)) mod 2
    [x(6) x(5) x(4) x(3) x(2) x(1) x(0)] = [1 1 1 0 1 1 0]

Mapping to resources (clause 7.4.3) is not followed: this waveform is "NR-numerology CP-OFDM
(not NR-compliant)"; see waveform/ofdm.py for where the PSS is placed.
"""

import numpy as np

PSS_LENGTH = 127
_X_INIT = (0, 1, 1, 0, 1, 1, 1)  # x(0) .. x(6), from [x(6) .. x(0)] = [1 1 1 0 1 1 0]


def _m_sequence():
    x = list(_X_INIT)
    while len(x) < PSS_LENGTH:
        i = len(x) - 7
        x.append((x[i + 4] + x[i]) % 2)
    return np.array(x, dtype=np.int64)


def pss_sequence(n_id2):
    """d_PSS(n), n = 0..126, as float64 +-1, for N_ID^(2) in {0, 1, 2}."""
    if n_id2 not in (0, 1, 2):
        raise ValueError("N_ID2 must be 0, 1 or 2")
    m = (np.arange(PSS_LENGTH) + 43 * n_id2) % PSS_LENGTH
    return (1 - 2 * _m_sequence()[m]).astype(np.float64)
