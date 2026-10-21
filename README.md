# NTN Channel Emulator

A channel emulator for the downlink service link of a LEO non-terrestrial network (NTN). It
targets S-band (~2 GHz) and follows 3GPP Rel-17 (TR 38.811). A satellite from a pinned TLE and a
fixed ground UE define one range track rho(t). Delay, Doppler and carrier phase are all derived
from that track. The plan is a C++20 engine (MSVC x64) with a Python/SGP4 control plane and a
float64 Python reference ("golden") model.

This is a portfolio project. See [SPEC.md](SPEC.md) for the full specification.

## Status
| Milestone | Scope | State |
|---|---|---|
| M1a | Python golden model, orbit geometry, segments, waveform, SigMF | tests green, numbers in results.md; awaiting gate tag |
| M1b | C++ scalar engine + pybind11 | not started |
| M2 | Accuracy, AVX2, benchmarks | not started |
| M3 | NTN-TDL-C/D fading | not started |
| M4 | Receiver harness | not started |
| M5 | Cross-checks, final docs | not started |

## Scope labels
- **Waveform**: CP-OFDM with NR numerology. It is **not NR-compliant**: it has the NR PSS, comb
  pilots and random QPSK/16QAM data, with no SSS, PBCH or LDPC.
- **Link**: service-link downlink only. There is no uplink, feeder link or timing advance.
- **Physics modelled**: geometric delay, Doppler and carrier phase from SGP4 orbits, with
  light-time correction. TDL fading and AWGN come in later milestones.
- **Not modelled**: path loss, link budget, antennas and beams, atmospheric and ionospheric
  effects.
- **Implementation so far**: Python golden model only. There is no C++ engine yet.

## Design note: Doppler on a wideband signal
Doppler is applied by delaying the baseband signal by the time-varying range delay and rotating
the carrier phase. Each frequency component therefore gets a Doppler shift proportional to its
own RF frequency, f_c + f, not to the carrier alone. A test tone at +2.8 MHz in the M1a demo is
shifted by -(f_c + f_tone)·rho_dot/c, about 55 Hz more than the carrier Doppler -f_c·rho_dot/c
(+54.9865 to +55.1002 Hz measured). That is why the demo predicts -(f_c + f_tone)·rho_dot/c.

## Known limitations
- **The SGP4 range track has µm-level steps.** The SGP4 library's Kepler solver stops at a
  fixed tolerance (1e-12), so the satellite position jumps slightly whenever the solver's
  iteration count changes. On the pinned passes the largest step in range is 6.29 µm. That is
  at most 2.6e-04 rad of carrier phase, and at most 1.4e-04 Hz of Doppler artefact for a
  fraction of a second. Skyfield's two-part Julian date adds range noise of rms 7.9e-08 m.
- **Doppler is not taken from SGP4 velocity.** SGP4's velocity differs from the derivative of
  its own position by up to 1.08e-02 m/s, so the range rate is derived from the position-based
  range track instead.

These numbers, and how they were measured, are in [docs/results.md](docs/results.md).

## Build
_Placeholder._ Build steps for the C++ engine and the Python extension will be added in M1b.
Until then, the Python tests run from the repo root with:

```powershell
.\.venv\Scripts\python.exe -m pytest            # fast suite
.\.venv\Scripts\python.exe -m pytest -m slow    # heavy variants
```

## Results
All measured numbers, with the commit, machine and command that produced them, are in
[docs/results.md](docs/results.md). M1a (Python golden model) is recorded there; later
milestones are pending. Any number in this README also appears in results.md.

Headline from the M1a demo: the measured Doppler of the test tone matches the prediction within
0.0015 Hz (limit 1.0 Hz) across all ten 0.1 s sub-windows of a 1 s window on the 84.54 deg pass.
Figure: [docs/figures/m1a_doppler.png](docs/figures/m1a_doppler.png).
