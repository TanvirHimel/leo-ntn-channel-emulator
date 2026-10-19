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
| M1a | Python golden model, orbit geometry, segments, waveform, SigMF | in progress |
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

## Build
_Placeholder._ Build steps for the C++ engine and the Python extension will be added in M1b.
Until then, the Python tests run from the repo root with:

```powershell
.\.venv\Scripts\python.exe -m pytest            # fast suite
.\.venv\Scripts\python.exe -m pytest -m slow    # heavy variants
```

## Results
**Results pending.** This README makes no performance or accuracy claims yet. Every number
that appears here will first be recorded, with how it was measured, in
[docs/results.md](docs/results.md).
