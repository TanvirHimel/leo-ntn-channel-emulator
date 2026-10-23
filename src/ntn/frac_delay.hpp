#pragma once

// Fractional delay: polyphase Kaiser-windowed sinc, 32 taps x 256 phases, linear interpolation
// between adjacent phases (SPEC §2). Provisional design for M1b; M2 finalises it.
//
// Output sample n reads x at position p = n - d. With d split like the golden model
// (python/ntn_emu/golden/fracdelay.py): d_int = floor(d), d_frac = d - d_int, snapped to an
// integer when within kSnap. An integer delay is a pure gather (bit-exact). Otherwise
// p = m0 + f with m0 = n - d_int - 1, f = 1 - d_frac in (0, 1), and
//     y[n] = sum_j x[m0 + j - (kHalf - 1)] * h_j(f),   j = 0 .. kTaps-1,
//     h_j(f) = sinc(f - off_j) * kaiser(f - off_j),     off_j = j - (kHalf - 1),
// so taps cover x[m0-15 .. m0+16]. Row r of the table is h(r / kPhases), r = 0 .. kPhases;
// row 0 is the exact unit impulse and row kPhases the impulse shifted by one.

#include <cstddef>
#include <cstdint>
#include <vector>

#include "ntn/types.hpp"

namespace ntn {

namespace fd {
inline constexpr int kTaps = 32;
inline constexpr int kHalf = kTaps / 2;
inline constexpr int kPhases = 256;
inline constexpr double kBeta = 9.7;   // Kaiser: ~97 dB for a 0.4 -> 0.6 fs transition, 32 taps
inline constexpr double kSnap = 1e-9;  // samples; same as the golden model
}  // namespace fd

class FracDelayTable {
 public:
    FracDelayTable();
    const float* row(int r) const { return &taps_[static_cast<std::size_t>(r) * fd::kTaps]; }

 private:
    std::vector<float> taps_;  // (kPhases + 1) x kTaps, row-major
};

// Process-wide table, built once.
const FracDelayTable& frac_delay_table();

// out[i] = x(n_first + i - d[i]) for i < count, where x[m] = hist[m - hist_first].
// The caller guarantees every read index m satisfies hist_first <= m < hist_first + len(hist);
// reads span [n - d - kHalf, n - d + kHalf] at most.
void frac_delay_scalar(const FracDelayTable& table, const cf32* hist, int64_t hist_first,
                       int64_t n_first, const double* d, cf32* out, std::size_t count);

}  // namespace ntn
