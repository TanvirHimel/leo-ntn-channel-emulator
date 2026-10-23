#pragma once

// Fractional delay: polyphase Kaiser-windowed sinc with linear interpolation between adjacent
// phases (SPEC §2). The tap count T and phase count P are compile-time parameters; the
// production design is 32 taps x 256 phases (fd::kTaps, fd::kPhases); the M2 design sweep
// (test FracDelayDesign.Sweep) compares it with 16/24 taps and 128 phases.
//
// Output sample n reads x at position p = n - d. With d split like the golden model
// (python/ntn_emu/golden/fracdelay.py): d_int = floor(d), d_frac = d - d_int, snapped to an
// integer when within kSnap. An integer delay is a pure gather (bit-exact). Otherwise
// p = m0 + f with m0 = n - d_int - 1, f = 1 - d_frac in (0, 1), and
//     y[n] = sum_j x[m0 + j - (H - 1)] * h_j(f),   j = 0 .. T-1,  H = T/2,
//     h_j(f) = sinc(f - off_j) * kaiser(f - off_j),  off_j = j - (H - 1),
// so taps cover x[m0-(H-1) .. m0+H]. Row r of the table is h(r / P), r = 0 .. P; row 0 is the
// exact unit impulse and row P the impulse shifted by one.
//
// Kaiser beta comes from the window design formula (not tuned): transition 0.4 -> 0.6 fs,
// dw = 0.4 pi, order M = T - 1, A = 2.285 M dw + 8 dB, beta = 0.1102 (A - 8.7) for A > 50,
// 0.5842 (A - 21)^0.4 + 0.07886 (A - 21) for 21 <= A <= 50, else 0.
//
// Templates are defined in frac_delay.cpp (target ntn_scalar, no /arch:AVX2) and explicitly
// instantiated there for T in {16, 24, 32} x P in {128, 256}; no other design links.

#include <cstddef>
#include <cstdint>
#include <vector>

#include "ntn/types.hpp"

namespace ntn {

namespace fd {
inline constexpr int kTaps = 32;  // production design
inline constexpr int kHalf = kTaps / 2;
inline constexpr int kPhases = 256;
inline constexpr double kSnap = 1e-9;  // samples; same as the golden model

double kaiser_design_atten_db(int taps);  // A from the design formula above, dB
double kaiser_beta(int taps);             // beta from A
}  // namespace fd

template <int T, int P>
class FracDelayTableT {
    static_assert(T >= 4 && T % 2 == 0, "tap count must be even");
    static_assert(P >= 2 && P % 2 == 0, "phase count must be even");

 public:
    static constexpr int kTaps = T;
    static constexpr int kHalf = T / 2;
    static constexpr int kPhases = P;

    FracDelayTableT();
    const float* row(int r) const { return &taps_[static_cast<std::size_t>(r) * T]; }
    double beta() const { return beta_; }

 private:
    double beta_;
    std::vector<float> taps_;  // (P + 1) x T, row-major
};

// Process-wide table per design, built once.
template <int T, int P>
const FracDelayTableT<T, P>& frac_delay_table_t();

// out[i] = x(n_first + i - d[i]) for i < count, where x[m] = hist[m - hist_first].
// The caller guarantees every read index m satisfies hist_first <= m < hist_first + len(hist);
// reads span [n - d - T/2, n - d + T/2] at most.
template <int T, int P>
void frac_delay_scalar_t(const FracDelayTableT<T, P>& table, const cf32* hist, int64_t hist_first,
                         int64_t n_first, const double* d, cf32* out, std::size_t count);

#define NTN_FD_DESIGN(T, P)                                                                           \
    extern template class FracDelayTableT<T, P>;                                                      \
    extern template const FracDelayTableT<T, P>& frac_delay_table_t<T, P>();                          \
    extern template void frac_delay_scalar_t<T, P>(const FracDelayTableT<T, P>&, const cf32*, int64_t, \
                                                   int64_t, const double*, cf32*, std::size_t);
NTN_FD_DESIGN(16, 128)
NTN_FD_DESIGN(16, 256)
NTN_FD_DESIGN(24, 128)
NTN_FD_DESIGN(24, 256)
NTN_FD_DESIGN(32, 128)
NTN_FD_DESIGN(32, 256)
#undef NTN_FD_DESIGN

// Production design. Out of line (frac_delay.cpp), so callers compiled with /arch:AVX2
// (ntn_core) get no inline copy of them.
using FracDelayTable = FracDelayTableT<fd::kTaps, fd::kPhases>;

const FracDelayTable& frac_delay_table();  // same object as frac_delay_table_t<32, 256>()

void frac_delay_scalar(const FracDelayTable& table, const cf32* hist, int64_t hist_first, int64_t n_first,
                       const double* d, cf32* out, std::size_t count);

}  // namespace ntn
