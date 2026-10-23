// AVX2 frac-delay kernel (xsimd). Compiled only in target ntn_avx2, with /arch:AVX2.
//
// Same contract and read span as frac_delay_scalar (frac_delay.hpp). Vectorised across the
// taps: the input is interleaved re/im, so T complex taps are 2T floats = 2T/8 AVX2 lanes of 8,
// multiplied by the duplicated coefficient rows FracDelayTable::row2/diff2. The delay split
// (floor, snap, phase row r, weight w) is the scalar kernel's float64 code, so both kernels pick
// the same row and weight; they differ only in the float32 summation (8 lanes + FMA).
//
// Helpers have internal linkage and the file avoids std:: inline templates (std::min etc.), so
// no function compiled here with VEX code is shared with the non-AVX2 objects
// (scripts/check_isa.py, check 4).
#include <xsimd/xsimd.hpp>

#include <cmath>
#include <cstdint>

#include "ntn/frac_delay.hpp"

// /arch:AVX2 must be active for this file (MSVC defines __AVX2__ under it). xsimd then also
// enables FMA3 (XSIMD_WITH_FMA3_AVX2), which the kernel uses.
#if !defined(__AVX2__) || !XSIMD_WITH_AVX2 || !XSIMD_WITH_FMA3_AVX2
#error "frac_delay_avx2.cpp must be compiled with /arch:AVX2 (target ntn_avx2)"
#endif

namespace ntn {

namespace {

using Arch = xsimd::fma3<xsimd::avx2>;
using BatchF = xsimd::batch<float, Arch>;
constexpr int kLanes = static_cast<int>(BatchF::size);  // 8
constexpr int kFloats = 2 * fd::kTaps;                  // interleaved re/im
static_assert(kFloats % (2 * kLanes) == 0, "two accumulators need 2T to be a multiple of 16 floats");

// y = sum_j x[j] * (a[j] + w (b[j] - a[j])) over interleaved complex x; returns re/im sums.
inline cf32 dot_interleaved(const float* x, const float* a2, const float* d2, float w) {
    const BatchF wv(w);
    BatchF acc0(0.0f), acc1(0.0f);
    for (int k = 0; k < kFloats; k += 2 * kLanes) {
        const BatchF c0 = xsimd::fma(wv, BatchF::load_unaligned(d2 + k), BatchF::load_unaligned(a2 + k));
        const BatchF c1 = xsimd::fma(wv, BatchF::load_unaligned(d2 + k + kLanes), BatchF::load_unaligned(a2 + k + kLanes));
        acc0 = xsimd::fma(BatchF::load_unaligned(x + k), c0, acc0);
        acc1 = xsimd::fma(BatchF::load_unaligned(x + k + kLanes), c1, acc1);
    }
    alignas(32) float t[kLanes];
    (acc0 + acc1).store_aligned(t);
    return {(t[0] + t[2]) + (t[4] + t[6]), (t[1] + t[3]) + (t[5] + t[7])};
}

}  // namespace

void frac_delay_avx2(const FracDelayTable& table, const cf32* hist, int64_t hist_first, int64_t n_first,
                     const double* d, cf32* out, std::size_t count) {
    for (std::size_t i = 0; i < count; ++i) {
        const int64_t n = n_first + static_cast<int64_t>(i);
        // Delay split: identical to frac_delay_scalar_t (float64).
        double d_int = std::floor(d[i]);
        double d_frac = d[i] - d_int;
        if (d_frac > 1.0 - fd::kSnap) {
            d_int += 1.0;
            d_frac = 0.0;
        } else if (d_frac < fd::kSnap) {
            d_frac = 0.0;
        }
        int64_t m0 = n - static_cast<int64_t>(d_int);
        if (d_frac == 0.0) {  // integer delay: pure gather, bit-exact
            out[i] = hist[m0 - hist_first];
            continue;
        }
        m0 -= 1;
        const double f = 1.0 - d_frac;  // read position m0 + f, f in (0, 1)
        const double p = f * fd::kPhases;
        int r = static_cast<int>(p);
        if (r > fd::kPhases - 1) r = fd::kPhases - 1;
        const auto w = static_cast<float>(p - r);
        // std::complex<float> arrays may be accessed as interleaved float pairs ([complex.numbers]).
        const auto* x = reinterpret_cast<const float*>(hist + (m0 - (fd::kHalf - 1) - hist_first));
        out[i] = dot_interleaved(x, table.row2(r), table.diff2(r), w);
    }
}

}  // namespace ntn
