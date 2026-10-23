// AVX2 rotator (xsimd). Compiled only in target ntn_avx2, with /arch:AVX2.
//
// Per block of 8 samples:
//  1. carrier cycles fc*rho/c in float64 (two 4-lane batches), the same operations in the same
//     order as rotate_scalar: (fc * rho) / c, frac = cycles - nearbyint(cycles) (round half to
//     even, like std::nearbyint), angle = (-2 pi) * frac. No multiply-add pair, so nothing can
//     be contracted to FMA; the angles are bit-identical to the scalar path.
//  2. cast to float32 per lane (the same rounding as the scalar static_cast).
//  3. xsimd::sincos on one 8-lane float batch: the only step that differs from the scalar path.
//  4. complex multiply on interleaved re/im: x * [c, c] + swap(x) * [-s, s], which rounds exactly
//     like the scalar xr*c - xi*s, xr*s + xi*c (a sign flip is exact; addition commutes).
// The last partial block is copied into a zero-padded block of 8 and run through the same code,
// so each output depends only on its own rho and input, never on the call's length.
//
// Helpers have internal linkage and the file avoids std:: inline templates, so no function with
// VEX code here is shared with the non-AVX2 objects (scripts/check_isa.py, check 4).
#include <xsimd/xsimd.hpp>

#include <cstddef>
#include <cstdint>

#include "ntn/constants.hpp"
#include "ntn/rotator.hpp"

#if !defined(__AVX2__) || !XSIMD_WITH_AVX2 || !XSIMD_WITH_FMA3_AVX2
#error "rotator_avx2.cpp must be compiled with /arch:AVX2 (target ntn_avx2)"
#endif

namespace ntn {

namespace {

using Arch = xsimd::fma3<xsimd::avx2>;
using BatchF = xsimd::batch<float, Arch>;   // 8 lanes
using BatchD = xsimd::batch<double, Arch>;  // 4 lanes
template <uint32_t... V>
using Mask = xsimd::batch_constant<uint32_t, Arch, V...>;
constexpr std::size_t kBlock = 8;
static_assert(BatchF::size == kBlock && BatchD::size * 2 == kBlock);

// Wrapped carrier angles of 8 samples, float32.
inline BatchF angles8(const double* rho, double fc_hz) {
    const BatchD fc(fc_hz), c(kC), minus_two_pi(-kTwoPi);
    alignas(32) double ad[kBlock];
    for (std::size_t h = 0; h < kBlock; h += BatchD::size) {
        const BatchD cycles = (fc * BatchD::load_unaligned(rho + h)) / c;
        const BatchD frac = cycles - xsimd::nearbyint(cycles);
        (minus_two_pi * frac).store_aligned(ad + h);
    }
    alignas(32) float af[kBlock];
    for (std::size_t k = 0; k < kBlock; ++k) af[k] = static_cast<float>(ad[k]);
    return BatchF::load_aligned(af);
}

// io: 8 interleaved complex samples (16 floats).
inline void rotate8(const double* rho, double fc_hz, float* io) {
    const auto [s, c] = xsimd::sincos(angles8(rho, fc_hz));
    const BatchF sign(-1.0f, 1.0f, -1.0f, 1.0f, -1.0f, 1.0f, -1.0f, 1.0f);
    const BatchF c_lo = xsimd::swizzle(c, Mask<0, 0, 1, 1, 2, 2, 3, 3>());
    const BatchF c_hi = xsimd::swizzle(c, Mask<4, 4, 5, 5, 6, 6, 7, 7>());
    const BatchF s_lo = xsimd::swizzle(s, Mask<0, 0, 1, 1, 2, 2, 3, 3>()) * sign;
    const BatchF s_hi = xsimd::swizzle(s, Mask<4, 4, 5, 5, 6, 6, 7, 7>()) * sign;
    const BatchF x_lo = BatchF::load_unaligned(io);
    const BatchF x_hi = BatchF::load_unaligned(io + kBlock);
    using Swap = Mask<1, 0, 3, 2, 5, 4, 7, 6>;  // (re, im) -> (im, re)
    (x_lo * c_lo + xsimd::swizzle(x_lo, Swap()) * s_lo).store_unaligned(io);
    (x_hi * c_hi + xsimd::swizzle(x_hi, Swap()) * s_hi).store_unaligned(io + kBlock);
}

}  // namespace

void rotate_avx2(const double* rho, double fc_hz, cf32* io, std::size_t count) {
    // std::complex<float> arrays may be accessed as interleaved float pairs ([complex.numbers]).
    auto* f = reinterpret_cast<float*>(io);
    std::size_t i = 0;
    for (; i + kBlock <= count; i += kBlock) rotate8(rho + i, fc_hz, f + 2 * i);
    if (i == count) return;
    alignas(32) double r[kBlock] = {};
    alignas(32) float t[2 * kBlock] = {};
    const std::size_t rem = count - i;
    for (std::size_t k = 0; k < rem; ++k) {
        r[k] = rho[i + k];
        t[2 * k] = f[2 * (i + k)];
        t[2 * k + 1] = f[2 * (i + k) + 1];
    }
    rotate8(r, fc_hz, t);
    for (std::size_t k = 0; k < 2 * rem; ++k) f[2 * i + k] = t[k];
}

void sincos_avx2(const float* ang, float* s, float* c, std::size_t count) {
    std::size_t i = 0;
    for (; i + kBlock <= count; i += kBlock) {
        const auto [vs, vc] = xsimd::sincos(BatchF::load_unaligned(ang + i));
        vs.store_unaligned(s + i);
        vc.store_unaligned(c + i);
    }
    if (i == count) return;
    alignas(32) float a[kBlock] = {}, ts[kBlock], tc[kBlock];
    const std::size_t rem = count - i;
    for (std::size_t k = 0; k < rem; ++k) a[k] = ang[i + k];
    const auto [vs, vc] = xsimd::sincos(BatchF::load_aligned(a));
    vs.store_aligned(ts);
    vc.store_aligned(tc);
    for (std::size_t k = 0; k < rem; ++k) {
        s[i + k] = ts[k];
        c[i + k] = tc[k];
    }
}

}  // namespace ntn
