// AVX2 kernels against their scalar oracles (SPEC §8 M2 test 3: relative error < -100 dB).
// Covers odd lengths and tails, buffers misaligned by 8 bytes, a mix of integer and fractional
// delays, and NaN sentinels just outside the allowed read range (an out-of-bounds read turns
// the output into NaN). Output guard cells catch out-of-bounds writes.
#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <limits>
#include <random>
#include <span>
#include <vector>

#include "ntn/cpu_check.hpp"
#include "ntn/frac_delay.hpp"
#include "test_util.hpp"

using ntn::cf32;
using ntn_test::cd;
namespace fd = ntn::fd;

namespace {

#define NTN_REQUIRE_AVX2()                                                  \
    do {                                                                    \
        if (!ntn::cpu_has_avx2()) GTEST_SKIP() << "CPU without AVX2";       \
    } while (0)

constexpr float kNaN = std::numeric_limits<float>::quiet_NaN();
constexpr int kGuard = 8;           // sentinel cells on each side
const cf32 kOutGuard(12345.0f, -54321.0f);

// n elements of T starting `offset` bytes after a 64-byte boundary.
template <class T>
struct Buf {
    std::vector<std::byte> raw;
    T* p;
    Buf(std::size_t n, std::size_t offset) : raw(n * sizeof(T) + 128 + offset) {
        const auto base = (reinterpret_cast<std::uintptr_t>(raw.data()) + 63) & ~std::uintptr_t{63};
        p = reinterpret_cast<T*>(base + offset);
    }
};

bool bit_equal(cf32 a, cf32 b) { return std::memcmp(&a, &b, sizeof(cf32)) == 0; }
bool finite(cf32 v) { return std::isfinite(v.real()) && std::isfinite(v.imag()); }

double rel_error_db(std::span<const cf32> y, std::span<const cf32> ref) {
    return ntn_test::error_db<cf32, cf32>(y, ref);
}

// Delay mix, cycling: exact integer, integer +1e-10 and -3e-10 (snap to the integer: gather),
// integer + 1e-6 (fractional, near the phase-0 row), and two uniform fractional delays.
std::vector<double> delay_mix(std::size_t n, uint64_t seed) {
    std::mt19937_64 rng(seed);
    std::uniform_int_distribution<int> k(-8, 8);
    std::uniform_real_distribution<double> u(-8.0, 8.0);
    std::vector<double> d(n);
    for (std::size_t i = 0; i < n; ++i) {
        switch (i % 6) {
            case 0: d[i] = k(rng); break;
            case 1: d[i] = k(rng) + 1e-10; break;
            case 2: d[i] = k(rng) - 3e-10; break;
            case 3: d[i] = k(rng) + 1e-6; break;
            default: d[i] = u(rng); break;
        }
    }
    return d;
}

struct CaseResult {
    double err_db = -400.0;  // AVX2 vs scalar
    std::size_t gather_samples = 0, differing_bits = 0;
};

// Runs both kernels on one case: `count` outputs from n_first, delays `d`, input history
// covering exactly the documented read envelope [n - d - kHalf, n - d + kHalf] (frac_delay.hpp),
// with NaN sentinels on both sides. Buffers start `offset` bytes past a 64-byte boundary.
CaseResult run_case(const std::vector<double>& dvals, int64_t n_first, std::size_t offset, const char* label) {
    const std::size_t count = dvals.size();
    CaseResult r;
    // Read envelope from the header contract.
    int64_t lo = 0, hi = -1;
    for (std::size_t i = 0; i < count; ++i) {
        const double p = static_cast<double>(n_first + static_cast<int64_t>(i)) - dvals[i];
        const auto a = static_cast<int64_t>(std::floor(p - fd::kHalf));
        const auto b = static_cast<int64_t>(std::ceil(p + fd::kHalf));
        lo = (i == 0) ? a : std::min(lo, a);
        hi = (i == 0) ? b : std::max(hi, b);
    }
    const std::size_t valid = count ? static_cast<std::size_t>(hi - lo + 1) : 0;
    Buf<cf32> hist(valid + 2 * kGuard, offset);
    std::mt19937_64 rng(static_cast<uint64_t>(count) * 7919u + offset);
    std::normal_distribution<float> g(0.0f, 1.0f);
    for (std::size_t i = 0; i < valid + 2 * kGuard; ++i) {
        const bool inside = i >= static_cast<std::size_t>(kGuard) && i < valid + kGuard;
        hist.p[i] = inside ? cf32(g(rng), g(rng)) : cf32(kNaN, kNaN);
    }
    const cf32* h = hist.p + kGuard;  // h[0] is x[lo]

    Buf<double> d(count, offset);
    std::copy(dvals.begin(), dvals.end(), d.p);
    Buf<cf32> ys(count + 2 * kGuard, offset), ya(count + 2 * kGuard, offset);
    std::fill(ys.p, ys.p + count + 2 * kGuard, kOutGuard);
    std::fill(ya.p, ya.p + count + 2 * kGuard, kOutGuard);

    ntn::frac_delay_scalar(ntn::frac_delay_table(), h, lo, n_first, d.p, ys.p + kGuard, count);
    ntn::frac_delay_avx2(ntn::frac_delay_table(), h, lo, n_first, d.p, ya.p + kGuard, count);

    for (int i = 0; i < kGuard; ++i) {
        EXPECT_TRUE(bit_equal(ya.p[i], kOutGuard) && bit_equal(ya.p[kGuard + count + i], kOutGuard))
            << label << ": AVX2 wrote outside the output, guard " << i;
        EXPECT_TRUE(bit_equal(ys.p[i], kOutGuard) && bit_equal(ys.p[kGuard + count + i], kOutGuard))
            << label << ": scalar wrote outside the output, guard " << i;
    }
    const std::span<const cf32> sa(ya.p + kGuard, count), ss(ys.p + kGuard, count);
    for (std::size_t i = 0; i < count; ++i) {
        EXPECT_TRUE(finite(sa[i])) << label << ": AVX2 output " << i << " not finite (read outside the envelope?)";
        EXPECT_TRUE(finite(ss[i])) << label << ": scalar output " << i << " not finite";
        const double frac = dvals[i] - std::floor(dvals[i]);
        if (frac < fd::kSnap || frac > 1.0 - fd::kSnap) {  // gather path: both bit-exact
            ++r.gather_samples;
            EXPECT_TRUE(bit_equal(sa[i], ss[i])) << label << ": integer delay " << dvals[i] << " at " << i;
        }
        r.differing_bits += !bit_equal(sa[i], ss[i]);
    }
    if (count) r.err_db = rel_error_db(sa, ss);
    return r;
}

}  // namespace

// ---- AVX2 table layout ------------------------------------------------------------------------

TEST(FracDelayAvx2Table, DuplicatedRowsAndDiffs) {
    // row2(r): each coefficient duplicated for interleaved re/im; diff2(r) = row(r+1) - row(r) in
    // float32 (the same rounded value the scalar kernel forms), duplicated.
    const auto& t = ntn::frac_delay_table();
    for (int r = 0; r <= fd::kPhases; ++r) {
        for (int j = 0; j < fd::kTaps; ++j) {
            ASSERT_EQ(t.row2(r)[2 * j], t.row(r)[j]) << r << "," << j;
            ASSERT_EQ(t.row2(r)[2 * j + 1], t.row(r)[j]) << r << "," << j;
            if (r < fd::kPhases) {
                const float diff = t.row(r + 1)[j] - t.row(r)[j];
                ASSERT_EQ(t.diff2(r)[2 * j], diff) << r << "," << j;
                ASSERT_EQ(t.diff2(r)[2 * j + 1], diff) << r << "," << j;
            }
        }
    }
}

// ---- AVX2 frac-delay kernel -------------------------------------------------------------------

TEST(FracDelayAvx2, IsNotTheScalarKernel) {
    // The AVX2 kernel sums in a different order (8 lanes, FMA), so on fractional delays it must
    // differ from the scalar kernel in some bits while agreeing to < -100 dB. A dispatch that
    // just forwards to the scalar kernel is bit-identical and fails here.
    NTN_REQUIRE_AVX2();
    std::mt19937_64 rng(5);
    std::uniform_real_distribution<double> u(-8.0, 8.0);
    std::vector<double> d(4096);
    for (auto& v : d) v = u(rng);
    const CaseResult r = run_case(d, 0, 0, "not-scalar");
    std::printf("[ measure  ] AVX2 vs scalar, 4096 fractional delays: %zu outputs differ in bits, error %.1f dB\n",
                r.differing_bits, r.err_db);
    EXPECT_GT(r.differing_bits, 0u);
    EXPECT_LT(r.err_db, -100.0);
}

TEST(FracDelayAvx2, MatchesScalarLengthsAlignmentDelayMix) {
    NTN_REQUIRE_AVX2();
    double worst = -400.0;
    for (std::size_t offset : {std::size_t{0}, std::size_t{8}}) {
        for (std::size_t len : {0, 1, 7, 8, 9, 31, 33, 4095, 4097}) {
            for (int64_t n_first : {int64_t{0}, int64_t{1'000'003}}) {
                char label[96];
                std::snprintf(label, sizeof label, "len %zu, offset %zu B, n_first %lld", len, offset,
                              static_cast<long long>(n_first));
                const CaseResult r = run_case(delay_mix(len, 100 + len), n_first, offset, label);
                if (len) {
                    EXPECT_LT(r.err_db, -100.0) << label;
                    worst = std::max(worst, r.err_db);
                }
                if (len >= 6) EXPECT_GT(r.gather_samples, 0u) << label;
            }
        }
    }
    std::printf("[ measure  ] AVX2 vs scalar frac delay: worst %.1f dB over lengths 0..4097, offsets 0/8 B, "
                "integer+fractional delay mix (limit -100 dB)\n",
                worst);
}

TEST(FracDelayAvx2, SingleSampleExactTapSpan) {
    // One output at a time with the input history holding exactly the documented taps
    // (frac_delay.hpp): x[m0 - (kHalf-1) .. m0 + kHalf], m0 = n - floor(d) - 1, for fractional
    // d (fraction in [0.01, 0.99], no snapping), and only x[n - d] for integer d. NaN sentinels
    // sit right next to the span, so an off-by-one read in either kernel gives NaN.
    NTN_REQUIRE_AVX2();
    std::mt19937_64 rng(77);
    std::uniform_int_distribution<int> k(-8, 8);
    std::uniform_real_distribution<double> fr(0.01, 0.99);
    std::normal_distribution<float> g(0.0f, 1.0f);
    double worst = -400.0;
    for (int trial = 0; trial < 300; ++trial) {
        const bool integer = trial % 4 == 0;
        const double d = k(rng) + (integer ? 0.0 : fr(rng));
        const int64_t n = 1000 + trial;
        const auto di = static_cast<int64_t>(std::floor(d));
        const int64_t first = integer ? n - di : n - di - fd::kHalf;
        const std::size_t span = integer ? 1 : static_cast<std::size_t>(fd::kTaps);
        for (std::size_t offset : {std::size_t{0}, std::size_t{8}}) {
            Buf<cf32> hist(span + 2, offset);
            hist.p[0] = hist.p[span + 1] = cf32(kNaN, kNaN);
            for (std::size_t i = 1; i <= span; ++i) hist.p[i] = {g(rng), g(rng)};
            cf32 ys, ya;
            ntn::frac_delay_scalar(ntn::frac_delay_table(), hist.p + 1, first, n, &d, &ys, 1);
            ntn::frac_delay_avx2(ntn::frac_delay_table(), hist.p + 1, first, n, &d, &ya, 1);
            ASSERT_TRUE(finite(ys)) << "scalar, d = " << d;
            ASSERT_TRUE(finite(ya)) << "AVX2, d = " << d;
            if (integer) {
                EXPECT_TRUE(bit_equal(ya, ys)) << "d = " << d;
            } else {
                const double e = rel_error_db(std::span<const cf32>(&ya, 1), std::span<const cf32>(&ys, 1));
                worst = std::max(worst, e);
                EXPECT_LT(e, -100.0) << "d = " << d;
            }
        }
    }
    std::printf("[ measure  ] AVX2 vs scalar, single samples on the exact tap span: worst %.1f dB\n", worst);
}
