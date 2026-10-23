// AVX2 kernels (frac delay, rotator, sincos) against their scalar oracles and closed forms
// (SPEC §8 M2 test 3: AVX2 vs scalar relative error < -100 dB).
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
#include <numbers>
#include <random>
#include <span>
#include <utility>
#include <vector>

#include "ntn/cpu_check.hpp"
#include "ntn/frac_delay.hpp"
#include "ntn/range_track.hpp"
#include "ntn/rotator.hpp"
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

// ---- AVX2 rotator -----------------------------------------------------------------------------
//
// rotate_avx2: carrier cycles fc*rho/c in float64 (same operations as rotate_scalar), wrapped by
// nearbyint, cast to float32, vectorised float32 sincos (sincos_avx2), complex multiply. Only the
// sincos differs from the scalar path; the tail is padded and run through the same vector code.

namespace {

constexpr double kCLight = 299'792'458.0;
constexpr double kFsRot = 7'680'000.0;
constexpr double kFcRot = 2'000'000'000.0;
constexpr double kPiD = std::numbers::pi;

double wrap_pi(double a) { return a - 2.0 * kPiD * std::nearbyint(a / (2.0 * kPiD)); }

std::vector<double> realistic_rho(std::size_t n, uint64_t seed) {
    // Pass-like ranges (0.5e6 .. 2.5e6 m) with sub-mm detail, so every phase is reachable.
    std::mt19937_64 rng(seed);
    std::uniform_real_distribution<double> u(0.5e6, 2.5e6);
    std::vector<double> r(n);
    for (auto& v : r) v = u(rng);
    return r;
}

struct RotCase {
    double err_db = -400.0;
    std::size_t differing_bits = 0;
};

// Both rotators on the same rho and input; io buffers start `offset` bytes past a 64-byte
// boundary with guard cells on both sides (out-of-bounds writes).
RotCase run_rot_case(std::size_t count, std::size_t offset, const char* label) {
    const auto rv = realistic_rho(count, 300 + count);
    Buf<double> rho(count, offset);
    std::copy(rv.begin(), rv.end(), rho.p);
    Buf<cf32> ys(count + 2 * kGuard, offset), ya(count + 2 * kGuard, offset);
    std::mt19937_64 rng(count + 17);
    std::normal_distribution<float> g(0.0f, 1.0f);
    for (std::size_t i = 0; i < count + 2 * kGuard; ++i) {
        const bool inside = i >= static_cast<std::size_t>(kGuard) && i < count + kGuard;
        ys.p[i] = ya.p[i] = inside ? cf32(g(rng), g(rng)) : kOutGuard;
    }
    ntn::rotate_scalar(rho.p, kFcRot, ys.p + kGuard, count);
    ntn::rotate_avx2(rho.p, kFcRot, ya.p + kGuard, count);
    for (int i = 0; i < kGuard; ++i) {
        EXPECT_TRUE(bit_equal(ya.p[i], kOutGuard) && bit_equal(ya.p[kGuard + count + i], kOutGuard))
            << label << ": AVX2 rotator wrote outside the buffer, guard " << i;
    }
    RotCase r;
    const std::span<const cf32> sa(ya.p + kGuard, count), ss(ys.p + kGuard, count);
    for (std::size_t i = 0; i < count; ++i) r.differing_bits += !bit_equal(sa[i], ss[i]);
    if (count) r.err_db = rel_error_db(sa, ss);
    return r;
}

// M1a test 1 on the AVX2 rotator: DC input, linear rho = rho0 + v t from a RangeTrack
// (rho_block, as the channel evaluates it). Reference phase, float64:
// -2 pi (frac(fc rho0 / c) + frac(B n)), B = fc v / (c fs); see test_geometric_channel.cpp.
struct RotPhaseStats {
    double max_err = 0.0, max_err_head = 0.0, max_err_tail = 0.0, max_mag_err = 0.0;
};

RotPhaseStats rot_tone_stats(double rho0, double v, int64_t count) {
    const ntn::RangeTrack track(kFsRot, {{0, count, {rho0, v, 0.0, 0.0}}});
    const double a = kFcRot * rho0 / kCLight;
    const double fa = a - std::nearbyint(a);
    const double b = kFcRot * v / (kCLight * kFsRot);
    RotPhaseStats st;
    constexpr int64_t kBlock = 1 << 20;
    std::vector<double> rho(kBlock);
    std::vector<cf32> y(kBlock);
    for (int64_t s = 0; s < count; s += kBlock) {
        const auto m = static_cast<std::size_t>(std::min(kBlock, count - s));
        track.rho_block(s, std::span<double>(rho.data(), m));
        std::fill(y.begin(), y.begin() + static_cast<std::ptrdiff_t>(m), cf32(1.0f, 0.0f));
        ntn::rotate_avx2(rho.data(), kFcRot, y.data(), m);
        for (std::size_t i = 0; i < m; ++i) {
            const int64_t n = s + static_cast<int64_t>(i);
            const double bn = b * static_cast<double>(n);
            const double ref = -2.0 * kPiD * (fa + (bn - std::nearbyint(bn)));
            const double err = std::abs(wrap_pi(std::arg(cd(y[i].real(), y[i].imag())) - ref));
            st.max_err = std::max(st.max_err, err);
            if (n < 100'000) st.max_err_head = std::max(st.max_err_head, err);
            if (n >= count - 100'000) st.max_err_tail = std::max(st.max_err_tail, err);
            st.max_mag_err = std::max(st.max_mag_err, std::abs(std::abs(cd(y[i].real(), y[i].imag())) - 1.0));
        }
    }
    return st;
}

void rot_tone_test(int64_t count) {
    // Same cases as M1a / M1b (python/tests/test_golden_geometric.py TONE_CASES).
    for (const auto [rho0, v] : {std::pair{1'000'000.0, -7000.0}, std::pair{1'234'567.891, 6853.37}}) {
        const RotPhaseStats st = rot_tone_stats(rho0, v, count);
        std::printf("[ measure  ] AVX2 rotator tone rho0=%.9g v=%.6g N=%lld: max phase err %.3e rad, first 1e5 %.3e, "
                    "last 1e5 %.3e, max ||y|-1| %.2e\n",
                    rho0, v, static_cast<long long>(count), st.max_err, st.max_err_head, st.max_err_tail,
                    st.max_mag_err);
        EXPECT_LT(st.max_err, 2e-6);
        EXPECT_LT(st.max_err_tail, 2e-6);  // no drift: the end of the run meets the same limit
        EXPECT_LT(st.max_mag_err, 1e-6);
    }
}

// Exact carrier-cycle cases (as tests/cpp/test_rotator.cpp): fc = c makes fc*rho/c = rho exactly.
void expect_rotated_avx2(double rho, cf32 factor) {
    const std::vector<cf32> x = {{1.0f, 0.0f}, {0.0f, 1.0f}, {0.6f, -0.8f}, {-3.25f, 1.5f}};
    auto y = x;
    std::vector<double> r(y.size(), rho);
    ntn::rotate_avx2(r.data(), kCLight, y.data(), y.size());
    for (std::size_t i = 0; i < x.size(); ++i) {
        EXPECT_LT(std::abs(y[i] - x[i] * factor), 2e-7f * std::abs(x[i])) << "rho = " << rho << ", i = " << i;
    }
}

// sincos_avx2 accuracy on float arguments `a`, against sin/cos of the same float in float64.
struct SincosErr {
    double sin_err = 0.0, cos_err = 0.0, cplx_err = 0.0, phase_err = 0.0, mag_err = 0.0;
    float worst_arg = 0.0f;
    void add(const SincosErr& o) {
        if (o.cplx_err > cplx_err) worst_arg = o.worst_arg;
        sin_err = std::max(sin_err, o.sin_err);
        cos_err = std::max(cos_err, o.cos_err);
        cplx_err = std::max(cplx_err, o.cplx_err);
        phase_err = std::max(phase_err, o.phase_err);
        mag_err = std::max(mag_err, o.mag_err);
    }
};

SincosErr sincos_errors(const std::vector<float>& a, bool use_avx2) {
    std::vector<float> s(a.size()), c(a.size());
    if (use_avx2) {
        ntn::sincos_avx2(a.data(), s.data(), c.data(), a.size());
    } else {
        for (std::size_t i = 0; i < a.size(); ++i) {
            s[i] = std::sin(a[i]);
            c[i] = std::cos(a[i]);
        }
    }
    SincosErr e;
    for (std::size_t i = 0; i < a.size(); ++i) {
        const double x = a[i];
        const double rs = std::sin(x), rc = std::cos(x);
        const double es = std::abs(s[i] - rs), ec = std::abs(c[i] - rc);
        const double ex = std::hypot(s[i] - rs, c[i] - rc);
        e.sin_err = std::max(e.sin_err, es);
        e.cos_err = std::max(e.cos_err, ec);
        if (ex > e.cplx_err) {
            e.cplx_err = ex;
            e.worst_arg = a[i];
        }
        e.phase_err = std::max(e.phase_err, std::abs(wrap_pi(std::atan2(double(s[i]), double(c[i])) - x)));
        e.mag_err = std::max(e.mag_err, std::abs(std::hypot(double(s[i]), double(c[i])) - 1.0));
    }
    return e;
}

void print_sincos(const char* label, const SincosErr& e, std::size_t n) {
    std::printf("[ measure  ] %s, %zu args in [-pi_f, pi_f]: max |sin err| %.3e, |cos err| %.3e, "
                "|e^ja err| %.3e (at a = %.9g), phase %.3e rad, ||.|-1| %.3e\n",
                label, n, e.sin_err, e.cos_err, e.cplx_err, double(e.worst_arg), e.phase_err, e.mag_err);
}

// Limit for the vectorised sincos alone: per-sample |(c + js) - e^{ja}| < 1e-6, half of the
// 2e-6 rad tone-phase budget (the float32 angle rounding and the complex multiply take the rest).
constexpr double kSincosLimit = 1e-6;

}  // namespace

TEST(RotatorAvx2, ExactCycles) {
    NTN_REQUIRE_AVX2();
    expect_rotated_avx2(1000000.25, {0.0f, -1.0f});
    expect_rotated_avx2(1000000.75, {0.0f, 1.0f});
    expect_rotated_avx2(1000000.5, {-1.0f, 0.0f});
    expect_rotated_avx2(1000001.5, {-1.0f, 0.0f});
    // Integer cycles and fc = 0: identity, bit-exact.
    const std::vector<cf32> x = {{1.0f, 0.0f}, {0.0f, 1.0f}, {0.6f, -0.8f}, {-3.25f, 1.5f}};
    for (const auto [rho, fc] : {std::pair{1234567.0, kCLight}, std::pair{987654.321, 0.0}}) {
        auto y = x;
        std::vector<double> r(y.size(), rho);
        ntn::rotate_avx2(r.data(), fc, y.data(), y.size());
        EXPECT_EQ(std::memcmp(x.data(), y.data(), x.size() * sizeof(cf32)), 0) << rho << ", fc " << fc;
    }
}

TEST(RotatorAvx2, IsNotTheScalarRotator) {
    // The vectorised sincos differs from the CRT's in some bits; forwarding to rotate_scalar is
    // bit-identical and fails here.
    NTN_REQUIRE_AVX2();
    const RotCase r = run_rot_case(4096, 0, "not-scalar");
    std::printf("[ measure  ] AVX2 vs scalar rotator, 4096 samples: %zu outputs differ in bits, error %.1f dB\n",
                r.differing_bits, r.err_db);
    EXPECT_GT(r.differing_bits, 0u);
    EXPECT_LT(r.err_db, -100.0);
}

TEST(RotatorAvx2, MatchesScalarLengthsAlignment) {
    NTN_REQUIRE_AVX2();
    double worst = -400.0;
    for (std::size_t offset : {std::size_t{0}, std::size_t{8}}) {
        for (std::size_t len : {0, 1, 7, 8, 9, 31, 33, 4095, 4097}) {
            char label[64];
            std::snprintf(label, sizeof label, "len %zu, offset %zu B", len, offset);
            const RotCase r = run_rot_case(len, offset, label);
            if (len) {
                EXPECT_LT(r.err_db, -100.0) << label;
                worst = std::max(worst, r.err_db);
            }
        }
    }
    std::printf("[ measure  ] AVX2 vs scalar rotator: worst %.1f dB over lengths 0..4097, offsets 0/8 B "
                "(limit -100 dB)\n",
                worst);
}

TEST(RotatorAvx2, BlockSizeInvariantBitExact) {
    // Each sample's result must not depend on where the call boundaries fall (padded tail).
    NTN_REQUIRE_AVX2();
    constexpr std::size_t kN = 10007;
    const auto rho = realistic_rho(kN, 9);
    std::mt19937_64 rng(10);
    std::normal_distribution<float> g(0.0f, 1.0f);
    std::vector<cf32> x(kN);
    for (auto& v : x) v = {g(rng), g(rng)};
    auto whole = x;
    ntn::rotate_avx2(rho.data(), kFcRot, whole.data(), kN);
    for (std::size_t block : {1, 7, 8, 9, 4096}) {
        auto y = x;
        for (std::size_t s = 0; s < kN; s += block) {
            ntn::rotate_avx2(rho.data() + s, kFcRot, y.data() + s, std::min(block, kN - s));
        }
        EXPECT_EQ(std::memcmp(whole.data(), y.data(), kN * sizeof(cf32)), 0) << "block " << block;
    }
}

TEST(RotatorAvx2, TonePhase) {
    NTN_REQUIRE_AVX2();
    rot_tone_test(1'000'000);
}

TEST(RotatorAvx2, TonePhaseSlow) {
    NTN_REQUIRE_SLOW();
    NTN_REQUIRE_AVX2();
    rot_tone_test(10'000'000);
}

TEST(SincosAvx2, AccuracyOverWrappedRange) {
    // The rotator's float32 argument is -2 pi frac with frac in [-0.5, 0.5], i.e. within
    // [-pi_f, pi_f], pi_f = float(pi). Fast: 2^22 + 1 evenly spaced arguments plus special
    // values; slow: every float in the range (AccuracyExhaustiveSlow).
    NTN_REQUIRE_AVX2();
    const float pi_f = static_cast<float>(kPiD);
    std::vector<float> a;
    constexpr int kGrid = 1 << 22;
    for (int i = 0; i <= kGrid; ++i) a.push_back(static_cast<float>(-kPiD + 2.0 * kPiD * i / kGrid));
    for (float v : {0.0f, -0.0f, pi_f, -pi_f, std::nextafter(pi_f, 0.0f), -std::nextafter(pi_f, 0.0f),
                    std::numeric_limits<float>::min(), std::numeric_limits<float>::denorm_min(), 1e-30f,
                    static_cast<float>(kPiD / 2), -static_cast<float>(kPiD / 2), static_cast<float>(kPiD / 4)}) {
        a.push_back(v);
    }
    for (auto& v : a) v = std::clamp(v, -pi_f, pi_f);
    const SincosErr ev = sincos_errors(a, true), ec = sincos_errors(a, false);
    print_sincos("sincos_avx2 (xsimd)", ev, a.size());
    print_sincos("CRT sinf/cosf (scalar rotator, for reference)", ec, a.size());
    EXPECT_LT(ev.cplx_err, kSincosLimit);
    EXPECT_LT(ev.phase_err, kSincosLimit);
}

TEST(SincosAvx2, AccuracyExhaustiveSlow) {
    NTN_REQUIRE_SLOW();
    NTN_REQUIRE_AVX2();
    // Every float32 in [-pi_f, pi_f]: bit patterns 0 .. bits(pi_f), both signs.
    const float pi_f = static_cast<float>(kPiD);
    uint32_t top;
    std::memcpy(&top, &pi_f, sizeof top);
    SincosErr total;
    std::size_t n = 0;
    std::vector<float> a;
    a.reserve(1 << 20);
    for (uint32_t sign : {0u, 0x8000'0000u}) {
        for (uint64_t u = 0; u <= top; u += (1u << 20)) {
            a.clear();
            const uint64_t end = std::min<uint64_t>(u + (1u << 20), uint64_t{top} + 1);
            for (uint64_t k = u; k < end; ++k) {
                const uint32_t bits = static_cast<uint32_t>(k) | sign;
                float f;
                std::memcpy(&f, &bits, sizeof f);
                a.push_back(f);
            }
            total.add(sincos_errors(a, true));
            n += a.size();
        }
    }
    print_sincos("sincos_avx2 (xsimd), exhaustive", total, n);
    EXPECT_LT(total.cplx_err, kSincosLimit);
    EXPECT_LT(total.phase_err, kSincosLimit);
}
