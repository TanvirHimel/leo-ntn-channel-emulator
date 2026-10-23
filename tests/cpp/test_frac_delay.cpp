// Scalar fractional-delay kernel: table structure, exact integer delays, and accuracy against
// the exact FFT-domain delay of band-limited noise.
#include "ntn/frac_delay.hpp"

#include <gtest/gtest.h>

#include <cmath>
#include <cstdio>
#include <cstring>
#include <numbers>
#include <random>
#include <vector>

#include "test_util.hpp"

using ntn::cf32;
using ntn_test::cd;
namespace fd = ntn::fd;

namespace {

std::vector<cf32> random_cf32(std::size_t n, uint64_t seed) {
    std::mt19937_64 rng(seed);
    std::normal_distribution<float> g(0.0f, 1.0f);
    std::vector<cf32> x(n);
    for (auto& v : x) v = {g(rng), g(rng)};
    return x;
}

bool bit_equal(cf32 a, cf32 b) { return std::memcmp(&a, &b, sizeof(cf32)) == 0; }

}  // namespace

// ---- Table structure ------------------------------------------------------------------------

TEST(FracDelayTable, PhaseZeroIsExactUnitImpulse) {
    const float* r0 = ntn::frac_delay_table().row(0);
    for (int j = 0; j < fd::kTaps; ++j) {
        EXPECT_EQ(r0[j], j == fd::kHalf - 1 ? 1.0f : 0.0f) << "tap " << j;
    }
}

TEST(FracDelayTable, LastRowIsImpulseShiftedByOne) {
    const float* r = ntn::frac_delay_table().row(fd::kPhases);
    for (int j = 0; j < fd::kTaps; ++j) {
        EXPECT_EQ(r[j], j == fd::kHalf ? 1.0f : 0.0f) << "tap " << j;
    }
}

TEST(FracDelayTable, HalfSampleRowIsSymmetric) {
    // f = 0.5: h(0.5 - off) is even about the centre, so tap j equals tap 31 - j.
    const float* r = ntn::frac_delay_table().row(fd::kPhases / 2);
    for (int j = 0; j < fd::kTaps / 2; ++j) EXPECT_EQ(r[j], r[fd::kTaps - 1 - j]) << "tap " << j;
}

TEST(FracDelayTable, DcGainNearOneEveryPhase) {
    // Sanity bound (-80 dB); the band-limited accuracy test below is the real criterion.
    for (int r = 0; r <= fd::kPhases; ++r) {
        const float* row = ntn::frac_delay_table().row(r);
        double s = 0.0;
        for (int j = 0; j < fd::kTaps; ++j) s += row[j];
        EXPECT_NEAR(s, 1.0, 1e-4) << "row " << r;
    }
}

// ---- Integer delays are bit-exact -------------------------------------------------------------

TEST(FracDelayKernel, IntegerDelayIsBitExact) {
    // Delays within kSnap (1e-9) of an integer snap to it (golden rule) and are a pure gather.
    constexpr int64_t kPad = 64;
    constexpr std::size_t kN = 1000;
    const auto hist = random_cf32(kN + 2 * kPad, 11);
    const int64_t hist_first = -kPad;
    for (double d : {0.0, 1.0, 3.0, -2.0, 37.0, 3.0 + 1e-10, 3.0 - 1e-10, -2.0 + 5e-10}) {
        const double k = std::round(d);
        std::vector<double> dv(kN, d);
        std::vector<cf32> out(kN);
        ntn::frac_delay_scalar(ntn::frac_delay_table(), hist.data(), hist_first, 0, dv.data(),
                               out.data(), kN);
        for (std::size_t n = 0; n < kN; ++n) {
            const int64_t m = static_cast<int64_t>(n) - static_cast<int64_t>(k);
            ASSERT_TRUE(bit_equal(out[n], hist[static_cast<std::size_t>(m - hist_first)]))
                << "d = " << d << ", n = " << n;
        }
    }
}

// ---- Accuracy: band-limited noise vs the exact FFT-domain delay -------------------------------

TEST(FracDelayKernel, FftSelfCheck) {
    // The helper FFT on a single complex tone at bin 5 gives N at bin 5, ~0 elsewhere.
    constexpr std::size_t kN = 4096;
    std::vector<cd> a(kN);
    for (std::size_t n = 0; n < kN; ++n) {
        const double ang = 2.0 * std::numbers::pi * 5.0 * static_cast<double>(n) / kN;
        a[n] = {std::cos(ang), std::sin(ang)};
    }
    ntn_test::fft_inplace(a, -1);
    for (std::size_t k = 0; k < kN; ++k) {
        ASSERT_NEAR(std::abs(a[k] - cd(k == 5 ? double(kN) : 0.0, 0.0)), 0.0, 1e-9) << "k = " << k;
    }
}

TEST(FracDelayKernel, BandLimitedNoiseMatchesFftDelay) {
    // Periodic band-limited noise: random complex Gaussian bins for |k| <= 0.4 N, zero above.
    // Reference for delay d: ifft(X[k] exp(-j 2 pi k d / N)), exact for a periodic signal.
    // The kernel reads a periodic extension, so there are no edges. Limit: -90 dB.
    constexpr std::size_t kN = 4096;
    constexpr int64_t kPad = 64;
    const int kMaxBin = static_cast<int>(0.4 * kN);  // 1638 -> |f| <= 0.39990 fs

    std::mt19937_64 rng(2024);
    std::normal_distribution<double> g(0.0, 1.0);
    std::vector<cd> X(kN, cd(0.0, 0.0));
    for (int k = -kMaxBin; k <= kMaxBin; ++k) X[static_cast<std::size_t>((k + kN) % kN)] = {g(rng), g(rng)};

    std::vector<cd> x = X;
    ntn_test::fft_inplace(x, +1);
    for (auto& v : x) v /= static_cast<double>(kN);

    std::vector<cf32> hist(kN + 2 * kPad);
    for (int64_t m = -kPad; m < static_cast<int64_t>(kN) + kPad; ++m) {
        const cd v = x[static_cast<std::size_t>((m % static_cast<int64_t>(kN) + kN) % kN)];
        hist[static_cast<std::size_t>(m + kPad)] = {static_cast<float>(v.real()), static_cast<float>(v.imag())};
    }

    double worst = -1e9;
    for (double d : {0.1, 0.25, 0.5, 0.73, 0.95, 5.3, -2.4, 0.999, 0.5 / 256.0, 0.0013}) {
        std::vector<cd> Y(kN);
        for (std::size_t i = 0; i < kN; ++i) {
            const double k = (i <= kN / 2) ? static_cast<double>(i) : static_cast<double>(i) - kN;
            const double ang = -2.0 * std::numbers::pi * k * d / static_cast<double>(kN);
            Y[i] = X[i] * cd(std::cos(ang), std::sin(ang));
        }
        ntn_test::fft_inplace(Y, +1);
        for (auto& v : Y) v /= static_cast<double>(kN);

        std::vector<double> dv(kN, d);
        std::vector<cf32> out(kN);
        ntn::frac_delay_scalar(ntn::frac_delay_table(), hist.data(), -kPad, 0, dv.data(), out.data(), kN);
        const double e = ntn_test::error_db<cf32, cd>(out, Y);
        std::printf("[ measure  ] frac delay d = %-9.6g error %.1f dB\n", d, e);
        worst = std::max(worst, e);
        EXPECT_LT(e, -90.0) << "d = " << d;
    }
    std::printf("[ measure  ] frac delay worst %.1f dB (limit -90 dB)\n", worst);
}
