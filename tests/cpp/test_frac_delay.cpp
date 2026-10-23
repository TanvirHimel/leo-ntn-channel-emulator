// Scalar fractional-delay kernel: table structure, exact integer delays, and accuracy against
// the exact FFT-domain delay of band-limited noise.
#include "ntn/frac_delay.hpp"

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <numbers>
#include <random>
#include <type_traits>
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

// ---- M2: final polyphase design ---------------------------------------------------------------
//
// The tap count and phase count are compile-time parameters; six designs (16/24/32 taps x
// 128/256 phases) are explicitly instantiated in frac_delay.cpp. Kaiser beta comes from the
// window design formula for a 0.4 -> 0.6 fs transition; nothing is tuned to these tests.

namespace {

// Periodic band-limited noise (random complex Gaussian bins for |k| <= 0.4 N, zero above) and
// the kernel's input as its periodic extension, so there are no edges. The reference for delay
// d is ifft(X[k] exp(-j 2 pi k d / N)), exact for a periodic signal.
struct PeriodicNoise {
    static constexpr std::size_t kN = 4096;
    static constexpr int64_t kPad = 64;
    std::vector<cd> X;
    std::vector<cf32> hist;  // x[m] for m in [-kPad, kN + kPad)

    explicit PeriodicNoise(uint64_t seed) : X(kN, cd(0.0, 0.0)), hist(kN + 2 * kPad) {
        const int max_bin = static_cast<int>(0.4 * kN);
        std::mt19937_64 rng(seed);
        std::normal_distribution<double> g(0.0, 1.0);
        for (int k = -max_bin; k <= max_bin; ++k) X[static_cast<std::size_t>((k + kN) % kN)] = {g(rng), g(rng)};
        std::vector<cd> x = X;
        ntn_test::fft_inplace(x, +1);
        for (int64_t m = -kPad; m < static_cast<int64_t>(kN) + kPad; ++m) {
            const cd v = x[static_cast<std::size_t>((m % static_cast<int64_t>(kN) + kN) % kN)] / double(kN);
            hist[static_cast<std::size_t>(m + kPad)] = {static_cast<float>(v.real()), static_cast<float>(v.imag())};
        }
    }

    std::vector<cd> delayed(double d) const {
        std::vector<cd> Y(kN);
        for (std::size_t i = 0; i < kN; ++i) {
            const double k = (i <= kN / 2) ? static_cast<double>(i) : static_cast<double>(i) - kN;
            const double ang = -2.0 * std::numbers::pi * k * d / static_cast<double>(kN);
            Y[i] = X[i] * cd(std::cos(ang), std::sin(ang));
        }
        ntn_test::fft_inplace(Y, +1);
        for (auto& v : Y) v /= static_cast<double>(kN);
        return Y;
    }
};

// Error of design <T, P> on the noise at constant delay d, dB relative to the reference.
template <int T, int P>
double noise_error_db(const PeriodicNoise& s, double d) {
    std::vector<double> dv(PeriodicNoise::kN, d);
    std::vector<cf32> out(PeriodicNoise::kN);
    ntn::frac_delay_scalar_t<T, P>(ntn::frac_delay_table_t<T, P>(), s.hist.data(), -PeriodicNoise::kPad, 0,
                                   dv.data(), out.data(), PeriodicNoise::kN);
    const auto ref = s.delayed(d);
    return ntn_test::error_db<cf32, cd>(out, ref);
}

// Peak over |f| <= 0.4 fs (grid step 0.001 fs) and over `delays` of |H(f) - exp(-j 2 pi f d)|,
// measured through the kernel: input x[m] = exp(j 2 pi f m) (unit magnitude), one output at
// n = 0, whose ideal value is exp(-j 2 pi f d).
template <int T, int P>
double peak_response_error(const std::vector<double>& delays) {
    constexpr int64_t kPad = 40;  // covers |d| <= 8 with kHalf <= 16
    std::vector<cf32> hist(2 * kPad + 1);
    double peak = 0.0;
    for (int i = -400; i <= 400; ++i) {
        const double f = i * 0.001;
        for (int64_t m = -kPad; m <= kPad; ++m) {
            const double a = 2.0 * std::numbers::pi * f * static_cast<double>(m);
            hist[static_cast<std::size_t>(m + kPad)] = {static_cast<float>(std::cos(a)), static_cast<float>(std::sin(a))};
        }
        for (double d : delays) {
            cf32 y;
            ntn::frac_delay_scalar_t<T, P>(ntn::frac_delay_table_t<T, P>(), hist.data(), -kPad, 0, &d, &y, 1);
            const double a = -2.0 * std::numbers::pi * f * d;
            peak = std::max(peak, std::abs(cd(y.real(), y.imag()) - cd(std::cos(a), std::sin(a))));
        }
    }
    return peak;
}

// ~200 delays: a grid of 100 fractions (i + 0.5)/100 and 100 seeded random delays in [-8, 8).
std::vector<double> sweep_delays() {
    std::vector<double> d;
    for (int i = 0; i < 100; ++i) d.push_back((i + 0.5) / 100.0);
    std::mt19937_64 rng(7);
    std::uniform_real_distribution<double> u(-8.0, 8.0);
    for (int i = 0; i < 100; ++i) d.push_back(u(rng));
    return d;
}

template <int T, int P>
void table_structure_checks() {
    const auto& t = ntn::frac_delay_table_t<T, P>();
    constexpr int half = T / 2;
    for (int j = 0; j < T; ++j) {
        EXPECT_EQ(t.row(0)[j], j == half - 1 ? 1.0f : 0.0f) << T << "x" << P << " row 0 tap " << j;
        EXPECT_EQ(t.row(P)[j], j == half ? 1.0f : 0.0f) << T << "x" << P << " last row tap " << j;
    }
    for (int j = 0; j < T / 2; ++j) {
        EXPECT_EQ(t.row(P / 2)[j], t.row(P / 2)[T - 1 - j]) << T << "x" << P << " half row tap " << j;
    }
}

struct SweepRow {
    int taps, phases;
    double beta, atten_db, worst_db, median_db, peak;
};

template <int T, int P>
SweepRow sweep_row(const PeriodicNoise& s, const std::vector<double>& delays) {
    std::vector<double> e;
    for (double d : delays) e.push_back(noise_error_db<T, P>(s, d));
    std::sort(e.begin(), e.end());
    const double median = (e[e.size() / 2 - 1] + e[e.size() / 2]) / 2.0;
    return {T, P, ntn::fd::kaiser_beta(T), ntn::fd::kaiser_design_atten_db(T), e.back(), median,
            peak_response_error<T, P>(delays)};
}

}  // namespace

TEST(FracDelayDesign, ProductionDesignIs32x256) {
    static_assert(fd::kTaps == 32 && fd::kPhases == 256 && fd::kHalf == 16);
    static_assert(std::is_same_v<ntn::FracDelayTable, ntn::FracDelayTableT<32, 256>>);
    const auto* production = &ntn::frac_delay_table_t<32, 256>();
    EXPECT_EQ(&ntn::frac_delay_table(), production);
}

TEST(FracDelayDesign, KaiserBetaFromDesignFormula) {
    // Kaiser window design formula, transition 0.4 -> 0.6 fs (dw = 0.4 pi), order M = taps - 1:
    //   A = 2.285 M dw + 8 dB;  beta = 0.1102 (A - 8.7) for A > 50,
    //   0.5842 (A - 21)^0.4 + 0.07886 (A - 21) for 21 <= A <= 50, 0 below.
    const double dw = 0.4 * std::numbers::pi;
    for (int taps : {16, 24, 32}) {
        const double a = 2.285 * (taps - 1) * dw + 8.0;
        const double beta = a > 50.0   ? 0.1102 * (a - 8.7)
                            : a >= 21.0 ? 0.5842 * std::pow(a - 21.0, 0.4) + 0.07886 * (a - 21.0)
                                        : 0.0;
        EXPECT_NEAR(ntn::fd::kaiser_design_atten_db(taps), a, 1e-12) << taps;
        EXPECT_NEAR(ntn::fd::kaiser_beta(taps), beta, 1e-12) << taps;
    }
    EXPECT_EQ(ntn::frac_delay_table().beta(), ntn::fd::kaiser_beta(32));
}

TEST(FracDelayDesign, TableStructureAllDesigns) {
    table_structure_checks<16, 128>();
    table_structure_checks<16, 256>();
    table_structure_checks<24, 128>();
    table_structure_checks<24, 256>();
    table_structure_checks<32, 128>();
    table_structure_checks<32, 256>();
}

TEST(FracDelayDesign, Sweep) {
    // Informational design table (M2 step 2). Asserts only that the production design meets the
    // SPEC target (-60 dB over |f| < 0.4 fs); the choice of design is made from the printed table.
    const PeriodicNoise s(2025);
    const auto delays = sweep_delays();
    const SweepRow rows[] = {
        sweep_row<16, 128>(s, delays), sweep_row<16, 256>(s, delays), sweep_row<24, 128>(s, delays),
        sweep_row<24, 256>(s, delays), sweep_row<32, 128>(s, delays), sweep_row<32, 256>(s, delays),
    };
    std::printf("[ measure  ] design sweep: %zu delays, band-limited periodic noise |f| <= 0.4 fs (N = %zu); "
                "peak |H(f) - exp(-j2pi f d)| over |f| <= 0.4 fs, step 0.001 fs\n",
                delays.size(), PeriodicNoise::kN);
    std::printf("[ measure  ] | taps | phases | Kaiser beta | design A (dB) | worst (dB) | median (dB) | peak |H-ideal| | peak (dB) |\n");
    for (const auto& r : rows) {
        std::printf("[ measure  ] | %d | %d | %.4f | %.2f | %.1f | %.1f | %.3e | %.1f |\n", r.taps, r.phases, r.beta,
                    r.atten_db, r.worst_db, r.median_db, r.peak, 20.0 * std::log10(r.peak));
    }
    EXPECT_LT(rows[5].worst_db, -60.0) << "production design 32 x 256";
}

TEST(FracDelayM2, RandomDelaysBandLimitedNoise) {
    // M2 acceptance test 1: 50 seeded random fractional delays in [-8, 8), band-limited noise,
    // production design; error against the FFT-domain reference < -60 dB for each.
    const PeriodicNoise s(31337);
    std::mt19937_64 rng(20261023);
    std::uniform_real_distribution<double> u(-8.0, 8.0);
    double worst = -1e9;
    for (int i = 0; i < 50; ++i) {
        const double d = u(rng);
        const double e = noise_error_db<fd::kTaps, fd::kPhases>(s, d);
        worst = std::max(worst, e);
        EXPECT_LT(e, -60.0) << "d = " << d;
    }
    std::printf("[ measure  ] M2 test 1: 50 random delays, worst %.1f dB (limit -60 dB)\n", worst);
}
