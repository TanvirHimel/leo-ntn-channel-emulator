// M2 acceptance test 2: slip continuity. A linear delay ramp crosses >= 20 integer-sample
// boundaries ("slips"); the output must match the ideal delayed signal to < -60 dB overall and
// within +-64 samples of every slip.
//
// Stimulus: a sum of K random complex tones with |f| < 0.4 fs, so the ideal output
// x(n - D(n)) has a closed form, evaluated per sample in float64 at the exact delay.
#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <numbers>
#include <random>
#include <span>
#include <vector>

#include "ntn/frac_delay.hpp"
#include "ntn/geometric_channel.hpp"
#include "test_util.hpp"

using ntn::cf32;
using ntn_test::cd;

namespace {

constexpr double kC = 299'792'458.0;
constexpr double kFs = 7'680'000.0;
constexpr int64_t kWin = 64;  // samples either side of a slip

// x(t) = sum_k a_k exp(j 2 pi f_k t), f_k uniform in [-0.4, 0.4) cycles/sample, total power ~1.
struct ToneSum {
    std::vector<double> f;
    std::vector<cd> a;
    ToneSum(int k, uint64_t seed) {
        std::mt19937_64 rng(seed);
        std::uniform_real_distribution<double> uf(-0.4, 0.4);
        std::normal_distribution<double> g(0.0, 1.0 / std::sqrt(2.0 * k));
        for (int i = 0; i < k; ++i) {
            f.push_back(uf(rng));
            a.emplace_back(g(rng), g(rng));
        }
    }
    cd operator()(double t) const {
        cd s(0.0, 0.0);
        for (std::size_t i = 0; i < f.size(); ++i) {
            const double ph = 2.0 * std::numbers::pi * f[i] * t;
            s += a[i] * cd(std::cos(ph), std::sin(ph));
        }
        return s;
    }
};

struct SlipResult {
    double overall_db;
    double worst_window_db;
    std::size_t slips;
};

// Compares y[i] with ref[i] for i in [0, y.size()); delay[i] is D at that sample. Slips are the
// i where floor(D) changes; each needs a full +-kWin window inside the range.
SlipResult evaluate(std::span<const cf32> y, std::span<const cd> ref, std::span<const double> delay,
                    const char* label) {
    SlipResult r{ntn_test::error_db<cf32, cd>(y, ref), -1e9, 0};
    for (std::size_t i = 1; i < delay.size(); ++i) {
        if (std::floor(delay[i]) == std::floor(delay[i - 1])) continue;
        ++r.slips;
        EXPECT_GE(i, static_cast<std::size_t>(kWin)) << label << ": slip too close to the start";
        EXPECT_LE(i + kWin + 1, delay.size()) << label << ": slip too close to the end";
        if (i < static_cast<std::size_t>(kWin) || i + kWin + 1 > delay.size()) continue;
        const std::size_t a = i - kWin, len = 2 * kWin + 1;
        const double e = ntn_test::error_db<cf32, cd>(y.subspan(a, len), ref.subspan(a, len));
        r.worst_window_db = std::max(r.worst_window_db, e);
        EXPECT_LT(e, -60.0) << label << ": window around slip at " << i << " (D " << delay[i - 1] << " -> "
                            << delay[i] << ")";
    }
    return r;
}

// Kernel-level ramp D(n) = d0 + slope n, n in [0, n_out), through design <T, P>.
template <int T, int P>
SlipResult kernel_ramp(double d0, double slope, int64_t n_out, const char* label, bool assert_limits = true) {
    constexpr int64_t kPad = 64;  // covers |D| <= ~40 with kHalf <= 16
    const ToneSum x(32, 99);
    std::vector<cf32> hist(static_cast<std::size_t>(n_out + 2 * kPad));
    for (int64_t m = -kPad; m < n_out + kPad; ++m) {
        const cd v = x(static_cast<double>(m));
        hist[static_cast<std::size_t>(m + kPad)] = {static_cast<float>(v.real()), static_cast<float>(v.imag())};
    }
    std::vector<double> d(static_cast<std::size_t>(n_out));
    std::vector<cd> ref(d.size());
    for (int64_t n = 0; n < n_out; ++n) {
        d[static_cast<std::size_t>(n)] = d0 + slope * static_cast<double>(n);
        ref[static_cast<std::size_t>(n)] = x(static_cast<double>(n) - d[static_cast<std::size_t>(n)]);
    }
    EXPECT_LE(std::max(std::abs(d.front()), std::abs(d.back())) + T / 2 + 1, kPad) << label;
    std::vector<cf32> y(d.size());
    ntn::frac_delay_scalar_t<T, P>(ntn::frac_delay_table_t<T, P>(), hist.data(), -kPad, 0, d.data(), y.data(),
                                   y.size());
    if (!assert_limits) {
        return {ntn_test::error_db<cf32, cd>(y, ref), 0.0, 0};
    }
    const SlipResult r = evaluate(y, ref, d, label);
    std::printf("[ measure  ] %s: %zu slips, overall %.1f dB, worst +-%lld window %.1f dB (limit -60 dB)\n", label,
                r.slips, r.overall_db, static_cast<long long>(kWin), r.worst_window_db);
    EXPECT_GE(r.slips, 20u) << label;
    EXPECT_LT(r.overall_db, -60.0) << label;
    return r;
}

}  // namespace

TEST(SlipContinuity, FastRampKernelIncreasingDelay) {
    // D: 0.3 -> 25.3 samples over 1e5 samples (25 slips), production design.
    kernel_ramp<ntn::fd::kTaps, ntn::fd::kPhases>(0.3, 2.5e-4, 100'000, "kernel ramp +2.5e-4");
}

TEST(SlipContinuity, FastRampKernelDecreasingDelay) {
    // D: 0.3 -> -24.7 samples over 1e5 samples (25 slips), production design.
    kernel_ramp<ntn::fd::kTaps, ntn::fd::kPhases>(0.3, -2.5e-4, 100'000, "kernel ramp -2.5e-4");
}

TEST(SlipContinuity, NegativeControlShortDesignFails) {
    // The metric must be able to fail: the 16-tap design's Kaiser formula attenuation is ~51 dB,
    // so the same ramp cannot reach -60 dB with it.
    const SlipResult r = kernel_ramp<16, 128>(0.3, 2.5e-4, 100'000, "control 16x128", false);
    std::printf("[ measure  ] negative control, 16 x 128 design, same ramp: overall %.1f dB (must be > -60 dB)\n",
                r.overall_db);
    EXPECT_GT(r.overall_db, -60.0);
}

TEST(SlipContinuity, RealisticRampThroughGeometricChannel) {
    // Linear rho = rho0 + v t with v = -7000 m/s: D drifts at v/c = -23.3 ppm, one slip per
    // ~42.8k samples, so 1e6 samples give >= 20 slips. Doppler off, delay on, Kernel::Scalar.
    constexpr int64_t kN = 1'000'000;
    constexpr double kRho0 = 1'000'000.0, kV = -7000.0;
    constexpr int64_t kTau0 = 25'600;  // rho0 fs / c = 25617.6 -> D from 17.6 down to -5.7
    ntn::GeoConfig g;
    g.fc_hz = 2e9;
    g.fs_hz = kFs;
    g.tau0_samples = kTau0;
    g.delay = true;
    g.doppler = false;
    ntn::GeometricChannel ch(ntn::RangeTrack(kFs, {{0, kN, {kRho0, kV, 0.0, 0.0}}}), g, ntn::Kernel::Scalar);
    const int64_t L = ch.latency();

    const ToneSum x(32, 4242);
    std::vector<cf32> in(static_cast<std::size_t>(kN));
    for (int64_t m = 0; m < kN; ++m) {
        const cd v = x(static_cast<double>(m));
        in[static_cast<std::size_t>(m)] = {static_cast<float>(v.real()), static_cast<float>(v.imag())};
    }
    std::vector<cf32> out(in.size());
    for (std::size_t s = 0; s < in.size(); s += 65536) {
        const std::size_t e = std::min(in.size(), s + 65536);
        ch.process(std::span<const cf32>(in.data() + s, e - s), std::span<cf32>(out.data() + s, e - s));
    }

    // out[i] is channel time n = i - L. Start where every read index is >= 0 (input before the
    // track start is zero): n - floor(D_max) - kHalf >= 0 with D_max < 18.
    const int64_t n_lo = 64, n_hi = kN - L;
    std::vector<double> d;
    std::vector<cd> ref;
    for (int64_t n = n_lo; n < n_hi; ++n) {
        const double rho = kRho0 + kV * (static_cast<double>(n) / kFs);
        d.push_back(rho * kFs / kC - static_cast<double>(kTau0));
        ref.push_back(x(static_cast<double>(n) - d.back()));
    }
    const std::span<const cf32> y(out.data() + (n_lo + L), d.size());
    const SlipResult r = evaluate(y, ref, d, "GeometricChannel ramp");
    std::printf("[ measure  ] GeometricChannel ramp v=-7000 m/s (%.2f ppm), L = %lld: %zu slips, overall %.1f dB, "
                "worst +-%lld window %.1f dB (limit -60 dB)\n",
                kV / kC * 1e6, static_cast<long long>(L), r.slips, r.overall_db, static_cast<long long>(kWin),
                r.worst_window_db);
    EXPECT_GE(r.slips, 20u);
    EXPECT_LT(r.overall_db, -60.0);
}
