// Geometric channel: M1a acceptance tests 1, 2 and 4 re-run against the same closed-form
// expectations (SPEC §8 M1b test 1), plus latency and a time-varying delay against its closed
// form. Since M2 the closed-form suites run on both kernels (value-parameterized, Scalar and
// Avx2) with unchanged expectations and limits.
#include "ntn/geometric_channel.hpp"

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <numbers>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

#include "ntn/cpu_check.hpp"
#include "test_util.hpp"

using ntn::cf32;
using ntn::GeoConfig;
using ntn::GeometricChannel;
using ntn::RangeSegment;
using ntn::RangeTrack;
using ntn_test::cd;
using ntn::Kernel;

class KernelParam : public ::testing::TestWithParam<Kernel> {};
class GeoTone : public KernelParam {};
class GeoShift : public KernelParam {};
class GeoQuadratic : public KernelParam {};
class GeoLatency : public KernelParam {};
class GeoDelay : public KernelParam {};

namespace {

std::string kernel_name(const ::testing::TestParamInfo<Kernel>& info) {
    return info.param == Kernel::Scalar ? "Scalar" : "Avx2";
}

constexpr double kC = 299'792'458.0;
constexpr double kFs = 7'680'000.0;
constexpr double kFc = 2'000'000'000.0;
constexpr double kPi = std::numbers::pi;

struct ToneCase {
    double rho0, v;  // m, m/s
};
// Same cases as M1a (python/tests/test_golden_geometric.py TONE_CASES).
constexpr ToneCase kToneCases[] = {{1'000'000.0, -7000.0}, {1'234'567.891, 6853.37}};

GeoConfig geo(double fc, int64_t tau0, bool delay, bool doppler) {
    GeoConfig g;
    g.fc_hz = fc;
    g.fs_hz = kFs;
    g.tau0_samples = tau0;
    g.delay = delay;
    g.doppler = doppler;
    return g;
}

// Streams `in` through the channel in fixed blocks and returns the output.
std::vector<cf32> run(GeometricChannel& ch, const std::vector<cf32>& in, std::size_t block = 65536) {
    std::vector<cf32> out(in.size());
    for (std::size_t s = 0; s < in.size(); s += block) {
        const std::size_t e = std::min(in.size(), s + block);
        ch.process(std::span<const cf32>(in.data() + s, e - s), std::span<cf32>(out.data() + s, e - s));
    }
    return out;
}

double wrap_pi(double a) { return a - 2.0 * kPi * std::nearbyint(a / (2.0 * kPi)); }

// ---- Test 1: tone / Doppler phase -------------------------------------------------------------
//
// DC input, linear rho = rho0 + v t, delay off. Exact output phase: -2 pi frac(fc (rho0 + v n/fs) / c).
// Reference, float64, split as  frac(fc rho0 / c) + frac(B n),  B = fc v / (c fs)  [cycles/sample].
// Accuracy of this reference: fc rho0 / c is ~6.7e6..8.2e6 cycles, so its float64 rounding is
// ~1e-9 cycles, i.e. ~6e-9 rad; B n reaches ~6e4 cycles at n = 1e7 with ~2e-11 cycles error.
// The total (~6e-9 rad) is ~300x below the 2e-6 rad limit, so the test measures the C++
// float32 output (float32 angle + sincos + multiply), not the reference.
struct PhaseStats {
    double max_err = 0.0, max_err_tail = 0.0, max_mag_err = 0.0;
};

PhaseStats tone_phase_stats(const ToneCase& tc, int64_t count, Kernel kernel, std::vector<cf32>* keep = nullptr) {
    GeometricChannel ch(RangeTrack(kFs, {{0, count, {tc.rho0, tc.v, 0.0, 0.0}}}), geo(kFc, 0, false, true), kernel);
    EXPECT_EQ(ch.latency(), 0);
    const double a = kFc * tc.rho0 / kC;
    const double fa = a - std::nearbyint(a);
    const double b = kFc * tc.v / (kC * kFs);
    const int64_t tail_from = count - 100'000;

    PhaseStats st;
    constexpr int64_t kBlock = 1 << 20;
    std::vector<cf32> ones(kBlock, cf32(1.0f, 0.0f)), y(kBlock);
    if (keep) keep->resize(static_cast<std::size_t>(count));
    for (int64_t s = 0; s < count; s += kBlock) {
        const auto m = static_cast<std::size_t>(std::min(kBlock, count - s));
        ch.process(std::span<const cf32>(ones.data(), m), std::span<cf32>(y.data(), m));
        for (std::size_t i = 0; i < m; ++i) {
            const int64_t n = s + static_cast<int64_t>(i);
            const double bn = b * static_cast<double>(n);
            const double ref = -2.0 * kPi * (fa + (bn - std::nearbyint(bn)));
            const double err = std::abs(wrap_pi(std::arg(cd(y[i].real(), y[i].imag())) - ref));
            st.max_err = std::max(st.max_err, err);
            if (n >= tail_from) st.max_err_tail = std::max(st.max_err_tail, err);
            st.max_mag_err = std::max(st.max_mag_err, std::abs(std::abs(cd(y[i].real(), y[i].imag())) - 1.0));
            if (keep) (*keep)[static_cast<std::size_t>(n)] = y[i];
        }
    }
    return st;
}

// |DTFT| of y at frequency f (Hz); phase by n*f/fs reduced per term in float64.
double dtft_mag(const std::vector<cf32>& y, double f) {
    double re = 0.0, im = 0.0;
    const double step = f / kFs;
    for (std::size_t n = 0; n < y.size(); ++n) {
        double cyc = step * static_cast<double>(n);
        cyc -= std::nearbyint(cyc);
        const double ang = -2.0 * kPi * cyc;
        const double c = std::cos(ang), s = std::sin(ang);
        re += y[n].real() * c - y[n].imag() * s;
        im += y[n].real() * s + y[n].imag() * c;
    }
    return std::hypot(re, im);
}

// Peak of |DTFT| by golden-section search in [lo, hi] (unimodal within +-1 bin of a single tone).
double dtft_peak(const std::vector<cf32>& y, double lo, double hi) {
    const double g = (std::sqrt(5.0) - 1.0) / 2.0;
    double a = lo, b = hi;
    double x1 = b - g * (b - a), x2 = a + g * (b - a);
    double f1 = dtft_mag(y, x1), f2 = dtft_mag(y, x2);
    while (b - a > 1e-5) {
        if (f1 < f2) {
            a = x1; x1 = x2; f1 = f2; x2 = a + g * (b - a); f2 = dtft_mag(y, x2);
        } else {
            b = x2; x2 = x1; f2 = f1; x1 = b - g * (b - a); f1 = dtft_mag(y, x1);
        }
    }
    return 0.5 * (a + b);
}

void tone_phase_test(int64_t count, Kernel kernel) {
    for (const auto& tc : kToneCases) {
        const PhaseStats st = tone_phase_stats(tc, count, kernel);
        std::printf("[ measure  ] tone rho0=%.9g v=%.6g N=%lld: max phase err %.3e rad, last 1e5 %.3e rad, "
                    "max ||y|-1| %.2e\n",
                    tc.rho0, tc.v, static_cast<long long>(count), st.max_err, st.max_err_tail, st.max_mag_err);
        EXPECT_LT(st.max_err, 2e-6);
        EXPECT_LT(st.max_err_tail, 2e-6);  // no drift: the end of the run meets the same limit
        EXPECT_LT(st.max_mag_err, 1e-6);
    }
}

}  // namespace

TEST_P(GeoTone, Phase) { tone_phase_test(1'000'000, GetParam()); }

TEST_P(GeoTone, PhaseSlow) {
    NTN_REQUIRE_SLOW();
    tone_phase_test(10'000'000, GetParam());
}

TEST_P(GeoTone, FftPeakAtDoppler) {
    // M1a test 1, second criterion: the tone sits at f_d = -fc v / c within 0.1 Hz.
    // Peak refined by |DTFT| maximization in +-1 bin (fs/N = 7.68 Hz) around f_d; a peak at the
    // bracket edge would mean the true peak lies outside, so that is rejected too.
    constexpr int64_t kN = 1'000'000;
    const double df = kFs / static_cast<double>(kN);
    for (const auto& tc : kToneCases) {
        std::vector<cf32> y;
        (void)tone_phase_stats(tc, kN, GetParam(), &y);
        const double fd = -kFc * tc.v / kC;
        const double peak = dtft_peak(y, fd - df, fd + df);
        std::printf("[ measure  ] tone v=%.6g: peak %.6f Hz, f_d %.6f Hz, error %.2e Hz\n", tc.v, peak, fd,
                    std::abs(peak - fd));
        EXPECT_LT(std::abs(peak - fd), 0.1);
        EXPECT_LT(std::abs(peak - fd), 0.9 * df);
    }
}

// ---- Test 2: integer shift, bit-exact ----------------------------------------------------------
//
// fc = 0, tau - tau0 = k/fs (rho = (tau0 + k) c / fs; rho fs / c lands within the 1e-9 snap of
// the integer). Expected: y(n) = x(n - k), bit-exact, zero for n - k < 0. Output index i maps to
// channel time i - L; the run feeds L extra samples so every channel time in [0, N) is produced.
TEST_P(GeoShift, IntegerShiftBitExact) {
    constexpr int64_t kN = 200'000;
    constexpr int64_t kTau0 = 15'360;  // 2 ms at 7.68 Msps
    std::mt19937_64 rng(2);
    std::normal_distribution<float> g(0.0f, 1.0f);
    for (int64_t k : {int64_t{0}, int64_t{1}, int64_t{37}, int64_t{70'001}}) {
        const double rho = static_cast<double>(kTau0 + k) * kC / kFs;
        GeometricChannel ch(RangeTrack(kFs, {{0, kN, {rho, 0.0, 0.0, 0.0}}}), geo(0.0, kTau0, true, true), GetParam());
        const int64_t L = ch.latency();
        // D = k (to float64 rounding); with the 1e-6 guard floor gives k - 1, so L = 16 - k, >= 0.
        const int64_t expected_L = std::max<int64_t>(0, ntn::fd::kHalf - k);
        std::printf("[ measure  ] integer shift k=%lld: latency L = %lld\n", static_cast<long long>(k),
                    static_cast<long long>(L));
        EXPECT_EQ(L, expected_L);

        std::vector<cf32> x(static_cast<std::size_t>(kN + L));
        for (int64_t i = 0; i < kN; ++i) x[static_cast<std::size_t>(i)] = {g(rng), g(rng)};
        const auto y = run(ch, x);
        for (int64_t i = 0; i < L; ++i) {
            ASSERT_EQ(y[static_cast<std::size_t>(i)], cf32(0.0f, 0.0f)) << "k=" << k << " i=" << i;
        }
        for (int64_t n = 0; n < kN; ++n) {
            const cf32 want = n - k >= 0 ? x[static_cast<std::size_t>(n - k)] : cf32(0.0f, 0.0f);
            const cf32 got = y[static_cast<std::size_t>(n + L)];
            ASSERT_EQ(std::memcmp(&got, &want, sizeof(cf32)), 0) << "k=" << k << " n=" << n;
        }
    }
}

// ---- Test 4: quadratic phase -------------------------------------------------------------------
//
// Cubic segment with known coefficients, window starting 5 s into it, delay off. The
// instantaneous frequency (+-100-sample central difference of the unwrapped float64 phase of
// the float32 output) equals -fc rho_dot(t) / c within 0.01 Hz. The central-difference bias,
// fc/c * rho''' * (M/fs)^2 / 6, is ~1e-9 Hz; float32 phase noise (~1e-7 rad) gives
// ~1e-7 * fs / (2 pi * 200) ~ 6e-4 Hz.
TEST_P(GeoQuadratic, InstantaneousFrequency) {
    constexpr double c0 = 1'200'000.0, c1 = -5000.0, c2 = 40.0, c3 = -0.5;  // M1a QUAD_COEFFS
    const int64_t n_start = 5 * static_cast<int64_t>(kFs);
    constexpr int64_t kN = 200'000;
    constexpr int64_t kM = 100;
    GeometricChannel ch(RangeTrack(kFs, {{0, 10 * static_cast<int64_t>(kFs), {c0, c1, c2, c3}}}),
                        geo(kFc, 0, false, true), GetParam());
    ch.reset(n_start);
    EXPECT_EQ(ch.sample_index(), n_start);
    const auto y = run(ch, std::vector<cf32>(kN, cf32(1.0f, 0.0f)));

    std::vector<double> phi(kN);
    phi[0] = std::arg(cd(y[0].real(), y[0].imag()));
    for (std::size_t i = 1; i < y.size(); ++i) {
        const double a = std::arg(cd(y[i].real(), y[i].imag()));
        phi[i] = phi[i - 1] + wrap_pi(a - phi[i - 1]);
    }
    double max_err = 0.0;
    for (int64_t i = kM; i < kN - kM; ++i) {
        const double f_meas = (phi[static_cast<std::size_t>(i + kM)] - phi[static_cast<std::size_t>(i - kM)]) *
                              kFs / (2.0 * kPi * 2.0 * kM);
        const double u = static_cast<double>(n_start + i) / kFs;
        const double f_true = -kFc * (c1 + 2.0 * c2 * u + 3.0 * c3 * u * u) / kC;
        max_err = std::max(max_err, std::abs(f_meas - f_true));
    }
    std::printf("[ measure  ] quadratic phase: max |f_inst - (-fc rho_dot / c)| = %.3e Hz\n", max_err);
    EXPECT_LT(max_err, 0.01);
}

// ---- Latency (closed form) ---------------------------------------------------------------------

// Read window: for non-integer D the taps read x[m0-15 .. m0+16], m0 = n - floor(D) - 1, so the
// highest index read is n + 15 - floor(D); an integer D reads x[n - D] only. Hence
//     L = max(0, 15 - floor(min D))
// (the implementation applies floor to min D - 1e-6 to absorb float64 rounding of the bound).
TEST_P(GeoLatency, ApproachingTrack) {
    // rho = (tau0 + 3) c / fs + v u, v = -7000 m/s, 1 s: D(n) = 3 + v n / c, decreasing.
    // min D at n = N - 1: 3 - 7000 * 7679999 / c = -176.324..., so L = 15 + 177 = 192.
    constexpr int64_t kTau0 = 26'244;
    constexpr int64_t kN = 7'680'000;
    constexpr double v = -7000.0;
    const double rho0 = (static_cast<double>(kTau0) + 3.0) * kC / kFs;
    GeometricChannel ch(RangeTrack(kFs, {{0, kN, {rho0, v, 0.0, 0.0}}}), geo(kFc, kTau0, true, true), GetParam());
    const double d_min = 3.0 + v * static_cast<double>(kN - 1) / kC;
    std::printf("[ measure  ] approaching track (v=-7000 m/s, 1 s): min D = %.4f, latency L = %lld\n", d_min,
                static_cast<long long>(ch.latency()));
    EXPECT_EQ(static_cast<int64_t>(15.0 - std::floor(d_min)), 192);
    EXPECT_EQ(ch.latency(), 192);
}

TEST_P(GeoLatency, RecedingTrackNeedsNone) {
    // D(n) = 20 + v n / c with v > 0: min D = 20 >= 16, so L = 0.
    constexpr int64_t kTau0 = 26'244;
    const double rho0 = (static_cast<double>(kTau0) + 20.0) * kC / kFs;
    GeometricChannel ch(RangeTrack(kFs, {{0, 7'680'000, {rho0, 7000.0, 0.0, 0.0}}}), geo(kFc, kTau0, true, true),
                        GetParam());
    EXPECT_EQ(ch.latency(), 0);
}

TEST_P(GeoLatency, DelayOffNeedsNone) {
    GeometricChannel ch(RangeTrack(kFs, {{0, 1000, {1e6, -7000.0, 0.0, 0.0}}}), geo(kFc, 0, false, true), GetParam());
    EXPECT_EQ(ch.latency(), 0);
}

// L is minimal (non-integer min D). Output for channel time n is emitted at output index n + L,
// in the same process() position as input sample n + L. The read-window arithmetic says that
// output reads x up to n + 15 - floor(D(n)), which equals n + L where D is minimal. An input
// impulse at m* = n + L must therefore change that output (the kernel really reads the newest
// sample, so with L - 1 the output would be emitted before in[m*] arrives), while an impulse
// at m* + 1 must leave it exactly zero (L is sufficient). fc = 0 isolates the delay.
namespace {

std::vector<cf32> impulse_response(const RangeTrack& track, int64_t tau0, std::size_t len, std::size_t at,
                                   std::size_t block, Kernel kernel) {
    GeometricChannel ch(track, geo(0.0, tau0, true, true), kernel);
    std::vector<cf32> in(len, cf32(0.0f, 0.0f));
    in[at] = cf32(1.0f, 0.0f);
    return run(ch, in, block);
}

}  // namespace

TEST_P(GeoLatency, MinimalForConstantFractionalDelay) {
    // D = 2.3 everywhere: L = 15 - floor(2.3) = 13; every output reads up to x[n + 13].
    constexpr int64_t kTau0 = 26'244;
    constexpr int64_t kN = 10'000;
    const RangeTrack track(kFs, {{0, kN, {(static_cast<double>(kTau0) + 2.3) * kC / kFs, 0.0, 0.0, 0.0}}});
    const int64_t L = GeometricChannel(track, geo(0.0, kTau0, true, true), GetParam()).latency();
    ASSERT_EQ(L, 13);

    const int64_t n = 5'000;
    const auto m_star = static_cast<std::size_t>(n + 15 - 2);  // read-window top = n + L
    ASSERT_EQ(static_cast<int64_t>(m_star), n + L);
    for (std::size_t block : {std::size_t{1}, std::size_t{4096}}) {
        const auto y_needed = impulse_response(track, kTau0, kN, m_star, block, GetParam());
        const auto y_beyond = impulse_response(track, kTau0, kN, m_star + 1, block, GetParam());
        std::printf("[ measure  ] D=2.3, L=13, block %zu: |out[n+L]| with impulse at n+L = %.3e, at n+L+1 = %.3e\n",
                    block, std::abs(y_needed[m_star]), std::abs(y_beyond[m_star]));
        EXPECT_NE(y_needed[m_star], cf32(0.0f, 0.0f));  // needs in[n + L]: L - 1 would be too small
        EXPECT_EQ(y_beyond[m_star], cf32(0.0f, 0.0f));  // does not need in[n + L + 1]: L suffices
    }
}

TEST_P(GeoLatency, MinimalForApproachingTrack) {
    // D(n) = 0.4 + v n / c, v = -7000 m/s, 1e6 samples: minimal at n = N - 1, D = -22.949...,
    // L = 15 + 23 = 38. Output for channel time N - 1 is the last output (index N - 1 + L) and
    // must depend on the last input sample in[N - 1 + L].
    constexpr int64_t kTau0 = 26'244;
    constexpr int64_t kN = 1'000'000;
    constexpr double v = -7000.0;
    const RangeTrack track(kFs, {{0, kN, {(static_cast<double>(kTau0) + 0.4) * kC / kFs, v, 0.0, 0.0}}});
    const double d_last = 0.4 + v * static_cast<double>(kN - 1) / kC;
    const int64_t L = GeometricChannel(track, geo(0.0, kTau0, true, true), GetParam()).latency();
    ASSERT_EQ(L, 15 - static_cast<int64_t>(std::floor(d_last)));
    ASSERT_EQ(L, 38);

    const auto last = static_cast<std::size_t>(kN - 1 + L);
    const auto y = impulse_response(track, kTau0, last + 1, last, 4096, GetParam());
    std::printf("[ measure  ] approaching, D(N-1)=%.4f, L=%lld: |last output| with impulse at last input = %.3e\n",
                d_last, static_cast<long long>(L), std::abs(y[last]));
    EXPECT_NE(y[last], cf32(0.0f, 0.0f));
}

// ---- Time-varying delay against its closed form ------------------------------------------------
//
// fc = 0, approaching satellite (v = -7000 m/s) over 1e6 samples: D(n) = 0.37 + v n / c falls by
// ~23 samples (23 integer slips, negative D for most of the run). Input: 32 complex tones at
// random frequencies in |f| < 0.4 fs, evaluated exactly, so the reference is x(n - D(n)) in
// closed form. Limit: -90 dB, the same as the kernel's band-limited test.
TEST_P(GeoDelay, TimeVaryingDelayMatchesClosedForm) {
    constexpr int64_t kTau0 = 26'244;
    constexpr int64_t kN = 1'000'000;
    constexpr int kTones = 32;
    constexpr double v = -7000.0;
    const double rho0 = (static_cast<double>(kTau0) + 0.37) * kC / kFs;
    GeometricChannel ch(RangeTrack(kFs, {{0, kN, {rho0, v, 0.0, 0.0}}}), geo(0.0, kTau0, true, true), GetParam());
    const int64_t L = ch.latency();
    const double d_min = 0.37 + v * static_cast<double>(kN - 1) / kC;  // -22.979...
    std::printf("[ measure  ] time-varying delay: min D = %.4f, latency L = %lld\n", d_min,
                static_cast<long long>(L));
    EXPECT_EQ(L, 15 - static_cast<int64_t>(std::floor(d_min)));  // 38

    std::mt19937_64 rng(7);
    std::uniform_real_distribution<double> uf(-0.4, 0.4), up(0.0, 1.0);
    double freq[kTones], ph[kTones];
    for (int t = 0; t < kTones; ++t) {
        freq[t] = uf(rng);
        ph[t] = up(rng);
    }
    auto x_at = [&](double t) {  // x(t) = sum exp(j 2 pi (f t + ph)) / sqrt(K), t in samples
        cd s(0.0, 0.0);
        for (int k = 0; k < kTones; ++k) {
            double cyc = freq[k] * t + ph[k];
            cyc -= std::nearbyint(cyc);
            s += cd(std::cos(2.0 * kPi * cyc), std::sin(2.0 * kPi * cyc));
        }
        return s / std::sqrt(static_cast<double>(kTones));
    };

    std::vector<cf32> in(static_cast<std::size_t>(kN + L));
    for (std::size_t m = 0; m < in.size(); ++m) {
        const cd v_ = x_at(static_cast<double>(m));
        in[m] = {static_cast<float>(v_.real()), static_cast<float>(v_.imag())};
    }
    const auto y = run(ch, in, 4096);

    constexpr int64_t kSkip = 64;  // first outputs read x before t = 0 (input is zero there)
    std::vector<cd> ref;
    std::vector<cf32> got;
    for (int64_t n = kSkip; n < kN; ++n) {
        const double d = 0.37 + v * static_cast<double>(n) / kC;
        ref.push_back(x_at(static_cast<double>(n) - d));
        got.push_back(y[static_cast<std::size_t>(n + L)]);
    }
    const double e = ntn_test::error_db<cf32, cd>(got, ref);
    std::printf("[ measure  ] time-varying delay vs closed form: %.1f dB\n", e);
    EXPECT_LT(e, -90.0);
}

// ---- Validation --------------------------------------------------------------------------------

TEST(GeoValidation, RejectsBadConfig) {
    const RangeTrack track(kFs, {{0, 1000, {1e6, 0.0, 0.0, 0.0}}});
    GeoConfig g = geo(kFc, 0, true, true);
    g.fs_hz = kFs + 1.0;
    EXPECT_THROW(GeometricChannel(track, g), std::invalid_argument);
    g = geo(kFc, 0, true, true);
    g.residual_ref_doppler_hz = 100.0;
    EXPECT_THROW(GeometricChannel(track, g), std::invalid_argument);
    g = geo(std::nan(""), 0, true, true);
    EXPECT_THROW(GeometricChannel(track, g), std::invalid_argument);
    // M2: Kernel::Avx2 works and is the default (M1b rejected it with std::invalid_argument).
    EXPECT_NO_THROW(GeometricChannel(track, geo(kFc, 0, true, true), ntn::Kernel::Avx2));
    EXPECT_EQ(GeometricChannel(track, geo(kFc, 0, true, true)).kernel(), ntn::Kernel::Avx2);
}

TEST(GeoValidation, BothKernelsRequireCpuFeatures) {
    // ntn_core is compiled with /arch:AVX2 (and BMI), so construction checks the CPU for either
    // kernel. A hidden feature (test-only switch, cpu_check.hpp) must make both throw.
    if (!ntn::cpu_has_avx2()) GTEST_SKIP() << "CPU without the AVX2 feature set";
    const RangeTrack track(kFs, {{0, 1000, {1e6, 0.0, 0.0, 0.0}}});
    const ntn_test::EnvGuard hide("NTN_TEST_HIDE_AVX2", "BMI2");
    for (const Kernel k : {Kernel::Scalar, Kernel::Avx2}) {
        try {
            GeometricChannel ch(track, geo(kFc, 0, true, true), k);
            ADD_FAILURE() << "constructed without BMI2, kernel " << (k == Kernel::Scalar ? "Scalar" : "Avx2");
        } catch (const std::runtime_error& e) {
            EXPECT_NE(std::string(e.what()).find("missing: BMI2"), std::string::npos) << e.what();
        }
    }
}

TEST(GeoValidation, SizeMismatchThrows) {
    GeometricChannel ch(RangeTrack(kFs, {{0, 1000, {1e6, 0.0, 0.0, 0.0}}}), geo(kFc, 0, false, true));
    std::vector<cf32> in(10), out(9);
    EXPECT_THROW(ch.process(in, out), std::invalid_argument);
}

TEST(GeoValidation, PastTrackEndThrowsWithoutStateChange) {
    GeometricChannel ch(RangeTrack(kFs, {{0, 1000, {1e6, 0.0, 0.0, 0.0}}}), geo(kFc, 0, false, true));
    std::vector<cf32> in(1000, cf32(1.0f, 0.0f)), out(1000);
    ch.process(in, out);
    EXPECT_EQ(ch.sample_index(), 1000);
    std::vector<cf32> one(1), o1(1);
    EXPECT_THROW(ch.process(one, o1), std::out_of_range);
    EXPECT_EQ(ch.sample_index(), 1000);
}

INSTANTIATE_TEST_SUITE_P(Kernels, GeoTone, ::testing::Values(Kernel::Scalar, Kernel::Avx2), kernel_name);
INSTANTIATE_TEST_SUITE_P(Kernels, GeoShift, ::testing::Values(Kernel::Scalar, Kernel::Avx2), kernel_name);
INSTANTIATE_TEST_SUITE_P(Kernels, GeoQuadratic, ::testing::Values(Kernel::Scalar, Kernel::Avx2), kernel_name);
INSTANTIATE_TEST_SUITE_P(Kernels, GeoLatency, ::testing::Values(Kernel::Scalar, Kernel::Avx2), kernel_name);
INSTANTIATE_TEST_SUITE_P(Kernels, GeoDelay, ::testing::Values(Kernel::Scalar, Kernel::Avx2), kernel_name);
