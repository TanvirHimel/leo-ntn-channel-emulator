// Pipeline (geo only): block-size invariance (scalar: bit-identical; AVX2: < -120 dB, SPEC §5.1),
// reset() and sample_index(), default kernel and the CPU check. The M1b state tests keep the
// scalar kernel explicitly.
#include "ntn/pipeline.hpp"

#include <gtest/gtest.h>

#include <cmath>
#include <cstdio>
#include <cstring>
#include <random>
#include <span>
#include <stdexcept>
#include <string>
#include <vector>

#include "ntn/cpu_check.hpp"
#include "test_util.hpp"

using ntn::cf32;
using ntn::GeoConfig;
using ntn::Kernel;
using ntn::Pipeline;
using ntn::RangeSegment;
using ntn::RangeTrack;

namespace {

constexpr double kC = 299'792'458.0;
constexpr double kFs = 7'680'000.0;
constexpr double kFc = 2'000'000'000.0;
constexpr int64_t kSeg = 200'000;  // samples per segment; knot at kSeg

// Mid-pass-like cubic (rho ~ 1024 km, rho_dot ~ -5.9 km/s), split into two segments at the
// knot by Taylor shift; c0', c1' are the Horner end values of segment 1, so C1 holds exactly.
// D falls by ~5899.5 * 400000 / c ~ 7.9 samples: several integer slips, one segment knot.
RangeTrack slip_track() {
    const std::array<double, 4> c{1'024'452.123, -5899.5, 13.4, -0.02};
    const double h = static_cast<double>(kSeg) / kFs;
    const std::array<double, 4> c2{ntn::horner_rho(c, h), ntn::horner_rho_dot(c, h), c[2] + 3.0 * c[3] * h, c[3]};
    return RangeTrack(kFs, {{0, kSeg, c}, {kSeg, 2 * kSeg, c2}});
}

GeoConfig slip_geo(const RangeTrack& t) {
    GeoConfig g;
    g.fc_hz = kFc;
    g.fs_hz = kFs;
    g.tau0_samples = static_cast<int64_t>(std::floor(t.rho(0) * kFs / kC));  // as in demo_m1a
    return g;
}

std::vector<cf32> random_input(std::size_t n, uint64_t seed) {
    std::mt19937_64 rng(seed);
    std::normal_distribution<float> g(0.0f, 1.0f);
    std::vector<cf32> x(n);
    for (auto& v : x) v = {g(rng), g(rng)};
    return x;
}

// Feeds `in` in blocks following `pattern` (cycled).
std::vector<cf32> run_blocks(Pipeline& p, const std::vector<cf32>& in, const std::vector<std::size_t>& pattern) {
    std::vector<cf32> out(in.size());
    std::size_t s = 0, k = 0;
    while (s < in.size()) {
        const std::size_t m = std::min(pattern[k++ % pattern.size()], in.size() - s);
        p.process(std::span<const cf32>(in.data() + s, m), std::span<cf32>(out.data() + s, m));
        s += m;
    }
    return out;
}

bool bit_identical(const std::vector<cf32>& a, const std::vector<cf32>& b) {
    return a.size() == b.size() && std::memcmp(a.data(), b.data(), a.size() * sizeof(cf32)) == 0;
}

}  // namespace

TEST(PipelineBlocks, ScalarBlockSizeInvarianceBitIdentical) {
    const RangeTrack track = slip_track();
    const GeoConfig g = slip_geo(track);
    const double d_first = track.rho(0) * kFs / kC - static_cast<double>(g.tau0_samples);
    const double d_last = track.rho(2 * kSeg - 1) * kFs / kC - static_cast<double>(g.tau0_samples);
    const auto slips = static_cast<int64_t>(std::floor(d_first) - std::floor(d_last));
    const int64_t L = Pipeline(track, g, Kernel::Scalar).latency();
    std::printf("[ measure  ] block invariance track: D %.4f -> %.4f (%lld integer slips), knot at %lld, L = %lld\n",
                d_first, d_last, static_cast<long long>(slips), static_cast<long long>(kSeg),
                static_cast<long long>(L));
    ASSERT_GE(slips, 5);

    const auto in = random_input(static_cast<std::size_t>(2 * kSeg + L), 99);
    Pipeline ref_p(track, g, Kernel::Scalar);
    const auto ref = run_blocks(ref_p, in, {4096});
    double energy = 0.0;
    for (std::size_t i = static_cast<std::size_t>(L); i < ref.size(); ++i) energy += std::norm(ref[i]);
    ASSERT_GT(energy, 0.0);

    for (const auto& pattern : std::vector<std::vector<std::size_t>>{{1}, {7}, {1, 7, 4096, 3, 1000, 13}}) {
        Pipeline p(track, g, Kernel::Scalar);
        const auto y = run_blocks(p, in, pattern);
        EXPECT_TRUE(bit_identical(y, ref)) << "first block size " << pattern[0];
        EXPECT_EQ(p.sample_index(), static_cast<int64_t>(in.size()));
    }
}

TEST(PipelineState, SampleIndexAdvancesByBlockSize) {
    const RangeTrack track = slip_track();
    Pipeline p(track, slip_geo(track), Kernel::Scalar);
    EXPECT_EQ(p.sample_index(), track.n_begin());
    std::vector<cf32> in(7), out(7);
    p.process(in, out);
    EXPECT_EQ(p.sample_index(), 7);
    std::vector<cf32> in1(1), out1(1);
    p.process(in1, out1);
    EXPECT_EQ(p.sample_index(), 8);
    p.reset(12'345);
    EXPECT_EQ(p.sample_index(), 12'345);
}

TEST(PipelineState, ResetToStartRepeatsOutputExactly) {
    const RangeTrack track = slip_track();
    const GeoConfig g = slip_geo(track);
    const auto in = random_input(50'000, 5);
    Pipeline fresh(track, g, Kernel::Scalar);
    const auto ref = run_blocks(fresh, in, {4096});

    Pipeline p(track, g, Kernel::Scalar);
    (void)run_blocks(p, random_input(30'000, 6), {777});  // unrelated history
    p.reset(track.n_begin());
    EXPECT_TRUE(bit_identical(run_blocks(p, in, {7}), ref));
}

TEST(PipelineState, ResetMidTrackEqualsZeroHistory) {
    // reset(n_start) means "input before n_start is zero". Reference: a fresh pipeline fed zeros
    // for [0, n_start), then the same input. Its outputs from index n_start on cover channel
    // times n_start - L onward, exactly the outputs of the reset pipeline.
    const RangeTrack track = slip_track();
    const GeoConfig g = slip_geo(track);
    constexpr int64_t kStart = 150'000;  // before the knot; the run crosses it
    const auto in = random_input(100'000, 8);

    std::vector<cf32> padded(static_cast<std::size_t>(kStart), cf32(0.0f, 0.0f));
    padded.insert(padded.end(), in.begin(), in.end());
    Pipeline zero_fed(track, g, Kernel::Scalar);
    const auto ref_full = run_blocks(zero_fed, padded, {4096});
    const std::vector<cf32> ref(ref_full.begin() + kStart, ref_full.end());

    Pipeline p(track, g, Kernel::Scalar);
    (void)run_blocks(p, random_input(20'000, 9), {4096});  // unrelated history
    p.reset(kStart);
    EXPECT_EQ(p.sample_index(), kStart);
    const auto y = run_blocks(p, in, {1, 7, 4096});
    EXPECT_TRUE(bit_identical(y, ref));
    EXPECT_EQ(p.sample_index(), kStart + static_cast<int64_t>(in.size()));
}

TEST(PipelineState, LatencyMatchesGeometricStage) {
    const RangeTrack track = slip_track();
    const GeoConfig g = slip_geo(track);
    EXPECT_EQ(Pipeline(track, g).latency(), ntn::GeometricChannel(track, g).latency());
}

// ---- M2: AVX2 path ----------------------------------------------------------------------------

TEST(PipelineBlocks, Avx2BlockSizeInvariance) {
    // SPEC §5.1: on the AVX2 path, block sizes 1, 7 and 4096 agree to < -120 dB (the same track,
    // input and patterns as the scalar test).
    const RangeTrack track = slip_track();
    const GeoConfig g = slip_geo(track);
    const int64_t L = Pipeline(track, g, Kernel::Avx2).latency();
    const auto in = random_input(static_cast<std::size_t>(2 * kSeg + L), 99);
    Pipeline ref_p(track, g, Kernel::Avx2);
    const auto ref = run_blocks(ref_p, in, {4096});
    const std::vector<cf32> ref_tail(ref.begin() + L, ref.end());
    for (const auto& pattern : std::vector<std::vector<std::size_t>>{{1}, {7}, {1, 7, 4096, 3, 1000, 13}}) {
        Pipeline p(track, g, Kernel::Avx2);
        const auto y = run_blocks(p, in, pattern);
        const std::vector<cf32> y_tail(y.begin() + L, y.end());
        std::size_t differ = 0;
        for (std::size_t i = 0; i < y.size(); ++i) differ += std::memcmp(&y[i], &ref[i], sizeof(cf32)) != 0;
        const double e = ntn_test::error_db<cf32, cf32>(y_tail, ref_tail);
        std::string label;
        for (std::size_t b : pattern) label += (label.empty() ? "" : ",") + std::to_string(b);
        std::printf("[ measure  ] AVX2 block invariance, blocks {%s} vs {4096}: %.1f dB, %zu of %zu samples "
                    "differ in bits (limit -120 dB)\n",
                    label.c_str(), e, differ, y.size());
        EXPECT_TRUE(differ == 0 || e < -120.0) << "first block size " << pattern[0];
        EXPECT_EQ(p.sample_index(), static_cast<int64_t>(in.size()));
    }
}

TEST(PipelineDefaults, Avx2IsTheDefaultKernel) {
    const RangeTrack track = slip_track();
    EXPECT_EQ(Pipeline(track, slip_geo(track)).kernel(), Kernel::Avx2);
    EXPECT_EQ(Pipeline(track, slip_geo(track), Kernel::Scalar).kernel(), Kernel::Scalar);
}

TEST(PipelineDefaults, BothKernelsRequireCpuFeatures) {
    // ntn_core is compiled with /arch:AVX2 and BMI: construction checks the CPU for either kernel.
    if (!ntn::cpu_has_avx2()) GTEST_SKIP() << "CPU without the AVX2 feature set";
    const RangeTrack track = slip_track();
    const ntn_test::EnvGuard hide("NTN_TEST_HIDE_AVX2", "FMA3");
    for (const Kernel k : {Kernel::Scalar, Kernel::Avx2}) {
        try {
            Pipeline p(track, slip_geo(track), k);
            ADD_FAILURE() << "constructed without FMA3, kernel " << (k == Kernel::Scalar ? "Scalar" : "Avx2");
        } catch (const std::runtime_error& e) {
            EXPECT_NE(std::string(e.what()).find("missing: FMA3"), std::string::npos) << e.what();
        }
    }
}
