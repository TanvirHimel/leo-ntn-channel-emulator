// RangeTrack: validation (contiguity, C1) and evaluation against closed-form expectations.
#include "ntn/range_track.hpp"

#include <gtest/gtest.h>

#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <vector>

using ntn::RangeSegment;
using ntn::RangeTrack;

namespace {

// ---- Exact closed form ------------------------------------------------------------------
// fs = 1024 and dyadic coefficients: u = k/1024 and every Horner intermediate are exactly
// representable (rho < 2^21 m with resolution 2^-31 m needs <= 52 bits), so evaluation must
// be bit-exact. The reference is computed in int64, scaled by 2^32, independently of Horner:
//   P(u) = 1048576 - 7000.5 u + 40.25 u^2 - 0.5 u^3,  u = k/1024
//   rho     * 2^32 = 2^52 - 14001 k 2^21 + 161 k^2 2^10 - 2 k^3
//   rho_dot * 2^32 = -14001 * 2^31 + 161 k 2^21 - 3 k^2 2^11
// Segment 2 is the same polynomial Taylor-shifted to h = 4 s (k = 4096):
//   c0' = P(4) = 1021186, c1' = P'(4) = -6702.5, c2' = c2 + 3 c3 h = 34.25, c3' = c3.
constexpr double kFsExact = 1024.0;
constexpr int64_t kN0 = 5000;  // global index of the first sample (u != n / fs)
constexpr int64_t kSegLen = 4096;

double exact_rho(int64_t k) {
    const int64_t v = (int64_t{1} << 52) - 14001 * k * (int64_t{1} << 21) +
                      161 * k * k * (int64_t{1} << 10) - 2 * k * k * k;
    return std::ldexp(static_cast<double>(v), -32);
}

double exact_rho_dot(int64_t k) {
    const int64_t v = -14001 * (int64_t{1} << 31) + 161 * k * (int64_t{1} << 21) -
                      3 * k * k * (int64_t{1} << 11);
    return std::ldexp(static_cast<double>(v), -32);
}

std::vector<RangeSegment> exact_segments() {
    return {{kN0, kN0 + kSegLen, {1048576.0, -7000.5, 40.25, -0.5}},
            {kN0 + kSegLen, kN0 + 2 * kSegLen, {1021186.0, -6702.5, 34.25, -0.5}}};
}

// A valid one-second segment at fs = 1024; its end values are exact:
// rho(1) = 1041615.25 m, rho_dot(1) = -6921.5 m/s.
RangeSegment seg_a() { return {0, 1024, {1048576.0, -7000.5, 40.25, -0.5}}; }
constexpr double kEndRhoA = 1041615.25;
constexpr double kEndRateA = -6921.5;

RangeSegment seg_after_a(double d_rho, double d_rate) {
    return {1024, 2048, {kEndRhoA + d_rho, kEndRateA + d_rate, 0.0, 0.0}};
}

}  // namespace

// ---- Validation ---------------------------------------------------------------------------

TEST(RangeTrackValidation, RejectsEmptySegmentList) {
    EXPECT_THROW(RangeTrack(kFsExact, {}), std::invalid_argument);
}

TEST(RangeTrackValidation, RejectsNonPositiveOrNonFiniteFs) {
    for (double fs : {0.0, -1024.0, std::numeric_limits<double>::quiet_NaN(),
                      std::numeric_limits<double>::infinity()}) {
        EXPECT_THROW(RangeTrack(fs, {seg_a()}), std::invalid_argument) << "fs = " << fs;
    }
}

TEST(RangeTrackValidation, RejectsEmptyOrReversedSegment) {
    EXPECT_THROW(RangeTrack(kFsExact, {{10, 10, {1.0, 0.0, 0.0, 0.0}}}), std::invalid_argument);
    EXPECT_THROW(RangeTrack(kFsExact, {{10, 5, {1.0, 0.0, 0.0, 0.0}}}), std::invalid_argument);
}

TEST(RangeTrackValidation, RejectsGap) {
    RangeSegment b = seg_after_a(0.0, 0.0);
    b.n0 += 1;
    b.n1 += 1;
    EXPECT_THROW(RangeTrack(kFsExact, {seg_a(), b}), std::invalid_argument);
}

TEST(RangeTrackValidation, RejectsOverlap) {
    RangeSegment b = seg_after_a(0.0, 0.0);
    b.n0 -= 1;
    EXPECT_THROW(RangeTrack(kFsExact, {seg_a(), b}), std::invalid_argument);
}

TEST(RangeTrackValidation, RejectsNonFiniteCoefficient) {
    RangeSegment a = seg_a();
    a.c[2] = std::numeric_limits<double>::quiet_NaN();
    EXPECT_THROW(RangeTrack(kFsExact, {a}), std::invalid_argument);
    a.c[2] = std::numeric_limits<double>::infinity();
    EXPECT_THROW(RangeTrack(kFsExact, {a}), std::invalid_argument);
}

TEST(RangeTrackValidation, AcceptsContiguousC1Track) {
    EXPECT_NO_THROW(RangeTrack(kFsExact, {seg_a(), seg_after_a(0.0, 0.0)}));
    EXPECT_NO_THROW(RangeTrack(kFsExact, exact_segments()));
}

// C1 tolerances are 1e-6 m and 1e-6 m/s: 2e-6 is rejected, 0.5e-6 accepted.
TEST(RangeTrackValidation, RejectsRhoJumpAboveTolerance) {
    EXPECT_THROW(RangeTrack(kFsExact, {seg_a(), seg_after_a(2e-6, 0.0)}), std::invalid_argument);
    EXPECT_THROW(RangeTrack(kFsExact, {seg_a(), seg_after_a(-2e-6, 0.0)}), std::invalid_argument);
}

TEST(RangeTrackValidation, AcceptsRhoJumpBelowTolerance) {
    EXPECT_NO_THROW(RangeTrack(kFsExact, {seg_a(), seg_after_a(0.5e-6, 0.0)}));
    EXPECT_NO_THROW(RangeTrack(kFsExact, {seg_a(), seg_after_a(-0.5e-6, 0.0)}));
}

TEST(RangeTrackValidation, RejectsRateJumpAboveTolerance) {
    EXPECT_THROW(RangeTrack(kFsExact, {seg_a(), seg_after_a(0.0, 2e-6)}), std::invalid_argument);
    EXPECT_THROW(RangeTrack(kFsExact, {seg_a(), seg_after_a(0.0, -2e-6)}), std::invalid_argument);
}

TEST(RangeTrackValidation, AcceptsRateJumpBelowTolerance) {
    EXPECT_NO_THROW(RangeTrack(kFsExact, {seg_a(), seg_after_a(0.0, 0.5e-6)}));
    EXPECT_NO_THROW(RangeTrack(kFsExact, {seg_a(), seg_after_a(0.0, -0.5e-6)}));
}

TEST(RangeTrackValidation, ChecksEveryInteriorKnot) {
    // Knot 1 is exact; knot 2 (end of the constant-slope segment B) has a 2e-6 m jump.
    const RangeSegment b = seg_after_a(0.0, 0.0);  // rho = kEndRhoA + kEndRateA * u
    const double end_b = kEndRhoA + kEndRateA * 1.0;
    const RangeSegment c_bad{2048, 3072, {end_b + 2e-6, kEndRateA, 0.0, 0.0}};
    const RangeSegment c_ok{2048, 3072, {end_b, kEndRateA, 0.0, 0.0}};
    EXPECT_NO_THROW(RangeTrack(kFsExact, {seg_a(), b, c_ok}));
    EXPECT_THROW(RangeTrack(kFsExact, {seg_a(), b, c_bad}), std::invalid_argument);
}

// ---- Evaluation: exact closed form --------------------------------------------------------

TEST(RangeTrackEval, ExactClosedFormBothSegments) {
    const RangeTrack track(kFsExact, exact_segments());
    EXPECT_EQ(track.n_begin(), kN0);
    EXPECT_EQ(track.n_end(), kN0 + 2 * kSegLen);
    for (int64_t k = 0; k < 2 * kSegLen; ++k) {
        const int64_t n = kN0 + k;
        ASSERT_EQ(track.rho(n), exact_rho(k)) << "k = " << k;
        ASSERT_EQ(track.rho_dot(n), exact_rho_dot(k)) << "k = " << k;
        ASSERT_EQ(track.tau(n), exact_rho(k) / 299'792'458.0) << "k = " << k;
    }
}

TEST(RangeTrackEval, SegmentSelectionAtKnot) {
    // Global k = 4095 belongs to segment 1 (u = 4095/1024), k = 4096 to segment 2 (u = 0).
    const RangeTrack track(kFsExact, exact_segments());
    const int64_t knot = kN0 + kSegLen;
    EXPECT_EQ(track.rho(knot - 1), exact_rho(kSegLen - 1));
    EXPECT_EQ(track.rho(knot), 1021186.0);  // c0' exactly: u = 0 in segment 2
    EXPECT_EQ(track.rho_dot(knot), -6702.5);
    EXPECT_EQ(track.rho(knot + 1), exact_rho(kSegLen + 1));
}

TEST(RangeTrackEval, OutOfRangeThrows) {
    const RangeTrack track(kFsExact, exact_segments());
    for (int64_t n : {kN0 - 1, kN0 + 2 * kSegLen, int64_t{0}, int64_t{-1}}) {
        EXPECT_THROW((void)track.rho(n), std::out_of_range) << "n = " << n;
        EXPECT_THROW((void)track.rho_dot(n), std::out_of_range) << "n = " << n;
        EXPECT_THROW((void)track.tau(n), std::out_of_range) << "n = " << n;
    }
    std::vector<double> buf(4);
    EXPECT_THROW(track.rho_block(kN0 - 1, buf), std::out_of_range);
    EXPECT_THROW(track.rho_block(kN0 + 2 * kSegLen - 3, buf), std::out_of_range);
    EXPECT_NO_THROW(track.rho_block(kN0 + 2 * kSegLen - 4, buf));
}

TEST(RangeTrackEval, BlockMatchesPointBitwiseAcrossKnot) {
    const RangeTrack track(kFsExact, exact_segments());
    std::vector<double> buf(301);
    const int64_t first = kN0 + kSegLen - 150;
    track.rho_block(first, buf);
    for (std::size_t i = 0; i < buf.size(); ++i) {
        ASSERT_EQ(buf[i], exact_rho(kSegLen - 150 + static_cast<int64_t>(i))) << "i = " << i;
    }
}

// ---- Evaluation: realistic rate, power form ------------------------------------------------

TEST(RangeTrackEval, PowerFormAt7p68Msps) {
    // Non-dyadic coefficients at 7.68 Msps: Horner vs. the power form, both float64, agree to
    // a few ulp; the bound is 1e-15 relative (ulp of 1.2e6 m is 2.3e-10 m, i.e. 1.9e-16).
    constexpr double fs = 7'680'000.0;
    constexpr int64_t n0 = 123'456'789;
    const std::array<double, 4> c{1234567.891, 6853.37, -12.3, 0.0456};
    const RangeTrack track(fs, {{n0, n0 + 2 * 7'680'000, c}});
    for (int64_t k = 0; k < 2 * 7'680'000; k += 9973) {
        const double u = static_cast<double>(k) / fs;
        const double rho = c[0] + c[1] * u + c[2] * u * u + c[3] * u * u * u;
        const double rate = c[1] + 2.0 * c[2] * u + 3.0 * c[3] * u * u;
        ASSERT_NEAR(track.rho(n0 + k), rho, 1e-15 * std::abs(rho)) << "k = " << k;
        ASSERT_NEAR(track.rho_dot(n0 + k), rate, 1e-15 * std::abs(rate)) << "k = " << k;
    }
}

// ---- Bounds (used for the pipeline latency) -----------------------------------------------

TEST(RangeTrackBounds, CubicWithInteriorExtrema) {
    // rho = 2^20 + 9u - 6u^2 + u^3 on [0, 4) s: rho' = 3(u-1)(u-3).
    // Local max 2^20 + 4 at u = 1, local min 2^20 at u = 3 (and at u = 0).
    const RangeTrack track(kFsExact, {{0, 4096, {1048576.0, 9.0, -6.0, 1.0}}});
    const auto b = track.rho_bounds(0, 4096);
    EXPECT_EQ(b.max, 1048580.0);
    EXPECT_EQ(b.min, 1048576.0);
}

TEST(RangeTrackBounds, SubrangeWithoutCriticalPointUsesEndpoints) {
    // Same cubic, samples [1536, 2560): u in [1.5, 2.5 - 1/1024], strictly decreasing.
    // Closed form at the endpoints (exact dyadics): rho(1.5) = 2^20 + 3.375,
    // rho(2559/1024) = 2^20 + 9u - 6u^2 + u^3 evaluated exactly below.
    const RangeTrack track(kFsExact, {{0, 4096, {1048576.0, 9.0, -6.0, 1.0}}});
    const auto b = track.rho_bounds(1536, 2560);
    const double u = 2559.0 / 1024.0;
    const double rho_end = 1048576.0 + (9.0 * u - 6.0 * u * u + u * u * u);  // exact
    EXPECT_EQ(b.max, 1048579.375);
    EXPECT_EQ(b.min, rho_end);
}

TEST(RangeTrackBounds, SpansSegments) {
    // Exact track: P decreasing on [0, 8) s (P' < 0 throughout), so the bounds are
    // P(0) = 2^20 and P at the last sample.
    const RangeTrack track(kFsExact, exact_segments());
    const auto b = track.rho_bounds(kN0, kN0 + 2 * kSegLen);
    EXPECT_EQ(b.max, 1048576.0);
    EXPECT_EQ(b.min, exact_rho(2 * kSegLen - 1));
}

TEST(RangeTrackBounds, InvalidRangeThrows) {
    const RangeTrack track(kFsExact, exact_segments());
    EXPECT_THROW((void)track.rho_bounds(kN0, kN0), std::out_of_range);
    EXPECT_THROW((void)track.rho_bounds(kN0 - 1, kN0 + 10), std::out_of_range);
    EXPECT_THROW((void)track.rho_bounds(kN0, kN0 + 2 * kSegLen + 1), std::out_of_range);
}
