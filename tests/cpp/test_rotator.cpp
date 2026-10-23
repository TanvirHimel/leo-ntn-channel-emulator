// Scalar rotator: exact carrier-cycle cases with known output.
#include "ntn/rotator.hpp"

#include <gtest/gtest.h>

#include <complex>
#include <cstring>
#include <vector>

using ntn::cf32;

namespace {

// fc = c makes fc*rho/c = rho exactly for these rho (the product fc*rho is an exact double,
// < 2^53, and the quotient is then exact), so the carrier cycles are known exactly.
constexpr double kFcEqC = 299'792'458.0;

std::vector<cf32> inputs() { return {{1.0f, 0.0f}, {0.0f, 1.0f}, {0.6f, -0.8f}, {-3.25f, 1.5f}}; }

void expect_rotated(double rho, cf32 factor) {
    auto x = inputs();
    auto y = x;
    std::vector<double> r(y.size(), rho);
    ntn::rotate_scalar(r.data(), kFcEqC, y.data(), y.size());
    for (std::size_t i = 0; i < x.size(); ++i) {
        const cf32 want = x[i] * factor;
        // float32 sincos of the float32-rounded angle: error <= ~1e-7 relative.
        EXPECT_LT(std::abs(y[i] - want), 2e-7f * std::abs(x[i])) << "rho = " << rho << ", i = " << i;
    }
}

}  // namespace

TEST(Rotator, QuarterCycleGivesMinusJ) { expect_rotated(1000000.25, {0.0f, -1.0f}); }   // e^{-j pi/2}
TEST(Rotator, ThreeQuarterCycleGivesPlusJ) { expect_rotated(1000000.75, {0.0f, 1.0f}); } // frac -0.25
TEST(Rotator, HalfCycleGivesMinusOne) {
    expect_rotated(1000000.5, {-1.0f, 0.0f});  // nearbyint -> 1000000, frac +0.5
    expect_rotated(1000001.5, {-1.0f, 0.0f});  // nearbyint -> 1000002, frac -0.5
}
TEST(Rotator, IntegerCyclesIsIdentityBitExact) {
    auto x = inputs();
    auto y = x;
    std::vector<double> r(y.size(), 1234567.0);
    ntn::rotate_scalar(r.data(), kFcEqC, y.data(), y.size());
    EXPECT_EQ(std::memcmp(x.data(), y.data(), x.size() * sizeof(cf32)), 0);
}
TEST(Rotator, ZeroCarrierIsIdentityBitExact) {
    auto x = inputs();
    auto y = x;
    std::vector<double> r(y.size(), 987654.321);
    ntn::rotate_scalar(r.data(), 0.0, y.data(), y.size());
    EXPECT_EQ(std::memcmp(x.data(), y.data(), x.size() * sizeof(cf32)), 0);
}
