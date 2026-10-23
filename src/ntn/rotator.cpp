#include "ntn/rotator.hpp"

#include <cmath>

#include "ntn/constants.hpp"

namespace ntn {

void rotate_scalar(const double* rho, double fc_hz, cf32* io, std::size_t count) {
    for (std::size_t i = 0; i < count; ++i) {
        const double cycles = fc_hz * rho[i] / kC;              // float64, ~1e7 cycles
        const double frac = cycles - std::nearbyint(cycles);   // [-0.5, 0.5], half-even
        const auto ang = static_cast<float>(-kTwoPi * frac);   // wrapped, then float32
        const float c = std::cos(ang), s = std::sin(ang);
        const float xr = io[i].real(), xi = io[i].imag();
        io[i] = {xr * c - xi * s, xr * s + xi * c};
    }
}

}  // namespace ntn
