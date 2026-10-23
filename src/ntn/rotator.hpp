#pragma once

// Carrier phase from the range track (SPEC §2 "Phase"):
//     io[i] *= exp(-j 2 pi frac(fc * rho[i] / c))
// Carrier cycles fc*rho/c are computed in float64 and reduced to [-0.5, 0.5] by
// nearbyint (half-even, like np.round), and only then cast to float32 for sincos.
// Phase is never accumulated; Doppler follows from rho alone.

#include <cstddef>

#include "ntn/types.hpp"

namespace ntn {

void rotate_scalar(const double* rho, double fc_hz, cf32* io, std::size_t count);

}  // namespace ntn
