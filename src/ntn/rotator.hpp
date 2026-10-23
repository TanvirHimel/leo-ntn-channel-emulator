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

// AVX2 rotator (rotator_avx2.cpp, target ntn_avx2). Same carrier cycles, wrap and float32 cast
// as rotate_scalar (float64, same operation order) and the same complex multiply; only the
// float32 sincos differs (sincos_avx2 instead of the CRT). The last partial block of 8 is padded
// and run through the same vector code, so each output is independent of the block size.
// The caller must have checked the CPU (cpu_check.hpp).
void rotate_avx2(const double* rho, double fc_hz, cf32* io, std::size_t count);

// The vectorised float32 sincos rotate_avx2 uses (xsimd), exposed for its accuracy test:
// s[i] = sin(ang[i]), c[i] = cos(ang[i]). Same padded-tail rule.
void sincos_avx2(const float* ang, float* s, float* c, std::size_t count);

}  // namespace ntn
