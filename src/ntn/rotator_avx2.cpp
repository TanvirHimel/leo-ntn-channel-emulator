// AVX2 rotator kernel (xsimd). Compiled only in target ntn_avx2, with /arch:AVX2.
// M2 step 1: build plumbing only; the kernel arrives in a later M2 step.
#include <xsimd/xsimd.hpp>

// /arch:AVX2 must be active for this file (MSVC defines __AVX2__ under it).
#if !defined(__AVX2__) || !XSIMD_WITH_AVX2
#error "rotator_avx2.cpp must be compiled with /arch:AVX2 (target ntn_avx2)"
#endif
