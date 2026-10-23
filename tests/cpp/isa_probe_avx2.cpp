// ISA probe for scripts/check_isa.py. Compiled with the same flags as ntn_avx2 (/arch:AVX2,
// /fp:precise) but never linked into anything.
//  - ntn_isa_probe_ymm: explicit 8-lane xsimd code, so its disassembly must contain ymm
//    registers (positive control for the "no ymm in the scalar objects" check).
//  - ntn_isa_probe_contract_*: plain a*b + c expressions, so the disassembly shows whether
//    MSVC contracts them to FMA under /fp:precise with /arch:AVX2 (M2 step 5).
#include <cstddef>

#include <xsimd/xsimd.hpp>

#if !defined(__AVX2__) || !XSIMD_WITH_AVX2
#error "isa_probe_avx2.cpp must be compiled with /arch:AVX2"
#endif

extern "C" {

void ntn_isa_probe_ymm(const float* a, const float* b, float* out, std::size_t n) {
    using batch = xsimd::batch<float, xsimd::avx2>;
    for (std::size_t i = 0; i + batch::size <= n; i += batch::size) {
        (batch::load_unaligned(a + i) * batch::load_unaligned(b + i)).store_unaligned(out + i);
    }
}

float ntn_isa_probe_contract_f32(float a, float b, float c) { return a * b + c; }

double ntn_isa_probe_contract_f64(double a, double b, double c) { return a * b + c; }

void ntn_isa_probe_contract_loop(const float* a, const float* b, float* acc, std::size_t n) {
    for (std::size_t i = 0; i < n; ++i) acc[i] = a[i] * b[i] + acc[i];
}

}  // extern "C"
