#include "ntn/cpu_check.hpp"

#include <intrin.h>  // __cpuidex, _xgetbv: CPU feature queries, not SIMD code

#include <stdexcept>

namespace ntn {

namespace {

bool detect_avx2() {
    int r[4];  // eax, ebx, ecx, edx
    __cpuidex(r, 0, 0);
    if (r[0] < 7) return false;  // leaf 7 (extended features) not available

    __cpuidex(r, 1, 0);
    const bool osxsave = (r[2] >> 27) & 1;  // OS uses XSAVE/XGETBV
    const bool avx = (r[2] >> 28) & 1;
    if (!osxsave || !avx) return false;
    // XCR0 bits 1 (SSE) and 2 (AVX/YMM): the OS saves YMM state on context switches.
    if ((_xgetbv(0) & 0x6) != 0x6) return false;

    __cpuidex(r, 7, 0);
    return (r[1] >> 5) & 1;  // leaf 7, ebx bit 5: AVX2
}

}  // namespace

bool cpu_has_avx2() {
    static const bool has = detect_avx2();
    return has;
}

void require_avx2() {
    if (!cpu_has_avx2()) {
        throw std::runtime_error(
            "ntn: this CPU (or OS) does not support AVX2, which the engine requires; "
            "the scalar kernels remain available via Kernel::Scalar");
    }
}

}  // namespace ntn
