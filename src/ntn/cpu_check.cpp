#include "ntn/cpu_check.hpp"

#include <intrin.h>  // __cpuidex, _xgetbv: CPU feature queries, not SIMD code

#include <cstdlib>
#include <stdexcept>
#include <string_view>

namespace ntn {

namespace {

CpuFeatures detect() {
    CpuFeatures f;
    int r[4];  // eax, ebx, ecx, edx
    __cpuidex(r, 0, 0);
    const int max_leaf = r[0];

    __cpuidex(r, 1, 0);
    const bool avx = (r[2] >> 28) & 1;
    const bool osxsave = (r[2] >> 27) & 1;
    f.fma3 = (r[2] >> 12) & 1;
    // XCR0 bits 1 (SSE) and 2 (AVX/YMM): the OS saves YMM state on context switches.
    f.os_avx = osxsave && (_xgetbv(0) & 0x6) == 0x6;

    if (max_leaf >= 7) {
        __cpuidex(r, 7, 0);
        f.avx2 = avx && ((r[1] >> 5) & 1);
        f.bmi1 = (r[1] >> 3) & 1;
        f.bmi2 = (r[1] >> 8) & 1;
    }
    return f;
}

// Test-only switch (cpu_check.hpp): clears one feature; any other value is ignored.
void apply_test_hide(CpuFeatures& f) {
    char* v = nullptr;
    std::size_t len = 0;
    if (_dupenv_s(&v, &len, "NTN_TEST_HIDE_AVX2") != 0 || v == nullptr) return;
    const std::string_view s(v);
    if (s == "1" || s == "AVX2") f.avx2 = false;
    else if (s == "FMA3") f.fma3 = false;
    else if (s == "BMI1") f.bmi1 = false;
    else if (s == "BMI2") f.bmi2 = false;
    else if (s == "OSAVX") f.os_avx = false;
    std::free(v);
}

}  // namespace

CpuFeatures cpu_features() {
    static const CpuFeatures detected = detect();
    CpuFeatures f = detected;
    apply_test_hide(f);
    return f;
}

std::vector<std::string> cpu_missing_features() {
    const CpuFeatures f = cpu_features();
    std::vector<std::string> m;
    if (!f.avx2) m.emplace_back("AVX2");
    if (!f.fma3) m.emplace_back("FMA3");
    if (!f.bmi1) m.emplace_back("BMI1");
    if (!f.bmi2) m.emplace_back("BMI2");
    if (!f.os_avx) m.emplace_back("OS AVX state (XGETBV)");
    return m;
}

bool cpu_has_avx2() { return cpu_missing_features().empty(); }

void require_avx2() {
    const std::vector<std::string> missing = cpu_missing_features();
    if (missing.empty()) return;
    std::string list;
    for (const auto& m : missing) list += (list.empty() ? "" : ", ") + m;
    throw std::runtime_error("ntn: this CPU/OS lacks instruction-set features the engine is built for "
                             "(missing: " + list + "; required: AVX2, FMA3, BMI1, BMI2, OS AVX state (XGETBV))");
}

}  // namespace ntn
