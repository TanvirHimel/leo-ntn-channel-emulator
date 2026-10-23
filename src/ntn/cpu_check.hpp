#pragma once

// CPU guard for the AVX2 build (SPEC §6). Compiled in the scalar target (no /arch:AVX2), so the
// check itself runs on any x64 CPU.
//
// "AVX2" here means everything the /arch:AVX2 objects (ntn_avx2, ntn_core) can execute:
//   AVX2   CPUID.1:ECX.AVX[bit 28] and CPUID.(7,0):EBX.AVX2[bit 5]
//   FMA3   CPUID.1:ECX.FMA[bit 12]      (xsimd enables fma3<avx2> under MSVC /arch:AVX2)
//   BMI1   CPUID.(7,0):EBX.BMI1[bit 3]  (MSVC emits andn etc. under /arch:AVX2)
//   BMI2   CPUID.(7,0):EBX.BMI2[bit 8]  (MSVC emits shlx/shrx/bzhi under /arch:AVX2)
//   OS AVX state (XGETBV): CPUID.1:ECX.OSXSAVE[bit 27] and XCR0 bits 1-2 set (the OS saves
//          XMM/YMM state); XGETBV is executed only when OSXSAVE is set.
// scripts/check_isa.py fails if an instruction family outside this set appears in the AVX2
// objects. ntn_core itself is compiled with /arch:AVX2, so without these features no part of
// the engine can run, the scalar kernels included.
//
// Test-only switch: environment variable NTN_TEST_HIDE_AVX2 hides one detected feature, so the
// error path can be tested on a capable machine. Values (exact, case-sensitive): "1" or "AVX2",
// "FMA3", "BMI1", "BMI2", "OSAVX". Any other value is ignored. Never set it outside tests.

#include <string>
#include <vector>

namespace ntn {

struct CpuFeatures {
    bool avx2 = false;
    bool fma3 = false;
    bool bmi1 = false;
    bool bmi2 = false;
    bool os_avx = false;
};

// Detected features (CPUID/XGETBV, cached), with the test-only switch applied.
CpuFeatures cpu_features();

// Names of the missing features, in the order AVX2, FMA3, BMI1, BMI2, "OS AVX state (XGETBV)".
std::vector<std::string> cpu_missing_features();

// True if every feature the AVX2 build needs is present (see above).
bool cpu_has_avx2();

// Throws std::runtime_error naming the missing features if cpu_has_avx2() is false.
void require_avx2();

}  // namespace ntn
