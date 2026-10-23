#pragma once

// AVX2 CPUID guard (SPEC §6). Compiled in the scalar target (no /arch:AVX2), so the check
// itself runs on any x64 CPU. M1b: detection only; M2 enforces it at import/startup.

namespace ntn {

// True if the CPU supports AVX2 and the OS saves/restores YMM state (CPUID + XGETBV).
bool cpu_has_avx2();

// Throws std::runtime_error with a clear message if cpu_has_avx2() is false.
void require_avx2();

}  // namespace ntn
