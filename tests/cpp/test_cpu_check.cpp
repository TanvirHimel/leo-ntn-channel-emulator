// cpu_check: CPUID/XGETBV detection against an independent oracle, the Windows API
// IsProcessorFeaturePresent(PF_AVX2_INSTRUCTIONS_AVAILABLE), which also reports OS support.
#include "ntn/cpu_check.hpp"

#include <gtest/gtest.h>

#include <cstdio>
#include <stdexcept>

#define NOMINMAX
#define WIN32_LEAN_AND_MEAN
#include <windows.h>

#ifndef PF_AVX2_INSTRUCTIONS_AVAILABLE
#define PF_AVX2_INSTRUCTIONS_AVAILABLE 40
#endif

TEST(CpuCheck, MatchesWindowsFeatureQuery) {
    const bool os = IsProcessorFeaturePresent(PF_AVX2_INSTRUCTIONS_AVAILABLE) != 0;
    const bool ours = ntn::cpu_has_avx2();
    std::printf("[ measure  ] cpu_has_avx2() = %d, IsProcessorFeaturePresent(AVX2) = %d\n", ours, os);
    EXPECT_EQ(ours, os);
    EXPECT_EQ(ntn::cpu_has_avx2(), ours);  // stable across calls
}

TEST(CpuCheck, RequireThrowsIffMissing) {
    if (ntn::cpu_has_avx2()) {
        EXPECT_NO_THROW(ntn::require_avx2());
    } else {
        EXPECT_THROW(ntn::require_avx2(), std::runtime_error);
    }
}
