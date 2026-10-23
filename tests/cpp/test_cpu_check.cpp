// cpu_check: CPUID/XGETBV detection against an independent oracle, the Windows API
// IsProcessorFeaturePresent(PF_AVX2_INSTRUCTIONS_AVAILABLE), which also reports OS support.
#include "ntn/cpu_check.hpp"

#include <gtest/gtest.h>

#include <cstdio>
#include <cstdlib>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

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

// ---- Exact feature set (M2 step 4) ------------------------------------------------------------
//
// cpu_has_avx2() means "everything the AVX2 build needs": AVX2, FMA3, BMI1, BMI2 and OS support
// for AVX state. The test-only switch NTN_TEST_HIDE_AVX2 hides one feature: "1" or "AVX2" hides
// AVX2; "FMA3", "BMI1", "BMI2", "OSAVX" hide that feature; any other value is ignored.


namespace {

class EnvGuard {
 public:
    EnvGuard(const char* name, const char* value) : name_(name) {
        char* v = nullptr;
        std::size_t len = 0;
        if (_dupenv_s(&v, &len, name) == 0 && v != nullptr) old_ = std::string(v);
        std::free(v);
        _putenv_s(name, value);
    }
    ~EnvGuard() { _putenv_s(name_.c_str(), old_ ? old_->c_str() : ""); }
    EnvGuard(const EnvGuard&) = delete;
    EnvGuard& operator=(const EnvGuard&) = delete;

 private:
    std::string name_;
    std::optional<std::string> old_;
};

constexpr const char* kSwitch = "NTN_TEST_HIDE_AVX2";

std::string require_message() {
    try {
        ntn::require_avx2();
    } catch (const std::runtime_error& e) {
        return e.what();
    }
    return {};
}

#define NTN_REQUIRE_REAL_FEATURES()                                                          \
    do {                                                                                     \
        if (!ntn::cpu_has_avx2()) GTEST_SKIP() << "this CPU lacks a feature; hiding is moot"; \
    } while (0)

}  // namespace

TEST(CpuCheck, Avx2AndOsBitsMatchWindowsQuery) {
    // PF_AVX2_INSTRUCTIONS_AVAILABLE covers the AVX2 instructions and OS support for their state.
    const bool os = IsProcessorFeaturePresent(PF_AVX2_INSTRUCTIONS_AVAILABLE) != 0;
    const ntn::CpuFeatures f = ntn::cpu_features();
    std::printf("[ measure  ] cpu features: AVX2 %d, FMA3 %d, BMI1 %d, BMI2 %d, OS AVX state %d; "
                "IsProcessorFeaturePresent(AVX2) %d\n",
                f.avx2, f.fma3, f.bmi1, f.bmi2, f.os_avx, os);
    EXPECT_EQ(f.avx2 && f.os_avx, os);
    EXPECT_EQ(ntn::cpu_has_avx2(), f.avx2 && f.fma3 && f.bmi1 && f.bmi2 && f.os_avx);
}

TEST(CpuCheck, HideSwitchEachFeatureIsReportedByName) {
    NTN_REQUIRE_REAL_FEATURES();
    struct Case {
        const char* value;
        const char* name;  // as it must appear in the error message
    };
    const Case cases[] = {{"1", "AVX2"},    {"AVX2", "AVX2"}, {"FMA3", "FMA3"},
                          {"BMI1", "BMI1"}, {"BMI2", "BMI2"}, {"OSAVX", "OS AVX state (XGETBV)"}};
    for (const auto& c : cases) {
        const EnvGuard g(kSwitch, c.value);
        EXPECT_FALSE(ntn::cpu_has_avx2()) << c.value;
        EXPECT_EQ(ntn::cpu_missing_features(), std::vector<std::string>{c.name}) << c.value;
        const std::string msg = require_message();
        EXPECT_NE(msg.find(std::string("missing: ") + c.name), std::string::npos) << c.value << ": " << msg;
        // The message also states the full requirement.
        EXPECT_NE(msg.find("AVX2, FMA3, BMI1, BMI2, OS AVX state (XGETBV)"), std::string::npos) << msg;
        std::printf("[ measure  ] %s=%s -> \"%s\"\n", kSwitch, c.value, msg.c_str());
    }
    EXPECT_TRUE(ntn::cpu_has_avx2());  // switch cleared again
}

TEST(CpuCheck, HideSwitchOtherValuesAreIgnored) {
    NTN_REQUIRE_REAL_FEATURES();
    for (const char* v : {"0", "2", "true", "yes", "fma3", "Fma3", " FMA3", "FMA3 ", "AVX512", "1,FMA3",
                          "FMA3,BMI2", "FMA", "BMI", "OS", "AVX"}) {
        const EnvGuard g(kSwitch, v);
        EXPECT_TRUE(ntn::cpu_has_avx2()) << "'" << v << "'";
        EXPECT_TRUE(ntn::cpu_missing_features().empty()) << "'" << v << "'";
        EXPECT_NO_THROW(ntn::require_avx2()) << "'" << v << "'";
    }
}

TEST(CpuCheck, MessageDoesNotOfferAScalarFallback) {
    // ntn_core is compiled with /arch:AVX2, so without these features nothing in the engine can
    // run, Kernel::Scalar included; the message must not suggest otherwise.
    NTN_REQUIRE_REAL_FEATURES();
    const EnvGuard g(kSwitch, "BMI2");
    const std::string msg = require_message();
    ASSERT_FALSE(msg.empty()) << "require_avx2() did not throw with BMI2 hidden";
    EXPECT_EQ(msg.find("Scalar"), std::string::npos) << msg;
}
