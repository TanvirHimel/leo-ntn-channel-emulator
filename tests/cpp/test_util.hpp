#pragma once

// Test helpers: slow-test gating, a float64 radix-2 FFT for exact FFT-domain references,
// and relative error power.

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdlib>
#include <numbers>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace ntn_test {

using cd = std::complex<double>;

// Slow variants run only with NTN_SLOW=1 in the environment (ctest label `slow`).
inline bool slow_enabled() {
    char* v = nullptr;
    std::size_t len = 0;
    const bool on = _dupenv_s(&v, &len, "NTN_SLOW") == 0 && v != nullptr && std::string_view(v) == "1";
    std::free(v);
    return on;
}

// Sets an environment variable for the guard's lifetime, then restores (or removes) it.
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

#define NTN_REQUIRE_SLOW()                                                   \
    do {                                                                     \
        if (!ntn_test::slow_enabled()) GTEST_SKIP() << "slow: set NTN_SLOW=1"; \
    } while (0)

// In-place iterative radix-2 FFT, X[k] = sum_n x[n] exp(-+j 2 pi k n / N) (sign = -1 forward).
// Inverse (sign = +1) is unscaled; callers divide by N.
inline void fft_inplace(std::vector<cd>& a, int sign) {
    const std::size_t n = a.size();
    if (n == 0 || (n & (n - 1)) != 0) throw std::invalid_argument("fft: size must be a power of 2");
    for (std::size_t i = 1, j = 0; i < n; ++i) {
        std::size_t bit = n >> 1;
        for (; j & bit; bit >>= 1) j ^= bit;
        j ^= bit;
        if (i < j) std::swap(a[i], a[j]);
    }
    for (std::size_t len = 2; len <= n; len <<= 1) {
        for (std::size_t i = 0; i < n; i += len) {
            for (std::size_t k = 0; k < len / 2; ++k) {
                // Twiddles computed directly (no recurrence) to keep float64 accuracy.
                const double ang = sign * 2.0 * std::numbers::pi * static_cast<double>(k) /
                                   static_cast<double>(len);
                const cd w(std::cos(ang), std::sin(ang));
                const cd u = a[i + k];
                const cd v = a[i + k + len / 2] * w;
                a[i + k] = u + v;
                a[i + k + len / 2] = u - v;
            }
        }
    }
}

// 10 log10( sum |y - ref|^2 / sum |ref|^2 ).
template <class A, class B>
double error_db(std::span<const A> y, std::span<const B> ref) {
    double num = 0.0, den = 0.0;
    for (std::size_t i = 0; i < ref.size(); ++i) {
        const cd r(ref[i].real(), ref[i].imag());
        const cd e = cd(y[i].real(), y[i].imag()) - r;
        num += std::norm(e);
        den += std::norm(r);
    }
    return 10.0 * std::log10(num / den);
}

}  // namespace ntn_test
