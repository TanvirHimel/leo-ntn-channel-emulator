#include "ntn/frac_delay.hpp"

#include <algorithm>
#include <cmath>
#include <numbers>

namespace ntn {

namespace fd {

double kaiser_design_atten_db(int taps) {
    constexpr double dw = 0.4 * std::numbers::pi;  // transition 0.4 -> 0.6 fs
    return 2.285 * (taps - 1) * dw + 8.0;
}

double kaiser_beta(int taps) {
    const double a = kaiser_design_atten_db(taps);
    if (a > 50.0) return 0.1102 * (a - 8.7);
    if (a >= 21.0) return 0.5842 * std::pow(a - 21.0, 0.4) + 0.07886 * (a - 21.0);
    return 0.0;
}

}  // namespace fd

template <int T, int P>
FracDelayTableT<T, P>::FracDelayTableT()
    : beta_(fd::kaiser_beta(T)), taps_(static_cast<std::size_t>(P + 1) * T, 0.0f) {
    constexpr double pi = std::numbers::pi;
    constexpr int half = T / 2;
    const double i0_beta = std::cyl_bessel_i(0.0, beta_);
    // Kaiser window of half-width `half` at offset t (samples); zero beyond.
    const auto kaiser = [&](double t) {
        const double r = 1.0 - (t / half) * (t / half);
        return std::cyl_bessel_i(0.0, beta_ * std::sqrt(std::max(r, 0.0))) / i0_beta;
    };
    for (int r = 0; r <= P; ++r) {
        float* row = &taps_[static_cast<std::size_t>(r) * T];
        if (r == 0 || r == P) {
            // Exact unit impulse at f = 0 (tap off = 0) and f = 1 (tap off = 1).
            row[r == 0 ? half - 1 : half] = 1.0f;
            continue;
        }
        const double f = static_cast<double>(r) / P;
        // sin(pi (f - off)) = (-1)^off sin(pi f); sin(pi f) = sin(pi (1 - f)) keeps rows f and
        // 1 - f mirror images (and row 1/2 exactly symmetric).
        const double sin_pf = std::sin(pi * std::min(f, 1.0 - f));
        for (int j = 0; j < T; ++j) {
            const int off = j - (half - 1);
            const double t = f - off;  // never 0 here: 0 < f < 1
            const double num = (off % 2 == 0) ? sin_pf : -sin_pf;
            row[j] = static_cast<float>(num / (pi * t) * kaiser(t));
        }
    }
}

template <int T, int P>
const FracDelayTableT<T, P>& frac_delay_table_t() {
    static const FracDelayTableT<T, P> table;
    return table;
}

template <int T, int P>
void frac_delay_scalar_t(const FracDelayTableT<T, P>& table, const cf32* hist, int64_t hist_first, int64_t n_first,
                         const double* d, cf32* out, std::size_t count) {
    constexpr int half = T / 2;
    for (std::size_t i = 0; i < count; ++i) {
        const int64_t n = n_first + static_cast<int64_t>(i);
        // Split before forming n - d, so precision is set by |d|, not by n (golden rule).
        double d_int = std::floor(d[i]);
        double d_frac = d[i] - d_int;
        if (d_frac > 1.0 - fd::kSnap) {
            d_int += 1.0;
            d_frac = 0.0;
        } else if (d_frac < fd::kSnap) {
            d_frac = 0.0;
        }
        int64_t m0 = n - static_cast<int64_t>(d_int);
        if (d_frac == 0.0) {  // integer delay: pure gather, bit-exact
            out[i] = hist[m0 - hist_first];
            continue;
        }
        m0 -= 1;
        const double f = 1.0 - d_frac;  // read position m0 + f, f in (0, 1)
        const double p = f * P;
        const int r = std::min(static_cast<int>(p), P - 1);
        const float w = static_cast<float>(p - r);
        const float* a = table.row(r);
        const float* b = table.row(r + 1);
        const cf32* x = hist + (m0 - (half - 1) - hist_first);
        float acc_re = 0.0f, acc_im = 0.0f;
        for (int j = 0; j < T; ++j) {
            const float c = a[j] + w * (b[j] - a[j]);
            acc_re += x[j].real() * c;
            acc_im += x[j].imag() * c;
        }
        out[i] = {acc_re, acc_im};
    }
}

const FracDelayTable& frac_delay_table() { return frac_delay_table_t<fd::kTaps, fd::kPhases>(); }

void frac_delay_scalar(const FracDelayTable& table, const cf32* hist, int64_t hist_first, int64_t n_first,
                       const double* d, cf32* out, std::size_t count) {
    frac_delay_scalar_t<fd::kTaps, fd::kPhases>(table, hist, hist_first, n_first, d, out, count);
}

#define NTN_FD_DESIGN(T, P)                                                                    \
    template class FracDelayTableT<T, P>;                                                      \
    template const FracDelayTableT<T, P>& frac_delay_table_t<T, P>();                          \
    template void frac_delay_scalar_t<T, P>(const FracDelayTableT<T, P>&, const cf32*, int64_t, \
                                            int64_t, const double*, cf32*, std::size_t);
NTN_FD_DESIGN(16, 128)
NTN_FD_DESIGN(16, 256)
NTN_FD_DESIGN(24, 128)
NTN_FD_DESIGN(24, 256)
NTN_FD_DESIGN(32, 128)
NTN_FD_DESIGN(32, 256)
#undef NTN_FD_DESIGN

}  // namespace ntn
