#include "ntn/frac_delay.hpp"

#include <algorithm>
#include <cmath>
#include <numbers>

namespace ntn {

namespace {

// Kaiser window of half-width kHalf at offset t (samples); zero beyond.
double kaiser(double t) {
    const double r = 1.0 - (t / fd::kHalf) * (t / fd::kHalf);
    static const double i0_beta = std::cyl_bessel_i(0.0, fd::kBeta);
    return std::cyl_bessel_i(0.0, fd::kBeta * std::sqrt(std::max(r, 0.0))) / i0_beta;
}

}  // namespace

FracDelayTable::FracDelayTable() : taps_(static_cast<std::size_t>(fd::kPhases + 1) * fd::kTaps, 0.0f) {
    constexpr double pi = std::numbers::pi;
    for (int r = 0; r <= fd::kPhases; ++r) {
        float* row = &taps_[static_cast<std::size_t>(r) * fd::kTaps];
        if (r == 0 || r == fd::kPhases) {
            // Exact unit impulse at f = 0 (tap off = 0) and f = 1 (tap off = 1).
            row[r == 0 ? fd::kHalf - 1 : fd::kHalf] = 1.0f;
            continue;
        }
        const double f = static_cast<double>(r) / fd::kPhases;
        // sin(pi (f - off)) = (-1)^off sin(pi f); sin(pi f) = sin(pi (1 - f)) keeps rows f and
        // 1 - f mirror images (and row 1/2 exactly symmetric).
        const double sin_pf = std::sin(pi * std::min(f, 1.0 - f));
        for (int j = 0; j < fd::kTaps; ++j) {
            const int off = j - (fd::kHalf - 1);
            const double t = f - off;  // never 0 here: 0 < f < 1
            const double num = (off % 2 == 0) ? sin_pf : -sin_pf;
            row[j] = static_cast<float>(num / (pi * t) * kaiser(t));
        }
    }
}

const FracDelayTable& frac_delay_table() {
    static const FracDelayTable table;
    return table;
}

void frac_delay_scalar(const FracDelayTable& table, const cf32* hist, int64_t hist_first, int64_t n_first,
                       const double* d, cf32* out, std::size_t count) {
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
        const double p = f * fd::kPhases;
        const int r = std::min(static_cast<int>(p), fd::kPhases - 1);
        const float w = static_cast<float>(p - r);
        const float* a = table.row(r);
        const float* b = table.row(r + 1);
        const cf32* x = hist + (m0 - (fd::kHalf - 1) - hist_first);
        float acc_re = 0.0f, acc_im = 0.0f;
        for (int j = 0; j < fd::kTaps; ++j) {
            const float c = a[j] + w * (b[j] - a[j]);
            acc_re += x[j].real() * c;
            acc_im += x[j].imag() * c;
        }
        out[i] = {acc_re, acc_im};
    }
}

}  // namespace ntn
