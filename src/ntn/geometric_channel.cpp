#include "ntn/geometric_channel.hpp"

#include <algorithm>
#include <cmath>
#include <format>
#include <stdexcept>
#include <utility>

#include "ntn/constants.hpp"
#include "ntn/rotator.hpp"

namespace ntn {

namespace {
constexpr std::size_t kChunk = 4096;  // internal scratch size; does not affect any output bit
constexpr double kBoundGuard = 1e-6;  // samples; absorbs float64 rounding of the D bounds
}  // namespace

GeometricChannel::GeometricChannel(RangeTrack track, GeoConfig cfg, Kernel kernel)
    : track_(std::move(track)), cfg_(cfg), kernel_(kernel), table_(frac_delay_table()),
      rho_(kChunk), d_(kChunk) {
    if (cfg_.fs_hz != track_.fs()) {
        throw std::invalid_argument(
            std::format("GeometricChannel: GeoConfig.fs_hz ({}) differs from the track fs ({})", cfg_.fs_hz, track_.fs()));
    }
    if (!std::isfinite(cfg_.fc_hz)) throw std::invalid_argument("GeometricChannel: fc_hz must be finite");
    if (cfg_.residual_ref_doppler_hz != 0.0) {
        throw std::invalid_argument("GeometricChannel: residual-CFO mode is not implemented before M4");
    }
    if (kernel_ != Kernel::Scalar && kernel_ != Kernel::Avx2) {
        throw std::invalid_argument("GeometricChannel: unknown Kernel value");
    }

    if (cfg_.delay) {
        // D(n) = rho(n) fs / c - tau0. Output n reads x[n - floor(D) - kHalf .. n - floor(D) + kHalf - 1]
        // (frac_delay.hpp), so the latency covers the highest read and the history the lowest.
        const RangeBounds b = track_.rho_bounds(track_.n_begin(), track_.n_end());
        const double tau0 = static_cast<double>(cfg_.tau0_samples);
        const double d_min = b.min * cfg_.fs_hz / kC - tau0;
        const double d_max = b.max * cfg_.fs_hz / kC - tau0;
        const auto fl_min = static_cast<int64_t>(std::floor(d_min - kBoundGuard));
        const auto fl_max = static_cast<int64_t>(std::floor(d_max + kBoundGuard));
        latency_ = std::max<int64_t>(0, fd::kHalf - 1 - fl_min);
        keep_behind_ = std::max<int64_t>(0, fl_max + fd::kHalf);
    }
    reset(track_.n_begin());
}

void GeometricChannel::reset(int64_t n_start) {
    n_in_ = n_start;
    if (cfg_.delay) {
        // Zero history for every read of the first output (channel time n_start - L).
        hist_first_ = n_start - latency_ - keep_behind_;
        hist_.assign(static_cast<std::size_t>(latency_ + keep_behind_), cf32(0.0f, 0.0f));
    }
}

void GeometricChannel::process(std::span<const cf32> in, std::span<cf32> out) {
    if (in.size() != out.size()) {
        throw std::invalid_argument(
            std::format("GeometricChannel: input size {} != output size {}", in.size(), out.size()));
    }
    if (in.empty()) return;
    const auto count = static_cast<int64_t>(in.size());
    const int64_t n_out_first = n_in_ - latency_;
    if (n_out_first + count > track_.n_end()) {
        throw std::out_of_range(std::format("GeometricChannel: channel times up to {} exceed the track end {}",
                                            n_out_first + count - 1, track_.n_end()));
    }

    if (cfg_.delay) {
        hist_.insert(hist_.end(), in.begin(), in.end());
    } else {
        std::copy(in.begin(), in.end(), out.begin());
    }
    compute(n_out_first, out);
    n_in_ += count;

    if (cfg_.delay) {
        // Drop history no later output can read; compact only when it is a large share.
        const int64_t drop = (n_in_ - latency_ - keep_behind_) - hist_first_;
        if (drop > 0 && static_cast<std::size_t>(drop) >= std::max<std::size_t>(kChunk, hist_.size() / 2)) {
            hist_.erase(hist_.begin(), hist_.begin() + drop);
            hist_first_ += drop;
        }
    }
}

// Channel outputs for times n_first .. n_first + out.size() - 1. With delay off, `out` holds the
// input on entry.
void GeometricChannel::compute(int64_t n_first, std::span<cf32> out) {
    std::size_t i = 0;
    int64_t n = n_first;
    for (; i < out.size() && n < track_.n_begin(); ++i, ++n) out[i] = cf32(0.0f, 0.0f);

    const double tau0 = static_cast<double>(cfg_.tau0_samples);
    while (i < out.size()) {
        const std::size_t m = std::min(kChunk, out.size() - i);
        const std::span<double> rho(rho_.data(), m);
        track_.rho_block(n, rho);  // one rho per sample: delay and phase both derive from it
        if (cfg_.delay) {
            for (std::size_t k = 0; k < m; ++k) d_[k] = rho[k] * cfg_.fs_hz / kC - tau0;
            if (kernel_ == Kernel::Avx2) {
                frac_delay_avx2(table_, hist_.data(), hist_first_, n, d_.data(), out.data() + i, m);
            } else {
                frac_delay_scalar(table_, hist_.data(), hist_first_, n, d_.data(), out.data() + i, m);
            }
        }
        if (cfg_.doppler) {
            if (kernel_ == Kernel::Avx2) {
                rotate_avx2(rho.data(), cfg_.fc_hz, out.data() + i, m);
            } else {
                rotate_scalar(rho.data(), cfg_.fc_hz, out.data() + i, m);
            }
        }
        i += m;
        n += static_cast<int64_t>(m);
    }
}

}  // namespace ntn
