#pragma once

// Geometric channel: delay and carrier phase from one RangeTrack (SPEC §2, §5.1).
//     y(n) = x(n - D(n)) * exp(-j 2 pi frac(fc rho(n) / c)),   D(n) = rho(n) fs / c - tau0
// Delay is relative to tau0_samples; phase uses the full rho. Delay first, then phase.
//
// Streaming with fixed latency L: process() maps input samples n_in .. n_in+M-1 to
//     out[i] = y(n_in + i - L),
//     L = max(0, 15 - floor(min D - 1e-6))   (0 with delay off; 15 = kHalf - 1)
// Derivation: the kernel reads x[n + 15 - floor(D)] at most (taps x[m0-15 .. m0+16],
// m0 = n - floor(D) - 1; an integer D reads only x[n - D]), and output n is emitted when input
// n + L arrives, so L >= 15 - floor(D(n)) for every n. The 1e-6-sample guard absorbs float64
// rounding of the bound. The bound is tight for non-integer min D (test GeoLatency.*Minimal*).
// Outputs at channel times before the track start are zero; channel times at or past the track
// end throw. Input before the reset index is zero. Each output depends only on its channel time
// and the input, so any block partition gives identical output.

#include <cstdint>
#include <span>
#include <vector>

#include "ntn/cpu_check.hpp"
#include "ntn/frac_delay.hpp"
#include "ntn/range_track.hpp"
#include "ntn/types.hpp"

namespace ntn {

struct GeoConfig {
    double fc_hz = 0.0;
    double fs_hz = 0.0;  // must equal the track's fs
    int64_t tau0_samples = 0;
    bool delay = true;
    bool doppler = true;
    double residual_ref_doppler_hz = 0.0;  // M4 residual-CFO mode; must be 0 until then
};

// Kernel selection: Avx2 (default) uses frac_delay_avx2 / rotate_avx2; Scalar uses the scalar
// kernels, kept as the test oracle and benchmark baseline (no runtime fallback: both need the
// CPU features in cpu_check.hpp, because this class is compiled with /arch:AVX2).
enum class Kernel { Scalar, Avx2 };

class GeometricChannel : private RequireAvx2 {
 public:
    // Throws std::invalid_argument on fs mismatch, non-finite fc, nonzero residual Doppler
    // (M4); std::runtime_error (from RequireAvx2, before anything else) if the CPU lacks a
    // required feature, for either kernel. Starts at the track's first sample.
    GeometricChannel(RangeTrack track, GeoConfig cfg, Kernel kernel = Kernel::Avx2);

    int64_t latency() const { return latency_; }
    Kernel kernel() const { return kernel_; }
    int64_t sample_index() const { return n_in_; }  // global index of the next input sample
    const RangeTrack& track() const { return track_; }
    const GeoConfig& config() const { return cfg_; }

    // Restart with an empty (zero) input history; the next input sample has index n_start.
    void reset(int64_t n_start);

    // in.size() must equal out.size(). Throws std::invalid_argument on size mismatch and
    // std::out_of_range (state unchanged) if a channel time reaches the track end.
    void process(std::span<const cf32> in, std::span<cf32> out);

 private:
    void compute(int64_t n_out_first, std::span<cf32> out);

    RangeTrack track_;
    GeoConfig cfg_;
    Kernel kernel_;
    const FracDelayTable& table_;
    int64_t latency_ = 0;
    int64_t keep_behind_ = 0;  // history kept behind the next channel time, samples
    int64_t n_in_ = 0;
    std::vector<cf32> hist_;
    int64_t hist_first_ = 0;   // global index of hist_[0]
    std::vector<double> rho_, d_;  // per-chunk scratch
};

}  // namespace ntn
