#pragma once

// Channel pipeline (SPEC §5.1). Order: TDL (M3) -> geometric delay/Doppler -> AWGN (M4).
// M1b: geometric stage only.
//
// Streaming: process() accepts any block size and gives the same output for any partition
// (bit-identical on the scalar path). out[i] corresponds to channel time n_in + i - latency(),
// where n_in = sample_index() before the call.

#include <cstdint>
#include <span>

#include "ntn/cpu_check.hpp"
#include "ntn/geometric_channel.hpp"
#include "ntn/range_track.hpp"
#include "ntn/types.hpp"

namespace ntn {

class Pipeline : private RequireAvx2 {
 public:
    // Kernel::Avx2 is the default (SPEC §5.1). Throws std::runtime_error (RequireAvx2, before
    // anything else) if the CPU lacks a feature the build needs, for either kernel.
    Pipeline(RangeTrack track, GeoConfig geo, Kernel kernel = Kernel::Avx2);

    void process(std::span<const cf32> in, std::span<cf32> out);
    int64_t sample_index() const;  // global index of the next input sample
    void reset(int64_t n_start);   // empty (zero) history; next input sample has index n_start
    int64_t latency() const;       // samples, see geometric_channel.hpp
    Kernel kernel() const { return geo_.kernel(); }

 private:
    GeometricChannel geo_;
};

}  // namespace ntn
