#pragma once

// Channel pipeline (SPEC §5.1). Order: TDL (M3) -> geometric delay/Doppler -> AWGN (M4).
// M1b: geometric stage only.
//
// Streaming: process() accepts any block size and gives the same output for any partition
// (bit-identical on the scalar path). out[i] corresponds to channel time n_in + i - latency(),
// where n_in = sample_index() before the call.

#include <cstdint>
#include <span>

#include "ntn/geometric_channel.hpp"
#include "ntn/range_track.hpp"
#include "ntn/types.hpp"

namespace ntn {

class Pipeline {
 public:
    // Kernel::Scalar is the M1b default; SPEC's default becomes Kernel::Avx2 in M2.
    Pipeline(RangeTrack track, GeoConfig geo, Kernel kernel = Kernel::Scalar);

    void process(std::span<const cf32> in, std::span<cf32> out);
    int64_t sample_index() const;  // global index of the next input sample
    void reset(int64_t n_start);   // empty (zero) history; next input sample has index n_start
    int64_t latency() const;       // samples, see geometric_channel.hpp

 private:
    GeometricChannel geo_;
};

}  // namespace ntn
