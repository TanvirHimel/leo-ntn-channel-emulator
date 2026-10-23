#include "ntn/pipeline.hpp"

#include <utility>

namespace ntn {

Pipeline::Pipeline(RangeTrack track, GeoConfig geo, Kernel kernel) : geo_(std::move(track), geo, kernel) {}

void Pipeline::process(std::span<const cf32> in, std::span<cf32> out) {
    // M1b: geometric stage only (TDL before it in M3, AWGN after it in M4).
    geo_.process(in, out);
}

int64_t Pipeline::sample_index() const { return geo_.sample_index(); }

void Pipeline::reset(int64_t n_start) { geo_.reset(n_start); }

int64_t Pipeline::latency() const { return geo_.latency(); }

}  // namespace ntn
