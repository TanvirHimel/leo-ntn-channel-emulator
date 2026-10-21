#pragma once

// Range track rho(n): the single source of truth for delay and carrier phase (SPEC §2).
//
// Segment k covers global samples [n0, n1) with local time u = (n - n0)/fs [s] and
//     rho(u) = c0 + c1 u + c2 u^2 + c3 u^3   [m],
// evaluated by Horner in float64. Segments must be contiguous and C1 at every interior knot
// (|d rho| < 1e-6 m, |d rho_dot| < 1e-6 m/s, same rule as python/ntn_emu/segments.check_c1).

#include <array>
#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

namespace ntn {

struct RangeSegment {
    int64_t n0, n1;            // [n0, n1) in global samples
    std::array<double, 4> c;   // m, m/s, m/s^2, m/s^3
};

struct RangeBounds {
    double min, max;  // m
};

inline double horner_rho(const std::array<double, 4>& c, double u) {
    return ((c[3] * u + c[2]) * u + c[1]) * u + c[0];
}

inline double horner_rho_dot(const std::array<double, 4>& c, double u) {
    return (3.0 * c[3] * u + 2.0 * c[2]) * u + c[1];
}

class RangeTrack {
 public:
    static constexpr double kC1TolRho = 1e-6;   // m
    static constexpr double kC1TolRate = 1e-6;  // m/s

    // Throws std::invalid_argument on fs <= 0, no segments, empty/non-contiguous segments,
    // non-finite coefficients or a C1 violation.
    RangeTrack(double fs, std::vector<RangeSegment> segments);

    double fs() const { return fs_; }
    int64_t n_begin() const { return segs_.front().n0; }
    int64_t n_end() const { return segs_.back().n1; }
    const std::vector<RangeSegment>& segments() const { return segs_; }

    // Point evaluation at global sample n; throws std::out_of_range outside [n_begin, n_end).
    double rho(int64_t n) const;      // m
    double rho_dot(int64_t n) const;  // m/s
    double tau(int64_t n) const;      // s, rho / c

    // out[i] = rho(n_first + i), bit-identical to rho(); throws std::out_of_range if any
    // index falls outside the track.
    void rho_block(int64_t n_first, std::span<double> out) const;

    // Bounds on rho over the continuous local-time span covering samples [n_a, n_b)
    // (endpoints n_a and n_b - 1, plus interior critical points). They contain every sampled
    // value to within float64 rounding; callers needing a hard bound add a margin.
    // Throws std::out_of_range unless n_begin <= n_a < n_b <= n_end.
    RangeBounds rho_bounds(int64_t n_a, int64_t n_b) const;

 private:
    std::size_t index_of(int64_t n) const;
    double local_u(const RangeSegment& s, int64_t n) const {
        return static_cast<double>(n - s.n0) / fs_;
    }

    double fs_;
    std::vector<RangeSegment> segs_;
};

}  // namespace ntn
