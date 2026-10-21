#include "ntn/range_track.hpp"

#include <algorithm>
#include <cmath>
#include <format>
#include <stdexcept>
#include <utility>

#include "ntn/constants.hpp"

namespace ntn {

RangeTrack::RangeTrack(double fs, std::vector<RangeSegment> segments)
    : fs_(fs), segs_(std::move(segments)) {
    if (!std::isfinite(fs_) || fs_ <= 0.0) {
        throw std::invalid_argument(std::format("RangeTrack: fs must be finite and > 0 (got {})", fs_));
    }
    if (segs_.empty()) throw std::invalid_argument("RangeTrack: no segments");

    for (std::size_t k = 0; k < segs_.size(); ++k) {
        const RangeSegment& s = segs_[k];
        if (s.n1 <= s.n0) {
            throw std::invalid_argument(
                std::format("RangeTrack: segment {} is empty or reversed [{}, {})", k, s.n0, s.n1));
        }
        for (double c : s.c) {
            if (!std::isfinite(c)) {
                throw std::invalid_argument(std::format("RangeTrack: segment {} has a non-finite coefficient", k));
            }
        }
        if (k > 0 && s.n0 != segs_[k - 1].n1) {
            throw std::invalid_argument(std::format(
                "RangeTrack: segments are not contiguous at {} (previous ends at {}, next starts at {})",
                k, segs_[k - 1].n1, s.n0));
        }
    }

    // C1 at every interior knot: end of segment k (u = h) against the start of segment k+1.
    for (std::size_t k = 0; k + 1 < segs_.size(); ++k) {
        const RangeSegment& a = segs_[k];
        const RangeSegment& b = segs_[k + 1];
        const double h = local_u(a, a.n1);
        const double d_rho = std::abs(horner_rho(a.c, h) - b.c[0]);
        const double d_rate = std::abs(horner_rho_dot(a.c, h) - b.c[1]);
        if (d_rho >= kC1TolRho || d_rate >= kC1TolRate) {
            throw std::invalid_argument(std::format(
                "RangeTrack: C1 violated at knot n = {}: |d rho| = {:.3e} m, |d rho_dot| = {:.3e} m/s "
                "(limits {:.0e} m, {:.0e} m/s)",
                a.n1, d_rho, d_rate, kC1TolRho, kC1TolRate));
        }
    }
}

std::size_t RangeTrack::index_of(int64_t n) const {
    if (n < n_begin() || n >= n_end()) {
        throw std::out_of_range(
            std::format("RangeTrack: sample {} outside the track [{}, {})", n, n_begin(), n_end()));
    }
    const auto it = std::upper_bound(segs_.begin(), segs_.end(), n,
                                     [](int64_t v, const RangeSegment& s) { return v < s.n0; });
    return static_cast<std::size_t>(it - segs_.begin()) - 1;
}

double RangeTrack::rho(int64_t n) const {
    const RangeSegment& s = segs_[index_of(n)];
    return horner_rho(s.c, local_u(s, n));
}

double RangeTrack::rho_dot(int64_t n) const {
    const RangeSegment& s = segs_[index_of(n)];
    return horner_rho_dot(s.c, local_u(s, n));
}

double RangeTrack::tau(int64_t n) const { return rho(n) / kC; }

void RangeTrack::rho_block(int64_t n_first, std::span<double> out) const {
    if (out.empty()) return;
    const auto count = static_cast<int64_t>(out.size());
    if (n_first < n_begin() || n_first > n_end() - count) {
        throw std::out_of_range(std::format("RangeTrack: block [{}, {}) outside the track [{}, {})",
                                            n_first, n_first + count, n_begin(), n_end()));
    }
    std::size_t k = index_of(n_first);
    int64_t n = n_first;
    std::size_t i = 0;
    while (i < out.size()) {
        const RangeSegment& s = segs_[k];
        const int64_t run_end = std::min(s.n1, n_first + count);
        for (; n < run_end; ++n, ++i) out[i] = horner_rho(s.c, local_u(s, n));
        ++k;
    }
}

RangeBounds RangeTrack::rho_bounds(int64_t n_a, int64_t n_b) const {
    if (n_a >= n_b || n_a < n_begin() || n_b > n_end()) {
        throw std::out_of_range(std::format("RangeTrack: bounds range [{}, {}) invalid for track [{}, {})",
                                            n_a, n_b, n_begin(), n_end()));
    }
    RangeBounds b{rho(n_a), rho(n_a)};
    auto take = [&b](double v) {
        b.min = std::min(b.min, v);
        b.max = std::max(b.max, v);
    };
    for (std::size_t k = index_of(n_a); k < segs_.size() && segs_[k].n0 < n_b; ++k) {
        const RangeSegment& s = segs_[k];
        const double ua = local_u(s, std::max(s.n0, n_a));
        const double ub = local_u(s, std::min(s.n1, n_b) - 1);
        take(horner_rho(s.c, ua));
        take(horner_rho(s.c, ub));

        // Critical points: roots of rho' = A u^2 + B u + C0 strictly inside (ua, ub).
        const double A = 3.0 * s.c[3], B = 2.0 * s.c[2], C0 = s.c[1];
        auto take_if_inside = [&](double u) {
            if (u > ua && u < ub) take(horner_rho(s.c, u));
        };
        if (A == 0.0) {
            if (B != 0.0) take_if_inside(-C0 / B);
        } else {
            const double disc = B * B - 4.0 * A * C0;
            if (disc >= 0.0) {
                // Numerically stable roots: q = -(B + sign(B) sqrt(disc)) / 2, u = q/A, C0/q.
                const double q = -0.5 * (B + std::copysign(std::sqrt(disc), B));
                take_if_inside(q / A);
                if (q != 0.0) take_if_inside(C0 / q);
            }
        }
    }
    return b;
}

}  // namespace ntn
