// M2 throughput benchmarks (Google Benchmark; Release only; single-threaded).
// Run through scripts/bench_m2.py, which adds repetitions, random interleaving, machine and power
// state, and the summary (median, spread, real-time factor). Throughput is reported, never asserted.
//
// Benchmarks (items = complex samples; items_per_second / 1e6 = Msps):
//   BM_FracDelay/kernel:k/block:B      fractional-delay kernel alone (delays from a real track)
//   BM_Rotator/kernel:k/block:B        carrier rotator alone (rho from the same track)
//   BM_GeoChannel/kernel:k/block:B     GeometricChannel::process (range eval, delay, Doppler, history)
//   BM_Pipeline/kernel:k/block:B/fs_kHz:F   Pipeline::process: the full chain, which in M2 is the
//                                      geometric stage only (TDL arrives in M3, AWGN in M4)
// Timing: wall clock (UseRealTime). Google Benchmark's default rate basis is process CPU time,
// which on Windows advances in 15.625 ms ticks (~1.1% of a 1.4 s repetition): the first M2 run
// showed medians and spreads snapped to that grid (out/bench_m2_cputime_quantized.*).
// kernel: 0 = Scalar, 1 = Avx2. Every benchmark processes the same block repeatedly, so the data
// stay in cache as far as the working set allows; counter ws_KiB is an estimate of the bytes the
// block touches (buffers + filter table), reported for comparison with the cache sizes.
#include <benchmark/benchmark.h>

#include <array>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <random>
#include <span>
#include <stdexcept>
#include <string>
#include <vector>

#include "ntn/cpu_check.hpp"
#include "ntn/frac_delay.hpp"
#include "ntn/geometric_channel.hpp"
#include "ntn/pipeline.hpp"
#include "ntn/range_track.hpp"
#include "ntn/rotator.hpp"

namespace {

constexpr double kC = 299'792'458.0;
constexpr double kFc = 2'000'000'000.0;
constexpr double kFsBase = 7'680'000.0;
constexpr double kTrackSeconds = 2.0;
// Mid-pass-like range (rho ~ 1024 km, rho_dot ~ -5.9 km/s), as in tests/cpp/test_pipeline.cpp.
constexpr std::array<double, 4> kCoeffs{1'024'452.123, -5899.5, 13.4, -0.02};
constexpr int64_t kPad = 64;  // input margin around the frac-delay reads

ntn::Kernel kernel_of(int64_t arg) { return arg != 0 ? ntn::Kernel::Avx2 : ntn::Kernel::Scalar; }
const char* kernel_name(ntn::Kernel k) { return k == ntn::Kernel::Avx2 ? "Avx2" : "Scalar"; }

ntn::RangeTrack make_track(double fs) {
    const auto n = static_cast<int64_t>(kTrackSeconds * fs);
    return ntn::RangeTrack(fs, {{0, n, kCoeffs}});
}

ntn::GeoConfig make_geo(const ntn::RangeTrack& t) {
    ntn::GeoConfig g;
    g.fc_hz = kFc;
    g.fs_hz = t.fs();
    g.tau0_samples = static_cast<int64_t>(std::floor(t.rho(t.n_begin()) * t.fs() / kC));
    return g;
}

std::vector<ntn::cf32> random_block(std::size_t n, uint64_t seed) {
    std::mt19937_64 rng(seed);
    std::normal_distribution<float> g(0.0f, 1.0f);
    std::vector<ntn::cf32> x(n);
    for (auto& v : x) v = {g(rng), g(rng)};
    return x;
}

// Filter-table bytes each kernel reads (scalar: rows; AVX2: duplicated rows + row differences).
double table_bytes(ntn::Kernel k) {
    constexpr double taps = ntn::fd::kTaps, phases = ntn::fd::kPhases;
    return k == ntn::Kernel::Avx2 ? ((phases + 1) + phases) * 2 * taps * 4 : (phases + 1) * taps * 4;
}

void set_common(benchmark::State& st, ntn::Kernel k, std::size_t block, double ws_bytes) {
    st.SetItemsProcessed(static_cast<int64_t>(st.iterations()) * static_cast<int64_t>(block));
    st.SetLabel(kernel_name(k));
    st.counters["block"] = static_cast<double>(block);
    st.counters["ws_KiB"] = ws_bytes / 1024.0;
}

void BM_FracDelay(benchmark::State& st) {
    const ntn::Kernel k = kernel_of(st.range(0));
    const auto block = static_cast<std::size_t>(st.range(1));
    const ntn::RangeTrack track = make_track(kFsBase);
    const ntn::GeoConfig g = make_geo(track);
    std::vector<double> rho(block), d(block);
    track.rho_block(0, rho);
    for (std::size_t i = 0; i < block; ++i) d[i] = rho[i] * kFsBase / kC - static_cast<double>(g.tau0_samples);
    const auto hist = random_block(block + 2 * kPad, 1);
    std::vector<ntn::cf32> out(block);
    const auto& table = ntn::frac_delay_table();
    for (auto _ : st) {
        if (k == ntn::Kernel::Avx2) {
            ntn::frac_delay_avx2(table, hist.data(), -kPad, 0, d.data(), out.data(), block);
        } else {
            ntn::frac_delay_scalar(table, hist.data(), -kPad, 0, d.data(), out.data(), block);
        }
        benchmark::DoNotOptimize(out.data());
        benchmark::ClobberMemory();
    }
    set_common(st, k, block, (block + 2 * kPad) * 8.0 + block * 8.0 + block * 8.0 + table_bytes(k));
}

void BM_Rotator(benchmark::State& st) {
    const ntn::Kernel k = kernel_of(st.range(0));
    const auto block = static_cast<std::size_t>(st.range(1));
    const ntn::RangeTrack track = make_track(kFsBase);
    std::vector<double> rho(block);
    track.rho_block(0, rho);
    std::vector<ntn::cf32> io(block, ntn::cf32(1.0f, 0.0f));  // rotated in place, magnitude stays ~1
    for (auto _ : st) {
        if (k == ntn::Kernel::Avx2) {
            ntn::rotate_avx2(rho.data(), kFc, io.data(), block);
        } else {
            ntn::rotate_scalar(rho.data(), kFc, io.data(), block);
        }
        benchmark::DoNotOptimize(io.data());
        benchmark::ClobberMemory();
    }
    set_common(st, k, block, block * 8.0 + block * 8.0);
}

// Streams blocks through `stage` (GeometricChannel or Pipeline); restarts at the track start,
// outside the timed region, before the channel time would pass the track end.
template <class Stage>
void stream(benchmark::State& st, Stage& stage, const ntn::RangeTrack& track, std::size_t block) {
    const auto in = random_block(block, 2);
    std::vector<ntn::cf32> out(block);
    const int64_t lat = stage.latency();
    for (auto _ : st) {
        if (stage.sample_index() - lat + static_cast<int64_t>(block) > track.n_end()) {
            st.PauseTiming();
            stage.reset(track.n_begin());
            st.ResumeTiming();
        }
        stage.process(in, out);
        benchmark::DoNotOptimize(out.data());
        benchmark::ClobberMemory();
    }
}

// Working-set estimate of the streaming stages: input + output blocks, the input history (about
// one block plus the latency, up to two blocks before compaction), the per-4096-sample rho and
// delay scratch, and the filter table.
double stream_ws_bytes(ntn::Kernel k, std::size_t block, int64_t latency) {
    const double hist = (2.0 * static_cast<double>(block) + static_cast<double>(latency)) * 8.0;
    return block * 8.0 * 2.0 + hist + 2.0 * 4096 * 8.0 + table_bytes(k);
}

void BM_GeoChannel(benchmark::State& st) {
    const ntn::Kernel k = kernel_of(st.range(0));
    const auto block = static_cast<std::size_t>(st.range(1));
    const ntn::RangeTrack track = make_track(kFsBase);
    ntn::GeometricChannel ch(track, make_geo(track), k);
    stream(st, ch, track, block);
    set_common(st, k, block, stream_ws_bytes(k, block, ch.latency()));
    st.counters["fs_Msps"] = kFsBase / 1e6;
}

void BM_Pipeline(benchmark::State& st) {
    const ntn::Kernel k = kernel_of(st.range(0));
    const auto block = static_cast<std::size_t>(st.range(1));
    const double fs = static_cast<double>(st.range(2)) * 1e3;
    const ntn::RangeTrack track = make_track(fs);
    ntn::Pipeline p(track, make_geo(track), k);
    stream(st, p, track, block);
    set_common(st, k, block, stream_ws_bytes(k, block, p.latency()));
    st.counters["fs_Msps"] = fs / 1e6;
}

}  // namespace

BENCHMARK(BM_FracDelay)->ArgNames({"kernel", "block"})->Args({0, 4096})->Args({1, 4096})->UseRealTime();
BENCHMARK(BM_Rotator)->ArgNames({"kernel", "block"})->Args({0, 4096})->Args({1, 4096})->UseRealTime();
BENCHMARK(BM_GeoChannel)->ArgNames({"kernel", "block"})->Args({0, 4096})->Args({1, 4096})->UseRealTime();
BENCHMARK(BM_Pipeline)
    ->ArgNames({"kernel", "block", "fs_kHz"})
    ->ArgsProduct({{0, 1}, {4096, 65536}, {7680, 30720}})
    ->UseRealTime();

int main(int argc, char** argv) {
    try {
        ntn::require_avx2();
    } catch (const std::runtime_error& e) {
        std::fprintf(stderr, "%s\n", e.what());
        return 1;
    }
    const ntn::RangeTrack t = make_track(kFsBase);
    benchmark::AddCustomContext("cpu_has_avx2", ntn::cpu_has_avx2() ? "true" : "false");
    benchmark::AddCustomContext("default_kernel", kernel_name(ntn::Pipeline(t, make_geo(t)).kernel()));
    benchmark::AddCustomContext("ntn_build_type", NTN_BUILD_TYPE);
    benchmark::AddCustomContext("frac_delay_design", std::to_string(ntn::fd::kTaps) + " taps x " +
                                                         std::to_string(ntn::fd::kPhases) + " phases");
    benchmark::AddCustomContext("full_chain", "geometric stage only (TDL in M3, AWGN in M4)");
    benchmark::Initialize(&argc, argv);
    if (benchmark::ReportUnrecognizedArguments(argc, argv)) return 1;
    benchmark::RunSpecifiedBenchmarks();
    benchmark::Shutdown();
    return 0;
}
