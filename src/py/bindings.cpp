// ntn_emu._core: pybind11 bindings for the C++ engine (SPEC §5.2).
//
// The control plane passes segment coefficients only, never sampled rho/tau/phase vectors.
// Pipeline.process() takes a 1-D C-contiguous complex64 array without conversion (zero-copy;
// any other dtype or layout raises TypeError instead of being copied silently), returns a new
// complex64 array, and runs with the GIL released.
// C++ exceptions map to Python: std::invalid_argument -> ValueError, std::out_of_range -> IndexError.
//
// CPU check (SPEC §6): module init calls ntn::require_avx2() before anything else and raises
// ImportError naming the missing features. This file is compiled without /arch:AVX2 (target _core,
// no LTCG), so the check runs before any AVX2 instruction of the engine.

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include <cstdint>
#include <format>
#include <memory>
#include <mutex>
#include <span>
#include <stdexcept>
#include <vector>

#include "ntn/cpu_check.hpp"
#include "ntn/geometric_channel.hpp"
#include "ntn/pipeline.hpp"
#include "ntn/range_track.hpp"

namespace py = pybind11;
using ntn::cf32;
using ntn::GeoConfig;
using ntn::Kernel;
using ntn::RangeTrack;

namespace {

using I64Array = py::array_t<int64_t, py::array::c_style | py::array::forcecast>;
using F64Array = py::array_t<double, py::array::c_style | py::array::forcecast>;
using CF32Array = py::array_t<cf32, py::array::c_style>;

RangeTrack make_track(double fs, const I64Array& n0, const I64Array& n1, const F64Array& coeffs) {
    if (n0.ndim() != 1 || n1.ndim() != 1) throw py::value_error("RangeTrack: n0 and n1 must be 1-D");
    if (coeffs.ndim() != 2 || coeffs.shape(1) != 4) throw py::value_error("RangeTrack: coeffs must have shape [N, 4]");
    const py::ssize_t n = n0.shape(0);
    if (n1.shape(0) != n || coeffs.shape(0) != n) {
        throw py::value_error(std::format("RangeTrack: n0, n1, coeffs lengths differ ({}, {}, {})", n,
                                          n1.shape(0), coeffs.shape(0)));
    }
    std::vector<ntn::RangeSegment> segs(static_cast<std::size_t>(n));
    const auto a0 = n0.unchecked<1>();
    const auto a1 = n1.unchecked<1>();
    const auto c = coeffs.unchecked<2>();
    for (py::ssize_t k = 0; k < n; ++k) {
        segs[static_cast<std::size_t>(k)] = {a0(k), a1(k), {c(k, 0), c(k, 1), c(k, 2), c(k, 3)}};
    }
    return RangeTrack(fs, std::move(segs));  // validates contiguity and C1
}

// process() runs without the GIL, so calls on one pipeline from several threads are serialised.
struct PyPipeline {
    PyPipeline(const RangeTrack& t, const GeoConfig& g, Kernel k) : p(t, g, k) {}
    ntn::Pipeline p;
    std::mutex m;
};

}  // namespace

PYBIND11_MODULE(_core, m) {
    try {
        ntn::require_avx2();  // first: nothing compiled with /arch:AVX2 has run yet
    } catch (const std::runtime_error& e) {
        throw py::import_error(e.what());
    }

    m.doc() = "NTN channel emulator C++ engine (AVX2 kernels by default; scalar kernels as oracle)";

    py::class_<RangeTrack>(m, "RangeTrack",
                           "Contiguous, C1 cubic range segments: rho(u) = c0 + c1 u + c2 u^2 + c3 u^3, "
                           "u = (n - n0)/fs, n in [n0, n1).")
        .def(py::init(&make_track), py::arg("fs"), py::arg("n0"), py::arg("n1"), py::arg("coeffs"))
        .def("rho", &RangeTrack::rho, py::arg("n"), "Range at global sample n [m].")
        .def("rho_dot", &RangeTrack::rho_dot, py::arg("n"), "Range rate at global sample n [m/s].")
        .def("tau", &RangeTrack::tau, py::arg("n"), "Delay rho/c at global sample n [s].")
        .def_property_readonly("fs", &RangeTrack::fs)
        .def_property_readonly("n_begin", &RangeTrack::n_begin)
        .def_property_readonly("n_end", &RangeTrack::n_end);

    py::enum_<Kernel>(m, "Kernel", "Avx2 (default) or Scalar (test oracle and benchmark baseline)")
        .value("Scalar", Kernel::Scalar)
        .value("Avx2", Kernel::Avx2);

    py::class_<GeoConfig>(m, "GeoConfig")
        .def(py::init([](double fc_hz, double fs_hz, int64_t tau0_samples, bool delay, bool doppler,
                         double residual_ref_doppler_hz) {
                 GeoConfig g;
                 g.fc_hz = fc_hz;
                 g.fs_hz = fs_hz;
                 g.tau0_samples = tau0_samples;
                 g.delay = delay;
                 g.doppler = doppler;
                 g.residual_ref_doppler_hz = residual_ref_doppler_hz;
                 return g;
             }),
             py::kw_only(), py::arg("fc_hz"), py::arg("fs_hz"), py::arg("tau0_samples") = 0,
             py::arg("delay") = true, py::arg("doppler") = true, py::arg("residual_ref_doppler_hz") = 0.0)
        .def_readwrite("fc_hz", &GeoConfig::fc_hz)
        .def_readwrite("fs_hz", &GeoConfig::fs_hz)
        .def_readwrite("tau0_samples", &GeoConfig::tau0_samples)
        .def_readwrite("delay", &GeoConfig::delay)
        .def_readwrite("doppler", &GeoConfig::doppler)
        .def_readwrite("residual_ref_doppler_hz", &GeoConfig::residual_ref_doppler_hz)
        .def("__repr__", [](const GeoConfig& g) {
            return std::format("GeoConfig(fc_hz={}, fs_hz={}, tau0_samples={}, delay={}, doppler={}, "
                               "residual_ref_doppler_hz={})",
                               g.fc_hz, g.fs_hz, g.tau0_samples, g.delay ? "True" : "False",
                               g.doppler ? "True" : "False", g.residual_ref_doppler_hz);
        });

    py::class_<PyPipeline>(m, "Pipeline",
                           "Geometric channel pipeline. out[i] is channel time n_in + i - latency(), "
                           "n_in = sample_index() before the call.")
        .def(py::init([](const RangeTrack& t, const GeoConfig& g, Kernel k) {
                 return std::make_unique<PyPipeline>(t, g, k);
             }),
             py::arg("track"), py::arg("geo"), py::arg("kernel") = Kernel::Avx2)
        .def_property_readonly("kernel", [](PyPipeline& self) { return self.p.kernel(); })
        .def(
            "process",
            [](PyPipeline& self, const CF32Array& x) {
                if (x.ndim() != 1) throw py::value_error("Pipeline.process: input must be 1-D");
                const py::ssize_t n = x.shape(0);
                CF32Array y(n);
                const cf32* in = x.data();  // no copy: x is complex64 and C-contiguous
                cf32* out = y.mutable_data();
                {
                    py::gil_scoped_release nogil;
                    std::lock_guard lock(self.m);
                    self.p.process(std::span<const cf32>(in, static_cast<std::size_t>(n)),
                                   std::span<cf32>(out, static_cast<std::size_t>(n)));
                }
                return y;
            },
            py::arg("x").noconvert(),
            "Process a 1-D C-contiguous complex64 block (no conversion); returns complex64 of the same length.")
        .def("sample_index", [](PyPipeline& self) {
            std::lock_guard lock(self.m);
            return self.p.sample_index();
        })
        .def("latency", [](PyPipeline& self) {
            std::lock_guard lock(self.m);
            return self.p.latency();
        })
        .def(
            "reset",
            [](PyPipeline& self, int64_t n_start) {
                std::lock_guard lock(self.m);
                self.p.reset(n_start);
            },
            py::arg("n_start"));

    m.def("cpu_has_avx2", &ntn::cpu_has_avx2,
          "True if the CPU and OS support everything the AVX2 build needs (AVX2, FMA3, BMI1, BMI2, OS AVX state).");
}
