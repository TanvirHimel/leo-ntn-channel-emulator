cmake_minimum_required(VERSION 3.25)
project(ntn_emulator LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(NOT MSVC)
  message(FATAL_ERROR "This project targets MSVC x64 only (SPEC §6).")
endif()

# /MD everywhere, matching python311.dll; /fp:precise everywhere, deps included (SPEC §6).
set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>DLL")
add_compile_options(/fp:precise)

option(NTN_BUILD_TESTS "Build the GoogleTest suite" ON)
option(NTN_BUILD_PYTHON "Build the ntn_emu._core pybind11 module" OFF)
option(NTN_BUILD_BENCH "Build the Google Benchmark suite (Release only)" OFF)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
include(Warnings)
include(Deps)

# Toolchain sanity target; not part of the Python wheel build (SKBUILD is set by scikit-build-core).
if(NOT SKBUILD)
  add_executable(hello src/hello.cpp)
  ntn_target_warnings(hello)
endif()

# Scalar kernels: their own target, never compiled with /arch:AVX2 (SPEC §6), so they stay
# the test oracle and an honest benchmark baseline.
add_library(ntn_scalar STATIC
  src/ntn/frac_delay.cpp
  src/ntn/rotator.cpp
  src/ntn/cpu_check.cpp
)
target_include_directories(ntn_scalar PUBLIC src)
ntn_target_warnings(ntn_scalar)

# AVX2 kernels (xsimd): the only kernel target compiled with /arch:AVX2. Selected at runtime via
# Kernel::Avx2; cpu_check (in ntn_scalar) guards it.
add_library(ntn_avx2 STATIC
  src/ntn/frac_delay_avx2.cpp
  src/ntn/rotator_avx2.cpp
)
target_include_directories(ntn_avx2 PUBLIC src)
target_link_libraries(ntn_avx2 PUBLIC xsimd)
target_compile_options(ntn_avx2 PRIVATE /arch:AVX2)
ntn_target_warnings(ntn_avx2)

# Core engine (range track, channel orchestration), /arch:AVX2. The bindings and ntn_scalar are
# built without it, so the import-time AVX2 check runs before any AVX2 instruction (SPEC §6 note).
add_library(ntn_core STATIC
  src/ntn/range_track.cpp
  src/ntn/geometric_channel.cpp
  src/ntn/pipeline.cpp
)
target_include_directories(ntn_core PUBLIC src)
target_link_libraries(ntn_core PUBLIC ntn_scalar ntn_avx2)
target_compile_options(ntn_core PRIVATE /arch:AVX2)
ntn_target_warnings(ntn_core)

# Python extension ntn_emu._core (pip install -e . turns this on). pybind11 comes from the venv;
# /MD matches python311.dll and every dependency is linked statically into the .pyd.
if(NTN_BUILD_PYTHON)
  find_package(Python 3.11 REQUIRED COMPONENTS Interpreter Development.Module)
  if(NOT pybind11_DIR)
    execute_process(
      COMMAND "${Python_EXECUTABLE}" -m pybind11 --cmakedir
      OUTPUT_VARIABLE pybind11_DIR OUTPUT_STRIP_TRAILING_WHITESPACE COMMAND_ERROR_IS_FATAL ANY)
  endif()
  set(PYBIND11_FINDPYTHON ON)
  # No LTCG (/GL): whole-program codegen would blur the per-object /arch:AVX2 boundary and leave
  # bindings.obj as IL that scripts/check_isa.py cannot inspect. Defining this stops
  # pybind11_add_module from adding pybind11::lto.
  set(CMAKE_INTERPROCEDURAL_OPTIMIZATION OFF)
  find_package(pybind11 CONFIG REQUIRED)
  pybind11_add_module(_core MODULE src/py/bindings.cpp)
  target_link_libraries(_core PRIVATE ntn_core)
  ntn_target_warnings(_core)
  install(TARGETS _core DESTINATION ntn_emu)
endif()

if(NTN_BUILD_TESTS)
  enable_testing()
  include(GoogleTest)
  add_executable(ntn_tests
    tests/cpp/test_smoke.cpp
    tests/cpp/test_range_track.cpp
    tests/cpp/test_frac_delay.cpp
    tests/cpp/test_rotator.cpp
    tests/cpp/test_geometric_channel.cpp
    tests/cpp/test_pipeline.cpp
    tests/cpp/test_cpu_check.cpp
  )
  ntn_target_warnings(ntn_tests)
  target_link_libraries(ntn_tests PRIVATE ntn_core GTest::gtest_main)
  # Fast suite: plain `ctest`. Slow variants (*Slow*) carry label `slow` and skip unless
  # NTN_SLOW=1 is set in the environment: `set NTN_SLOW=1 && ctest -L slow`.
  gtest_discover_tests(ntn_tests TEST_FILTER "-*Slow*")
  gtest_discover_tests(ntn_tests TEST_FILTER "*Slow*"
    PROPERTIES LABELS slow SKIP_REGULAR_EXPRESSION "\\[  SKIPPED \\]")

  # ISA probe: never linked. Compiled with the ntn_avx2 flags so scripts/check_isa.py has a
  # positive control (ymm present) and can inspect FMA contraction under /fp:precise.
  add_library(ntn_isa_probe OBJECT tests/cpp/isa_probe_avx2.cpp)
  target_link_libraries(ntn_isa_probe PRIVATE xsimd)
  target_compile_options(ntn_isa_probe PRIVATE /arch:AVX2)
  ntn_target_warnings(ntn_isa_probe)
endif()
