# Third-party dependencies via FetchContent: pinned tags, shallow clones, SYSTEM includes,
# dependency tests/examples off (SPEC §6). pybind11 comes from the venv, not from here.
# Tags verified with `git ls-remote` (M2): xsimd 13.2.0 -> 53341e94 (annotated; commit 1f8dd9c8),
# benchmark v1.9.4 -> 31c66f47.

include(FetchContent)

if(NTN_BUILD_TESTS)
  set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)  # /MD, matching our targets
  set(BUILD_GMOCK OFF CACHE BOOL "" FORCE)
  set(INSTALL_GTEST OFF CACHE BOOL "" FORCE)
  FetchContent_Declare(
    googletest
    GIT_REPOSITORY https://github.com/google/googletest.git
    GIT_TAG v1.17.0
    GIT_SHALLOW TRUE
    SYSTEM
  )
  FetchContent_MakeAvailable(googletest)
endif()

# xsimd (header-only): the AVX2 kernels and the .pyd need it, so it is always fetched.
FetchContent_Declare(
  xsimd
  GIT_REPOSITORY https://github.com/xtensor-stack/xsimd.git
  GIT_TAG 13.2.0
  GIT_SHALLOW TRUE
  SYSTEM
)
FetchContent_MakeAvailable(xsimd)

if(NTN_BUILD_BENCH)
  set(BENCHMARK_ENABLE_TESTING OFF CACHE BOOL "" FORCE)
  set(BENCHMARK_ENABLE_GTEST_TESTS OFF CACHE BOOL "" FORCE)
  set(BENCHMARK_ENABLE_INSTALL OFF CACHE BOOL "" FORCE)
  set(BENCHMARK_INSTALL_DOCS OFF CACHE BOOL "" FORCE)
  set(BENCHMARK_ENABLE_WERROR OFF CACHE BOOL "" FORCE)
  FetchContent_Declare(
    benchmark
    GIT_REPOSITORY https://github.com/google/benchmark.git
    GIT_TAG v1.9.4
    GIT_SHALLOW TRUE
    SYSTEM
  )
  FetchContent_MakeAvailable(benchmark)
endif()
