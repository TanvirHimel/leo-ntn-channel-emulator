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

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
include(Warnings)
include(Deps)

add_executable(hello src/hello.cpp)
ntn_target_warnings(hello)

if(NTN_BUILD_TESTS)
  enable_testing()
  include(GoogleTest)
  add_executable(ntn_tests tests/cpp/test_smoke.cpp)
  ntn_target_warnings(ntn_tests)
  target_link_libraries(ntn_tests PRIVATE GTest::gtest_main)
  gtest_discover_tests(ntn_tests)
endif()
