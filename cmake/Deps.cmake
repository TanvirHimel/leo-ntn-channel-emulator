# Third-party dependencies via FetchContent: pinned tags, shallow clones, SYSTEM includes,
# dependency tests/examples off (SPEC §6). pybind11 comes from the venv, not from here.

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
