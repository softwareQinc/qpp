message(STATUS "Detecting pybind11...")
find_package(pybind11 CONFIG)
if(NOT pybind11_FOUND)
  # Install pybind11 on demand
  include(FetchContent)
  set(FETCHCONTENT_QUIET FALSE)
  set(CMAKE_POLICY_DEFAULT_CMP0077 NEW)
  message(STATUS "Fetching pybind11...")
  FetchContent_Declare(
    pybind11
    GIT_REPOSITORY https://github.com/pybind/pybind11
    GIT_TAG master
    GIT_SHALLOW TRUE
    GIT_PROGRESS TRUE)
  FetchContent_MakeAvailable(pybind11)
  set(PYBIND11_INCLUDE_DIRS ${pybind11_SOURCE_DIR})
endif()
message(STATUS "Detected pybind11 in: ${PYBIND11_INCLUDE_DIRS}")
