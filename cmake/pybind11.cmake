# pybind11

message(STATUS "Detecting pybind11...")

find_package(pybind11 CONFIG)

if(NOT pybind11_FOUND)
  include(FetchContent)
  message(STATUS "pybind11 not detected, fetching it...")
  FetchContent_Declare(
    pybind11
    GIT_REPOSITORY https://github.com/pybind/pybind11
    GIT_TAG master
    GIT_PROGRESS TRUE
    GIT_SHALLOW TRUE)
  FetchContent_MakeAvailable(pybind11)
  set(PYBIND11_INCLUDE_DIRS ${pybind11_SOURCE_DIR})
  message(STATUS "Installed pybind11: ${PYBIND11_INCLUDE_DIRS}")
endif()
