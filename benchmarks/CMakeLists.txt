# Benchmarks with Catch2
cmake_minimum_required(VERSION 3.20)

project(benchmarks)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

include(FetchContent)
message(STATUS "Fetching Catch2...")
FetchContent_Declare(
  Catch2
  GIT_REPOSITORY https://github.com/catchorg/Catch2.git
  GIT_TAG v3.11.0)
FetchContent_MakeAvailable(Catch2)

# Quantum++ interface library and headers
add_library(libqpp INTERFACE)
set(QPP_PROJECT_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/../)
target_include_directories(
  libqpp INTERFACE $<BUILD_INTERFACE:${QPP_PROJECT_ROOT}/include>)

# qasmtools library headers (qasmtools/include/)
target_include_directories(
  libqpp INTERFACE $<BUILD_INTERFACE:${QPP_PROJECT_ROOT}/qasmtools/include>)

# Quantum++ dependencies
include(../cmake/qpp_types.cmake)
include(../cmake/qpp_compiler_flags.cmake)
include(../cmake/qpp_eigen3.cmake)

# Internal interface only, do not link against it, use libqpp instead. Used to
# build benchmarks.
add_library(libqpp_internal INTERFACE)
target_link_libraries(libqpp_internal INTERFACE libqpp)

include(../cmake/qpp_openmp.cmake)
include(../cmake/qpp_openqasm.cmake)
include(../cmake/qpp_sanitize.cmake)
include(../cmake/qpp_MATLAB.cmake)

# Common benchmark interface
add_library(benchmarks_interface INTERFACE)
target_include_directories(
  benchmarks_interface
  INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)
target_link_libraries(benchmarks_interface)

# Benchmark source files
set(SRC_DIR "src")
set(SRC_FILES "")
file(GLOB CURRENT_FILES "${SRC_DIR}/*.cpp" "${SRC_DIR}/*.cxx" "${SRC_DIR}/*.cc"
     "${SRC_DIR}/*.c++")
list(APPEND SRC_FILES ${CURRENT_FILES})

if(NOT SRC_FILES)
  message(WARNING "No C++ source files found in the ${SRC_DIR} directory.")
endif()

set(BENCHMARK_TARGETS "")

# cmake-lint: disable=C0103
foreach(SOURCE_FILE ${SRC_FILES})
  get_filename_component(TARGET_NAME ${SOURCE_FILE} NAME_WE)
  add_executable(${TARGET_NAME} ${SOURCE_FILE})
  set_target_properties(
    ${TARGET_NAME}
    PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/benchmarks"
               RUNTIME_OUTPUT_DIRECTORY_DEBUG "${CMAKE_BINARY_DIR}/benchmarks"
               RUNTIME_OUTPUT_DIRECTORY_RELEASE
               "${CMAKE_BINARY_DIR}/benchmarks")
  target_link_libraries(
    ${TARGET_NAME} PRIVATE benchmarks_interface libqpp_internal
                           Catch2::Catch2WithMain)
  list(APPEND BENCHMARK_TARGETS ${TARGET_NAME})
endforeach()
