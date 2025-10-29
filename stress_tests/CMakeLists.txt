cmake_minimum_required(VERSION 3.20)

project(stress_tests LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
  message(
    FATAL_ERROR
      "In-source builds not allowed. Please instruct CMake to use an\
      out-of-source build, e.g.,\n\
      cmake -B build && cmake --build build\n\
      You may need to remove CMakeCache.txt.")
endif()

# Use INTERFACE library for headers-only libraries
add_library(libqpp INTERFACE)

# Quantum++ headers
target_include_directories(
  libqpp INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../include/>)

# qasmtools library headers
target_include_directories(
  libqpp
  INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../qasmtools/include/>
)

# Dependencies
include(../cmake/qpp_types.cmake)
include(../cmake/qpp_compiler_flags.cmake)
include(../cmake/qpp_eigen3.cmake)
include(../cmake/qpp_openmp.cmake)
include(../cmake/qpp_openqasm.cmake)
include(../cmake/qpp_sanitize.cmake)
include(../cmake/qpp_MATLAB.cmake)

file(GLOB STRESS_TEST_SOURCES "${CMAKE_SOURCE_DIR}/src/*.cpp")

foreach(file ${STRESS_TEST_SOURCES})
  get_filename_component(TARGET_NAME ${file} NAME_WE)
  add_executable(${TARGET_NAME} ${file})
  target_link_libraries(${TARGET_NAME} PRIVATE libqpp)
endforeach()
