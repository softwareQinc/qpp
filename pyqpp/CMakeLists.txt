# pyqpp - Python3 bindings for Quantum++

cmake_minimum_required(VERSION 3.20)
project(pyqpp LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Define libqpp if not already defined
if(NOT TARGET libqpp)
  message(STATUS "Defining local libqpp interface target")
  add_library(libqpp INTERFACE)

  include("../cmake/qpp_openmp.cmake")
  include("../cmake/qpp_MATLAB.cmake")

  target_include_directories(
    libqpp INTERFACE ${CMAKE_SOURCE_DIR}/../include
                     ${CMAKE_SOURCE_DIR}/../qasmtools/include)

  target_compile_definitions(libqpp INTERFACE QPP_IDX_DEFAULT
                                              QPP_BIGINT_DEFAULT QPP_FP_DEFAULT)
endif()

# Detect Python and pybind11
find_package(Python3 REQUIRED COMPONENTS Interpreter Development)
include("../cmake/pybind11.cmake")

# Detect Eigen3
include("../cmake/qpp_eigen3.cmake")

# Add Python extension
pybind11_add_module(pyqpp src/pyqpp.cpp)

# Include directories
target_include_directories(
  pyqpp
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include ${CMAKE_SOURCE_DIR}/../include
          ${CMAKE_SOURCE_DIR}/../qasmtools/include)

# Link libraries
target_link_libraries(pyqpp PRIVATE pybind11::module libqpp)

# Install the module
install(TARGETS pyqpp LIBRARY DESTINATION pyqpp)

# Generate stubs by running pybind11-stubgen after the library is compiled
set(STUBGEN_OUTPUT_PATH "${CMAKE_CURRENT_BINARY_DIR}/stubs")

file(MAKE_DIRECTORY "${STUBGEN_OUTPUT_PATH}")

if(WIN32)
  set(PYTHONPATH_SEP ";")
else()
  set(PYTHONPATH_SEP ":")
endif()

add_custom_command(
  TARGET pyqpp
  POST_BUILD
  COMMAND
  COMMAND
    ${CMAKE_COMMAND} -E env
    "PYTHONPATH=$<TARGET_FILE_DIR:pyqpp>${PYTHONPATH_SEP}$ENV{PYTHONPATH}"
    ${Python3_EXECUTABLE} -m pybind11_stubgen pyqpp --output-dir
    "${STUBGEN_OUTPUT_PATH}" --ignore-all-errors
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  COMMENT "Generating Python stubs for pyqpp..."
  VERBATIM) # Install the generated stub directory alongside the module
install(DIRECTORY "${STUBGEN_OUTPUT_PATH}/pyqpp/" DESTINATION pyqpp)
