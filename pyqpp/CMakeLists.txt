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
add_custom_command(
  TARGET pyqpp
  POST_BUILD
  COMMAND
    ${CMAKE_COMMAND} -E env "PYTHONPATH=$<TARGET_FILE_DIR:pyqpp>"
    ${Python3_EXECUTABLE} -m pybind11_stubgen pyqpp --output-dir
    ${CMAKE_CURRENT_BINARY_DIR} --ignore-all-errors
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  COMMENT "Generating Python stubs for pyqpp...")

# Create an empty "py.typed" file in the build directory
file(TOUCH "${CMAKE_CURRENT_BINARY_DIR}/py.typed")

# Install the generated stub directory alongside the module
install(DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/pyqpp" DESTINATION .)
