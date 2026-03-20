# pyqpp - Python 3 bindings for Quantum++

cmake_minimum_required(VERSION 3.20)
project(pyqpp LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# -----------------------------------------------------------------------------
# Detect Quantum++
# -----------------------------------------------------------------------------
# In-tree building from the Quantum++ source tree has the highest priority.
if(TARGET libqpp_internal)
  set(QPP_TARGET libqpp_internal)
  message(STATUS "pyqpp: using in-tree Quantum++ target")
else()
  find_package(qpp CONFIG REQUIRED)
  set(QPP_TARGET qpp::qpp)
  message(STATUS "pyqpp: using installed Quantum++ package")
endif()

# -----------------------------------------------------------------------------
# Detect Python
# -----------------------------------------------------------------------------
# If a virtual environment is active, prefer its interpreter.
if(DEFINED ENV{VIRTUAL_ENV})
  set(Python3_EXECUTABLE "$ENV{VIRTUAL_ENV}/bin/python")
endif()

message(STATUS "Detecting Python 3...")
find_package(Python3 REQUIRED COMPONENTS Interpreter Development)

# -----------------------------------------------------------------------------
# Detect pybind11
# -----------------------------------------------------------------------------
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

# -----------------------------------------------------------------------------
# Build the extension module
# -----------------------------------------------------------------------------
pybind11_add_module(_pyqpp src/pyqpp.cpp)
target_include_directories(_pyqpp PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(_pyqpp PRIVATE pybind11::module ${QPP_TARGET})

# -----------------------------------------------------------------------------
# Build-tree package staging
# -----------------------------------------------------------------------------
set(PYQPP_PYTHON_PACKAGE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/python/pyqpp")
set(PYQPP_PACKAGE_DIR "${CMAKE_CURRENT_BINARY_DIR}/pyqpp")
file(MAKE_DIRECTORY "${PYQPP_PACKAGE_DIR}")
file(MAKE_DIRECTORY "${PYQPP_PYTHON_PACKAGE_DIR}")

if(NOT EXISTS "${PYQPP_PYTHON_PACKAGE_DIR}/py.typed")
  file(WRITE "${PYQPP_PYTHON_PACKAGE_DIR}/py.typed" "")
endif()

configure_file(${PYQPP_PYTHON_PACKAGE_DIR}/__init__.py
               ${PYQPP_PACKAGE_DIR}/__init__.py COPYONLY)

configure_file(${PYQPP_PYTHON_PACKAGE_DIR}/gates.py
               ${PYQPP_PACKAGE_DIR}/gates.py COPYONLY)

configure_file(${PYQPP_PYTHON_PACKAGE_DIR}/states.py
               ${PYQPP_PACKAGE_DIR}/states.py COPYONLY)

configure_file(${PYQPP_PYTHON_PACKAGE_DIR}/qasm.py ${PYQPP_PACKAGE_DIR}/qasm.py
               COPYONLY)

if(EXISTS "${PYQPP_PYTHON_PACKAGE_DIR}/__init__.pyi")
  configure_file(${PYQPP_PYTHON_PACKAGE_DIR}/__init__.pyi
                 ${PYQPP_PACKAGE_DIR}/__init__.pyi COPYONLY)
endif()

if(EXISTS "${PYQPP_PYTHON_PACKAGE_DIR}/py.typed")
  configure_file(${PYQPP_PYTHON_PACKAGE_DIR}/py.typed
                 ${PYQPP_PACKAGE_DIR}/py.typed COPYONLY)
endif()

# -----------------------------------------------------------------------------
# Stub generation
# -----------------------------------------------------------------------------
option(PYQPP_ENABLE_STUBGEN "Enable stub generation for pyqpp" OFF)

if(PYQPP_ENABLE_STUBGEN)
  set(STUBGEN_OUTPUT_PATH "${CMAKE_CURRENT_BINARY_DIR}/stubs")
  file(MAKE_DIRECTORY "${STUBGEN_OUTPUT_PATH}")

  if(WIN32)
    set(PYTHONPATH_SEP ";")
  else()
    set(PYTHONPATH_SEP ":")
  endif()

  add_custom_command(
    OUTPUT "${PYQPP_PYTHON_PACKAGE_DIR}/_pyqpp.pyi"
           "${PYQPP_PYTHON_PACKAGE_DIR}/__init__.pyi"
    COMMAND ${CMAKE_COMMAND} -E make_directory "${PYQPP_PACKAGE_DIR}"
    COMMAND ${CMAKE_COMMAND} -E make_directory "${STUBGEN_OUTPUT_PATH}"
    COMMAND ${CMAKE_COMMAND} -E copy_if_different $<TARGET_FILE:_pyqpp>
            "${PYQPP_PACKAGE_DIR}/$<TARGET_FILE_NAME:_pyqpp>"
    COMMAND
      ${CMAKE_COMMAND} -E env
      "PYTHONPATH=${CMAKE_CURRENT_BINARY_DIR}${PYTHONPATH_SEP}$ENV{PYTHONPATH}"
      ${Python3_EXECUTABLE} -m pybind11_stubgen pyqpp --output-dir
      "${STUBGEN_OUTPUT_PATH}" --ignore-all-errors
      --print-invalid-expressions-as-is
    COMMAND
      ${CMAKE_COMMAND} -E copy_if_different
      "${STUBGEN_OUTPUT_PATH}/pyqpp/_pyqpp/__init__.pyi"
      "${PYQPP_PYTHON_PACKAGE_DIR}/_pyqpp.pyi"
    COMMAND
      ${CMAKE_COMMAND} -E copy_if_different
      "${STUBGEN_OUTPUT_PATH}/pyqpp/__init__.pyi"
      "${PYQPP_PYTHON_PACKAGE_DIR}/__init__.pyi"
    DEPENDS _pyqpp
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Generating Python stubs for pyqpp..."
    VERBATIM)

  # cmake-lint: disable=C0113
  add_custom_target(
    pqpp_stubgen ALL DEPENDS "${PYQPP_PYTHON_PACKAGE_DIR}/_pyqpp.pyi"
                             "${PYQPP_PYTHON_PACKAGE_DIR}/__init__.pyi")
endif()

# -----------------------------------------------------------------------------
# Install
# -----------------------------------------------------------------------------
# Install the compiled extension.
install(
  TARGETS _pyqpp
  COMPONENT python
  LIBRARY DESTINATION pyqpp)

# Install the pure-Python package files.
install(
  FILES ${PYQPP_PYTHON_PACKAGE_DIR}/__init__.py
        ${PYQPP_PYTHON_PACKAGE_DIR}/gates.py
        ${PYQPP_PYTHON_PACKAGE_DIR}/states.py
        ${PYQPP_PYTHON_PACKAGE_DIR}/qasm.py
  DESTINATION pyqpp
  COMPONENT python)

if(EXISTS "${PYQPP_PYTHON_PACKAGE_DIR}/__init__.pyi")
  install(
    FILES "${PYQPP_PYTHON_PACKAGE_DIR}/__init__.pyi"
    DESTINATION pyqpp
    COMPONENT python)
endif()

if(EXISTS "${PYQPP_PYTHON_PACKAGE_DIR}/_pyqpp.pyi")
  install(
    FILES "${PYQPP_PYTHON_PACKAGE_DIR}/_pyqpp.pyi"
    DESTINATION pyqpp
    COMPONENT python)
endif()

if(EXISTS "${PYQPP_PYTHON_PACKAGE_DIR}/py.typed")
  install(
    FILES "${PYQPP_PYTHON_PACKAGE_DIR}/py.typed"
    DESTINATION pyqpp
    COMPONENT python)
endif()
