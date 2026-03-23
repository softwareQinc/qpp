# pyqpp - Python 3 bindings for Quantum++

cmake_minimum_required(VERSION 3.20)
project(pyqpp LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# -----------------------------------------------------------------------------
# Detect Quantum++
# -----------------------------------------------------------------------------
# In-tree building from the Quantum++ source tree has the highest priority
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
# If a virtual environment is active, prefer its interpreter
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

# Ensure directories and type marker exist
file(MAKE_DIRECTORY "${PYQPP_PACKAGE_DIR}")
file(TOUCH "${PYQPP_PYTHON_PACKAGE_DIR}/py.typed")

# Stage core Python files
set(PY_FILES
    __init__.py
    codes.py
    gates.py
    qasm.py
    random_devices.py
    states.py
    py.typed)
foreach(file ${PY_FILES})
  configure_file("${PYQPP_PYTHON_PACKAGE_DIR}/${file}"
                 "${PYQPP_PACKAGE_DIR}/${file}" COPYONLY)
endforeach()

# -----------------------------------------------------------------------------
# Stub generation
# -----------------------------------------------------------------------------
option(PYQPP_ENABLE_STUBGEN "Enable stub generation for pyqpp" OFF)

if(PYQPP_ENABLE_STUBGEN)
  set(STUBGEN_OUTPUT_PATH "${CMAKE_CURRENT_BINARY_DIR}/stubs")

  # Windows/UNIX path separator
  set(PYTHONPATH_SEP "$<IF:$<BOOL:${WIN32}>,;,:>")

  # List of sub-modules to copy stubs for
  set(PYI_FILES codes gates qasm random_devices states)

  # Prepare the individual copy commands
  set(COPY_STUB_COMMANDS)
  foreach(stub ${PYI_FILES})
    list(
      APPEND
      COPY_STUB_COMMANDS
      COMMAND
      ${CMAKE_COMMAND}
      -E
      copy_if_different
      "${STUBGEN_OUTPUT_PATH}/pyqpp/_pyqpp/${stub}.pyi"
      "${PYQPP_PYTHON_PACKAGE_DIR}/${stub}.pyi")
  endforeach()

  add_custom_command(
    OUTPUT "${PYQPP_PYTHON_PACKAGE_DIR}/_pyqpp.pyi"
           "${PYQPP_PYTHON_PACKAGE_DIR}/__init__.pyi"
    COMMAND ${CMAKE_COMMAND} -E make_directory "${PYQPP_PACKAGE_DIR}"
            "${STUBGEN_OUTPUT_PATH}"
    COMMAND ${CMAKE_COMMAND} -E copy_if_different $<TARGET_FILE:_pyqpp>
            "${PYQPP_PACKAGE_DIR}/"
    COMMAND
      ${CMAKE_COMMAND} -E env
      "PYTHONPATH=${CMAKE_CURRENT_BINARY_DIR}${PYTHONPATH_SEP}$ENV{PYTHONPATH}"
      ${Python3_EXECUTABLE} -m pybind11_stubgen pyqpp --output-dir
      "${STUBGEN_OUTPUT_PATH}" --ignore-all-errors
      --print-invalid-expressions-as-is
    # Copy main package stubs
    COMMAND
      ${CMAKE_COMMAND} -E copy_if_different
      "${STUBGEN_OUTPUT_PATH}/pyqpp/_pyqpp/__init__.pyi"
      "${PYQPP_PYTHON_PACKAGE_DIR}/_pyqpp.pyi"
      # Execute the generated list of copy commands for sub-modules
    COMMAND
      ${CMAKE_COMMAND} -E copy_if_different
      "${STUBGEN_OUTPUT_PATH}/pyqpp/__init__.pyi"
      "${PYQPP_PYTHON_PACKAGE_DIR}/__init__.pyi" ${COPY_STUB_COMMANDS}
    DEPENDS _pyqpp
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Generating Python stubs for pyqpp..."
    VERBATIM)

  # cmake-lint: disable=C0113
  add_custom_target(pqpp_stubgen ALL
                    DEPENDS "${PYQPP_PYTHON_PACKAGE_DIR}/_pyqpp.pyi")
endif()

# -----------------------------------------------------------------------------
# Install
# -----------------------------------------------------------------------------
# Install the compiled extension
install(
  TARGETS _pyqpp
  COMPONENT python
  DESTINATION pyqpp)

# Install all Python files, stubs, and type markers from the package directory
install(
  DIRECTORY "${PYQPP_PYTHON_PACKAGE_DIR}/"
  DESTINATION pyqpp
  COMPONENT python
  FILES_MATCHING
  PATTERN "*.py"
  PATTERN "*.pyi"
  PATTERN "py.typed")
