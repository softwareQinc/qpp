cmake_minimum_required(VERSION 3.20)

# ----------------------------------------------------------------------------
# Project version
# ----------------------------------------------------------------------------
project(
  qpp
  VERSION 7.0.2
  LANGUAGES CXX)

# ----------------------------------------------------------------------------
# C++ standard
# ----------------------------------------------------------------------------
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# ----------------------------------------------------------------------------
# Global install & packaging helpers
# ----------------------------------------------------------------------------
include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

# ---------------------------------------------------------------------------
# Version definitions
# ---------------------------------------------------------------------------
# Default patch to 0 if not provided
if("${PROJECT_VERSION_PATCH}" STREQUAL "")
  set(PROJECT_VERSION_PATCH 0)
endif()

# The project() command automatically creates PROJECT_VERSION_MAJOR, _MINOR, and
# _PATCH. Define QPP_VERSION_NUM as an integer (e.g., 70000 for 7.0.0, 70103 for
# 7.1.3)
math(EXPR QPP_VERSION_NUM
     "${PROJECT_VERSION_MAJOR}*10000 + ${PROJECT_VERSION_MINOR}*100 +
  ${PROJECT_VERSION_PATCH}")

# Define QPP_VERSION_STR as the unquoted string representation of
# ${PROJECT_VERSION}, e.g., 7.0.0
set(QPP_VERSION_STR ${PROJECT_VERSION})

# ----------------------------------------------------------------------------
# Guard against in-source builds
# ----------------------------------------------------------------------------
if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
  message(
    FATAL_ERROR
      "In-source builds not allowed. Please instruct CMake to use an\
      out-of-source build, e.g.,\n\
      cmake -B build && cmake --build build\n\
      You may need to remove CMakeCache.txt.")
endif()

# ----------------------------------------------------------------------------
# Quantum++ interface library and headers
# ----------------------------------------------------------------------------
# Use INTERFACE library for headers-only libraries
add_library(qpp INTERFACE)

# Interface definitions using derived version variables
target_compile_definitions(
  qpp
  INTERFACE "QPP_VERSION_NUM=${QPP_VERSION_NUM}"
            "QPP_VERSION_STR=\"${QPP_VERSION_STR}\""
            $<BUILD_INTERFACE:PROJECT_ROOT_DIR=\"${PROJECT_SOURCE_DIR}\">)

# Quantum++ & qasmtools library headers
target_include_directories(
  qpp
  INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
            $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/qasmtools/include>
            $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/qpp>)

# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------
include(cmake/qpp_eigen.cmake)

# Internal interface only, do not link against it, use qpp::qpp instead for
# downstream projects. Used to build examples and unit tests.
add_library(libqpp_internal INTERFACE)
target_link_libraries(libqpp_internal INTERFACE qpp)

# The dependencies below will link to libqpp_internal if it exists (defined
# above), otherwise they link to the public qpp target. This ensures internal
# targets are used only when building within the source tree and avoids
# propagating them globally. See cmake/qpp_select_target.cmake for
# implementation details.
include(cmake/qpp_compiler_flags.cmake)
include(cmake/qpp_openmp.cmake)
include(cmake/qpp_openqasm.cmake)
include(cmake/qpp_sanitize.cmake)
include(cmake/qpp_types.cmake)
include(cmake/qpp_MATLAB.cmake)

# ----------------------------------------------------------------------------
# Examples
# ----------------------------------------------------------------------------
add_custom_target(examples COMMENT "Builds and groups all Quantum++ examples")
add_subdirectory(examples _examples EXCLUDE_FROM_ALL)

# ----------------------------------------------------------------------------
# Unit tests, disabled by default
# ----------------------------------------------------------------------------
option(QPP_ENABLE_TESTING "Enable Quantum++ unit tests" OFF)
if(QPP_ENABLE_TESTING)
  enable_testing()
  add_subdirectory(unit_tests EXCLUDE_FROM_ALL)
endif()

# ----------------------------------------------------------------------------
# pyqpp C++ development, disabled by default
# ----------------------------------------------------------------------------
option(QPP_ENABLE_PYQPP "Enable pyqpp development" OFF)
if(QPP_ENABLE_PYQPP)
  add_subdirectory(pyqpp EXCLUDE_FROM_ALL)
endif()

# ----------------------------------------------------------------------------
# pyqpp PACKAGE BUILD: Bypass EXCLUDE_FROM_ALL
# ----------------------------------------------------------------------------
# If this flag is defined by an external package manager (e.g., scikit-build),
# the pyqpp target is explicitly included in the default 'ALL' build target,
# overriding any 'EXCLUDE_FROM_ALL' used for C++ development builds
if(__INTERNAL_PYQPP_SCIKIT_BUILD__)
  add_subdirectory(pyqpp)
endif()

# ----------------------------------------------------------------------------
# Stage auxiliary CMake files into the build tree
#
# These files are NOT generated by CMake; they are copied verbatim so that - the
# build tree can be used as a development prefix - qppConfig.cmake can include
# them via relative paths - consumers can use the build directory without
# installing first
#
# Note: Exported target files (e.g., qpp_targets.cmake) are intentionally
# generated separately via `export(EXPORT ...)` and must not be copied
# ----------------------------------------------------------------------------
message(STATUS "Staging auxiliary CMake files for development prefix...")
set(QPP_AUX_CMAKE_FILES
    "qpp_select_target.cmake"
    "qpp_types.cmake"
    "qpp_compiler_flags.cmake"
    "qpp_eigen.cmake"
    "qpp_openmp.cmake"
    "qpp_openqasm.cmake"
    "qpp_sanitize.cmake"
    "qpp_MATLAB.cmake")
foreach(file IN LISTS QPP_AUX_CMAKE_FILES)
  configure_file("${CMAKE_CURRENT_SOURCE_DIR}/cmake/${file}"
                 "${CMAKE_CURRENT_BINARY_DIR}/${file}" COPYONLY)
endforeach()

# ----------------------------------------------------------------------------
# Installation
# ----------------------------------------------------------------------------
# Install Quantum++ & qasmtools library headers
set(QPP_INSTALL_DIR
    ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME})
install(DIRECTORY include/
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME})
install(DIRECTORY qasmtools/include/
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME})

# Install target and config files
install(TARGETS qpp EXPORT qpp_targets)
install(
  EXPORT qpp_targets
  NAMESPACE qpp::
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME})

# Generate exported target definitions for the build-tree development prefix
export(
  EXPORT qpp_targets
  NAMESPACE qpp::
  FILE "${CMAKE_CURRENT_BINARY_DIR}/qpp_targets.cmake")

# Generate qppConfig.cmake
configure_package_config_file(
  "cmake/qppConfig.cmake.in" "${CMAKE_CURRENT_BINARY_DIR}/qppConfig.cmake"
  INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME})

# Generate qppConfigVersion.cmake
write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/qppConfigVersion.cmake"
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY SameMajorVersion)

# Install configuration files
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/qppConfig.cmake"
              "${CMAKE_CURRENT_BINARY_DIR}/qppConfigVersion.cmake"
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME})

# Install additional dependencies
install(
  FILES "cmake/qpp_select_target.cmake"
        "cmake/qpp_types.cmake"
        "cmake/qpp_compiler_flags.cmake"
        "cmake/qpp_eigen.cmake"
        "cmake/qpp_openmp.cmake"
        "cmake/qpp_openqasm.cmake"
        "cmake/qpp_sanitize.cmake"
        "cmake/qpp_MATLAB.cmake"
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME})

# ----------------------------------------------------------------------------
# Uninstall target
# ----------------------------------------------------------------------------
# https://gitlab.kitware.com/cmake/community/-/wikis/FAQ#can-i-do-make-uninstall-with-cmake

# UNIX/Linux: sudo cmake --build build --target uninstall

# Windows: cmake --build build --target uninstall

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/cmake/qpp_uninstall.cmake.in"
               "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake" @ONLY)
add_custom_target(
  uninstall
  COMMAND ${CMAKE_COMMAND} -P
          "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
  COMMENT "Uninstalling Quantum++"
  VERBATIM)
