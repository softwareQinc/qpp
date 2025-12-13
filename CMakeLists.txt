cmake_minimum_required(VERSION 3.20)

# -----------------------------------------------------------------------------
# Project version
# -----------------------------------------------------------------------------
project(
  qpp
  VERSION 7.0.0
  LANGUAGES CXX)

# -----------------------------------------------------------------------------
# C++ standard
# -----------------------------------------------------------------------------
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# --------------------------------------------------
# Global install & packaging helpers
# --------------------------------------------------
include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

# -----------------------------------------------------------------------------
# Version Definitions
# -----------------------------------------------------------------------------
# Default patch to 0 if not provided
if("${PROJECT_VERSION_PATCH}" STREQUAL "")
  set(PROJECT_VERSION_PATCH 0)
endif()

# The project() command automatically creates PROJECT_VERSION_MAJOR, _MINOR, and
# _PATCH. Define QPP_VERSION_NUM as an integer (e.g., 70000 for 7.0.0)
math(EXPR QPP_VERSION_NUM
     "${PROJECT_VERSION_MAJOR}*10000 + ${PROJECT_VERSION_MINOR}*100 +
  ${PROJECT_VERSION_PATCH}")

# Define QPP_VERSION_STR as the unquoted string representation of
# ${PROJECT_VERSION}, e.g. 7.0.0
set(QPP_VERSION_STR ${PROJECT_VERSION})

# -----------------------------------------------------------------------------
# Guard against in-source builds
# -----------------------------------------------------------------------------
if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
  message(
    FATAL_ERROR
      "In-source builds not allowed. Please instruct CMake to use an\
      out-of-source build, e.g.,\n\
      cmake -B build && cmake --build build\n\
      You may need to remove CMakeCache.txt.")
endif()

# -----------------------------------------------------------------------------
# Quantum++ interface library and headers
# -----------------------------------------------------------------------------
# Use INTERFACE library for headers-only libraries
add_library(libqpp INTERFACE)

# Interface definitions using derived version variables
target_compile_definitions(
  libqpp
  INTERFACE "QPP_VERSION_NUM=${QPP_VERSION_NUM}"
            "QPP_VERSION_STR=\"${QPP_VERSION_STR}\""
            $<BUILD_INTERFACE:PROJECT_ROOT_DIR=\"${PROJECT_SOURCE_DIR}\">)

# Quantum++ & qasmtools library headers
target_include_directories(
  libqpp
  INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
            $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/qasmtools/include>
            $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/qpp>)

# -----------------------------------------------------------------------------
# Dependencies
# -----------------------------------------------------------------------------
include(cmake/qpp_types.cmake)
include(cmake/qpp_compiler_flags.cmake)
include(cmake/qpp_eigen3.cmake)

# Internal interface only, do not link against it, use libqpp instead. Used to
# build examples and unit tests.
add_library(libqpp_internal INTERFACE)
target_link_libraries(libqpp_internal INTERFACE libqpp)

# The dependencies below will link to libqpp_internal if it exists (defined
# above), otherwise they link to the public libqpp target. This ensures internal
# targets are used only when building within the source tree and avoids
# propagating them globally. See cmake/qpp_select_target.cmake for
# implementation details.
include(cmake/qpp_openmp.cmake)
include(cmake/qpp_openqasm.cmake)
include(cmake/qpp_sanitize.cmake)
include(cmake/qpp_MATLAB.cmake)

# -----------------------------------------------------------------------------
# Optional pyqpp development
# -----------------------------------------------------------------------------
add_subdirectory(pyqpp EXCLUDE_FROM_ALL)

# -----------------------------------------------------------------------------
# Examples
# -----------------------------------------------------------------------------
add_custom_target(examples COMMENT "Builds and groups all Quantum++ examples")
add_subdirectory(examples _examples EXCLUDE_FROM_ALL)

# -----------------------------------------------------------------------------
# Unit tests
# -----------------------------------------------------------------------------
enable_testing()
add_subdirectory(unit_tests EXCLUDE_FROM_ALL)

# -----------------------------------------------------------------------------
# Installation
# -----------------------------------------------------------------------------
# Install Quantum++ & qasmtools library headers
install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/qpp)
install(DIRECTORY qasmtools/include/
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/qpp)

# Install target and config files
install(TARGETS libqpp EXPORT qpp_targets)
install(EXPORT qpp_targets
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME})

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
        "cmake/qpp_eigen3.cmake"
        "cmake/qpp_openmp.cmake"
        "cmake/qpp_openqasm.cmake"
        "cmake/qpp_sanitize.cmake"
        "cmake/qpp_MATLAB.cmake"
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME})

# -----------------------------------------------------------------------------
# Uninstall target
# -----------------------------------------------------------------------------
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
