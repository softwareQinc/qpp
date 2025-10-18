cmake_minimum_required(VERSION 3.20)

# -----------------------------------------------------------------------------
# Project version
# -----------------------------------------------------------------------------
project(
  qpp
  VERSION 6.0.0
  LANGUAGES CXX)

# -----------------------------------------------------------------------------
# C++ standard
# -----------------------------------------------------------------------------
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
enable_testing()

# -----------------------------------------------------------------------------
# Version Definitions
# -----------------------------------------------------------------------------
# The project() command automatically creates PROJECT_VERSION_MAJOR, _MINOR, and
# _PATCH. Define QPP_VERSION_NUM as an integer (e.g., 60000 for 6.0.0)
math(
  EXPR
  QPP_VERSION_NUM
  "${PROJECT_VERSION_MAJOR}*10000 + ${PROJECT_VERSION_MINOR}*100 + ${PROJECT_VERSION_PATCH}"
)

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
# Quantum++ Headers and Interface Library
# -----------------------------------------------------------------------------
# Use INTERFACE library for headers-only libraries
add_library(libqpp INTERFACE)

# Interface definitions using derived version variables
target_compile_definitions(
  libqpp
  INTERFACE "QPP_VERSION_NUM=${QPP_VERSION_NUM}"
            "QPP_VERSION_STR=\"${PROJECT_VERSION}\""
            $<BUILD_INTERFACE:PROJECT_ROOT_DIR=\"${PROJECT_SOURCE_DIR}\">)

# Quantum++ headers (include/)
target_include_directories(
  libqpp INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
                   $<INSTALL_INTERFACE:include/qpp>)

# qasmtools library headers (qasmtools/include/)
target_include_directories(
  libqpp
  INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/qasmtools/include>
            $<INSTALL_INTERFACE:include>)

# -----------------------------------------------------------------------------
# Optional Features
# -----------------------------------------------------------------------------
find_package(Python3 QUIET COMPONENTS Interpreter Development)
if(Python3_FOUND)
  include(cmake/pybind11.cmake)
  include(cmake/pyqpp.cmake)
endif()

# -----------------------------------------------------------------------------
# BEGIN LOCAL stuff, you don't need those in standalone projects
# -----------------------------------------------------------------------------
option(SANITIZE "Enable code sanitizing (only for GCC/Clang)" OFF)

# Dependencies
include(cmake/qpp_eigen3.cmake)
include(cmake/qpp_MATLAB.cmake)
include(cmake/qpp_dependencies.cmake)
target_link_libraries(libqpp INTERFACE ${QPP_LINK_DEPS})

# Unit testing
add_subdirectory(unit_tests EXCLUDE_FROM_ALL SYSTEM)

# Compiler Flags
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang" OR CMAKE_CXX_COMPILER_ID STREQUAL
                                            "GNU")
  target_compile_options(libqpp INTERFACE "-pedantic" "-Wall" "-Wextra"
                                          "-Weffc++")

  if(SANITIZE)
    if(NOT CMAKE_CXX_COMPILER_ID MATCHES "AppleClang")
      list(APPEND SANITIZE_OPTIONS -fsanitize=undefined)
      set_property(
        TARGET libqpp
        APPEND
        PROPERTY INTERFACE_COMPILE_OPTIONS ${SANITIZE_OPTIONS})
      # Apply to executable linker flags for the build targets
      set(CMAKE_EXE_LINKER_FLAGS
          "${CMAKE_EXE_LINKER_FLAGS} ${SANITIZE_OPTIONS}")
    endif()
  endif()
endif()

# Examples
include(cmake/examples.cmake)
# -----------------------------------------------------------------------------
# END LOCAL stuff
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Installation
# -----------------------------------------------------------------------------
include(GNUInstallDirs)
set(QPP_INSTALL_DIR
    "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}")

# Install headers
install(
  DIRECTORY include/
  DESTINATION ${QPP_INSTALL_DIR}
  FILES_MATCHING
  PATTERN "*.hpp"
  PATTERN "*.h")
install(
  DIRECTORY qasmtools/include/
  DESTINATION ${QPP_INSTALL_DIR}
  FILES_MATCHING
  PATTERN "*.hpp"
  PATTERN "*.h")

# Install target and config files
install(TARGETS libqpp EXPORT qpp_targets)
install(EXPORT qpp_targets
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME})

include(CMakePackageConfigHelpers)
# Generate qppConfig.cmake
configure_package_config_file(
  "cmake/qppConfig.cmake.in" "${CMAKE_CURRENT_BINARY_DIR}/qppConfig.cmake"
  INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME})
# Generate qppConfigVersion.cmake
write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/qppConfigVersion.cmake"
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY SameMajorVersion)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/qppConfig.cmake"
              "${CMAKE_CURRENT_BINARY_DIR}/qppConfigVersion.cmake"
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME})
install(FILES "cmake/qpp_dependencies.cmake" "cmake/qpp_eigen3.cmake"
              "cmake/qpp_MATLAB.cmake"
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME})

# -----------------------------------------------------------------------------
# Uninstall Target
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
  COMMENT "Uninstalling Quantum++")
