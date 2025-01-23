cmake_minimum_required(VERSION 3.15)
set(QPP_VERSION_NUM 6.0)
set(QPP_VERSION_STR "${QPP_VERSION_NUM}")
project(
  qpp
  VERSION ${QPP_VERSION_NUM}
  LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
enable_testing()

# Quantum++ version
add_compile_definitions(QPP_VERSION_NUM=${QPP_VERSION_NUM})
add_compile_definitions(QPP_VERSION_STR="${QPP_VERSION_STR}")

# Quantum++ root directory
add_compile_definitions(PROJECT_ROOT_DIR="${PROJECT_SOURCE_DIR}")

# Guard against in-source builds (snippet from Eigen's CMakeLists.txt)
if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
  message(
    FATAL_ERROR
      "In-source builds not allowed. Please instruct CMake to use an\
      out-of-source build, e.g.,
      cmake -B build && cmake --build build
You may need to remove CMakeCache.txt.")
endif()

# Quantum++ headers
add_library(libqpp INTERFACE)
target_compile_definitions(libqpp
                           INTERFACE -DQPP_VERSION_NUM=${QPP_VERSION_NUM})
target_compile_definitions(libqpp
                           INTERFACE -DQPP_VERSION_STR="${QPP_VERSION_STR}")
target_include_directories(
  libqpp INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/>
                   $<INSTALL_INTERFACE:include/qpp/>)

# qasmtools library
target_include_directories(
  libqpp
  INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/qasmtools/include/>
            $<INSTALL_INTERFACE:include/>)

# pyqpp and pybind11 (only if the Python development kit is detected)
find_package(Python3 QUIET COMPONENTS Interpreter Development)
if(${Python3_FOUND})
  include(cmake/pybind11.cmake)
  include(cmake/pyqpp.cmake)
endif()

# BEGIN LOCAL stuff, you don't need those in standalone projects
option(SANITIZE "Enable code sanitizing (only for GCC/Clang)" OFF)

# Dependencies, do not modify unless you know what you're doing
include(cmake/qpp_eigen3.cmake)

# Dependencies, do not modify unless you know what you're doing
include(cmake/qpp_MATLAB.cmake)

# Dependencies, do not modify unless you know what you're doing
include(cmake/qpp_dependencies.cmake)

# Unit testing
add_subdirectory(${CMAKE_SOURCE_DIR}/unit_tests/ EXCLUDE_FROM_ALL SYSTEM)

# Enable all warnings for GCC and Clang/AppleClang
if(${CMAKE_CXX_COMPILER_ID} MATCHES "Clang" OR ${CMAKE_CXX_COMPILER_ID}
                                               STREQUAL "GNU")
  add_compile_options("-pedantic" "-Wall" "-Wextra" "-Weffc++")
  if(${SANITIZE})
    if(NOT (${CMAKE_CXX_COMPILER_ID} MATCHES "AppleClang"))
      list(APPEND SANITIZE_OPTIONS -fsanitize=undefined)
      add_compile_options("${SANITIZE_OPTIONS}")
      set(CMAKE_EXE_LINKER_FLAGS
          "${CMAKE_EXE_LINKER_FLAGS} ${SANITIZE_OPTIONS}")
    endif()
  endif()
endif()

# Examples
include(cmake/examples.cmake)
# END LOCAL stuff

# Installation
set(QPP_INSTALL_DIR "${CMAKE_INSTALL_PREFIX}/include/${PROJECT_NAME}")
install(DIRECTORY include/ DESTINATION ${QPP_INSTALL_DIR})
install(DIRECTORY qasmtools/include/ DESTINATION ${QPP_INSTALL_DIR})
install(TARGETS libqpp EXPORT qpp_targets)
install(EXPORT qpp_targets DESTINATION "lib/cmake/${PROJECT_NAME}")
include(CMakePackageConfigHelpers)
configure_package_config_file(
  "cmake/qppConfig.cmake.in" "${CMAKE_CURRENT_BINARY_DIR}/qppConfig.cmake"
  INSTALL_DESTINATION "lib/cmake/${PROJECT_NAME}")
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/qppConfig.cmake"
        DESTINATION "lib/cmake/${PROJECT_NAME}")
install(FILES "${CMAKE_SOURCE_DIR}/cmake/qpp_dependencies.cmake"
        DESTINATION "lib/cmake/${PROJECT_NAME}")
install(FILES "${CMAKE_SOURCE_DIR}/cmake/qpp_eigen3.cmake"
        DESTINATION "lib/cmake/${PROJECT_NAME}")
install(FILES "${CMAKE_SOURCE_DIR}/cmake/qpp_MATLAB.cmake"
        DESTINATION "lib/cmake/${PROJECT_NAME}")

# Uninstall
# https://gitlab.kitware.com/cmake/community/-/wikis/FAQ#can-i-do-make-uninstall-with-cmake
#
# UNIX/Linux: sudo cmake --build build --target uninstall

# Windows: cmake --build build --target uninstall
if(NOT TARGET uninstall)
  configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake/qpp_uninstall.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake" IMMEDIATE @ONLY)
  if(NOT MSVC)
    add_custom_target(
      uninstall
      COMMAND ${CMAKE_COMMAND} -P
              "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
      COMMAND ${CMAKE_COMMAND} -E remove_directory
              "${CMAKE_INSTALL_PREFIX}/lib/cmake/${PROJECT_NAME}"
      COMMAND ${CMAKE_COMMAND} -E remove_directory "${QPP_INSTALL_DIR}"
      COMMENT "Uninstall Quantum++")
  else()
    add_custom_target(
      uninstall
      COMMAND ${CMAKE_COMMAND} -P
              "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
      COMMAND ${CMAKE_COMMAND} -E remove_directory "${CMAKE_INSTALL_PREFIX}"
      COMMENT "Uninstall Quantum++")
  endif()
endif()
