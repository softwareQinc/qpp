# Compiler flags/options

# Select the target
include(${CMAKE_CURRENT_LIST_DIR}/qpp_select_target.cmake)
qpp_select_target(QPP_TARGET "qpp_compiler_flags")

option(QPP_QUBIT_OPTIMIZATIONS "Enable qubit-specific optimizations" ON)
message(STATUS "Qubit optimizations - ${QPP_QUBIT_OPTIMIZATIONS}")

option(QPP_NO_THREAD_LOCAL "Disable thread_local storage duration" OFF)
if(QPP_NO_THREAD_LOCAL)
  message(STATUS "thread_local support - OFF")
endif()

target_compile_options(
  ${QPP_TARGET}
  INTERFACE
    # MSVC: Disable spurious Eigen warning and enable /bigobj
    $<$<CXX_COMPILER_ID:MSVC>:-D_SILENCE_CXX17_ADAPTOR_TYPEDEFS_DEPRECATION_WARNING;/bigobj>
    #
    # MinGW or Cygwin: Handle issues with large object files
    $<$<AND:$<OR:$<CXX_COMPILER_ID:GNU>,$<CXX_COMPILER_ID:Clang>>,$<PLATFORM_ID:Windows>>:-Wa,-mbig-obj>
    #
    # Cygwin: Use std=gnu++ standard for compatibility (Cygwin uses GNU tools)
    $<$<AND:$<CXX_COMPILER_ID:GNU>,$<PLATFORM_ID:CYGWIN>>:-std=gnu++17>
    #
    # Clang: Force the use of libc++
    $<$<CXX_COMPILER_ID:Clang>:-stdlib=libc++>)

# Linker options
target_link_libraries(
  ${QPP_TARGET}
  INTERFACE
    $<$<AND:$<PLATFORM_ID:Linux>,$<CXX_COMPILER_ID:Clang>,$<VERSION_LESS:$<CXX_COMPILER_VERSION>,4.0>>:supc++>
)

# Debug macros
target_compile_definitions(
  ${QPP_TARGET}
  INTERFACE # Set the DEBUG macro for the Debug configuration.
            $<$<CONFIG:Debug>:DEBUG>
            # Set the EIGEN_NO_DEBUG macro for all non-Debug configurations
            # (MinSizeRel, Release, RelWithDebInfo).
            $<$<NOT:$<CONFIG:Debug>>:EIGEN_NO_DEBUG>)

# Qubit optimization macros
target_compile_definitions(
  ${QPP_TARGET}
  INTERFACE $<$<BOOL:${QPP_QUBIT_OPTIMIZATIONS}>:QPP_QUBIT_OPTIMIZATIONS>)

# thread_local storage duration
target_compile_definitions(
  ${QPP_TARGET} INTERFACE $<$<BOOL:${QPP_NO_THREAD_LOCAL}>:QPP_NO_THREAD_LOCAL>)

# Default build type
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE
      "Release"
      CACHE
        STRING
        "Choose the type of build: None Debug Release MinSizeRel RelWithDebInfo."
        FORCE)
endif()
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}${CMAKE_CONFIGURATION_TYPES}")

# Internal compiler flags/options. This block should remain the final section of
# the file.
#
# These settings are for in-tree Quantum++ development targets only (examples,
# unit tests, optional components, pyqpp development, etc.)
if(TARGET libqpp_internal)
  target_compile_definitions(
    libqpp_internal
    INTERFACE # Debug configuration
              $<$<CONFIG:Debug>:DEBUG>
              $<$<AND:$<CXX_COMPILER_ID:GNU>,$<CONFIG:Debug>>:_GLIBCXX_DEBUG>
              # Non-Debug configurations
              $<$<NOT:$<CONFIG:Debug>>:EIGEN_NO_DEBUG>
              # Enable qubit-specific optimizations if requested
              $<$<BOOL:${QPP_QUBIT_OPTIMIZATIONS}>:QPP_QUBIT_OPTIMIZATIONS>)

  target_compile_options(
    libqpp_internal
    INTERFACE
      # GCC: use -O0 in Debug builds
      $<$<AND:$<CXX_COMPILER_ID:GNU>,$<CONFIG:Debug>>:-O0>
      #
      # GCC or Clang warnings
      $<$<OR:$<CXX_COMPILER_ID:GNU>,$<CXX_COMPILER_ID:Clang>>:
      -pedantic
      -Wall
      -Wextra
      -Weffc++
      >
      #
      # MSVC warnings
      $<$<CXX_COMPILER_ID:MSVC>:
      /W3
      >
      #
      # Use the "no-weak" debugging flag only when debugging under OS X, as gdb
      # cannot step in template functions when debugging code produced by g++,
      # see
      # https://stackoverflow.com/questions/23330641/gnu-gdb-can-not-step-into-template-functions-os-x-mavericks
      # $<$<AND:$<CXX_COMPILER_ID:GNU>,$<CONFIG:Debug>,$<PLATFORM_ID:Darwin>>:-fno-weak>
  )
endif()
