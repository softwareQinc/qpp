# Compiler options
target_compile_options(
  libqpp
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
    $<$<CXX_COMPILER_ID:Clang>:-stdlib=libc++>
    #
    # GCC: Additional debug settings (-Og and _GLIBCXX_DEBUG)
    $<$<AND:$<CXX_COMPILER_ID:GNU>,$<CONFIG:Debug>>:-Og;-D_GLIBCXX_DEBUG>
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
    /W4
    /WX
    >
    #
    # Use the "no-weak" debugging flag only when debugging under OS X, as gdb
    # cannot step in template functions when debugging code produced by g++, see
    # https://stackoverflow.com/questions/23330641/gnu-gdb-can-not-step-into-template-functions-os-x-mavericks
    # $<$<AND:$<CXX_COMPILER_ID:GNU>,$<CONFIG:Debug>,$<PLATFORM_ID:Darwin>>:-fno-weak>
)

# Linker options
target_link_libraries(
  libqpp
  INTERFACE
    $<$<AND:$<PLATFORM_ID:Linux>,$<CXX_COMPILER_ID:Clang>,$<VERSION_LESS:$<CXX_COMPILER_VERSION>,4.0>>:supc++>
)

# Debug macros
target_compile_definitions(
  libqpp
  INTERFACE # Set the DEBUG macro for the Debug configuration.
            $<$<CONFIG:Debug>:DEBUG>
            # Set the EIGEN_NO_DEBUG macro for all non-Debug configurations
            # (MinSizeRel, Release, RelWithDebInfo).
            $<$<NOT:$<CONFIG:Debug>>:EIGEN_NO_DEBUG>)

# Default build type
set(CMAKE_BUILD_TYPE
    Release
    CACHE STRING "Choose the type of build, options are:
        None Debug Release MinSizeRel RelWithDebInfo.")
