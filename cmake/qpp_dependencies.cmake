# Quantum++ additional dependencies. Do not modify unless you know what you're
# doing.

# Custom index type. If none selected, a default one is selected in
# <qpp/types.hpp> (usually std::size_t).
set(QPP_IDX
    "default"
    CACHE STRING "Default index type, see <qpp/types.hpp>")
set_property(
  CACHE QPP_IDX
  PROPERTY STRINGS
           "default"
           "short"
           "int"
           "long"
           "long long"
           "unsigned short"
           "unsigned int"
           "unsigned long"
           "unsigned long long")

# Determine the index type compile definition
if(QPP_IDX STREQUAL "default")
  message(STATUS "Index type - default (see <qpp/types.hpp>)")
  set(QPP_IDX_DEFINITION QPP_IDX_DEFAULT)
elseif(QPP_IDX STREQUAL "short")
  message(STATUS "Index type - short")
  set(QPP_IDX_DEFINITION QPP_IDX_SHORT)
elseif(QPP_IDX STREQUAL "int")
  message(STATUS "Index type - int")
  set(QPP_IDX_DEFINITION QPP_IDX_INT)
elseif(QPP_IDX STREQUAL "long")
  message(STATUS "Index type - long")
  set(QPP_IDX_DEFINITION QPP_IDX_LONG)
elseif(QPP_IDX STREQUAL "long long")
  message(STATUS "Index type - long long")
  set(QPP_IDX_DEFINITION QPP_IDX_LONG_LONG)
elseif(QPP_IDX STREQUAL "unsigned short")
  message(STATUS "Index type - unsigned short")
  set(QPP_IDX_DEFINITION QPP_IDX_USHORT)
elseif(QPP_IDX STREQUAL "unsigned int")
  message(STATUS "Index type - unsigned int")
  set(QPP_IDX_DEFINITION QPP_IDX_UINT)
elseif(QPP_IDX STREQUAL "unsigned long")
  message(STATUS "Index type - unsigned long")
  set(QPP_IDX_DEFINITION QPP_IDX_ULONG)
elseif(QPP_IDX STREQUAL "unsigned long long")
  message(STATUS "Index type - unsigned long long")
  set(QPP_IDX_DEFINITION QPP_IDX_ULONG_LONG)
endif()

# Custom signed big integer type. If none selected, a default one is selected in
# <qpp/types.hpp> (usually long long).
set(QPP_BIGINT
    "default"
    CACHE STRING "Default big integer type, see <qpp/types.hpp>")
set_property(CACHE QPP_BIGINT PROPERTY STRINGS "default" "short" "int" "long"
                                       "long long")

# Determine the big integer type compile definition
if(QPP_BIGINT STREQUAL "default")
  message(STATUS "Big integer type - default (see <qpp/types.hpp>)")
  set(QPP_BIGINT_DEFINITION QPP_BIGINT_DEFAULT)
elseif(QPP_BIGINT STREQUAL "short")
  message(STATUS "Big integer type - short")
  set(QPP_BIGINT_DEFINITION QPP_BIGINT_SHORT)
elseif(QPP_BIGINT STREQUAL "int")
  message(STATUS "Big integer type - int")
  set(QPP_BIGINT_DEFINITION QPP_BIGINT_INT)
elseif(QPP_BIGINT STREQUAL "long")
  message(STATUS "Big integer type - long")
  set(QPP_BIGINT_DEFINITION QPP_BIGINT_LONG)
elseif(QPP_BIGINT STREQUAL "long long")
  message(STATUS "Big integer type - long long")
  set(QPP_BIGINT_DEFINITION QPP_BIGINT_LONG_LONG)
endif()

# Custom floating-point type. If none selected, a default one is selected in
# <qpp/types.hpp> (usually double).
set(QPP_FP
    "default"
    CACHE STRING "Default floating-point type, see <qpp/types.hpp>")
set_property(CACHE QPP_FP PROPERTY STRINGS "default" "float" "double"
                                   "long double")

# Determine the floating-point type compile definition
if(QPP_FP STREQUAL "default")
  message(STATUS "Floating-point type - default (see <qpp/types.hpp>)")
  set(QPP_FP_DEFINITION QPP_FP_DEFAULT)
elseif(QPP_FP STREQUAL "float")
  message(STATUS "Floating-point type - float")
  set(QPP_FP_DEFINITION QPP_FP_FLOAT)
elseif(QPP_FP STREQUAL "double")
  message(STATUS "Floating-point type - double")
  set(QPP_FP_DEFINITION QPP_FP_DOUBLE)
elseif(QPP_FP STREQUAL "long double")
  message(STATUS "Floating-point type - long double")
  set(QPP_FP_DEFINITION QPP_FP_LONG_DOUBLE)
endif()

# OpenQASM 2.0 specs, see DISCREPANCIES.md for a comparison with Qiskit
option(QASMTOOLS_QASM2_SPECS
       "Use OpenQASM 2.0 standard instead of Qiskit gate specifications" OFF)

if(${QASMTOOLS_QASM2_SPECS})
  message(STATUS "OpenQASM2 specs - ON")
  set(QPP_QASM_DEFINITION QASMTOOLS_QASM2_SPECS=true)
else()
  message(STATUS "OpenQASM2 specs - OFF")
  set(QPP_QASM_DEFINITION QASMTOOLS_QASM2_SPECS=false)
endif()

# OpenMP support
option(QPP_OPENMP "OpenMP support" ON)
set(QPP_OPENMP_LINK_DEPS "") # Initialize OpenMP link dependency
set(QPP_OPENMP_DEFINITION "") # Initialize OpenMP compile definition
if(${QPP_OPENMP})
  find_package(OpenMP)
  if(OpenMP_CXX_FOUND)
    if(OpenMP_CXX_VERSION_MAJOR GREATER_EQUAL 3)
      # Inject definition (as #define) and link dependency
      set(QPP_OPENMP_DEFINITION QPP_OPENMP)
      set(QPP_OPENMP_LINK_DEPS OpenMP::OpenMP_CXX)
    else()
      message(NOTICE "-- Found OpenMP version ${OpenMP_CXX_VERSION}, \
Quantum++ requires OpenMP 3.0 or later")
    endif()
  endif()
endif()

# List of all feature definitions to be passed to a target later
set(QPP_FEATURE_DEFINITIONS
    ${QPP_IDX_DEFINITION} ${QPP_BIGINT_DEFINITION} ${QPP_FP_DEFINITION}
    ${QPP_QASM_DEFINITION} ${QPP_OPENMP_DEFINITION})

# Windows issues with Microsoft Visual Studio
set(MSVC_COMPILE_DEFS "")
set(MSVC_COMPILE_OPTIONS "")
if(MSVC)
  # Disable spurious Eigen warnings with MSVC (warning STL4007)
  set(MSVC_COMPILE_DEFS _SILENCE_CXX17_ADAPTOR_TYPEDEFS_DEPRECATION_WARNING)
  set(MSVC_COMPILE_OPTIONS /bigobj)
endif()

# MinGW or Cygwin have issues with object files that are too large
set(MINGW_CYGWIN_COMPILE_OPTIONS "")
if(MINGW OR CYGWIN)
  set(MINGW_CYGWIN_COMPILE_OPTIONS "-Wa,-mbig-obj")
endif()

# Cygwin used to have issues with std=c++, use std=gnu++ instead
set(CYGWIN_CXX_OPTIONS "")
if(CYGWIN)
  set(CYGWIN_CXX_OPTIONS -std=gnu++17)
endif()

# Force clang to use libc++
set(CLANG_CXX_OPTIONS "")
if(${CMAKE_CXX_COMPILER_ID} MATCHES "Clang")
  set(CLANG_CXX_OPTIONS -stdlib=libc++)
endif()

set(QPP_GCC_DEBUG_OPTIONS "")
# GCC additional debug settings
if(${CMAKE_CXX_COMPILER_ID} MATCHES "GNU")
  set(QPP_GCC_DEBUG_OPTIONS "-Og -D_GLIBCXX_DEBUG")
endif()

# Configurations: Move all CMAKE_CXX_FLAGS_<CONFIG> changes to variables
set(QPP_CONFIG_DEFS_DEBUG "-DDEBUG") # This one remains simple, or can be moved to target
set(QPP_CONFIG_DEFS_EIGEN_NO_DEBUG "-DEIGEN_NO_DEBUG")

# Default build type
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE
      Release
      CACHE STRING "Choose the type of build, options are:
        None Debug Release MinSizeRel RelWithDebInfo." FORCE)
endif()

# On Linux we get a linker error related to __cxa_thread_atexit()
set(QPP_SUPC_LINK_DEPS "")
if(UNIX
   AND NOT (APPLE OR CYGWIN)
   AND (${CMAKE_CXX_COMPILER_ID} MATCHES "Clang")
   AND (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 4.0))
  set(QPP_SUPC_LINK_DEPS supc++)
endif()

# Export all settings to outside
list(APPEND QPP_LINK_DEPS ${QPP_EIGEN3_LINK_DEPS} ${QPP_OPENMP_LINK_DEPS}
     ${QPP_MATLAB_LINK_DEPS} ${QPP_SUPC_LINK_DEPS})

# New variables created for easier use with target_ commands in main CMakeLists.txt:
# QPP_FEATURE_DEFINITIONS (for type selection)
# MSVC_COMPILE_DEFS / MSVC_COMPILE_OPTIONS
# MINGW_CYGWIN_COMPILE_OPTIONS
# CYGWIN_CXX_OPTIONS
# CLANG_CXX_OPTIONS
# QPP_GCC_DEBUG_OPTIONS
# QPP_CONFIG_DEFS_DEBUG
# QPP_CONFIG_DEFS_EIGEN_NO_DEBUG
# QPP_LINK_DEPS (Final list of libraries to link)
