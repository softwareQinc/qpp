# Quantum++ additional dependencies Do not modify unless you know what you're
# doing

# Custom index type. If none selected, a default one is selected by
# include/types.hpp (usually std::size_t).
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

if(QPP_IDX STREQUAL "default")
  message(STATUS "Index type - default (see <qpp/types.hpp>)")
  add_compile_definitions(QPP_IDX_DEFAULT)
elseif(QPP_IDX STREQUAL "short")
  message(STATUS "Index type - short")
  add_compile_definitions(QPP_IDX_SHORT)
elseif(QPP_IDX STREQUAL "int")
  message(STATUS "Index type - int")
  add_compile_definitions(QPP_IDX_INT)
elseif(QPP_IDX STREQUAL "long")
  message(STATUS "Index type - long")
  add_compile_definitions(QPP_IDX_LONG)
elseif(QPP_IDX STREQUAL "long long")
  message(STATUS "Index type - long long")
  add_compile_definitions(QPP_IDX_LONG_LONG)
elseif(QPP_IDX STREQUAL "unsigned short")
  message(STATUS "Index type - unsigned short")
  add_compile_definitions(QPP_IDX_USHORT)
elseif(QPP_IDX STREQUAL "unsigned int")
  message(STATUS "Index type - unsigned int")
  add_compile_definitions(QPP_IDX_UINT)
elseif(QPP_IDX STREQUAL "unsigned long")
  message(STATUS "Index type - unsigned long")
  add_compile_definitions(QPP_IDX_ULONG)
elseif(QPP_IDX STREQUAL "unsigned long long")
  message(STATUS "Index type - unsigned long long")
  add_compile_definitions(QPP_IDX_ULONG_LONG)
endif()

# Custom signed big integer type. If none selected, a default one is selected by
# include/types.hpp (usually long long).
set(QPP_BIGINT
    "default"
    CACHE STRING "Default big integer type, see <qpp/types.hpp>")
set_property(CACHE QPP_BIGINT PROPERTY STRINGS "default" "short" "int" "long"
                                       "long long")

if(QPP_BIGINT STREQUAL "default")
  message(STATUS "Big integer type - default (see <qpp/types.hpp>)")
  add_compile_definitions(QPP_BIGINT_DEFAULT)
elseif(QPP_BIGINT STREQUAL "short")
  message(STATUS "Big integer type - short")
  add_compile_definitions(QPP_BIGINT_SHORT)
elseif(QPP_BIGINT STREQUAL "int")
  message(STATUS "Big integer type - int")
  add_compile_definitions(QPP_BIGINT_INT)
elseif(QPP_BIGINT STREQUAL "long")
  message(STATUS "Big integer type - long")
  add_compile_definitions(QPP_BIGINT_LONG)
elseif(QPP_BIGINT STREQUAL "long long")
  message(STATUS "Big integer type - long long")
  add_compile_definitions(QPP_BIGINT_LONG_LONG)
endif()

# Custom floating-point type. If none selected, a default one is selected by
# include/types.hpp (usually double).
set(QPP_FP
    "default"
    CACHE STRING "Default floating-point type, see <qpp/types.hpp>")
set_property(CACHE QPP_FP PROPERTY STRINGS "default" "float" "double"
                                   "long double")

if(QPP_FP STREQUAL "default")
  message(STATUS "Floating-point type - default (see <qpp/types.hpp>)")
  add_compile_definitions(QPP_FP_DEFAULT)
elseif(QPP_FP STREQUAL "float")
  message(STATUS "Floating-point type - float")
  add_compile_definitions(QPP_FP_FLOAT)
elseif(QPP_FP STREQUAL "double")
  message(STATUS "Floating-point type - double")
  add_compile_definitions(QPP_FP_DOUBLE)
elseif(QPP_FP STREQUAL "long double")
  message(STATUS "Floating-point type - long double")
  add_compile_definitions(QPP_FP_LONG_DOUBLE)
endif()

# OpenQASM 2.0 specs, see DISCREPANCIES.md for a comparison with Qiskit
option(QASMTOOLS_QASM2_SPECS
       "Use OpenQASM 2.0 standard instead of Qiskit gate specifications" OFF)
if(${QASMTOOLS_QASM2_SPECS})
  add_compile_definitions(QASMTOOLS_QASM2_SPECS=true)
  message(STATUS "OpenQASM2 specs - ON")
else()
  add_compile_definitions(QASMTOOLS_QASM2_SPECS=false)
  message(STATUS "OpenQASM2 specs - OFF")
endif()

# OpenMP support
option(QPP_OPENMP "OpenMP support" ON)
if(${QPP_OPENMP})
  find_package(OpenMP)
  if(OpenMP_CXX_FOUND)
    if(OpenMP_CXX_VERSION_MAJOR GREATER_EQUAL 3)
      # inject definition (as #define) in the source files
      add_compile_definitions(QPP_OPENMP)
      # OpenMP linking dependencies to be injected in the main CMakeLists.txt
      set(QPP_OPENMP_LINK_DEPS OpenMP::OpenMP_CXX)
    else()
      message(NOTICE "Found OpenMP version ${OpenMP_CXX_VERSION}, Quantum++ \
requires OpenMP 3.0 or later")
    endif()
  endif()
endif()

# Disable support for thread_local storage duration specifier when using
# AppleClang as libc++ doesn't yet support it if (${CMAKE_CXX_COMPILER_ID}
# STREQUAL "AppleClang") #### inject definition (as #define) in the source files
# add_definitions(-DNO_THREAD_LOCAL_) message(WARNING "Detected compiler:
# ${CMAKE_CXX_COMPILER_ID} \ ${CMAKE_CXX_COMPILER_VERSION}. thread_local not
# supported.") endif ()

# Windows issues with Microsoft Visual Studio
if(MSVC)
  # Disable spurious Eigen warnings with MSVC (warning STL4007)
  add_compile_definitions(_SILENCE_CXX17_ADAPTOR_TYPEDEFS_DEPRECATION_WARNING)
  add_compile_options(-bigobj)
endif()

# MinGW or Cygwin have issues with object files that are too large
if(MINGW OR CYGWIN)
  add_compile_options("-Wa,-mbig-obj")
endif()

# Cygwin has issues with std=c++11, use std=gnu++11 instead
if(CYGWIN)
  add_compile_options(-std=gnu++11)
endif()

# Force clang to use libc++
if(${CMAKE_CXX_COMPILER_ID} MATCHES "Clang")
  add_compile_options(-stdlib=libc++)
endif()

# GCC additional debug settings
if(${CMAKE_CXX_COMPILER_ID} MATCHES "GNU")
  # if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin") ## use the "no-weak" debugging
  # flag only when debugging under OS X, ## as gdb cannot step in template
  # functions when debugging code ## produced by g++ ## see
  # https://stackoverflow.com/questions/23330641/gnu-gdb-can-not-step-into-template-functions-os-x-mavericks
  # set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fno-weak") endif ()
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Og -D_GLIBCXX_DEBUG")
endif()

# Configurations
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DDEBUG")
set(CMAKE_CXX_FLAGS_MINSIZEREL "${CMAKE_CXX_FLAGS_MINSIZEREL} -DEIGEN_NO_DEBUG")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -DEIGEN_NO_DEBUG")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO
    "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} -DEIGEN_NO_DEBUG")

# Default build type
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE
      Release
      CACHE STRING "Choose the type of build, options are: \
         None Debug Release MinSizeRel RelWithDebInfo." FORCE)
endif()

# On Linux we get a linker error related to __cxa_thread_atexit()
if(UNIX
   AND NOT (APPLE OR CYGWIN)
   AND (${CMAKE_CXX_COMPILER_ID} MATCHES "Clang")
   AND (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 4.0))
  list(APPEND QPP_LINK_DEPS supc++)
endif()

# Export all settings to outside
list(APPEND QPP_LINK_DEPS ${QPP_EIGEN3_LINK_DEPS} ${QPP_OPENMP_LINK_DEPS}
     ${QPP_MATLAB_LINK_DEPS})
