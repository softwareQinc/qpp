#### Quantum++ additional dependencies
#### Do not modify unless you know what you're doing

#### Custom index type. If none selected, a default one is selected by
#### include/types.hpp (usually std::size_t)
set(TYPE_IDX "default" CACHE
        STRING "Default index type, see include/types.hpp")
set_property(CACHE TYPE_IDX PROPERTY
        STRINGS "default" "short" "int" "long" "long long" "unsigned short"
        "unsigned int" "unsigned long" "unsigned long long")

if (TYPE_IDX STREQUAL "default")
    message(STATUS "Index type - default (see include/types.hpp)")
    add_compile_definitions(TYPE_IDX_DEFAULT)
elseif (TYPE_IDX STREQUAL "short")
    message(STATUS "Index type - short")
    add_compile_definitions(TYPE_IDX_SHORT)
elseif (TYPE_IDX STREQUAL "int")
    message(STATUS "Index type - int")
    add_compile_definitions(TYPE_IDX_INT)
elseif (TYPE_IDX STREQUAL "long")
    message(STATUS "Index type - long")
    add_compile_definitions(TYPE_IDX_LONG)
elseif (TYPE_IDX STREQUAL "long long")
    message(STATUS "Index type - long long")
    add_compile_definitions(TYPE_IDX_LONG_LONG)
elseif (TYPE_IDX STREQUAL "unsigned short")
    message(STATUS "Index type - unsigned short")
    add_compile_definitions(TYPE_IDX_USHORT)
elseif (TYPE_IDX STREQUAL "unsigned int")
    message(STATUS "Index type - unsigned int")
    add_compile_definitions(TYPE_IDX_UINT)
elseif (TYPE_IDX STREQUAL "unsigned long")
    message(STATUS "Index type - unsigned long")
    add_compile_definitions(TYPE_IDX_ULONG)
elseif (TYPE_IDX STREQUAL "unsigned long long")
    message(STATUS "Index type - unsigned long long")
    add_compile_definitions(TYPE_IDX_ULONG_LONG)
endif ()

#### Custom signed big integer type. If none selected, a default one is
#### selected by include/types.hpp (usually long long)
set(TYPE_BIGINT "default" CACHE
        STRING "Default big integer type, see include/types.hpp")
set_property(CACHE TYPE_BIGINT PROPERTY
        STRINGS "default" "short" "int" "long" "long long")

if (TYPE_BIGINT STREQUAL "default")
    message(STATUS "Big integer type - default (see include/types.hpp)")
    add_compile_definitions(TYPE_BIGINT_DEFAULT)
elseif (TYPE_BIGINT STREQUAL "short")
    message(STATUS "Big integer type - short")
    add_compile_definitions(TYPE_BIGINT_SHORT)
elseif (TYPE_BIGINT STREQUAL "int")
    message(STATUS "Big integer type - int")
    add_compile_definitions(TYPE_BIGINT_INT)
elseif (TYPE_BIGINT STREQUAL "long")
    message(STATUS "Big integer type - long")
    add_compile_definitions(TYPE_BIGINT_LONG)
elseif (TYPE_BIGINT STREQUAL "long long")
    message(STATUS "Big integer type - long long")
    add_compile_definitions(TYPE_BIGINT_LONG_LONG)
endif ()

#### Custom floating-point type. If none selected, a default one is selected by
#### include/types.hpp (usually double)
set(TYPE_FP "default" CACHE
        STRING "Default floating-point type, see include/types.hpp")
set_property(CACHE TYPE_FP PROPERTY
        STRINGS "default" "float" "double" "long double")

if (TYPE_FP STREQUAL "default")
    message(STATUS "Floating-point type - default (see include/types.hpp)")
    add_compile_definitions(TYPE_FP_DEFAULT)
elseif (TYPE_FP STREQUAL "float")
    message(STATUS "Floating-point type - float")
    add_compile_definitions(TYPE_FP_FLOAT)
elseif (TYPE_FP STREQUAL "double")
    message(STATUS "Floating-point type - double")
    add_compile_definitions(TYPE_FP_DOUBLE)
elseif (TYPE_FP STREQUAL "long double")
    message(STATUS "Floating-point type - long double")
    add_compile_definitions(TYPE_FP_LONG_DOUBLE)
endif ()

#### Enable OpenQASM 2.0 specs, see DISCREPANCIES.md for a comparison with Qiskit
option(USE_OPENQASM2_SPECS "Use OpenQASM 2.0 standard instead of Qiskit gate specifications" OFF)
if (${USE_OPENQASM2_SPECS})
    add_compile_definitions(USE_OPENQASM2_SPECS=true)
else ()
    add_compile_definitions(USE_OPENQASM2_SPECS=false)
endif ()

#### Eigen3
message(STATUS "Detecting Eigen3")
#### Location specified via an environment variable
set(LOCATION_SET_VIA_ENV FALSE)
if (DEFINED ENV{EIGEN3_INSTALL_DIR})
    set(EIGEN3_INSTALL_DIR_ENV $ENV{EIGEN3_INSTALL_DIR})
    set(LOCATION_SET_VIA_ENV TRUE)
endif ()
#### Location set via CMake variable, trumps all other settings
set(EIGEN3_INSTALL_DIR "" CACHE PATH "Path to Eigen3")
if (NOT ${EIGEN3_INSTALL_DIR} STREQUAL "")
    message(STATUS "Overriding automatic Eigen3 detection (EIGEN3_INSTALL_DIR CMake variable)")
    if (IS_DIRECTORY ${EIGEN3_INSTALL_DIR})
        message(STATUS "Detecting Eigen3 - done (in ${EIGEN3_INSTALL_DIR})")
        include_directories(SYSTEM "${EIGEN3_INSTALL_DIR}")
    else ()
        message(FATAL_ERROR "Invalid path to Eigen3 installation")
    endif ()
#### Location set via environment variable
elseif (LOCATION_SET_VIA_ENV)
    message(STATUS "Overriding automatic Eigen3 detection (EIGEN3_INSTALL_DIR environment variable)")
    if (IS_DIRECTORY ${EIGEN3_INSTALL_DIR_ENV})
        message(STATUS "Detecting Eigen3 - done (in ${EIGEN3_INSTALL_DIR_ENV})")
        include_directories(SYSTEM "${EIGEN3_INSTALL_DIR_ENV}")
    else ()
        message(FATAL_ERROR "Invalid path to Eigen3 installation")
    endif ()
#### Try to find the location automatically
else ()
    find_package(Eigen3 3.0 QUIET NO_MODULE)
    if (NOT TARGET Eigen3::Eigen) # did not find Eigen3 automatically
        message(FATAL_ERROR "Eigen3 not detected! Please point EIGEN3_INSTALL_DIR to your Eigen3 location when building with cmake, for example
    cmake --build build -DEIGEN3_INSTALL_DIR=$HOME/eigen3
or set the EIGEN3_INSTALL_DIR environment variable to point to your Eigen3 installation, for example (UNIX/Linux)
    export EIGEN3_INSTALL_DIR=$HOME/eigen3")
        #### Eigen3
        message(STATUS "Detecting Eigen3")
    endif ()
    message(STATUS "Detecting Eigen3 - done (in ${EIGEN3_INCLUDE_DIR})")
    # Eigen3 header-only dependencies to be injected in the main CMakeLists.txt
    set(QPP_EIGEN3_LINK_DEPS Eigen3::Eigen)
endif ()

#### OpenMP support
option(WITH_OPENMP "OpenMP support" ON)
if (${WITH_OPENMP})
    find_package(OpenMP)
    if (OpenMP_CXX_FOUND)
        if (OpenMP_CXX_VERSION_MAJOR GREATER_EQUAL 3)
            #### inject definition (as #define) in the source files
            add_compile_definitions(HAS_OPENMP)
            # OpenMP linking dependencies to be injected in the main CMakeLists.txt
            set(QPP_OPENMP_LINK_DEPS OpenMP::OpenMP_CXX)
        else ()
            message(NOTICE "Found OpenMP version ${OpenMP_CXX_VERSION}, Quantum++ \
requires OpenMP 3.0 or later")
        endif ()
    endif ()
endif ()

#### MATLAB support, disabled by default
option(WITH_MATLAB "MATLAB support" OFF)
if (${WITH_MATLAB})
    message(STATUS "Detecting MATLAB")
    #### Try to find it automatically
    find_package(Matlab OPTIONAL_COMPONENTS MX_LIBRARY MAT_LIBRARY QUIET)
    if (MATLAB_FOUND)
        message(STATUS "Detecting MATLAB - done (in ${Matlab_ROOT_DIR})")
        include_directories(SYSTEM ${Matlab_INCLUDE_DIRS})
        if (WIN32)
            if (MSVC)
                set(MATLAB_LIB_DIR
                        "${Matlab_ROOT_DIR}/extern/lib/win64/microsoft"
                        CACHE PATH "Custom path to MATLAB lib directory")
            elseif (MINGW)
                set(MATLAB_LIB_DIR "${Matlab_ROOT_DIR}/extern/lib/win64/mingw64"
                        CACHE PATH "Custom path to MATLAB lib directory")
            else ()
                message(FATAL_ERROR "Platform not supported, aborting.")
            endif ()
        elseif (UNIX AND NOT APPLE)
            set(MATLAB_LIB_DIR "${Matlab_ROOT_DIR}/bin/glnxa64" CACHE
                    PATH "Custom path to MATLAB lib directory")
        elseif (APPLE)
            if (${CMAKE_HOST_SYSTEM_PROCESSOR} MATCHES "arm64")
                message(FATAL_ERROR "arm64 architecture is not (yet) supported by MATLAB, aborting.")
            endif ()
            set(MATLAB_LIB_DIR "${Matlab_ROOT_DIR}/bin/maci64" CACHE
                    PATH "Custom path to MATLAB lib directory")
        else ()
            message(FATAL_ERROR "Platform not supported, aborting.")
        endif ()
        link_directories(${MATLAB_LIB_DIR})
        add_compile_definitions(HAS_MATLAB_ENABLED)
        set(BUILD_WITH_MATLAB TRUE)
    else () #### Location manually specified
        set(MATLAB_INSTALL_DIR "" CACHE PATH
                "Custom path to MATLAB installation")
        if (IS_DIRECTORY ${MATLAB_INSTALL_DIR})
            #### MATLAB include files
            set(MATLAB_INCLUDE_DIR "${MATLAB_INSTALL_DIR}/extern/include"
                    CACHE PATH "Custom path to MATLAB include directory")
            if (IS_DIRECTORY ${MATLAB_INCLUDE_DIR})
                include_directories(SYSTEM ${MATLAB_INCLUDE_DIR})
            else ()
                message(FATAL_ERROR
                        "Possibly corrupted MATLAB include headers")
            endif ()
            #### MATLAB linker files
            if (WIN32)
                if (MSVC)
                    set(MATLAB_LIB_DIR
                            "${MATLAB_INSTALL_DIR}/extern/lib/win64/microsoft"
                            CACHE PATH "Custom path to MATLAB lib directory")
                elseif (MINGW64)
                    set(MATLAB_LIB_DIR
                            "${MATLAB_INSTALL_DIR}/extern/lib/win64/mingw64"
                            CACHE PATH "Custom path to MATLAB lib directory")
                else ()
                    message(FATAL_ERROR "Platform not supported, aborting.")
                endif ()
            elseif (UNIX AND NOT APPLE)
                set(MATLAB_LIB_DIR "${MATLAB_INSTALL_DIR}/bin/glnxa64" CACHE
                        PATH "Custom path to MATLAB lib directory")
            elseif (APPLE)
                set(MATLAB_LIB_DIR "${MATLAB_INSTALL_DIR}/bin/maci64" CACHE
                        PATH "Custom path to MATLAB lib directory")
            else ()
                message(FATAL_ERROR "Platform not supported, aborting.")
            endif ()
            if (IS_DIRECTORY ${MATLAB_LIB_DIR})
                link_directories(${MATLAB_LIB_DIR})
            else ()
                message(FATAL_ERROR
                        "Possibly corrupted MATLAB compiler libraries")
            endif ()
            #### Everything is OK, inject definition (as #define) in the source
            message(STATUS
                    "Detecting MATLAB - done (in ${MATLAB_INSTALL_DIR})")
            add_compile_definitions(HAS_MATLAB_ENABLED)
            set(BUILD_WITH_MATLAB TRUE)
        else ()
            message(FATAL_ERROR "Could not detect MATLAB, aborting")
        endif ()
    endif ()
    # MATLAB linking dependencies to be injected in the main CMakeLists.txt
    if (${BUILD_WITH_MATLAB})
        if (WIN32)
            if (MSVC)
                set(QPP_MATLAB_LINK_DEPS libmx libmat)
            elseif (MINGW)
                set(QPP_MATLAB_LINK_DEPS mx mat)
            endif ()
        else ()
            set(QPP_MATLAB_LINK_DEPS mx mat)
        endif ()
    endif ()
endif ()

##### Disable support for thread_local storage duration specifier
##### when using AppleClang as libc++ doesn't yet support it
#if (${CMAKE_CXX_COMPILER_ID} STREQUAL "AppleClang")
#    #### inject definition (as #define) in the source files
#    add_definitions(-DNO_THREAD_LOCAL_)
#    message(WARNING "Detected compiler: ${CMAKE_CXX_COMPILER_ID} \
#    ${CMAKE_CXX_COMPILER_VERSION}. thread_local not supported.")
#endif ()

#### Windows issues with Microsoft Visual Studio
if (MSVC)
    # Disable spurious Eigen warnings with MSVC (warning STL4007)
    add_compile_definitions(_SILENCE_CXX17_ADAPTOR_TYPEDEFS_DEPRECATION_WARNING)
    add_compile_options(-bigobj)
    add_compile_definitions(NOMINMAX)
endif ()

#### MinGW or Cygwin have issues with object files that are too large
if (MINGW OR CYGWIN)
    add_compile_options("-Wa,-mbig-obj")
endif ()

#### Cygwin has issues with std=c++11, use std=gnu++11 instead
if (CYGWIN)
    add_compile_options(-std=gnu++11)
endif ()

#### Force clang to use libc++
if (${CMAKE_CXX_COMPILER_ID} MATCHES "Clang")
    add_compile_options(-stdlib=libc++)
    list(APPEND QPP_LINK_DEPS c++)
endif ()

#### GNU gcc additional debug settings
if (${CMAKE_CXX_COMPILER_ID} MATCHES "GNU")
    #if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    #    ## use the "no-weak" debugging flag only when debugging under OS X,
    #    ## as gdb cannot step in template functions when debugging code
    #    ## produced by g++
    #    ## see https://stackoverflow.com/questions/23330641/gnu-gdb-can-not-step-into-template-functions-os-x-mavericks
    #    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fno-weak")
    #endif ()
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Og -D_GLIBCXX_DEBUG")
endif ()

#### Configurations
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DDEBUG")
set(CMAKE_CXX_FLAGS_MINSIZEREL
        "${CMAKE_CXX_FLAGS_MINSIZEREL} -DEIGEN_NO_DEBUG")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -DEIGEN_NO_DEBUG")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO
        "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} -DEIGEN_NO_DEBUG")

#### Default build type
if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING
            "Choose the type of build, options are: \
         None Debug Release MinSizeRel RelWithDebInfo."
            FORCE)
endif ()

#### On Linux we get a linker error related to __cxa_thread_atexit()
if (UNIX AND NOT (APPLE OR CYGWIN)
        AND (${CMAKE_CXX_COMPILER_ID} MATCHES "Clang")
        AND (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 4.0))
    list(APPEND QPP_LINK_DEPS supc++)
endif ()

#### Export all settings to outside
list(APPEND QPP_LINK_DEPS ${QPP_EIGEN3_LINK_DEPS} ${QPP_OPENMP_LINK_DEPS} ${QPP_MATLAB_LINK_DEPS})
