if (APPLE)
    cmake_minimum_required(VERSION 3.12)
else ()
    cmake_minimum_required(VERSION 3.9)
endif ()
project(qpp_testing)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

#### Source directory
add_definitions(-DPATH="${CMAKE_SOURCE_DIR}")

#### Quantum++ headers
include_directories(../include)

#### Link Google Test dynamically if using MSVC
if (MSVC)
    set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
endif ()

#### Eigen3
message(STATUS "Detecting Eigen3")
set(EIGEN3_INSTALL_DIR "" CACHE PATH "Path to Eigen3")
#### Location manually specified
if (NOT ${EIGEN3_INSTALL_DIR} STREQUAL "")
    if (IS_DIRECTORY ${EIGEN3_INSTALL_DIR})
        message(STATUS "Detecting Eigen3 - done (in ${EIGEN3_INSTALL_DIR})")
        include_directories(SYSTEM "${EIGEN3_INSTALL_DIR}")
    else ()
        message(FATAL_ERROR "Invalid path to Eigen3 installation")
    endif ()
else () #### Try to find it automatically
    find_package(Eigen3 3.0 QUIET NO_MODULE)
    if (NOT TARGET Eigen3::Eigen) # did not find Eigen3 automatically
        message(FATAL_ERROR
                "Eigen3 not detected! Please point EIGEN3_INSTALL_DIR\
            to your Eigen3 location when building with cmake,\
            for example\

            cmake .. -DEIGEN3_INSTALL_DIR=$HOME/eigen")
    endif ()
    message(STATUS "Detecting Eigen3 - done (in ${EIGEN3_INCLUDE_DIR})")
endif ()

#### OpenMP support
find_package(OpenMP)
if (OpenMP_CXX_FOUND AND (OpenMP_CXX_VERSION_MAJOR GREATER_EQUAL 3))
    #### inject definition (as #define) in the source files
    add_definitions(-DHAS_OPENMP)
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
                message(FATAL "Platform not supported, aborting.")
            endif ()
        elseif (UNIX AND NOT APPLE)
            set(MATLAB_LIB_DIR "${Matlab_ROOT_DIR}/bin/glnxa64" CACHE
                    PATH "Custom path to MATLAB lib directory")
        elseif (APPLE)
            set(MATLAB_LIB_DIR "${Matlab_ROOT_DIR}/bin/maci64" CACHE
                    PATH "Custom path to MATLAB lib directory")
        else ()
            message(FATAL "Platform not supported, aborting.")
        endif ()
        link_directories(${MATLAB_LIB_DIR})
        add_definitions(-DHAS_MATLAB_ENABLED)
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
                    message(FATAL "Platform not supported, aborting.")
                endif ()
            elseif (UNIX AND NOT APPLE)
                set(MATLAB_LIB_DIR "${MATLAB_INSTALL_DIR}/bin/glnxa64" CACHE
                        PATH "Custom path to MATLAB lib directory")
            elseif (APPLE)
                set(MATLAB_LIB_DIR "${MATLAB_INSTALL_DIR}/bin/maci64" CACHE
                        PATH "Custom path to MATLAB lib directory")
            else ()
                message(FATAL "Platform not supported, aborting.")
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
            add_definitions(-DHAS_MATLAB_ENABLED)
            set(BUILD_WITH_MATLAB TRUE)
        else ()
            message(FATAL_ERROR "Could not detect MATLAB, aborting")
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

#### Force clang to use libc++
if (${CMAKE_CXX_COMPILER_ID} MATCHES "Clang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
endif ()

#### Cygwin has issues with std=c++11, use std=gnu++11 instead
if (CYGWIN)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++11")
endif ()

#### Additional flags for GNU gcc and Clang/AppleClang
if (${CMAKE_CXX_COMPILER_ID} MATCHES "Clang" OR ${CMAKE_CXX_COMPILER_ID}
        STREQUAL "GNU")
    set(CMAKE_INCLUDE_SYSTEM_FLAG_CXX "-isystem")
endif ()

#### Google Test
add_subdirectory(lib/googletest-release-1.10.0)

#### Enable all warnings only after compiling Google Test
#### Additional flags for GNU gcc and Clang/AppleClang
if (${CMAKE_CXX_COMPILER_ID} MATCHES "Clang" OR ${CMAKE_CXX_COMPILER_ID}
        STREQUAL "GNU")
    #### Enable all warnings
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pedantic -Wall -Wextra -Weffc++")
endif ()

#### Unit tests
add_subdirectory(tests)

#### GNU gcc additional debug settings
if (${CMAKE_CXX_COMPILER_ID} MATCHES "GNU")
    #    if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    #        ## use the "no-weak" debugging flag only when debugging under OS X,
    #        ## as gdb cannot step in template functions when debugging code
    #        ## produced by g++
    #        ## see http://stackoverflow.com/questions/
    #        ## 23330641/gnu-gdb-can-not-step-into-template-functions-os-x-mavericks
    #        set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fno-weak")
    #    endif ()
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
