#
# This file is part of Quantum++.
#
# MIT License
#
# Copyright (c) 2013 - 2019 Vlad Gheorghiu (vgheorgh@gmail.com)
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
#

CMAKE_MINIMUM_REQUIRED(VERSION 3.2)
PROJECT(qpp VERSION 1.2 LANGUAGES CXX)
SET(CMAKE_CXX_STANDARD 11)
SET(CMAKE_CXX_STANDARD_REQUIRED ON)

#### Quantum++ headers
INCLUDE_DIRECTORIES(include)

#### Source file(s) to be compiled, modify as needed
AUX_SOURCE_DIRECTORY(examples EXAMPLE_FILES)
AUX_SOURCE_DIRECTORY(examples/circuits EXAMPLE_FILES)

#### Eigen3
MESSAGE(STATUS "Detecting Eigen3")
SET(EIGEN3_INCLUDE_DIR "" CACHE PATH "Path to Eigen3")
#### Location manually specified
IF (NOT ${EIGEN3_INCLUDE_DIR} STREQUAL "")
    IF (IS_DIRECTORY ${EIGEN3_INCLUDE_DIR})
        MESSAGE(STATUS "Detecting Eigen3 - done (in ${EIGEN3_INCLUDE_DIR})")
        INCLUDE_DIRECTORIES(SYSTEM "${EIGEN3_INCLUDE_DIR}")
    ELSE ()
        MESSAGE(FATAL_ERROR "Invalid path to Eigen3 installation")
    ENDIF ()
ELSE () #### Try to find it automatically
    FIND_PACKAGE(Eigen3 3.0 QUIET NO_MODULE)
    IF (NOT TARGET Eigen3::Eigen) # did not find Eigen3 automatically
        MESSAGE(FATAL_ERROR
                "Eigen3 not detected! Please point EIGEN3_INCLUDE_DIR\
            to your Eigen3 location when building with cmake,\
            for example\

            cmake -DEIGEN3_INCLUDE_DIR=$HOME/eigen ..")
    ENDIF ()
    MESSAGE(STATUS "Detecting Eigen3 - done")
ENDIF ()

#### Cygwin has issues with std=c++11, use std=gnu++11 instead
IF (CYGWIN)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++11")
ENDIF ()

#### Additional flags for GNU gcc and Clang/AppleClang
IF (${CMAKE_CXX_COMPILER_ID} MATCHES "Clang" OR ${CMAKE_CXX_COMPILER_ID}
        STREQUAL "GNU")
    SET(CMAKE_INCLUDE_SYSTEM_FLAG_CXX "-isystem")
    # Enable all warnings
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pedantic -Wall -Wextra -Weffc++")
ENDIF ()

#### Disable support for thread_local storage duration specifier
#### when using AppleClang as libc++ doesn't yet support it
IF (${CMAKE_CXX_COMPILER_ID} STREQUAL "AppleClang")
    #### inject definition (as #define) in the source files
    ADD_DEFINITIONS(-DNO_THREAD_LOCAL_)
    MESSAGE(WARNING "Detected compiler: ${CMAKE_CXX_COMPILER_ID} \
    ${CMAKE_CXX_COMPILER_VERSION}. thread_local not supported.")
ENDIF ()

#### MATLAB support
SET(WITH_MATLAB "" CACHE PATH "Path to MATLAB (to enable MATLAB support)")
# WITH_MATLAB option was activated (set)
IF (NOT ${WITH_MATLAB} STREQUAL "")
    IF (IS_DIRECTORY ${WITH_MATLAB})
        MESSAGE(STATUS "MATLAB support enabled")
        # MATLAB include files
        SET(MATLAB_INCLUDE "${WITH_MATLAB}/extern/include" CACHE PATH
                "Path to MATLAB include directory")
        IF (IS_DIRECTORY ${MATLAB_INCLUDE})
            INCLUDE_DIRECTORIES(SYSTEM ${MATLAB_INCLUDE})
        ELSE ()
            MESSAGE(FATAL_ERROR "Possibly corrupted MATLAB include headers")
        ENDIF ()
        # MATLAB linker files
        SET(MATLAB_LINK "${WITH_MATLAB}/bin/maci64" CACHE PATH
                "Path to MATLAB lib directory")
        IF (IS_DIRECTORY ${MATLAB_LINK})
            LINK_DIRECTORIES(${MATLAB_LINK})
        ELSE ()
            MESSAGE(FATAL_ERROR "Possibly corrupted MATLAB linker files")
        ENDIF ()
        #### Everything is OK, inject definition (as #define) in the source
        ADD_DEFINITIONS(-DWITH_MATLAB_)
        SET(BUILD_WITH_MATLAB TRUE)
    ELSE ()
        MESSAGE(FATAL_ERROR "Invalid path to MATLAB installation")
    ENDIF ()
ENDIF ()

#### OpenMP support
OPTION(WITH_OPENMP "OpenMP support" ON)
IF (${WITH_OPENMP})
    #### Enable only for GNU gcc or clang >= 3.8
    IF (${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU" OR
    (${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang" AND
    (NOT CMAKE_CXX_COMPILER_VERSION VERSION_LESS "3.8")))
        SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
        #### inject definition (as #define) in the source files
        ADD_DEFINITIONS(-DWITH_OPENMP_)
    ELSE ()
        MESSAGE(WARNING "Detected compiler: ${CMAKE_CXX_COMPILER_ID} \
        ${CMAKE_CXX_COMPILER_VERSION}. Support for OpenMP is enabled only for \
        the GNU gcc compiler or the clang compiler version version 3.8 or \
        higher.")
    ENDIF ()
ENDIF ()

#### GNU gcc additional debug settings
IF (${CMAKE_CXX_COMPILER_ID} MATCHES "GNU")
    IF (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
        ## use the "no-weak" debugging flag only when debugging under OS X,
        ## as gdb cannot step in template functions when debugging code
        ## produced by g++
        ## see http://stackoverflow.com/questions/
        ## 23330641/gnu-gdb-can-not-step-into-template-functions-os-x-mavericks
        SET(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fno-weak")
    ENDIF ()
    SET(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Og -D_GLIBCXX_DEBUG")
ENDIF ()

#### Configurations
SET(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DDEBUG")
SET(CMAKE_CXX_FLAGS_MINSIZEREL
        "${CMAKE_CXX_FLAGS_MINSIZEREL} -DEIGEN_NO_DEBUG")
SET(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -DEIGEN_NO_DEBUG")
SET(CMAKE_CXX_FLAGS_RELWITHDEBINFO
        "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} -DEIGEN_NO_DEBUG")

#### Default build type
IF (NOT CMAKE_BUILD_TYPE)
    SET(CMAKE_BUILD_TYPE Release CACHE STRING
            "Choose the type of build, options are: \
         None Debug Release MinSizeRel RelWithDebInfo."
            FORCE)
ENDIF ()

#### Build all examples in ${EXAMPLE_FILES}
FOREACH (FILE ${EXAMPLE_FILES})
    GET_FILENAME_COMPONENT(TARGET_NAME ${FILE} NAME_WE)
    #### MATALB support
    IF (${BUILD_WITH_MATLAB})
        ADD_EXECUTABLE(${TARGET_NAME} ${FILE})
        TARGET_LINK_LIBRARIES(${TARGET_NAME} mx mat)
    ELSE () #### No MATLAB support
        #### Do not build examples/matlab.cpp if there's no MATLAB support
        IF (${FILE} STREQUAL "examples/matlab.cpp")
            CONTINUE()
        ELSE ()
            ADD_EXECUTABLE(${TARGET_NAME} ${FILE})
        ENDIF ()
    ENDIF ()
    #### Eigen3 was found automatically
    IF (TARGET Eigen3::Eigen)
        TARGET_LINK_LIBRARIES(${TARGET_NAME} Eigen3::Eigen)
    ENDIF ()
    #### OpenMP was enabled
    IF ($WITH_OPENMP$ AND ${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang"
            AND CMAKE_CXX_COMPILER_VERSION VERSION_GREATER "3.7")
        TARGET_LINK_LIBRARIES(TARGET_NAME omp)
    ENDIF ()
ENDFOREACH ()
