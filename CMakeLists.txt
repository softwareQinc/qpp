# Quantum++
#
# Copyright (c) 2013 - 2016 Vlad Gheorghiu (vgheorgh@gmail.com)
#
# This file is part of Quantum++.
#
# Quantum++ is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Quantum++ is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Quantum++.  If not, see <http://www.gnu.org/licenses/>.

CMAKE_MINIMUM_REQUIRED(VERSION 3.0.0)
PROJECT(qpp VERSION 1.0.0 LANGUAGES CXX)

#### Use C++11
IF(CYGWIN)
    #### cygwin has issues with std=c++11, use std=gnu++11 instead
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++11")
ELSE()
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
ENDIF()

SET(ADDITIONAL_FLAGS "-pedantic -Wall -Wextra -Weffc++")
SET(CMAKE_INCLUDE_SYSTEM_FLAG_CXX "-isystem")

#### CMake support only for GNU gcc and Clang/AppleClang
#### AppleClang does not support OpenMP; only Clang >= 3.7 supports OpenMP
IF(NOT ${CMAKE_CXX_COMPILER_ID} MATCHES "Clang"
        AND NOT ${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
    MESSAGE(FATAL_ERROR
            "This CMakeLists.txt works only with GNU gcc/Clang/AppleClang")
ENDIF()

#### Quantum++ headers
INCLUDE_DIRECTORIES(include)

#### Source file(s) to be compiled, modify as needed
SET(SOURCE_FILES
        examples/_test.cpp)

#### Eigen3, modify as needed
INCLUDE_DIRECTORIES(SYSTEM "$ENV{HOME}/eigen")

#### MATLAB support
OPTION(WITH_MATLAB "MATLAB support" OFF)
IF(${WITH_MATLAB})
    #### inject definition (as #define) in the source files
    ADD_DEFINITIONS(-DWITH_MATLAB_)
    #### modify as needed
    SET(MATLAB "/Applications/MATLAB_R2016a.app")
    INCLUDE_DIRECTORIES(SYSTEM "${MATLAB}/extern/include")
    SET(CMAKE_EXE_LINKER_FLAGS
            "${CMAKE_EXE_LINKER_FLAGS} -L${MATLAB}/bin/maci64")
ENDIF()

#### OpenMP support
OPTION(WITH_OPENMP "OpenMP support" ON)
IF(${WITH_OPENMP})
    #### inject definition (as #define) in the source files
    ADD_DEFINITIONS(-DWITH_OPENMP_)
    #### GNU gcc
    IF(${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
        SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
        SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS}")
    #### Clang but not AppleClang
    ELSEIF(${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang")
        IF(CMAKE_CXX_COMPILER_VERSION VERSION_GREATER "3.7")
            #### modify as needed
            SET(CLANG_LIBOMP "/opt/local/lib/libomp/")
            SET(CLANG_LIBOMP_INCLUDE "/opt/local/include/libomp/")
            INCLUDE_DIRECTORIES(SYSTEM "${CLANG_LIBOMP_INCLUDE}")
            SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp=libomp")
            SET(CMAKE_EXE_LINKER_FLAGS
                    "${CMAKE_EXE_LINKER_FLAGS} -L${CLANG_LIBOMP}")
        ENDIF()
    #### no OpenMP support for other compilers
    ELSE()
        MESSAGE(WARNING "Detected compiler: ${CMAKE_CXX_COMPILER_ID} \
        ${CMAKE_CXX_COMPILER_VERSION}. Support for OpenMP is enabled only for \
        the GNU gcc compiler or Clang >= 3.7.")
    ENDIF()
ENDIF()

#### GoogleTest
OPTION(WITH_UNIT_TESTS "Enable unit tests" ON)
IF(${WITH_UNIT_TESTS})
    ADD_SUBDIRECTORY(unit_tests)
ENDIF()

#### Enable all warnings
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${ADDITIONAL_FLAGS}")

#### Disable support for thread_local storage duration specifier
#### when using AppleClang as libc++ doesn't yet support it
IF(${CMAKE_CXX_COMPILER_ID} STREQUAL "AppleClang")
    #### inject definition (as #define) in the source files
    ADD_DEFINITIONS(-DNO_THREAD_LOCAL_)
    MESSAGE(WARNING "Detected compiler: ${CMAKE_CXX_COMPILER_ID} \
    ${CMAKE_CXX_COMPILER_VERSION}. thread_local not supported.")
ENDIF()

#### GNU gcc additional debug settings
IF(${CMAKE_CXX_COMPILER_ID} MATCHES "GNU")
    IF(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
        ## use the "no-weak" debugging flag only when debugging under OS X,
        ## as gdb cannot step in template functions when debugging code
        ## produced by g++
        ## see http://stackoverflow.com/questions/
        ## 23330641/gnu-gdb-can-not-step-into-template-functions-os-x-mavericks
        SET(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fno-weak")
    ENDIF()
    SET(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Og -D_GLIBCXX_DEBUG")
ENDIF()

#### Configurations
SET(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DDEBUG")
SET(CMAKE_CXX_FLAGS_MINSIZEREL
        "${CMAKE_CXX_FLAGS_MINSIZEREL} -DEIGEN_NO_DEBUG")
SET(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -DEIGEN_NO_DEBUG")
SET(CMAKE_CXX_FLAGS_RELWITHDEBINFO
        "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} -DEIGEN_NO_DEBUG")

IF(NOT CMAKE_BUILD_TYPE)
    SET(CMAKE_BUILD_TYPE Release CACHE STRING
        "Choose the type of build, options are: \
         None Debug Release MinSizeRel RelWithDebInfo."
    FORCE)
ENDIF()

ADD_EXECUTABLE(qpp ${SOURCE_FILES})

IF(${WITH_MATLAB})
    TARGET_LINK_LIBRARIES(qpp mx mat)
ENDIF()

IF($WITH_OPENMP$ AND ${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang"
        AND CMAKE_CXX_COMPILER_VERSION VERSION_GREATER "3.7")
    TARGET_LINK_LIBRARIES(qpp omp)
ENDIF()
