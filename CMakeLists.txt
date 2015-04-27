cmake_minimum_required(VERSION 2.8.4)
project(qpp)

# modify as needed
set(SOURCE_FILES
    examples/example.cpp)

# modify as needed
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -pedantic -Wall -Wextra -Weffc++")
set(EIGEN "$ENV{HOME}/eigen")

option(WITH_OPENMP "OpenMP support" ON)
if(${WITH_OPENMP})
    # no OpenMP support for clang++ yet
    if(NOT ${CMAKE_CXX_COMPILER_ID} MATCHES "Clang")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
	      set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fopenmp")
    endif()
endif(${WITH_OPENMP})

# modify as needed
option(WITH_MATLAB "MATLAB integration" OFF)
if(${WITH_MATLAB})
    set(MATLAB "/Applications/MATLAB_R2014b.app")
    include_directories(SYSTEM "${MATLAB}/extern/include")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -lmx -lmat -L${MATLAB}/bin/maci64")
endif()

# use the "no-weak" debugging flag only when debugging under OS X,
# as gdb cannot step in template functions when debugging code produced by g++
# see http://stackoverflow.com/questions/23330641/gnu-gdb-can-not-step-into-template-functions-os-x-mavericks
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    if(${CMAKE_CXX_COMPILER_ID} MATCHES "GNU")
        set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fno-weak")
        set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Og")
    endif()
endif()

# disable support for thread_local storage when using clang++
# as libc++ doesn't yet support it
if (${CMAKE_CXX_COMPILER_ID} MATCHES "Clang")
    add_definitions(-D_NO_THREAD_LOCAL)
endif()

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -D_GLIBCXX_DEBUG")
set(CMAKE_CXX_FLAGS_MINSIZEREL "${CMAKE_CXX_FLAGS_MINSIZEREL} -DEIGEN_NO_DEBUG")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -DEIGEN_NO_DEBUG")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} -DEIGEN_NO_DEBUG")

include_directories(include)
include_directories(SYSTEM ${EIGEN})

if( NOT CMAKE_BUILD_TYPE )
  set( CMAKE_BUILD_TYPE Release CACHE STRING
       "Choose the type of build, options are: None Debug Release MinSizeRel RelWithDebInfo."
       FORCE )
endif()

add_executable(qpp ${SOURCE_FILES})
