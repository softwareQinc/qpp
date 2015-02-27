cmake_minimum_required(VERSION 2.8.0)
project(qpp)

# modify as needed
set(SOURCE_FILES
    test.cpp)

# modify as needed
set(CMAKE_CXX_COMPILER "g++")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -pedantic -Wall -Wextra -Weffc++")
set(EIGEN "$ENV{HOME}/eigen")
set(GMOCK "$ENV{HOME}/gmock-1.7.0")

option(WITH_OPENMP "OpenMP support" ON)
if(${WITH_OPENMP})
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
	set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fopenmp")
endif(${WITH_OPENMP})

# modify as needed
option(WITH_MATLAB "MATLAB integration" OFF)
if(${WITH_MATLAB})
    include_directories("/Applications/MATLAB_R2014b.app/extern/include")
    set (CMAKE_EXE_LINKER_FLAGS " -lmx -lmat -L/Applications/MATLAB_R2014b.app/bin/maci64")
endif(${WITH_MATLAB})

# use the "no-weak" debugging flag only when debugging under OS X,
# as gdb cannot step in template functions when debugging code produced by g++
# see http://stackoverflow.com/questions/23330641/gnu-gdb-can-not-step-into-template-functions-os-x-mavericks
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fno-weak")
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

# set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Og") # only for g++>=4.8
set(CMAKE_CXX_FLAGS_MINSIZEREL "${CMAKE_CXX_FLAGS_MINSIZEREL} -DEIGEN_NO_DEBUG")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -DEIGEN_NO_DEBUG")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} -DEIGEN_NO_DEBUG")

include_directories(include)
include_directories(SYSTEM ${EIGEN} ${GMOCK})
set (CMAKE_EXE_LINKER_FLAGS " -lgmock -L/Users/vlad/gmock-1.7.0/build/ ")



if( NOT CMAKE_BUILD_TYPE )
  set( CMAKE_BUILD_TYPE Release CACHE STRING
       "Choose the type of build, options are: None Debug Release MinSizeRel RelWithDebInfo MinSizeRel."
       FORCE )
endif()

add_executable(qpp ${SOURCE_FILES})
