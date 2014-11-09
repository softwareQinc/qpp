cmake_minimum_required(VERSION 2.7.0)
project(qpp)

set(CMAKE_CXX_COMPILER "g++") # modify as needed

set(EIGEN "$ENV{HOME}/eigen_3.2.2") # modify as needed

option(WITH_MATLAB "MATLAB integration" OFF)
if(${WITH_MATLAB})
    include_directories("/Applications/MATLAB_R2014b.app/extern/include") # modify as needed
    set (CMAKE_EXE_LINKER_FLAGS " -lmx -lmat -L/Applications/MATLAB_R2014b.app/bin/maci64") # modify as needed
endif(${WITH_MATLAB})

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -pedantic -Wall -Wextra -Weffc++ -mtune=native -msse3 -fopenmp")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fopenmp")

set(CMAKE_CXX_FLAGS_MINSIZEREL "${CMAKE_CXX_FLAGS_MINSIZEREL} -DEIGEN_NO_DEBUG")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -DEIGEN_NO_DEBUG")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} -DEIGEN_NO_DEBUG")

include_directories(include)
include_directories(SYSTEM ${EIGEN})

set(SOURCE_FILES
    example.cpp)

if( NOT CMAKE_BUILD_TYPE )
  set( CMAKE_BUILD_TYPE Release CACHE STRING
       "Choose the type of build, options are: None Debug Release MinSizeRel RelWithDebInfo
MinSizeRel."
       FORCE )
endif()

add_executable(qpp ${SOURCE_FILES})
