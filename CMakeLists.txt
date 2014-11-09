cmake_minimum_required(VERSION 2.8.4)
project(qppclinon)

SET(CMAKE_CXX_COMPILER "g++")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -fopenmp")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fopenmp")

set(EIGEN "$ENV{HOME}/eigen_3.2.2") # modify as needed

include_directories(include)
include_directories(SYSTEM ${EIGEN})

option(WITH_MATLAB "MATLAB integration" OFF)
if(${WITH_MATLAB})
    include_directories("/Applications/MATLAB_R2014b.app/extern/include") # modify as needed
    set (CMAKE_EXE_LINKER_FLAGS " -lmx -lmat -L/Applications/MATLAB_R2014b.app/bin/maci64")
endif(${WITH_MATLAB})

set(SOURCE_FILES
    example.cpp)

add_executable(qppclinon ${SOURCE_FILES})
