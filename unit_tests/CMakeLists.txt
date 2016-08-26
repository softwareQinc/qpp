CMAKE_MINIMUM_REQUIRED(VERSION 2.8.4)
PROJECT(testing)

#### supports only g++ and clang++
IF (NOT ${CMAKE_CXX_COMPILER_ID} MATCHES "Clang" AND NOT ${CMAKE_CXX_COMPILER_ID} MATCHES "GNU")
    MESSAGE(FATAL_ERROR "This CMakeLists.txt works only with g++/clang++")
ENDIF ()

#### modify as needed
SET(SOURCE_FILES
        testing.cpp)

#### modify as needed
SET(EIGEN "$ENV{HOME}/eigen")
SET(GMOCK "$ENV{HOME}/gmock-1.7.0")
INCLUDE_DIRECTORIES(../include)
INCLUDE_DIRECTORIES(SYSTEM ${EIGEN} "${GMOCK}/gtest/include")

#### modify as needed
SET(ADDITIONAL_FLAGS "-std=c++11 -pedantic -Wall -Wextra -Weffc++")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${ADDITIONAL_FLAGS}")
SET(CMAKE_INCLUDE_SYSTEM_FLAG_CXX "-isystem")
SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -L${GMOCK}/build/ ")

OPTION(WITH_OPENMP "OpenMP support" ON)
IF (${WITH_OPENMP})
    IF (${CMAKE_CXX_COMPILER_ID} MATCHES "GNU")
        SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
        SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fopenmp")
        #### no OpenMP support for clang++ yet
    ELSEIF (NOT ${CMAKE_CXX_COMPILER_ID} MATCHES "GNU")
        MESSAGE(WARNING "Detected compiler: ${CMAKE_CXX_COMPILER_ID}. Support for OpenMP is enabled only for the GNU gcc compiler.")
    ENDIF ()
ENDIF ()

OPTION(WITH_MATLAB "MATLAB integration" OFF)
IF (${WITH_MATLAB})
    #### modify as needed 
    SET(MATLAB "/Applications/MATLAB_R2015b.app")
    INCLUDE_DIRECTORIES(SYSTEM "${MATLAB}/extern/include")
    SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -lmx -lmat -L${MATLAB}/bin/maci64")
ENDIF ()

#### disable support for thread_local storage when using clang++
#### as libc++ doesn't yet support it
IF (${CMAKE_CXX_COMPILER_ID} MATCHES "Clang")
    ADD_DEFINITIONS(-D_NO_THREAD_LOCAL)
ENDIF ()

#### g++ additional debug settings
IF (${CMAKE_CXX_COMPILER_ID} MATCHES "GNU")
    IF (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
        ## use the "no-weak" debugging flag only when debugging under OS X,
        ## as gdb cannot step in template functions when debugging code produced by g++
        ## see http://stackoverflow.com/questions/23330641/gnu-gdb-can-not-step-into-template-functions-os-x-mavericks
        set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fno-weak")
    ENDIF ()
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Og -D_GLIBCXX_DEBUG")
ENDIF ()

SET(CMAKE_CXX_FLAGS_MINSIZEREL "${CMAKE_CXX_FLAGS_MINSIZEREL} -DEIGEN_NO_DEBUG")
SET(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -DEIGEN_NO_DEBUG")
SET(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} -DEIGEN_NO_DEBUG")

IF (NOT CMAKE_BUILD_TYPE)
    SET(CMAKE_BUILD_TYPE Release CACHE STRING
            "Choose the type of build, options are: None Debug Release MinSizeRel RelWithDebInfo."
            FORCE)
ENDIF ()

ADD_EXECUTABLE(testing ${SOURCE_FILES})
# make sure -lgmock and -lpthread are added at the end of the compiling command
TARGET_LINK_LIBRARIES(testing gmock pthread) 
