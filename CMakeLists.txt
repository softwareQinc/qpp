CMAKE_MINIMUM_REQUIRED(VERSION 2.8.4)
PROJECT(qpp)

#### CMake support only for g++ and clang++
IF (NOT ${CMAKE_CXX_COMPILER_ID} MATCHES "Clang" AND NOT ${CMAKE_CXX_COMPILER_ID} MATCHES "GNU")
    MESSAGE(FATAL_ERROR "This CMakeLists.txt works only with g++/clang++")
ENDIF ()

#### Quantum++ headers
INCLUDE_DIRECTORIES(include)

#### Source file(s) to be compiled, modify as needed
SET(SOURCE_FILES
        examples/minimal.cpp)

#### Eigen3, modify as needed
SET(EIGEN3_INCLUDE_DIR "$ENV{HOME}/eigen")

#### Requires Eigen3
IF (NOT EIGEN3_INCLUDE_DIR)
    MESSAGE(FATAL_ERROR "Please point the environment variable EIGEN3_INCLUDE_DIR to the include directory of your Eigen3 installation.")
ENDIF ()
INCLUDE_DIRECTORIES(SYSTEM "${EIGEN3_INCLUDE_DIR}")

#### Compiler flags
SET(ADDITIONAL_FLAGS "-std=c++11 -pedantic -Wall -Wextra -Weffc++")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${ADDITIONAL_FLAGS}")
SET(CMAKE_INCLUDE_SYSTEM_FLAG_CXX "-isystem")

#### OpenMP support
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

#### MATLAB integration
OPTION(WITH_MATLAB "MATLAB integration" OFF)
IF (${WITH_MATLAB})
    #### modify as needed 
    SET(MATLAB "/Applications/MATLAB_R2015b.app")
    INCLUDE_DIRECTORIES(SYSTEM "${MATLAB}/extern/include")
    SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -lmx -lmat -L${MATLAB}/bin/maci64")
ENDIF ()

#### Disable support for thread_local storage when using clang++
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
        SET(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fno-weak")
    ENDIF ()
    SET(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Og -D_GLIBCXX_DEBUG")
ENDIF ()

#### Configurations
SET(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DDEBUG")
SET(CMAKE_CXX_FLAGS_MINSIZEREL "${CMAKE_CXX_FLAGS_MINSIZEREL} -DEIGEN_NO_DEBUG")
SET(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -DEIGEN_NO_DEBUG")
SET(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} -DEIGEN_NO_DEBUG")

IF (NOT CMAKE_BUILD_TYPE)
    SET(CMAKE_BUILD_TYPE Release CACHE STRING
            "Choose the type of build, options are: None Debug Release MinSizeRel RelWithDebInfo."
            FORCE)
ENDIF ()

ADD_EXECUTABLE(qpp ${SOURCE_FILES})
