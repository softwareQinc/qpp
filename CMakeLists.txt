CMAKE_MINIMUM_REQUIRED(VERSION 3.0)
PROJECT(qpp VERSION 2.2 LANGUAGES CXX)
SET(CMAKE_CXX_STANDARD 11)
SET(CMAKE_CXX_STANDARD_REQUIRED ON)

#### Source directory
ADD_DEFINITIONS(-DPATH="${CMAKE_SOURCE_DIR}")

#### Quantum++ headers
INCLUDE_DIRECTORIES(include)

#### Source file(s) to be compiled, modify as needed
AUX_SOURCE_DIRECTORY(examples EXAMPLE_FILES)
AUX_SOURCE_DIRECTORY(examples/circuits EXAMPLE_FILES)
AUX_SOURCE_DIRECTORY(examples/qasm EXAMPLE_FILES)

#### Eigen3
MESSAGE(STATUS "Detecting Eigen3")
SET(EIGEN3_INSTALL_DIR "" CACHE PATH "Path to Eigen3")
#### Location manually specified
IF (NOT ${EIGEN3_INSTALL_DIR} STREQUAL "")
    IF (IS_DIRECTORY ${EIGEN3_INSTALL_DIR})
        MESSAGE(STATUS "Detecting Eigen3 - done (in ${EIGEN3_INSTALL_DIR})")
        INCLUDE_DIRECTORIES(SYSTEM "${EIGEN3_INSTALL_DIR}")

    ELSE ()
        MESSAGE(FATAL_ERROR "Invalid path to Eigen3 installation")
    ENDIF ()
ELSE () #### Try to find it automatically
    FIND_PACKAGE(Eigen3 3.0 QUIET NO_MODULE)
    IF (NOT TARGET Eigen3::Eigen) # did not find Eigen3 automatically
        MESSAGE(FATAL_ERROR
                "Eigen3 not detected! Please point EIGEN3_INSTALL_DIR\
            to your Eigen3 location when building with cmake,\
            for example\

            cmake .. -DEIGEN3_INSTALL_DIR=$HOME/eigen")
    ENDIF ()
    MESSAGE(STATUS "Detecting Eigen3 - done (in ${EIGEN3_INCLUDE_DIR})")
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

#### MATLAB support, disabled by default
OPTION(WITH_MATLAB "MATLAB support" OFF)
IF (${WITH_MATLAB})
    MESSAGE(STATUS "Detecting MATLAB")
    #### Try to find it automatically
    FIND_PACKAGE(Matlab OPTIONAL_COMPONENTS MX_LIBRARY MAT_LIBRARY QUIET)
    IF (MATLAB_FOUND)
        MESSAGE(STATUS "Detecting MATLAB - done (in ${Matlab_ROOT_DIR})")
        INCLUDE_DIRECTORIES(SYSTEM ${Matlab_INCLUDE_DIRS})
        IF (WIN32)
            IF (MSVC)
                SET(MATLAB_LIB_DIR
                        "${Matlab_ROOT_DIR}/extern/lib/win64/microsoft"
                        CACHE PATH "Custom path to MATLAB lib directory")
            ELSEIF (MINGW)
                SET(MATLAB_LIB_DIR "${Matlab_ROOT_DIR}/extern/lib/win64/mingw64"
                        CACHE PATH "Custom path to MATLAB lib directory")
            ELSE ()
                MESSAGE(FATAL "Platform not supported, aborting.")
            ENDIF ()
        ELSEIF (UNIX AND NOT APPLE)
            SET(MATLAB_LIB_DIR "${Matlab_ROOT_DIR}/bin/glnxa64" CACHE
                    PATH "Custom path to MATLAB lib directory")
        ELSEIF (APPLE)
            SET(MATLAB_LIB_DIR "${Matlab_ROOT_DIR}/bin/maci64" CACHE
                    PATH "Custom path to MATLAB lib directory")
        ELSE ()
            MESSAGE(FATAL "Platform not supported, aborting.")
        ENDIF ()
        LINK_DIRECTORIES(${MATLAB_LIB_DIR})
        ADD_DEFINITIONS(-DHAS_MATLAB_ENABLED)
        SET(BUILD_WITH_MATLAB TRUE)
    ELSE () #### Location manually specified
        SET(MATLAB_INSTALL_DIR "" CACHE PATH
                "Custom path to MATLAB installation")
        IF (IS_DIRECTORY ${MATLAB_INSTALL_DIR})
            # MATLAB include files
            SET(MATLAB_INCLUDE_DIR "${MATLAB_INSTALL_DIR}/extern/include"
                    CACHE PATH "Custom path to MATLAB include directory")
            IF (IS_DIRECTORY ${MATLAB_INCLUDE_DIR})
                INCLUDE_DIRECTORIES(SYSTEM ${MATLAB_INCLUDE_DIR})
            ELSE ()
                MESSAGE(FATAL_ERROR
                        "Possibly corrupted MATLAB include headers")
            ENDIF ()
            # MATLAB linker files
            IF (WIN32)
                IF (MSVC)
                    SET(MATLAB_LIB_DIR
                            "${MATLAB_INSTALL_DIR}/extern/lib/win64/microsoft"
                            CACHE PATH "Custom path to MATLAB lib directory")
                ELSEIF (MINGW64)
                    SET(MATLAB_LIB_DIR
                            "${MATLAB_INSTALL_DIR}/extern/lib/win64/mingw64"
                            CACHE PATH "Custom path to MATLAB lib directory")
                ELSE ()
                    MESSAGE(FATAL "Platform not supported, aborting.")
                ENDIF ()
            ELSEIF (UNIX AND NOT APPLE)
                SET(MATLAB_LIB_DIR "${MATLAB_INSTALL_DIR}/bin/glnxa64" CACHE
                        PATH "Custom path to MATLAB lib directory")
            ELSEIF (APPLE)
                SET(MATLAB_LIB_DIR "${MATLAB_INSTALL_DIR}/bin/maci64" CACHE
                        PATH "Custom path to MATLAB lib directory")
            ELSE ()
                MESSAGE(FATAL "Platform not supported, aborting.")
            ENDIF ()
            IF (IS_DIRECTORY ${MATLAB_LIB_DIR})
                LINK_DIRECTORIES(${MATLAB_LIB_DIR})
            ELSE ()
                MESSAGE(FATAL_ERROR
                        "Possibly corrupted MATLAB compiler libraries")
            ENDIF ()
            #### Everything is OK, inject definition (as #define) in the source
            MESSAGE(STATUS
                    "Detecting MATLAB - done (in ${MATLAB_INSTALL_DIR})")
            ADD_DEFINITIONS(-DHAS_MATLAB_ENABLED)
            SET(BUILD_WITH_MATLAB TRUE)
        ELSE ()
            MESSAGE(FATAL_ERROR "Could not detect MATLAB, aborting")
        ENDIF ()
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
        ADD_DEFINITIONS(-DHAS_OPENMP)
    ELSE ()
        MESSAGE(WARNING "Detected compiler: ${CMAKE_CXX_COMPILER_ID} \
        ${CMAKE_CXX_COMPILER_VERSION}. Support for OpenMP is enabled only for \
        the GNU gcc compiler or the clang compiler version 3.8 or higher.")
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
    #### Do not build "examples/matlab_io.cpp" if there's no MATLAB support
    IF (${FILE} STREQUAL "examples/matlab_io.cpp" AND NOT BUILD_WITH_MATLAB)
        CONTINUE()
    ENDIF ()
    ADD_EXECUTABLE(${TARGET_NAME} ${FILE})
    IF (${BUILD_WITH_MATLAB})
        IF (WIN32)
            IF (MSVC)
                TARGET_LINK_LIBRARIES(${TARGET_NAME} libmx libmat)
            ELSEIF (MINGW)
                TARGET_LINK_LIBRARIES(${TARGET_NAME} mx mat)
            ENDIF ()
        ELSE ()
            TARGET_LINK_LIBRARIES(${TARGET_NAME} mx mat)
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
