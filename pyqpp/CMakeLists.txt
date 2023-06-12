cmake_minimum_required(VERSION 3.15)
project(pyqpp)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

#### Eigen3
message(STATUS "Detecting Eigen3")
#### Location specified via an environment variable
set(LOCATION_SET_VIA_ENV FALSE)
if (DEFINED ENV{EIGEN3_INSTALL_DIR})
    set(EIGEN3_INSTALL_DIR_ENV $ENV{EIGEN3_INSTALL_DIR})
    set(LOCATION_SET_VIA_ENV TRUE)
endif ()
#### Location set via CMake variable, trumps all other settings
set(EIGEN3_INSTALL_DIR "" CACHE PATH "Path to Eigen3")
if (NOT ${EIGEN3_INSTALL_DIR} STREQUAL "")
    message(STATUS "Overriding automatic Eigen3 detection (EIGEN3_INSTALL_DIR CMake variable)")
    if (IS_DIRECTORY ${EIGEN3_INSTALL_DIR})
        message(STATUS "Detecting Eigen3 - done (in ${EIGEN3_INSTALL_DIR})")
        include_directories(SYSTEM "${EIGEN3_INSTALL_DIR}")
    else ()
        message(FATAL_ERROR "Invalid path to Eigen3 installation")
    endif ()
#### Location set via environment variable
elseif (LOCATION_SET_VIA_ENV)
    message(STATUS "Overriding automatic Eigen3 detection (EIGEN3_INSTALL_DIR environment variable)")
    if (IS_DIRECTORY ${EIGEN3_INSTALL_DIR_ENV})
        message(STATUS "Detecting Eigen3 - done (in ${EIGEN3_INSTALL_DIR_ENV})")
        include_directories(SYSTEM "${EIGEN3_INSTALL_DIR_ENV}")
    else ()
        message(FATAL_ERROR "Invalid path to Eigen3 installation")
    endif ()
#### Try to find the location automatically
else ()
    find_package(Eigen3 3.0 QUIET NO_MODULE)
    if (NOT TARGET Eigen3::Eigen) # did not find Eigen3 automatically
        message(FATAL_ERROR "Eigen3 not detected! Please point EIGEN3_INSTALL_DIR to your Eigen3 location when building with cmake, for example
    cmake .. -DEIGEN3_INSTALL_DIR=$HOME/eigen3
or set the EIGEN3_INSTALL_DIR environment variable to point to your Eigen3 installation, for example (UNIX/Linux)
    export EIGEN3_INSTALL_DIR=$HOME/eigen3")
    endif ()
    message(STATUS "Detecting Eigen3 - done (in ${EIGEN3_INCLUDE_DIR})")
endif ()
