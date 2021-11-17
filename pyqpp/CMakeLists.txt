cmake_minimum_required(VERSION 3.12)

#### Eigen3
message(STATUS "Detecting Eigen3")
set(EIGEN3_INSTALL_DIR "" CACHE PATH "Path to Eigen3")
#### Location manually specified
if (NOT ${EIGEN3_INSTALL_DIR} STREQUAL "")
    if (IS_DIRECTORY ${EIGEN3_INSTALL_DIR})
        message("Detecting Eigen3 - done in ${EIGEN3_INSTALL_DIR}")
    else ()
        message(FATAL_ERROR "Invalid path to Eigen3 installation")
    endif ()
else () #### Try to find it automatically
    find_package(Eigen3 3.0 QUIET NO_MODULE)
    if (NOT TARGET Eigen3::Eigen) # did not find Eigen3 automatically
        message(FATAL_ERROR
                "Eigen3 not detected! Please point EIGEN3_INSTALL_DIR\
            to your Eigen3 location when building with cmake,\
            for example\

            cmake .. -DEIGEN3_INSTALL_DIR=$HOME/eigen3")
    endif ()
    message("Detecting Eigen3 - done in ${EIGEN3_INCLUDE_DIR}")
endif ()
