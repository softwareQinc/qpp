message(STATUS "Detecting Eigen3...")
find_package(Eigen3 5.0 QUIET NO_MODULE)

if(NOT TARGET Eigen3::Eigen)
  # Install Eigen3 if not found by find_package()
  include(FetchContent)
  message(STATUS "Eigen3 not detected, fetching Eigen3...")
  FetchContent_Declare(
    Eigen3
    SYSTEM
    GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
    GIT_TAG 5.0.0
    GIT_SHALLOW TRUE
    SOURCE_SUBDIR cmake)
  FetchContent_MakeAvailable(Eigen3)
endif()

# In case FetchContent does not make Eigen3::Eigen available
if(NOT TARGET Eigen3::Eigen)
  add_library(Eigen3::Eigen INTERFACE IMPORTED)
  set_target_properties(Eigen3::Eigen PROPERTIES INTERFACE_INCLUDE_DIRECTORIES
                                                 ${eigen3_SOURCE_DIR})
  set(EIGEN3_INCLUDE_DIR ${eigen3_SOURCE_DIR})
endif()

message(STATUS "Detected Eigen3 in: ${EIGEN3_INCLUDE_DIR}")
set(QPP_EIGEN3_LINK_DEPS Eigen3::Eigen)
