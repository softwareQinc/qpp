message(STATUS "Detecting Eigen3...")
find_package(Eigen3 QUIET NO_MODULE)
if(NOT TARGET Eigen3::Eigen)
  # Install Eigen3 on demand
  include(FetchContent)
  set(CMAKE_POLICY_DEFAULT_CMP0077 NEW)
  message(STATUS "Eigen3 not detected, fetching Eigen3...")
  FetchContent_Declare(
    Eigen3
    SYSTEM
    GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
    GIT_TAG 3.4.0 # 3.4.0
    GIT_SHALLOW TRUE
    # no CMakeLists.txt in cmake, so this turns off configure. Recommend also
    # to add `FIND_PACKAGE_ARGS CONFIG` so that FetchContent checks to see if
    # Eigen is installed on the system, via the OS, or a package manager
    SOURCE_SUBDIR cmake)
  FetchContent_MakeAvailable(Eigen3)
endif()

if(NOT TARGET Eigen3::Eigen)
  add_library(Eigen3::Eigen INTERFACE IMPORTED)
  set_target_properties(Eigen3::Eigen PROPERTIES INTERFACE_INCLUDE_DIRECTORIES
                                                 ${eigen3_SOURCE_DIR})
  set(EIGEN3_INCLUDE_DIR ${eigen3_SOURCE_DIR})
endif()

message(STATUS "Detected Eigen3 in: ${EIGEN3_INCLUDE_DIR}")
set(QPP_EIGEN3_LINK_DEPS Eigen3::Eigen)
