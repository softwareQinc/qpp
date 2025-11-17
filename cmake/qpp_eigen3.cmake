# Eigen3

message(STATUS "Detecting Eigen3 (>= 5.0.0)...")

find_package(Eigen3 5.0 QUIET NO_MODULE)

# Check if Eigen3 was found by find_package
if(TARGET Eigen3::Eigen)
  get_target_property(EIGEN3_INCLUDE_DIRS Eigen3::Eigen
                      INTERFACE_INCLUDE_DIRECTORIES)
  # Use the first one found, which should be the main include path
  list(GET EIGEN3_INCLUDE_DIRS 0 EIGEN3_INCLUDE_DIR)
  message(
    STATUS "Found Eigen3 version ${Eigen3_VERSION}: ${EIGEN3_INCLUDE_DIR}")
else() # if(NOT TARGET Eigen3::Eigen)
  # Install Eigen3 if not found by find_package()
  include(FetchContent)
  message(STATUS "Eigen3 (>= 5.0.0) not detected, fetching it...")
  # cmake-lint: disable=C0103
  set(Eigen3_VERSION "5.0.1") # version we want to fetch
  FetchContent_Declare(
    Eigen3
    SYSTEM
    GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
    GIT_TAG ${Eigen3_VERSION}
    GIT_PROGRESS TRUE
    GIT_SHALLOW TRUE
    SOURCE_SUBDIR cmake)
  FetchContent_MakeAvailable(Eigen3)
  # FetchContent does not make Eigen3::Eigen available
  add_library(Eigen3::Eigen INTERFACE IMPORTED)
  # cmake-lint: disable=C0307
  set_target_properties(Eigen3::Eigen PROPERTIES INTERFACE_INCLUDE_DIRECTORIES
                                                 ${eigen3_SOURCE_DIR})
  set(EIGEN3_INCLUDE_DIR ${eigen3_SOURCE_DIR})
  message(
    STATUS "Installed Eigen3 version ${Eigen3_VERSION}: ${EIGEN3_INCLUDE_DIR}")
endif()

target_link_libraries(libqpp INTERFACE Eigen3::Eigen)
