message(STATUS "Detecting Eigen3...")
find_package(Eigen3 3.0 QUIET NO_MODULE)
if(NOT TARGET Eigen3::Eigen)
  # Install Eigen3 on demand
  include(FetchContent)
  set(FETCHCONTENT_QUIET FALSE)
  set(CMAKE_POLICY_DEFAULT_CMP0077 NEW)
  message(STATUS "Fetching Eigen3...")
  FetchContent_Declare(
    Eigen
    GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
    GIT_TAG 3.4.0
    GIT_SHALLOW TRUE
    GIT_PROGRESS TRUE)
  # note: To disable eigen tests, you should put this code in a add_subdirectory
  # to avoid to change BUILD_TESTING for your own project too since variables
  # are directory scoped
  set(BUILD_TESTING OFF)
  set(EIGEN_BUILD_TESTING OFF)
  set(EIGEN_MPL2_ONLY ON)
  set(EIGEN_BUILD_PKGCONFIG OFF)
  set(EIGEN_BUILD_DOC OFF)
  FetchContent_MakeAvailable(Eigen)
  set(EIGEN3_INCLUDE_DIR ${eigen_SOURCE_DIR})
endif()
set(QPP_EIGEN3_LINK_DEPS Eigen3::Eigen)
message(STATUS "Detected Eigen3 in: ${EIGEN3_INCLUDE_DIR}")
