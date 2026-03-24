# MATLAB support

option(QPP_MATLAB "Enable MATLAB support" OFF)
message(STATUS "MATLAB support - ${QPP_MATLAB}")

if(QPP_MATLAB)
  # Select the target
  include(${CMAKE_CURRENT_LIST_DIR}/qpp_select_target.cmake)
  qpp_select_target(QPP_TARGET "MATLAB")

  message(STATUS "Detecting MATLAB...")
  find_package(Matlab REQUIRED COMPONENTS MX_LIBRARY MAT_LIBRARY)
  if(MATLAB_FOUND)
    target_link_libraries(${QPP_TARGET} INTERFACE Matlab::mat Matlab::mx)
  endif()
endif()
