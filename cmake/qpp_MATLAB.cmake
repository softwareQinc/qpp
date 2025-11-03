# MATLAB support (disabled by default)
option(QPP_MATLAB "Enable MATLAB support" OFF)

if(QPP_MATLAB)
  message(STATUS "MATLAB support - ON")

  # Determine which target exists
  include(${CMAKE_CURRENT_LIST_DIR}/qpp_detect_target.cmake)
  qpp_detect_target(QPP_TARGET "MATLAB")

  message(STATUS "Detecting MATLAB...")

  find_package(Matlab REQUIRED COMPONENTS MX_LIBRARY MAT_LIBRARY)
  if(MATLAB_FOUND)
    target_link_libraries(${QPP_TARGET} INTERFACE Matlab::mat Matlab::mx)
  endif()
else()
  message(STATUS "MATLAB support - OFF")
endif()
