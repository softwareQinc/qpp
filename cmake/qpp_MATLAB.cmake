# MATLAB support (disabled by default)
option(QPP_MATLAB "Enable MATLAB support" OFF)

if(${QPP_MATLAB})
  message(STATUS "Detecting MATLAB...")
  find_package(Matlab REQUIRED COMPONENTS MX_LIBRARY MAT_LIBRARY)
  if(MATLAB_FOUND)
    set(QPP_MATLAB_LINK_DEPS Matlab::mat Matlab::mx)
    message(STATUS "Detected MATLAB in: ${Matlab_ROOT_DIR}")
    include_directories(SYSTEM ${Matlab_INCLUDE_DIRS})
  else()
    message(FATAL_ERROR "Could not detect MATLAB, aborting")
    message(
      FATAL_ERROR
        "Could not detect MATLAB, aborting build. Please ensure MATLAB is in your system PATH or Matlab_ROOT_DIR is set."
    )
  endif()
endif()
