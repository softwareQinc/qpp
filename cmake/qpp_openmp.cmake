# OpenMP support
option(QPP_OPENMP "Enable OpenMP support" ON)

if(QPP_OPENMP)
  message(STATUS "OpenMP support - ON")

  # Determine which target exists
  include(${CMAKE_CURRENT_LIST_DIR}/qpp_detect_target.cmake)
  qpp_detect_target(QPP_TARGET "OpenMP")

  find_package(OpenMP 3.0 COMPONENTS CXX)
  if(OpenMP_CXX_FOUND)
    target_compile_definitions(${QPP_TARGET} INTERFACE QPP_OPENMP)
    target_link_libraries(${QPP_TARGET} INTERFACE OpenMP::OpenMP_CXX)
  else()
    message(
      WARNING
        "OpenMP support requested, but OpenMP (>=3.0) not found — building without parallel support."
    )
  endif()
else()
  message(STATUS "OpenMP support - OFF")
endif()
