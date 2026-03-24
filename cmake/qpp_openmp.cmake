# OpenMP support

option(QPP_OPENMP "Enable OpenMP support" ON)
message(STATUS "OpenMP support - ${QPP_OPENMP}")

if(QPP_OPENMP)
  # Select the target
  include(${CMAKE_CURRENT_LIST_DIR}/qpp_select_target.cmake)
  qpp_select_target(QPP_TARGET "OpenMP")

  find_package(OpenMP 3.0 COMPONENTS CXX)
  if(OpenMP_CXX_FOUND)
    target_compile_definitions(${QPP_TARGET} INTERFACE QPP_OPENMP)
    target_link_libraries(${QPP_TARGET} INTERFACE OpenMP::OpenMP_CXX)
    if(APPLE)
      execute_process(
        COMMAND brew --prefix libomp
        OUTPUT_VARIABLE BREW_LIBOMP_PREFIX
        OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
      if(BREW_LIBOMP_PREFIX)
        message(STATUS "Homebrew libomp found at: ${BREW_LIBOMP_PREFIX}")
        target_include_directories(${QPP_TARGET}
                                   INTERFACE "${BREW_LIBOMP_PREFIX}/include")
      endif()
    endif()
  else()
    message(
      WARNING
        "OpenMP support requested, but OpenMP (>=3.0) not found - building without parallel support"
    )
  endif()
endif()
