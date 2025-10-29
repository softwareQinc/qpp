# OpenMP support
option(QPP_OPENMP "Enable OpenMP support" ON)

if(QPP_OPENMP)
  find_package(OpenMP REQUIRED)
  if(OpenMP_CXX_FOUND)
    target_compile_definitions(libqpp INTERFACE QPP_OPENMP)
    target_link_libraries(libqpp INTERFACE OpenMP::OpenMP_CXX)
  else()
    message(
      FATAL_ERROR "OpenMP support requested but OpenMP_CXX was not found.")
  endif()
endif()
