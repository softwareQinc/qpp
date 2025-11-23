# Examples source file(s) to be compiled, modify as needed

set(EXAMPLE_DIRS
    ${CMAKE_SOURCE_DIR}/examples ${CMAKE_SOURCE_DIR}/examples/circuits
    ${CMAKE_SOURCE_DIR}/examples/qasm)
set(EXAMPLE_FILES "")
foreach(dir ${EXAMPLE_DIRS})
  file(GLOB CURRENT_FILES "${dir}/*.cpp" "${dir}/*.cxx" "${dir}/*.cc"
       "${dir}/*.c++")
  list(APPEND EXAMPLE_FILES ${CURRENT_FILES})
endforeach()

# Build all examples in ${EXAMPLE_FILES}
foreach(file ${EXAMPLE_FILES})
  get_filename_component(TARGET_NAME ${file} NAME_WE)
  # Do not build "examples/matlab_io.cpp" if there's no MATLAB support
  if(${TARGET_NAME} STREQUAL "matlab_io" AND NOT MATLAB_FOUND)
    message(STATUS "Skipping example: ${file} (MATLAB support not enabled)")
    continue()
  endif()
  add_executable(${TARGET_NAME} EXCLUDE_FROM_ALL ${file})
  set_target_properties(
    ${TARGET_NAME}
    PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/examples"
               RUNTIME_OUTPUT_DIRECTORY_DEBUG "${CMAKE_BINARY_DIR}/examples"
               RUNTIME_OUTPUT_DIRECTORY_RELEASE "${CMAKE_BINARY_DIR}/examples")
  add_dependencies(examples ${TARGET_NAME})
  target_link_libraries(${TARGET_NAME} PRIVATE libqpp_internal)
endforeach()
