# Examples Source file(s) to be compiled, modify as needed
aux_source_directory(${CMAKE_SOURCE_DIR}/examples EXAMPLE_FILES)
aux_source_directory(${CMAKE_SOURCE_DIR}/examples/circuits EXAMPLE_FILES)
aux_source_directory(${CMAKE_SOURCE_DIR}/examples/qasm EXAMPLE_FILES)

# Build all examples in ${EXAMPLE_FILES}
add_custom_target(examples COMMENT "Examples")
foreach(file ${EXAMPLE_FILES})
  get_filename_component(TARGET_NAME ${file} NAME_WE)
  # Do not build "examples/matlab_io.cpp" if there's no MATLAB support
  if(${TARGET_NAME} STREQUAL "matlab_io" AND NOT MATLAB_FOUND)
    continue()
  endif()
  add_executable(${TARGET_NAME} EXCLUDE_FROM_ALL ${file})
  add_dependencies(examples ${TARGET_NAME})
  target_link_libraries(${TARGET_NAME} PUBLIC ${QPP_LINK_DEPS} libqpp)
endforeach()
