include_directories(${CMAKE_SOURCE_DIR}/benchmarks/include)

include(FetchContent)
message(STATUS "Fetching Catch2...")
FetchContent_Declare(
  Catch2
  GIT_REPOSITORY https://github.com/catchorg/Catch2.git
  GIT_TAG v3.11.0)
FetchContent_MakeAvailable(Catch2)

set(EXECUTABLE_DIR "src")
set(SRC_FILES "")
file(GLOB CURRENT_FILES "${EXECUTABLE_DIR}/*.cpp" "${EXECUTABLE_DIR}/*.cxx"
     "${EXECUTABLE_DIR}/*.cc" "${EXECUTABLE_DIR}/*.c++")
list(APPEND SRC_FILES ${CURRENT_FILES})

if(NOT SRC_FILES)
  message(
    WARNING "No C++ source files found in the ${EXECUTABLE_DIR} directory.")
endif()

set(BENCHMARK_TARGETS "")

# cmake-lint: disable=C0103
foreach(SOURCE_FILE ${SRC_FILES})
  get_filename_component(TARGET_NAME ${SOURCE_FILE} NAME_WE)
  add_executable(${TARGET_NAME} ${SOURCE_FILE})
  set_target_properties(
    ${TARGET_NAME}
    PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/benchmarks"
               RUNTIME_OUTPUT_DIRECTORY_DEBUG "${CMAKE_BINARY_DIR}/benchmarks"
               RUNTIME_OUTPUT_DIRECTORY_RELEASE
               "${CMAKE_BINARY_DIR}/benchmarks")
  target_link_libraries(${TARGET_NAME} PRIVATE Catch2::Catch2WithMain libqpp)
  list(APPEND BENCHMARK_TARGETS ${TARGET_NAME})
endforeach()

if(TARGET benchmarks)
  add_dependencies(benchmarks ${BENCHMARK_TARGETS})
endif()
