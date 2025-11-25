# Benchmarks with Catch2

include(FetchContent)
message(STATUS "Fetching Catch2...")
FetchContent_Declare(
  Catch2
  GIT_REPOSITORY https://github.com/catchorg/Catch2.git
  GIT_TAG v3.11.0)
FetchContent_MakeAvailable(Catch2)

add_library(benchmarks_interface INTERFACE)
target_include_directories(
  benchmarks_interface
  INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)
target_link_libraries(benchmarks_interface INTERFACE Catch2::Catch2WithMain
                                                     libqpp_internal)

set(SRC_DIR "src")
set(SRC_FILES "")
file(GLOB CURRENT_FILES "${SRC_DIR}/*.cpp" "${SRC_DIR}/*.cxx" "${SRC_DIR}/*.cc"
     "${SRC_DIR}/*.c++")
list(APPEND SRC_FILES ${CURRENT_FILES})

if(NOT SRC_FILES)
  message(WARNING "No C++ source files found in the ${SRC_DIR} directory.")
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
  target_link_libraries(${TARGET_NAME} PRIVATE benchmarks_interface)
  list(APPEND BENCHMARK_TARGETS ${TARGET_NAME})
endforeach()

if(TARGET benchmarks)
  add_dependencies(benchmarks ${BENCHMARK_TARGETS})
endif()
