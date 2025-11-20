cmake_minimum_required(VERSION 3.20)

project(benchmarks)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

find_package(qpp REQUIRED)

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

# cmake-lint: disable=C0103
foreach(SOURCE_FILE ${SRC_FILES})
  get_filename_component(TARGET_NAME ${SOURCE_FILE} NAME_WE)
  add_executable(${TARGET_NAME} ${SOURCE_FILE})
  target_link_libraries(${TARGET_NAME} PRIVATE Catch2::Catch2WithMain libqpp)
  message(
    STATUS "Created executable target: ${TARGET_NAME} from ${SOURCE_FILE}")
endforeach()
