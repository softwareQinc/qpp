cmake_minimum_required(VERSION 3.15)
project(standalone)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

find_package(qpp REQUIRED)

enable_testing()

include(FetchContent)
FetchContent_Declare(
  Catch2
  GIT_REPOSITORY https://github.com/catchorg/Catch2.git
  GIT_TAG v3.8.0)
FetchContent_MakeAvailable(Catch2)

file(GLOB SRC_FILES src/*.cpp)
add_executable(bench ${SRC_FILES})

target_link_libraries(bench PRIVATE Catch2::Catch2WithMain)
target_link_libraries(bench PUBLIC ${QPP_LINK_DEPS} libqpp)
