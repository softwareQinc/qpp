cmake_minimum_required(VERSION 3.12)
project(standalone)
set(CMAKE_CXX_STANDARD 11)

find_package(qpp REQUIRED)

add_executable(standalone src/main.cpp)
target_link_libraries(standalone PUBLIC ${QPP_LINK_DEPS} libqpp)
