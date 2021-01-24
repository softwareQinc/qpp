cmake_minimum_required(VERSION 3.12)
project(test)
set(CMAKE_CXX_STANDARD 11)

find_package(qpp REQUIRED)

add_executable(test src/main.cpp)
target_link_libraries(test PUBLIC ${QPP_LINK_DEPS} libqpp)
