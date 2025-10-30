cmake_minimum_required(VERSION 3.20)
project(standalone)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# If Quantum++'s installation path was non-standard, i.e., specified by
#
# cmake -B build -DCMAKE_INSTALL_PREFIX=/path/to/installed/qpp
#
# uncomment the following line and replace the installation path with yours

# set(CMAKE_PREFIX_PATH "/path/to/installed/qpp")

find_package(qpp REQUIRED)
add_executable(standalone src/main.cpp)
target_link_libraries(standalone PRIVATE libqpp)
