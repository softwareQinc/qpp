cmake_minimum_required(VERSION 3.12)
project(standalone)
set(CMAKE_CXX_STANDARD 17)

#### Disable spurious Eigen warnings if using C++17 with MSVC (warning STL4007)
if (MSVC)
    add_compile_definitions(_SILENCE_CXX17_ADAPTOR_TYPEDEFS_DEPRECATION_WARNING)
endif ()

# If the Quantum++ installation path was non-standard, i.e., specified by
#
# cmake .. -DCMAKE_INSTALL_PREFIX=/path/to/installed/qpp
#
# then uncomment the following line and replace the installation path with yours

# set(CMAKE_PREFIX_PATH "/path/to/installed/qpp")

find_package(qpp REQUIRED)
add_executable(standalone src/main.cpp)
target_link_libraries(standalone PUBLIC ${QPP_LINK_DEPS} libqpp)
