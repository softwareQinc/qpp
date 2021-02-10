cmake_minimum_required(VERSION 3.12)
project(standalone)
set(CMAKE_CXX_STANDARD 17)

# If Quantum++ was installed in a non-standard location via
#
# 'cmake .. -DCMAKE_INSTALL_PREFIX=/path/to/qpp'
#
# then uncomment the following line and replace the installation path with yours

# set(CMAKE_PREFIX_PATH "/path/to/qpp")

find_package(qpp REQUIRED)
add_compile_definitions(_SILENCE_CXX17_ADAPTOR_TYPEDEFS_DEPRECATION_WARNING)
add_executable(standalone src/main.cpp)
target_link_libraries(standalone PUBLIC ${QPP_LINK_DEPS} libqpp)
