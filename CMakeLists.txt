cmake_minimum_required(VERSION 3.12)
project(qpp VERSION 2.6 LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
enable_testing()

#### Guard against in-source builds (snippet from Eigen's CMakeLists.txt)
if (${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
    message(FATAL_ERROR "In-source builds not allowed. Please make a new directory (called a build directory) and run CMake from there. You may need to remove CMakeCache.txt.")
endif ()

add_library(libqpp INTERFACE)
target_include_directories(libqpp INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/>
        $<INSTALL_INTERFACE:include/>)

######## BEGIN LOCAL stuff, you don't need it in your project
#### Dependencies, do not modify unless you know what you're doing
include(cmake/qpp_dependencies.cmake)

#### Quantum++ root directory
add_compile_definitions(QPP_ROOT_DIR="${CMAKE_SOURCE_DIR}")

#### Unit testing
include(cmake/qpp_unit_tests.cmake)

#### Enable all warnings for GNU gcc and Clang/AppleClang
if (${CMAKE_CXX_COMPILER_ID} MATCHES "Clang" OR ${CMAKE_CXX_COMPILER_ID}
        STREQUAL "GNU")
    add_compile_options("-pedantic" "-Wall" "-Wextra" "-Weffc++")
endif ()

#### Examples
include(cmake/qpp_examples.cmake)
######## END LOCAL stuff

#### Quantum++ installation instructions, which can then be found via CMake's
#### find_package(qpp REQUIRED)

#### For both UNIX/Windows, first run from an out-of-source directory
####
#### cmake .. [additional options]
####
#### where [additional options] are optional arguments, one or more of
#### [-DEIGEN3_INSTALL_DIR=/path/to/eigen3] [-DWITH_MATLAB=ON/OFF] [-DWITH_OPENMP=ON/OFF]
####
#### If you want to install Quantum++ in a non-standard location (e.g., due to
#### lack of admin. rights), then pass the additional flag
#### '-DCMAKE_INSTALL_PREFIX=/path/to/qpp' to cmake above, i.e.,
####
#### cmake .. -DCMAKE_INSTALL_PREFIX=/path/to/qpp [additional options]

#### Next, follow the platform-specific instructions below

#### UNIX/Linux: sudo make install (or sudo cmake --build . --target install)
#### Will install in: /usr/local/include/qpp and /usr/local/lib/cmake/qpp

#### Windows: cmake --build . --target INSTALL (from an Admin. Command Prompt)
#### Will install in: C:\Program Files (x86)\qpp

#### Minimal CMakeLists.txt of a standalone application that uses Quantum++
#### Always add the ${QPP_LINK_DEPS} to target_link_libraries()!

# cmake_minimum_required(VERSION 3.12)
# project(my_qpp_project)
# set(CMAKE_CXX_STANDARD 11)
#
# # The line below is optional (Quantum++ installed in a non-standard location)
# set(CMAKE_PREFIX_PATH "/path/to/installed/qpp")
#
# find_package(qpp REQUIRED)
# add_executable(my_qpp_project src/main.cpp)
# target_link_libraries(my_qpp_project PUBLIC ${QPP_LINK_DEPS} libqpp)

set(QPP_INSTALL_DIR "${CMAKE_INSTALL_PREFIX}/include/${PROJECT_NAME}")
install(DIRECTORY include/ DESTINATION ${QPP_INSTALL_DIR})
install(TARGETS libqpp EXPORT qpp_targets)
install(EXPORT qpp_targets DESTINATION "lib/cmake/${PROJECT_NAME}")
include(CMakePackageConfigHelpers)
configure_package_config_file(
        "cmake/qppConfig.cmake.in" "${CMAKE_CURRENT_BINARY_DIR}/qppConfig.cmake"
        INSTALL_DESTINATION "lib/cmake/${PROJECT_NAME}"
)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/qppConfig.cmake" DESTINATION "lib/cmake/${PROJECT_NAME}")
install(FILES "${CMAKE_SOURCE_DIR}/cmake/qpp_dependencies.cmake" DESTINATION "lib/cmake/${PROJECT_NAME}")

#### Uninstall
#### https://gitlab.kitware.com/cmake/community/-/wikis/FAQ#can-i-do-make-uninstall-with-cmake
#### UNIX/Linux: sudo make uninstall (or sudo cmake --build . --target uninstall)
#### Windows:    cmake --build . --target UNINSTALL
if (NOT TARGET uninstall)
    configure_file(
            "${CMAKE_CURRENT_SOURCE_DIR}/cmake/qpp_uninstall.cmake.in"
            "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
            IMMEDIATE @ONLY)
    if (NOT MSVC)
        add_custom_target(uninstall
                COMMAND ${CMAKE_COMMAND} -P "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
                COMMAND ${CMAKE_COMMAND} -E remove_directory "${CMAKE_INSTALL_PREFIX}/lib/cmake/${PROJECT_NAME}"
                COMMAND ${CMAKE_COMMAND} -E remove_directory "${QPP_INSTALL_DIR}"
                )
    else ()
        add_custom_target(uninstall
                COMMAND ${CMAKE_COMMAND} -P "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
                COMMAND ${CMAKE_COMMAND} -E remove_directory "${CMAKE_INSTALL_PREFIX}"
                )
    endif ()
endif ()
