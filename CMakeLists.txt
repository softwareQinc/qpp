cmake_minimum_required(VERSION 3.12)
project(qpp VERSION 2.6 LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
enable_testing()

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

#### Package installation, can then be found via CMake's find_package(qpp)

#### UNIX/Linux: sudo make install (or sudo cmake --build . --target install)
#### Destination(s): /usr/local/include/qpp && /usr/local/lib/cmake/qpp

#### Windows: cmake --build . --target INSTALL (from an Admin. Command Prompt)
#### Destination: C:\Program Files (x86)\qpp

#### Minimal example, always add the ${QPP_LINK_DEPS} to target_link_libraries()

# cmake_minimum_required(VERSION 3.9)
# project(test)
# set(CMAKE_CXX_STANDARD 11)
# find_package(qpp REQUIRED NO_MODULE)
# add_executable(test src/main.cpp)
# target_link_libraries(test PUBLIC ${QPP_LINK_DEPS} libqpp)

set(QPP_INSTALL_DIR "${CMAKE_INSTALL_PREFIX}/include/${PROJECT_NAME}")
install(DIRECTORY include/ DESTINATION ${QPP_INSTALL_DIR})
install(TARGETS libqpp EXPORT qpp_targets)
install(EXPORT qpp_targets DESTINATION lib/cmake/qpp)
include(CMakePackageConfigHelpers)
configure_package_config_file(
        "cmake/qppConfig.cmake.in" "${CMAKE_CURRENT_BINARY_DIR}/qppConfig.cmake"
        INSTALL_DESTINATION "lib/cmake/qpp"
)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/qppConfig.cmake" DESTINATION "lib/cmake/qpp")
install(FILES "${CMAKE_SOURCE_DIR}/cmake/qpp_dependencies.cmake" DESTINATION "lib/cmake/qpp")

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
