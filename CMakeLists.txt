cmake_minimum_required(VERSION 3.15)
set(QPP_VERSION_NUM 4.3.3)
set(QPP_VERSION_STR "${QPP_VERSION_NUM}")
project(qpp VERSION ${QPP_VERSION_NUM} LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
enable_testing()

#### Quantum++ version
add_compile_definitions(QPP_VERSION_NUM=${QPP_VERSION_NUM})
add_compile_definitions(QPP_VERSION_STR="${QPP_VERSION_STR}")

#### Quantum++ root directory
add_compile_definitions(PROJECT_ROOT_DIR="${PROJECT_SOURCE_DIR}")

#### Guard against in-source builds (snippet from Eigen's CMakeLists.txt)
if (${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
    message(FATAL_ERROR "In-source builds not allowed. Please instruct CMake to use an out-of-source build, e.g.,
    cmake --build build [...]
You may need to remove CMakeCache.txt.")
endif ()

add_library(libqpp INTERFACE)
target_include_directories(libqpp INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/>
        $<INSTALL_INTERFACE:include/>)
target_compile_definitions(libqpp INTERFACE -DQPP_VERSION_NUM=${QPP_VERSION_NUM})
target_compile_definitions(libqpp INTERFACE -DQPP_VERSION_STR="${QPP_VERSION_STR}")

#### qasmtools library
target_include_directories(libqpp INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/qasmtools/include/>
        $<INSTALL_INTERFACE:include/>)

######## BEGIN LOCAL stuff, you don't need those in standalone projects
option(WITH_EXAMPLES "Enable example compilation" OFF)
option(WITH_UNIT_TESTS "Enable unit testing" OFF)
option(SANITIZE "Enable code sanitizing (only for gcc/clang)" OFF)

#### Dependencies, do not modify unless you know what you're doing
if (${WITH_EXAMPLES} OR ${WITH_UNIT_TESTS})
    include(cmake/qpp_dependencies.cmake)
endif ()

#### Unit testing
if (${WITH_UNIT_TESTS})
    include(cmake/qpp_unit_tests.cmake)
endif ()

#### Enable all warnings for GNU gcc and Clang/AppleClang
if (${CMAKE_CXX_COMPILER_ID} MATCHES "Clang" OR ${CMAKE_CXX_COMPILER_ID}
        STREQUAL "GNU")
    add_compile_options("-pedantic" "-Wall" "-Wextra" "-Weffc++")
    if (${SANITIZE})
        if (NOT (${CMAKE_CXX_COMPILER_ID} MATCHES "AppleClang"))
            list(APPEND SANITIZE_OPTIONS -fsanitize=undefined)
            add_compile_options("${SANITIZE_OPTIONS}")
            set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${SANITIZE_OPTIONS}")
        endif ()
    endif ()
endif ()

#### Examples
if (${WITH_EXAMPLES})
    include(cmake/qpp_examples.cmake)
endif ()
######## END LOCAL stuff

set(QPP_INSTALL_DIR "${CMAKE_INSTALL_PREFIX}/include/${PROJECT_NAME}")
install(DIRECTORY include/ DESTINATION ${QPP_INSTALL_DIR})
install(DIRECTORY qasmtools/include/ DESTINATION ${QPP_INSTALL_DIR})
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
#### UNIX/Linux: sudo cmake --build build --target uninstall
#### Windows:    cmake --build build --target UNINSTALL
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
