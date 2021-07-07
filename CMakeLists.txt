cmake_minimum_required(VERSION 3.12)
project(qpp VERSION 2.6 LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
enable_testing()

#### Disable spurious Eigen warnings if using C++17 with MSVC (warning STL4007)
if (MSVC)
    add_compile_definitions(_SILENCE_CXX17_ADAPTOR_TYPEDEFS_DEPRECATION_WARNING)
endif ()

#### Guard against in-source builds (snippet from Eigen's CMakeLists.txt)
if (${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
    message(FATAL_ERROR "In-source builds not allowed. Please make a new directory (called a build directory) and run CMake from there. You may need to remove CMakeCache.txt.")
endif ()

add_library(libqpp INTERFACE)
target_include_directories(libqpp INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/>
        $<INSTALL_INTERFACE:include/>)

######## BEGIN LOCAL stuff, you don't need it in your project
option(WITH_EXAMPLES "Enable example compilation" OFF)
option(WITH_UNIT_TESTS "Enable unit testing" OFF)

#### Dependencies, do not modify unless you know what you're doing
if (${WITH_EXAMPLES} OR ${WITH_UNIT_TESTS})
    include(cmake/qpp_dependencies.cmake)
endif ()

#### Quantum++ root directory
add_compile_definitions(QPP_ROOT_DIR="${CMAKE_SOURCE_DIR}")

#### Unit testing
if (${WITH_UNIT_TESTS})
    include(cmake/qpp_unit_tests.cmake)
endif ()

#### Enable all warnings for GNU gcc and Clang/AppleClang
if (${CMAKE_CXX_COMPILER_ID} MATCHES "Clang" OR ${CMAKE_CXX_COMPILER_ID}
        STREQUAL "GNU")
    add_compile_options("-pedantic" "-Wall" "-Wextra" "-Weffc++")
endif ()

#### Examples
if (${WITH_EXAMPLES})
    include(cmake/qpp_examples.cmake)
endif ()
######## END LOCAL stuff

#[[
Quantum++ installation instructions

For both UNIX/Windows, first configure the system by running from an
out-of-source directory

    cmake .. [optional arguments]

where [optional arguments] are passed with -DOPTIONAL_ARGUMENT=VALUE. In most
cases `cmake ..` should suffice. The optional arguments are:

    DEIGEN3_INSTALL_DIR="/path/to/eigen3"   - Path to Eigen3 installation, if not automatically detected
    WITH_OPENMP=ON/OFF     [ON by default]  - Enables/disables OpenMP
    WITH_EXAMPLES=ON/OFF   [OFF by default] - Enables/disables building examples as a cmake target.
                                              If enabled, adds "examples" as a cmake target.
    WITH_UNIT_TESTS=ON/OFF [OFF by default] - Enables/disables building unit tests as a cmake target.
                                              If enabled, adds "unit_tests" as a cmake target.
    WITH_MATLAB=ON/OFF     [OFF by default] - Enables/disables interoperability with MATLAB, will detect MATLAB
                                              installation automatically. If enabled, allows saving/loading qpp
                                              matrices and vectors to/from MATLAB.
    CMAKE_INSTALL_PREFIX="/path/to/install" - Installs Quantum++ header files in a non-standard location (e.g., due to
                                              lack of admin. rights).

If WITH_MATLAB=ON and the system could not detect your MATLAB installation, you
can manually specify the path to MATLAB's installation directory via the
additional argument

    MATLAB_INSTALL_DIR="/path/to/MATLAB"

If you still receive errors, you can manually specify the path to MATLAB's
required libraries and header files via the additional arguments

    MATLAB_LIB_DIR="/path/to/MATLAB/libs"
    MATLAB_INCLUDE_DIR="/path/to/MATLAB/headers"

Finally, to install the library into a target directory, follow the
platform-specific instructions below.

UNIX/Linux: sudo make install (or sudo cmake --build . --target install)
Will install in: /usr/local/include/qpp and /usr/local/lib/cmake/qpp

Windows: cmake --build . --target INSTALL (from an Admin. Command Prompt)
Will install in: C:\Program Files (x86)\qpp

Below is a minimal CMakeLists.txt of a standalone application that uses
Quantum++. For simplicity, we assume that the whole application is located in a
single file `src/main.cpp'.

Do not forget to ALWAYS add ${QPP_LINK_DEPS} (verbatim) to target_link_libraries()!

cmake_minimum_required(VERSION 3.12)
project(my_qpp_project)
set(CMAKE_CXX_STANDARD 11) # can set the C++ standard to C++17 as well

# Disable spurious Eigen warnings if using C++17 with MSVC (warning STL4007)
if (MSVC)
    add_compile_definitions(_SILENCE_CXX17_ADAPTOR_TYPEDEFS_DEPRECATION_WARNING)
endif ()

# The line below is optional (uncomment if Quantum++ was installed in a
# non-standard location)

# set(CMAKE_PREFIX_PATH "/path/to/installed/qpp")

find_package(qpp REQUIRED)
add_executable(my_qpp_project src/main.cpp)
target_link_libraries(my_qpp_project PUBLIC ${QPP_LINK_DEPS} libqpp)
]]

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
