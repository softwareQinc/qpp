include(GoogleTest)
set(TARGET_NAME "qpp_testing")

#### Link Google Test dynamically if using MSVC
if (MSVC)
    set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
endif ()

add_subdirectory(lib/googletest-release-1.10.0)

aux_source_directory(tests TEST_FILES)
aux_source_directory(tests/MATLAB TEST_FILES)
aux_source_directory(tests/classes TEST_FILES)
aux_source_directory(tests/classes/circuits TEST_FILES)
aux_source_directory(tests/qasm TEST_FILES)

add_executable(${TARGET_NAME} tests/main.cpp)

if (NOT ${CMAKE_VERSION} VERSION_LESS "3.13")
    cmake_policy(SET CMP0076 NEW)
endif ()

#### Build all tests in ${TEST_FILES}
foreach (FILE ${TEST_FILES})
    #### Do not build "tests/MATLAB/matlab.cpp" if there's no MATLAB support
    if (${FILE} STREQUAL "tests/MATLAB/matlab.cpp" AND NOT BUILD_WITH_MATLAB)
        continue()
    endif ()
    target_sources(${TARGET_NAME} PUBLIC ${FILE})
    if (${BUILD_WITH_MATLAB})
        if (WIN32)
            if (MSVC)
                target_link_libraries(${TARGET_NAME} PUBLIC libmx libmat)
            elseif (MINGW)
                target_link_libraries(${TARGET_NAME} PUBLIC mx mat)
            endif ()
        else ()
            target_link_libraries(${TARGET_NAME} PUBLIC mx mat)
        endif ()
    endif ()
    #### Eigen3 was found automatically
    if (TARGET Eigen3::Eigen)
        target_link_libraries(${TARGET_NAME} PUBLIC Eigen3::Eigen)
    endif ()
    #### OpenMP was found automatically
    if (OpenMP_CXX_FOUND AND (OpenMP_CXX_VERSION_MAJOR GREATER_EQUAL 3))
        target_link_libraries(${TARGET_NAME} PUBLIC OpenMP::OpenMP_CXX)
    endif ()
    target_link_libraries(${TARGET_NAME} PUBLIC gmock)
    target_link_libraries(${TARGET_NAME} PUBLIC libqpp)
    #### On Linux we get a linker error related to __cxa_thread_atexit()
    if (UNIX AND NOT (APPLE OR CYGWIN)
            AND (${CMAKE_CXX_COMPILER_ID} MATCHES "Clang")
            AND (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 4.0))
        target_link_libraries(${TARGET_NAME} supc++)
    endif ()
endforeach ()

gtest_discover_tests(${TARGET_NAME})
