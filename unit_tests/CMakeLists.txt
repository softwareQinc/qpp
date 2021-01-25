include(GoogleTest)
set(TARGET_NAME "unit_tests")

#### Link Google Test dynamically if using MSVC
if (MSVC)
    set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
endif ()

add_subdirectory(lib/googletest-release-1.10.0 EXCLUDE_FROM_ALL)

aux_source_directory(tests TEST_FILES)
aux_source_directory(tests/MATLAB TEST_FILES)
aux_source_directory(tests/classes TEST_FILES)
aux_source_directory(tests/classes/circuits TEST_FILES)
aux_source_directory(tests/qasm TEST_FILES)

add_executable(${TARGET_NAME} EXCLUDE_FROM_ALL tests/main.cpp)
add_dependencies(unit_tests ${TARGET_NAME})

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
    target_link_libraries(${TARGET_NAME} PUBLIC ${QPP_LINK_DEPS} libqpp gmock)
endforeach ()

gtest_discover_tests(${TARGET_NAME})
