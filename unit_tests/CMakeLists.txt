set(TARGET_NAME "unit_tests")

include(FetchContent)
message(STATUS "Fetching GoogleTest...")
FetchContent_Declare(
  googletest
  GIT_REPOSITORY https://github.com/google/googletest.git
  GIT_TAG main
  GIT_PROGRESS TRUE
  GIT_SHALLOW TRUE)

if(MSVC)
  # Ensure Google Test uses the Shared C Runtime (/MD or /MDd) to prevent linker
  # and runtime errors due to mixed CRT usage between the main project and the
  # GTest libraries.

  # cmake-lint: disable=C0103
  set(gtest_force_shared_crt
      ON
      CACHE BOOL "Force Google Test to use shared C Runtime on MSVC" FORCE)
endif()

FetchContent_MakeAvailable(googletest)

file(GLOB_RECURSE TEST_SOURCES "tests/*.cpp" "tests/MATLAB/*.cpp"
     "tests/classes/*.cpp" "tests/qasm/*.cpp")

set(TEST_SOURCES_FILTERED "")
foreach(file ${TEST_SOURCES})
  # Filter: Do not build "tests/MATLAB/matlab.cpp" if there's no MATLAB support
  if(NOT ("${file}" MATCHES "matlab\\.cpp" AND NOT MATLAB_FOUND))
    list(APPEND TEST_SOURCES_FILTERED "${file}")
  else()
    message(STATUS "Skipping test file: ${file} (MATLAB support not enabled)")
  endif()
endforeach()

add_executable(${TARGET_NAME} EXCLUDE_FROM_ALL tests/main.cpp
                                               ${TEST_SOURCES_FILTERED})
set_target_properties(
  ${TARGET_NAME}
  PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/tests"
             RUNTIME_OUTPUT_DIRECTORY_DEBUG "${CMAKE_BINARY_DIR}/tests"
             RUNTIME_OUTPUT_DIRECTORY_RELEASE "${CMAKE_BINARY_DIR}/tests")
target_link_libraries(${TARGET_NAME} PRIVATE libqpp_internal GTest::gmock)

include(GoogleTest)
gtest_discover_tests(${TARGET_NAME})
