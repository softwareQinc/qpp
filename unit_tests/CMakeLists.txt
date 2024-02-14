include(GoogleTest)
set(TARGET_NAME "unit_tests")
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# GoogleTest
include(FetchContent)
message(STATUS "Fetching GoogleTest...")
FetchContent_Declare(
  googletest
  GIT_REPOSITORY https://github.com/google/googletest.git
  GIT_TAG main
  GIT_SHALLOW TRUE
  GIT_PROGRESS TRUE)
# For Windows: Prevent overriding the parent project's compiler/linker settings,
# and add bigobj option
if(MSVC)
  add_compile_options(-bigobj)
  set(gtest_force_shared_crt
      ON
      CACHE BOOL "" FORCE)
endif()
FetchContent_MakeAvailable(googletest)

aux_source_directory(tests TEST_FILES)
aux_source_directory(tests/MATLAB TEST_FILES)
aux_source_directory(tests/classes TEST_FILES)
aux_source_directory(tests/qasm TEST_FILES)

add_executable(${TARGET_NAME} EXCLUDE_FROM_ALL tests/main.cpp)
add_dependencies(unit_tests ${TARGET_NAME})

if(${CMAKE_VERSION} VERSION_GREATER_EQUAL "3.13")
  cmake_policy(SET CMP0076 NEW)
endif()

# Build all tests in ${TEST_FILES}
foreach(file ${TEST_FILES})
  # Do not build "tests/MATLAB/matlab.cpp" if there's no MATLAB support
  if(${file} STREQUAL "tests/MATLAB/matlab.cpp" AND NOT BUILD_WITH_MATLAB)
    continue()
  endif()
  target_sources(${TARGET_NAME} PUBLIC ${file})
endforeach()

target_link_libraries(${TARGET_NAME} PUBLIC ${QPP_LINK_DEPS} gmock libqpp)
gtest_discover_tests(${TARGET_NAME})
