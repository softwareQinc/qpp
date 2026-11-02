# Code coverage support

message(STATUS "Code coverage - ${QPP_ENABLE_COVERAGE}")
if(QPP_ENABLE_COVERAGE)
  if(NOT CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    message(FATAL_ERROR "LLVM coverage requires Clang")
  endif()

  # Enable LLVM coverage instrumentation and UBSan
  target_compile_options(
    unit_tests PRIVATE -O0 -g -fprofile-instr-generate -fcoverage-mapping
                       -fsanitize=undefined -fno-omit-frame-pointer)

  target_link_options(unit_tests PRIVATE -fprofile-instr-generate
                      -fsanitize=undefined)

  # Coverage output paths
  set(COVERAGE_DIR "${CMAKE_BINARY_DIR}/coverage")
  set(COVERAGE_RAW "${COVERAGE_DIR}/coverage.profraw")
  set(COVERAGE_DATA "${COVERAGE_DIR}/coverage.profdata")

  # Exclude third-party dependencies and unit tests
  set(COVERAGE_EXCLUDE "(/Eigen/|/googletest/|/qasmtools/|/unit_tests/)")

  # Additionally exclude Quantum++ internal implementations
  set(COVERAGE_PUBLIC_EXCLUDE
      "(/Eigen/|/googletest/|/qasmtools/|/unit_tests/|/qpp/internal/)")

  add_custom_target(
    coverage
    # Create coverage directory
    COMMAND ${CMAKE_COMMAND} -E make_directory "${COVERAGE_DIR}"
    # Remove previous profiling data
    COMMAND ${CMAKE_COMMAND} -E rm -f "${COVERAGE_RAW}"
    # Execute all unit tests
    COMMAND ${CMAKE_COMMAND} -E env "LLVM_PROFILE_FILE=${COVERAGE_RAW}"
            $<TARGET_FILE:unit_tests>
    # Merge LLVM profiling data
    COMMAND llvm-profdata merge -sparse "${COVERAGE_RAW}" -o "${COVERAGE_DATA}"
    # Generate public API HTML report
    COMMAND
      llvm-cov show $<TARGET_FILE:unit_tests> -instr-profile=${COVERAGE_DATA}
      -format=html -output-dir=${COVERAGE_DIR}/public
      -show-line-counts-or-regions -show-branches=count
      -ignore-filename-regex=${COVERAGE_PUBLIC_EXCLUDE}
    # Generate full Quantum++ HTML report
    COMMAND
      llvm-cov show $<TARGET_FILE:unit_tests> -instr-profile=${COVERAGE_DATA}
      -format=html -output-dir=${COVERAGE_DIR}/full -show-line-counts-or-regions
      -show-branches=count -ignore-filename-regex=${COVERAGE_EXCLUDE}
      # Export public API coverage in LCOV format (for Codecov)
    COMMAND
      /bin/sh -c
      "llvm-cov export \"$1\" -instr-profile=\"$2\" -format=lcov -ignore-filename-regex=\"$3\" > \"$4\""
      _ "$<TARGET_FILE:unit_tests>" "${COVERAGE_DATA}"
      "${COVERAGE_PUBLIC_EXCLUDE}" "${COVERAGE_DIR}/coverage-public.info"
    # Export full Quantum++ coverage in LCOV format
    COMMAND
      /bin/sh -c
      "llvm-cov export \"$1\" -instr-profile=\"$2\" -format=lcov -ignore-filename-regex=\"$3\" > \"$4\""
      _ "$<TARGET_FILE:unit_tests>" "${COVERAGE_DATA}" "${COVERAGE_EXCLUDE}"
      "${COVERAGE_DIR}/coverage-full.info"
      # Print public API coverage summary
    COMMAND ${CMAKE_COMMAND} -E echo ""
    COMMAND ${CMAKE_COMMAND} -E echo "========== Public API Coverage =========="
    COMMAND
      llvm-cov report $<TARGET_FILE:unit_tests> -instr-profile=${COVERAGE_DATA}
      -ignore-filename-regex=${COVERAGE_PUBLIC_EXCLUDE}
    # Print full coverage summary
    COMMAND ${CMAKE_COMMAND} -E echo ""
    COMMAND ${CMAKE_COMMAND} -E echo
            "========== Full Quantum++ Coverage =========="
    COMMAND
      llvm-cov report $<TARGET_FILE:unit_tests> -instr-profile=${COVERAGE_DATA}
      -ignore-filename-regex=${COVERAGE_EXCLUDE}
    DEPENDS unit_tests
    WORKING_DIRECTORY "${CMAKE_BINARY_DIR}"
    COMMENT "Generating Quantum++ coverage reports"
    VERBATIM USES_TERMINAL)
endif()
