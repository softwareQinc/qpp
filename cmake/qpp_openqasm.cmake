# OpenQASM 2.0 specs, see DISCREPANCIES.md for a comparison with Qiskit

option(QASMTOOLS_QASM2_SPECS
       "Use OpenQASM 2.0 standard instead of Qiskit gate specifications" OFF)

if(QASMTOOLS_QASM2_SPECS)
  message(STATUS "OpenQASM2 specs - ON")

  # Select the target
  include(${CMAKE_CURRENT_LIST_DIR}/qpp_select_target.cmake)
  qpp_select_target(QPP_TARGET "QASMTOOLS_QASM2_SPECS")

  target_compile_definitions(${QPP_TARGET} INTERFACE QASMTOOLS_QASM2_SPECS)
else()
  message(STATUS "OpenQASM2 specs - OFF")
endif()
