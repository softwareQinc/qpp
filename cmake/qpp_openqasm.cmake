# OpenQASM 2.0 specs, see DISCREPANCIES.md for a comparison with Qiskit
option(QASMTOOLS_QASM2_SPECS
       "Use OpenQASM 2.0 standard instead of Qiskit gate specifications" OFF)
target_compile_definitions(
  libqpp
  INTERFACE $<$<BOOL:${QASMTOOLS_QASM2_SPECS}>:QASMTOOLS_QASM2_SPECS=true>)
if(QASMTOOLS_QASM2_SPECS)
  message(STATUS "OpenQASM2 specs - ON")
else()
  message(STATUS "OpenQASM2 specs - OFF")
endif()
