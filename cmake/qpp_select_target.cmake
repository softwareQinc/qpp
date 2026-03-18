# -----------------------------------------------------------------------------
# Function: qpp_select_target
#
# Selects the appropriate Quantum++ target depending on context
#
# If the internal development target `libqpp_internal` exists (i.e., when
# building inside the Quantum++ source tree), it is selected. Otherwise, the
# public installed target `qpp::qpp` is used.
#
# This allows the same configuration logic to work correctly in both
#   - build-tree usage (examples, tests, optional components)
#   - install-tree usage (find_package(qpp))
#
# Usage: qpp_select_target(<output_var> [feature_name])
#
# Arguments: <output_var> Name of the variable that will receive the selected
# target
#
# [feature_name] Optional human-readable feature name used only for status
# messages (e.g., "OpenMP", "MATLAB", "OpenQASM")
#
# Example: qpp_select_target(QPP_TARGET "OpenMP")
# target_link_libraries(my_target PRIVATE ${QPP_TARGET})
#
# message(STATUS "Configuring OpenMP support using target ${QPP_TARGET}")
# ----------------------------------------------------------------------------
function(qpp_select_target result_var)
  # Optional feature name for status/warning messages
  set(feature_name "${ARGV1}")
  if(NOT feature_name)
    set(feature_name "QPP feature")
  endif()

  # Determine which target to use
  if(TARGET libqpp_internal)
    set(target_name libqpp_internal)
  elseif(TARGET qpp)
    set(target_name qpp)
  elseif(TARGET qpp::qpp)
    set(target_name qpp::qpp)
  else()
    message(
      FATAL_ERROR "Quantum++: No target found for ${feature_name} configuration"
    )
  endif()

  # Return the detected target to the caller
  set(${result_var}
      "${target_name}"
      PARENT_SCOPE)
endfunction()
