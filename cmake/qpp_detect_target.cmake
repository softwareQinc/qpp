# Function: qpp_detect_target
#
# Detects the main QPP target in both build and install contexts.
#
# Usage: qpp_detect_target(<output_var> <feature_name>)
#
# <output_var>  : Name of the variable to store the detected target
# <feature_name>: Optional string for messages (e.g., "OpenMP", "MATLAB")
#
# Example:
#
# qpp_detect_target(QPP_TARGET "OpenMP")
#
# message(STATUS "Using target ${QPP_TARGET} for OpenMP configuration")
function(qpp_detect_target result_var)
  # Optional second argument for feature name
  set(feature_name "${ARGV1}")
  if(NOT feature_name)
    set(feature_name "QPP feature")
  endif()

  if(TARGET libqpp_internal)
    set(target_name libqpp_internal)
  elseif(TARGET qpp::libqpp)
    set(target_name qpp::libqpp)
  elseif(TARGET libqpp)
    set(target_name libqpp)
  else()
    message(
      WARNING "Quantum++: No target found for ${feature_name} configuration")
    set(target_name "")
  endif()

  # Return detected target name to the caller
  set(${result_var}
      "${target_name}"
      PARENT_SCOPE)
endfunction()
