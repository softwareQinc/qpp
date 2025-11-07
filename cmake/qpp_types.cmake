# ============================================================================
# Quantum++ type configuration
# ============================================================================
# This file defines the fundamental type aliases used by Quantum++. Each type
# can be customized via CMake cache options. If none is provided, defaults are
# selected in <qpp/types.hpp>.
# ============================================================================

# ----------------------------------------------------------------------------
# Index type
# ----------------------------------------------------------------------------
set(QPP_IDX
    "default"
    CACHE STRING "Default index type, see <qpp/types.hpp>")

set_property(
  CACHE QPP_IDX
  PROPERTY STRINGS
           "default"
           "short"
           "int"
           "long"
           "long long"
           "unsigned short"
           "unsigned int"
           "unsigned long"
           "unsigned long long")

# ----------------------------------------------------------------------------
# Big integer type
# ----------------------------------------------------------------------------
set(QPP_BIGINT
    "default"
    CACHE STRING "Default big integer type, see <qpp/types.hpp>")

set_property(CACHE QPP_BIGINT PROPERTY STRINGS "default" "short" "int" "long"
                                       "long long")

# ----------------------------------------------------------------------------
# Floating-point type
# ----------------------------------------------------------------------------
set(QPP_FP
    "default"
    CACHE STRING "Default floating-point type, see <qpp/types.hpp>")

set_property(CACHE QPP_FP PROPERTY STRINGS "default" "float" "double"
                                   "long double")

# ----------------------------------------------------------------------------
# Print configuration summary
# ----------------------------------------------------------------------------
message(STATUS "Quantum++ type configuration")
message(STATUS "Index type        : ${QPP_IDX}")
message(STATUS "Big integer type  : ${QPP_BIGINT}")
message(STATUS "Floating-point    : ${QPP_FP}")

# ----------------------------------------------------------------------------
# Resolve compile-time definitions
# ----------------------------------------------------------------------------

# ---- Index type ------------------------------------------------------------
if(QPP_IDX STREQUAL "default")
  set(QPP_IDX_DEFINITION QPP_IDX_DEFAULT)
elseif(QPP_IDX STREQUAL "short")
  set(QPP_IDX_DEFINITION QPP_IDX_SHORT)
elseif(QPP_IDX STREQUAL "int")
  set(QPP_IDX_DEFINITION QPP_IDX_INT)
elseif(QPP_IDX STREQUAL "long")
  set(QPP_IDX_DEFINITION QPP_IDX_LONG)
elseif(QPP_IDX STREQUAL "long long")
  set(QPP_IDX_DEFINITION QPP_IDX_LONG_LONG)
elseif(QPP_IDX STREQUAL "unsigned short")
  set(QPP_IDX_DEFINITION QPP_IDX_USHORT)
elseif(QPP_IDX STREQUAL "unsigned int")
  set(QPP_IDX_DEFINITION QPP_IDX_UINT)
elseif(QPP_IDX STREQUAL "unsigned long")
  set(QPP_IDX_DEFINITION QPP_IDX_ULONG)
elseif(QPP_IDX STREQUAL "unsigned long long")
  set(QPP_IDX_DEFINITION QPP_IDX_ULONG_LONG)
else()
  message(FATAL_ERROR "Invalid QPP_IDX value: ${QPP_IDX}")
endif()

# ---- Big integer type ------------------------------------------------------
if(QPP_BIGINT STREQUAL "default")
  set(QPP_BIGINT_DEFINITION QPP_BIGINT_DEFAULT)
elseif(QPP_BIGINT STREQUAL "short")
  set(QPP_BIGINT_DEFINITION QPP_BIGINT_SHORT)
elseif(QPP_BIGINT STREQUAL "int")
  set(QPP_BIGINT_DEFINITION QPP_BIGINT_INT)
elseif(QPP_BIGINT STREQUAL "long")
  set(QPP_BIGINT_DEFINITION QPP_BIGINT_LONG)
elseif(QPP_BIGINT STREQUAL "long long")
  set(QPP_BIGINT_DEFINITION QPP_BIGINT_LONG_LONG)
else()
  message(FATAL_ERROR "Invalid QPP_BIGINT value: ${QPP_BIGINT}")
endif()

# ---- Floating-point type ---------------------------------------------------
if(QPP_FP STREQUAL "default")
  set(QPP_FP_DEFINITION QPP_FP_DEFAULT)
elseif(QPP_FP STREQUAL "float")
  set(QPP_FP_DEFINITION QPP_FP_FLOAT)
elseif(QPP_FP STREQUAL "double")
  set(QPP_FP_DEFINITION QPP_FP_DOUBLE)
elseif(QPP_FP STREQUAL "long double")
  set(QPP_FP_DEFINITION QPP_FP_LONG_DOUBLE)
else()
  message(FATAL_ERROR "Invalid QPP_FP value: ${QPP_FP}")
endif()

# ----------------------------------------------------------------------------
# Apply to Quantum++ target
# ----------------------------------------------------------------------------
target_compile_definitions(
  libqpp INTERFACE ${QPP_IDX_DEFINITION} ${QPP_BIGINT_DEFINITION}
                   ${QPP_FP_DEFINITION})

# ----------------------------------------------------------------------------
# Final summary (for clarity in build logs)
# ----------------------------------------------------------------------------
message(STATUS "Configured compile definitions:")
message(STATUS "  ${QPP_IDX_DEFINITION}")
message(STATUS "  ${QPP_BIGINT_DEFINITION}")
message(STATUS "  ${QPP_FP_DEFINITION}")
