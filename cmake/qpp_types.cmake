# Custom index type. If none selected, a default one is selected in
# <qpp/types.hpp> (usually std::size_t).
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

# Determine the index type compile definition
if(QPP_IDX STREQUAL "default")
  message(STATUS "Index type - default (see <qpp/types.hpp>)")
  set(QPP_IDX_DEFINITION QPP_IDX_DEFAULT)
elseif(QPP_IDX STREQUAL "short")
  message(STATUS "Index type - short")
  set(QPP_IDX_DEFINITION QPP_IDX_SHORT)
elseif(QPP_IDX STREQUAL "int")
  message(STATUS "Index type - int")
  set(QPP_IDX_DEFINITION QPP_IDX_INT)
elseif(QPP_IDX STREQUAL "long")
  message(STATUS "Index type - long")
  set(QPP_IDX_DEFINITION QPP_IDX_LONG)
elseif(QPP_IDX STREQUAL "long long")
  message(STATUS "Index type - long long")
  set(QPP_IDX_DEFINITION QPP_IDX_LONG_LONG)
elseif(QPP_IDX STREQUAL "unsigned short")
  message(STATUS "Index type - unsigned short")
  set(QPP_IDX_DEFINITION QPP_IDX_USHORT)
elseif(QPP_IDX STREQUAL "unsigned int")
  message(STATUS "Index type - unsigned int")
  set(QPP_IDX_DEFINITION QPP_IDX_UINT)
elseif(QPP_IDX STREQUAL "unsigned long")
  message(STATUS "Index type - unsigned long")
  set(QPP_IDX_DEFINITION QPP_IDX_ULONG)
elseif(QPP_IDX STREQUAL "unsigned long long")
  message(STATUS "Index type - unsigned long long")
  set(QPP_IDX_DEFINITION QPP_IDX_ULONG_LONG)
endif()

# Custom signed big integer type. If none selected, a default one is selected in
# <qpp/types.hpp> (usually long long).
set(QPP_BIGINT
        "default"    
    CACHE STRING "Default big integer type, see <qpp/types.hpp>")
set_property(CACHE QPP_BIGINT PROPERTY STRINGS "default" "short" "int" "long"
                                        "long long")

# Determine the big integer type compile definition
if(QPP_BIGINT STREQUAL "default")
  message(STATUS "Big integer type - default (see <qpp/types.hpp>)")
  set(QPP_BIGINT_DEFINITION QPP_BIGINT_DEFAULT)
elseif(QPP_BIGINT STREQUAL "short")
  message(STATUS "Big integer type - short")
  set(QPP_BIGINT_DEFINITION QPP_BIGINT_SHORT)
elseif(QPP_BIGINT STREQUAL "int")
  message(STATUS "Big integer type - int")
  set(QPP_BIGINT_DEFINITION QPP_BIGINT_INT)
elseif(QPP_BIGINT STREQUAL "long")
  message(STATUS "Big integer type - long")
  set(QPP_BIGINT_DEFINITION QPP_BIGINT_LONG)
elseif(QPP_BIGINT STREQUAL "long long")
  message(STATUS "Big integer type - long long")
  set(QPP_BIGINT_DEFINITION QPP_BIGINT_LONG_LONG)
endif()

# Custom floating-point type. If none selected, a default one is selected in
# <qpp/types.hpp> (usually double).
set(QPP_FP
        "default"    
    CACHE STRING "Default floating-point type, see <qpp/types.hpp>")
set_property(CACHE QPP_FP PROPERTY STRINGS "default" "float" "double"
                                    "long double")

# Determine the floating-point type compile definition
if(QPP_FP STREQUAL "default")
  message(STATUS "Floating-point type - default (see <qpp/types.hpp>)")
  set(QPP_FP_DEFINITION QPP_FP_DEFAULT)
elseif(QPP_FP STREQUAL "float")
  message(STATUS "Floating-point type - float")
  set(QPP_FP_DEFINITION QPP_FP_FLOAT)
elseif(QPP_FP STREQUAL "double")
  message(STATUS "Floating-point type - double")
  set(QPP_FP_DEFINITION QPP_FP_DOUBLE)
elseif(QPP_FP STREQUAL "long double")
  message(STATUS "Floating-point type - long double")
  set(QPP_FP_DEFINITION QPP_FP_LONG_DOUBLE)
endif()

target_compile_definitions(
  libqpp INTERFACE ${QPP_IDX_DEFINITION} ${QPP_BIGINT_DEFINITION}
                   ${QPP_FP_DEFINITION})
