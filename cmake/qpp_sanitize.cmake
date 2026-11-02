# Code sanitizing support

option(QPP_ENABLE_SANITIZE "Enable code sanitizing" OFF)
message(STATUS "Code sanitizing - ${QPP_ENABLE_SANITIZE}")

# Sanitizing (only AddressSanitizer is available on MSVC)
if(QPP_ENABLE_SANITIZE)
  # Select the target
  include(${CMAKE_CURRENT_LIST_DIR}/qpp_select_target.cmake)
  qpp_select_target(QPP_TARGET "QPP_ENABLE_SANITIZE")

  set(QPP_ENABLE_SANITIZE_FLAGS)
  list(
    APPEND
    QPP_ENABLE_SANITIZE_FLAGS
    $<$<OR:$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:GNU>>:-fsanitize=undefined>
  )
  list(APPEND QPP_ENABLE_SANITIZE_FLAGS $<$<CXX_COMPILER_ID:MSVC>:/fsanitize=address>)
  target_compile_options(${QPP_TARGET} INTERFACE ${QPP_ENABLE_SANITIZE_FLAGS})
  target_link_options(
    ${QPP_TARGET}
    INTERFACE
    $<$<OR:$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:GNU>>:-fsanitize=undefined>
  )
endif()
