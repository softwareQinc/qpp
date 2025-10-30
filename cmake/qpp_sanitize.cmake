# Code sanitizing
option(QPP_SANITIZE "Enable code sanitizing" OFF)

# Sanitizing (only AddressSanitizer is available on MSVC)
if(QPP_SANITIZE)
  set(QPP_SANITIZE_FLAGS)
  list(
    APPEND
    QPP_SANITIZE_FLAGS
    $<$<OR:$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:GNU>>:-fsanitize=undefined>
  )
  list(APPEND QPP_SANITIZE_FLAGS $<$<CXX_COMPILER_ID:MSVC>:/fsanitize=address>)
  target_compile_options(libqpp INTERFACE ${QPP_SANITIZE_FLAGS})
  target_link_options(
    libqpp
    INTERFACE
    $<$<OR:$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:GNU>>:-fsanitize=undefined>
  )
endif()
