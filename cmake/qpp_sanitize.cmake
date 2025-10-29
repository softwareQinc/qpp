# Code sanitizing
option(SANITIZE "Enable code sanitizing (only for GCC/Clang)" OFF)

# Sanitizing (only AddressSanitizer is available on MSVC)
if(SANITIZE)
  set(SANITIZE_FLAGS)
  list(
    APPEND
    SANITIZE_FLAGS
    $<$<OR:$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:GNU>>:-fsanitize=undefined>
  )
  list(APPEND SANITIZE_FLAGS $<$<CXX_COMPILER_ID:MSVC>:/fsanitize=address>)
  target_compile_options(libqpp INTERFACE ${SANITIZE_FLAGS})
  target_link_options(
    libqpp
    INTERFACE
    $<$<OR:$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:GNU>>:-fsanitize=undefined>
  )
endif()
