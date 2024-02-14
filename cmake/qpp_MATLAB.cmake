# MATLAB support, disabled by default
option(QPP_MATLAB "MATLAB support" OFF)
if(${QPP_MATLAB})
  message(STATUS "Detecting MATLAB")
  # Try to find it automatically
  find_package(
    Matlab
    OPTIONAL_COMPONENTS MX_LIBRARY MAT_LIBRARY
    QUIET)
  if(MATLAB_FOUND)
    message(STATUS "Detecting MATLAB - done (in ${Matlab_ROOT_DIR})")
    include_directories(SYSTEM ${Matlab_INCLUDE_DIRS})
    if(WIN32)
      if(MSVC)
        set(MATLAB_LIB_DIR
            "${Matlab_ROOT_DIR}/extern/lib/win64/microsoft"
            CACHE PATH "Custom path to MATLAB lib directory")
      elseif(MINGW)
        set(MATLAB_LIB_DIR
            "${Matlab_ROOT_DIR}/extern/lib/win64/mingw64"
            CACHE PATH "Custom path to MATLAB lib directory")
      else()
        message(FATAL_ERROR "Platform not supported, aborting.")
      endif()
    elseif(UNIX AND NOT APPLE)
      set(MATLAB_LIB_DIR
          "${Matlab_ROOT_DIR}/bin/glnxa64"
          CACHE PATH "Custom path to MATLAB lib directory")
    elseif(APPLE)
      if(${CMAKE_HOST_SYSTEM_PROCESSOR} MATCHES "arm64")
        message(
          FATAL_ERROR
            "arm64 architecture is not (yet) supported by MATLAB, aborting.")
      endif()
      set(MATLAB_LIB_DIR
          "${Matlab_ROOT_DIR}/bin/maci64"
          CACHE PATH "Custom path to MATLAB lib directory")
    else()
      message(FATAL_ERROR "Platform not supported, aborting.")
    endif()
    link_directories(${MATLAB_LIB_DIR})
    add_compile_definitions(QPP_MATLAB)
    set(BUILD_WITH_MATLAB TRUE)
  else() # Location manually specified
    set(MATLAB_INSTALL_DIR
        ""
        CACHE PATH "Custom path to MATLAB installation")
    if(IS_DIRECTORY ${MATLAB_INSTALL_DIR})
      # MATLAB include files
      set(MATLAB_INCLUDE_DIR
          "${MATLAB_INSTALL_DIR}/extern/include"
          CACHE PATH "Custom path to MATLAB include directory")
      if(IS_DIRECTORY ${MATLAB_INCLUDE_DIR})
        include_directories(SYSTEM ${MATLAB_INCLUDE_DIR})
      else()
        message(FATAL_ERROR "Possibly corrupted MATLAB include headers")
      endif()
      # MATLAB linker files
      if(WIN32)
        if(MSVC)
          set(MATLAB_LIB_DIR
              "${MATLAB_INSTALL_DIR}/extern/lib/win64/microsoft"
              CACHE PATH "Custom path to MATLAB lib directory")
        elseif(MINGW64)
          set(MATLAB_LIB_DIR
              "${MATLAB_INSTALL_DIR}/extern/lib/win64/mingw64"
              CACHE PATH "Custom path to MATLAB lib directory")
        else()
          message(FATAL_ERROR "Platform not supported, aborting.")
        endif()
      elseif(UNIX AND NOT APPLE)
        set(MATLAB_LIB_DIR
            "${MATLAB_INSTALL_DIR}/bin/glnxa64"
            CACHE PATH "Custom path to MATLAB lib directory")
      elseif(APPLE)
        set(MATLAB_LIB_DIR
            "${MATLAB_INSTALL_DIR}/bin/maci64"
            CACHE PATH "Custom path to MATLAB lib directory")
      else()
        message(FATAL_ERROR "Platform not supported, aborting.")
      endif()
      if(IS_DIRECTORY ${MATLAB_LIB_DIR})
        link_directories(${MATLAB_LIB_DIR})
      else()
        message(FATAL_ERROR "Possibly corrupted MATLAB compiler libraries")
      endif()
      # Everything is OK, inject definition (as #define) in the source
      message(STATUS "Detecting MATLAB - done (in ${MATLAB_INSTALL_DIR})")
      add_compile_definitions(QPP_MATLAB)
      set(BUILD_WITH_MATLAB TRUE)
    else()
      message(FATAL_ERROR "Could not detect MATLAB, aborting")
    endif()
  endif()
  # MATLAB linking dependencies to be injected in the main CMakeLists.txt
  if(${BUILD_WITH_MATLAB})
    if(WIN32)
      if(MSVC)
        set(QPP_MATLAB_LINK_DEPS libmx libmat)
      elseif(MINGW)
        set(QPP_MATLAB_LINK_DEPS mx mat)
      endif()
    else()
      set(QPP_MATLAB_LINK_DEPS mx mat)
    endif()
  endif()
endif()
