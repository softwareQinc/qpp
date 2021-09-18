#### Quantum++ additional dependencies
#### Do not modify unless you know what you're doing

#### Eigen3
message(STATUS "Detecting Eigen3")
set(EIGEN3_INSTALL_DIR "" CACHE PATH "Path to Eigen3")
#### Location manually specified
if (NOT ${EIGEN3_INSTALL_DIR} STREQUAL "")
    if (IS_DIRECTORY ${EIGEN3_INSTALL_DIR})
        message(STATUS "Detecting Eigen3 - done (in ${EIGEN3_INSTALL_DIR})")
        include_directories(SYSTEM "${EIGEN3_INSTALL_DIR}")
    else ()
        message(FATAL_ERROR "Invalid path to Eigen3 installation")
    endif ()
else () #### Try to find it automatically
    find_package(Eigen3 3.0 QUIET NO_MODULE)
    if (NOT TARGET Eigen3::Eigen) # did not find Eigen3 automatically
        message(FATAL_ERROR
                "Eigen3 not detected! Please point EIGEN3_INSTALL_DIR\
            to your Eigen3 location when building with cmake,\
            for example\

            cmake .. -DEIGEN3_INSTALL_DIR=$HOME/eigen3")
    endif ()
    message(STATUS "Detecting Eigen3 - done (in ${EIGEN3_INCLUDE_DIR})")
    # Eigen3 header-only dependencies to be injected in the main CMakeLists.txt
    set(QPP_EIGEN3_LINK_DEPS Eigen3::Eigen)
endif ()

#### OpenMP support
option(WITH_OPENMP "OpenMP support" ON)
if (${WITH_OPENMP})
    find_package(OpenMP)
    if (OpenMP_CXX_FOUND)
        if (OpenMP_CXX_VERSION_MAJOR GREATER_EQUAL 3)
            #### inject definition (as #define) in the source files
            add_compile_definitions(HAS_OPENMP)
            # OpenMP linking dependencies to be injected in the main CMakeLists.txt
            set(QPP_OPENMP_LINK_DEPS OpenMP::OpenMP_CXX)
        else ()
            message(NOTICE "Found OpenMP version ${OpenMP_CXX_VERSION}, Quantum++ \
requires OpenMP 3.0 or later")
        endif ()
    endif ()
endif ()

#### MATLAB support, disabled by default
option(WITH_MATLAB "MATLAB support" OFF)
if (${WITH_MATLAB})
    message(STATUS "Detecting MATLAB")
    #### Try to find it automatically
    find_package(Matlab OPTIONAL_COMPONENTS MX_LIBRARY MAT_LIBRARY QUIET)
    if (MATLAB_FOUND)
        message(STATUS "Detecting MATLAB - done (in ${Matlab_ROOT_DIR})")
        include_directories(SYSTEM ${Matlab_INCLUDE_DIRS})
        if (WIN32)
            if (MSVC)
                set(MATLAB_LIB_DIR
                        "${Matlab_ROOT_DIR}/extern/lib/win64/microsoft"
                        CACHE PATH "Custom path to MATLAB lib directory")
            elseif (MINGW)
                set(MATLAB_LIB_DIR "${Matlab_ROOT_DIR}/extern/lib/win64/mingw64"
                        CACHE PATH "Custom path to MATLAB lib directory")
            else ()
                message(FATAL "Platform not supported, aborting.")
            endif ()
        elseif (UNIX AND NOT APPLE)
            set(MATLAB_LIB_DIR "${Matlab_ROOT_DIR}/bin/glnxa64" CACHE
                    PATH "Custom path to MATLAB lib directory")
        elseif (APPLE)
            set(MATLAB_LIB_DIR "${Matlab_ROOT_DIR}/bin/maci64" CACHE
                    PATH "Custom path to MATLAB lib directory")
        else ()
            message(FATAL "Platform not supported, aborting.")
        endif ()
        link_directories(${MATLAB_LIB_DIR})
        add_compile_definitions(HAS_MATLAB)
        set(BUILD_WITH_MATLAB TRUE)
    else () #### Location manually specified
        set(MATLAB_INSTALL_DIR "" CACHE PATH
                "Custom path to MATLAB installation")
        if (IS_DIRECTORY ${MATLAB_INSTALL_DIR})
            #### MATLAB include files
            set(MATLAB_INCLUDE_DIR "${MATLAB_INSTALL_DIR}/extern/include"
                    CACHE PATH "Custom path to MATLAB include directory")
            if (IS_DIRECTORY ${MATLAB_INCLUDE_DIR})
                include_directories(SYSTEM ${MATLAB_INCLUDE_DIR})
            else ()
                message(FATAL_ERROR
                        "Possibly corrupted MATLAB include headers")
            endif ()
            #### MATLAB linker files
            if (WIN32)
                if (MSVC)
                    set(MATLAB_LIB_DIR
                            "${MATLAB_INSTALL_DIR}/extern/lib/win64/microsoft"
                            CACHE PATH "Custom path to MATLAB lib directory")
                elseif (MINGW64)
                    set(MATLAB_LIB_DIR
                            "${MATLAB_INSTALL_DIR}/extern/lib/win64/mingw64"
                            CACHE PATH "Custom path to MATLAB lib directory")
                else ()
                    message(FATAL "Platform not supported, aborting.")
                endif ()
            elseif (UNIX AND NOT APPLE)
                set(MATLAB_LIB_DIR "${MATLAB_INSTALL_DIR}/bin/glnxa64" CACHE
                        PATH "Custom path to MATLAB lib directory")
            elseif (APPLE)
                set(MATLAB_LIB_DIR "${MATLAB_INSTALL_DIR}/bin/maci64" CACHE
                        PATH "Custom path to MATLAB lib directory")
            else ()
                message(FATAL "Platform not supported, aborting.")
            endif ()
            if (IS_DIRECTORY ${MATLAB_LIB_DIR})
                link_directories(${MATLAB_LIB_DIR})
            else ()
                message(FATAL_ERROR
                        "Possibly corrupted MATLAB compiler libraries")
            endif ()
            #### Everything is OK, inject definition (as #define) in the source
            message(STATUS
                    "Detecting MATLAB - done (in ${MATLAB_INSTALL_DIR})")
            add_compile_definitions(HAS_MATLAB_ENABLED)
            set(BUILD_WITH_MATLAB TRUE)
        else ()
            message(FATAL_ERROR "Could not detect MATLAB, aborting")
        endif ()
    endif ()
    # MATLAB linking dependencies to be injected in the main CMakeLists.txt
    if (${BUILD_WITH_MATLAB})
        if (WIN32)
            if (MSVC)
                set(QPP_MATLAB_LINK_DEPS libmx libmat)
            elseif (MINGW)
                set(QPP_MATLAB_LINK_DEPS mx mat)
            endif ()
        else ()
            set(QPP_MATLAB_LINK_DEPS mx mat)
        endif ()
    endif ()
endif ()

##### Disable support for thread_local storage duration specifier
##### when using AppleClang as libc++ doesn't yet support it
#if (${CMAKE_CXX_COMPILER_ID} STREQUAL "AppleClang")
#    #### inject definition (as #define) in the source files
#    add_definitions(-DNO_THREAD_LOCAL_)
#    message(WARNING "Detected compiler: ${CMAKE_CXX_COMPILER_ID} \
#    ${CMAKE_CXX_COMPILER_VERSION}. thread_local not supported.")
#endif ()

#### Windows issues with Microsoft Visual Studio
if (MSVC)
    # Disable spurious Eigen warnings with MSVC (warning STL4007)
    add_compile_definitions(_SILENCE_CXX17_ADAPTOR_TYPEDEFS_DEPRECATION_WARNING)
    add_compile_options(-bigobj)
    add_definitions(-DNOMINMAX)
endif ()

#### MinGW or Cygwin have issues with object files that are too large
if (MINGW OR CYGWIN)
    add_compile_options("-Wa,-mbig-obj")
endif ()

#### Cygwin has issues with std=c++11, use std=gnu++11 instead
if (CYGWIN)
    add_compile_options(-std=gnu++11)
endif ()

#### Force clang to use libc++
if (${CMAKE_CXX_COMPILER_ID} MATCHES "Clang")
    add_compile_options(-stdlib=libc++)
    list(APPEND QPP_LINK_DEPS c++)
endif ()

#### GNU gcc additional debug settings
if (${CMAKE_CXX_COMPILER_ID} MATCHES "GNU")
    #if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    #    ## use the "no-weak" debugging flag only when debugging under OS X,
    #    ## as gdb cannot step in template functions when debugging code
    #    ## produced by g++
    #    ## see http://stackoverflow.com/questions/
    #    ## 23330641/gnu-gdb-can-not-step-into-template-functions-os-x-mavericks
    #    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fno-weak")
    #endif ()
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Og -D_GLIBCXX_DEBUG")
endif ()

#### Configurations
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DDEBUG")
set(CMAKE_CXX_FLAGS_MINSIZEREL
        "${CMAKE_CXX_FLAGS_MINSIZEREL} -DEIGEN_NO_DEBUG")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -DEIGEN_NO_DEBUG")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO
        "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} -DEIGEN_NO_DEBUG")

#### Default build type
if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING
            "Choose the type of build, options are: \
         None Debug Release MinSizeRel RelWithDebInfo."
            FORCE)
endif ()

#### On Linux we get a linker error related to __cxa_thread_atexit()
if (UNIX AND NOT (APPLE OR CYGWIN)
        AND (${CMAKE_CXX_COMPILER_ID} MATCHES "Clang")
        AND (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 4.0))
    list(APPEND QPP_LINK_DEPS supc++)
endif ()

#### Export all settings to outside
list(APPEND QPP_LINK_DEPS ${QPP_EIGEN3_LINK_DEPS} ${QPP_OPENMP_LINK_DEPS} ${QPP_MATLAB_LINK_DEPS})
