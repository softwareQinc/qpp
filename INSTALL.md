# Quantum++ installation instructions

Quantum++ is a header-only library that uses [CMake](https://www.cmake.org/) as its build/install system. Quantum++ is platform-independent, supporting both [UNIX](https://www.opengroup.org/membership/forums/platform/unix) (including [macOS](https://www.apple.com/macos/)) and UNIX-like operating systems (e.g., [Linux](https://www.linux.org)), as well as [Windows](https://www.microsoft.com/en-us/windows). 

## Pre-requisites
- [CMake](http://www.cmake.org/) version 3.12 or later
- [Eigen 3](http://eigen.tuxfamily.org) linear algebra library
- A C++ 11 compliant compiler, e.g., [gcc](https://gcc.gnu.org/)/[clang](https://clang.llvm.org)/[MSVC](https://visualstudio.microsoft.com/vs/) etc.

## Optional
- [MATLAB](http://www.mathworks.com/products/matlab/) compiler shared libraries and include header files, in case you want to enable interoperability with MATLAB. If enabled, allows applications build with Quantum++ to save/load Quantum++ matrices and vectors to/from MATLAB. The MATLAB compiler shared libaries and header files locations are platform-specific, e.g., under Linux, they may be located under `/usr/local/MATLAB/R2021a/bin/glnxa64` and `/usr/local/MATLAB/R2021a/extern/include`, respectively. On your platform the locations may of course differ.

## Configuring the system

For UNIX/UNIX-like/Windows, first create an out-of-source build directory, e.g., from the project's root directory type in a terminal/console/command prompt 

```
mkdir build && cd build
```

Then, configure the system by typing (from inside the newly created `build` directory)

```bash
cmake ..
```

If you plan to compile the [examples](https://github.com/softwareQinc/qpp/tree/main/examples) or the [unit tests](https://github.com/softwareQinc/qpp/tree/main/unit_tests) you may need to pass additional arguments to CMake,

```bash
cmake .. [optional arguments]
```

where [optional arguments] are passed as `-DOPTIONAL_ARGUMENT=VALUE`. The optional arguments are:

Optional argument | Value | Description
| --- | --- | --- |
`EIGEN3_INSTALL_DIR` | `/path/to/eigen3` | Path to Eigen3 installation, if not automatically detected
`WITH_OPENMP` | `ON/OFF` [`ON` by default] | Enables (if available)/disables OpenMP multi-processing library
`WITH_EXAMPLES` | `ON/OFF` [`OFF` by default] | Enables/disables building examples as a CMake target
`WITH_UNIT_TESTS` | `ON/OFF` [`OFF` by default] |  Enables/disables building unit tests as a CMake target
`WITH_MATLAB=ON/OFF` | `ON/OFF` [`OFF` by default] | Enables (if available)/disables interoperability with MATLAB, allowing to detect MATLAB installation automatically. If enabled, allows applications to save/load Quantum++ matrices and vectors to/from MATLAB.
`CMAKE_INSTALL_PREFIX` | `/path/to/install` | Installs Quantum++ header files in a non-standard location (e.g., due to lack of admin. rights)

If `WITH_MATLAB=ON` and the system could not detect your MATLAB installation, you
can manually specify the path to MATLAB's installation directory via the
additional argument

    MATLAB_INSTALL_DIR=/path/to/MATLAB

If you still receive errors, you can manually specify the path to MATLAB's
required libraries and header files via the additional arguments

    MATLAB_LIB_DIR=/path/to/MATLAB/libs
    MATLAB_INCLUDE_DIR=/path/to/MATLAB/headers

## Installing Quantum++ into the default/custom target directory

To install Quantum++ into a target directory (after [Configuring the system](#Configuring-the-system)) follow the platform-specific instructions below.

### UNIX/UNIX-like
From the same out-of-source build directory execute

```bash
sudo make install 
```	

or 

```bash
sudo cmake --build . --target install
```

The commands above will install Quantum++ in: `/usr/local/include/qpp` and `/usr/local/lib/cmake/qpp`.

### Windows 
From the same out-of-source build directory execute under an Administrator Command Prompt

	cmake --build . --target INSTALL 

The commands above will install Quantum++ in: `C:\Program Files (x86)\qpp`

### FreeBSD

We are proud to be part of the [FreeBSD](https://www.freebsd.org/) operating system as an official package. If you are running FreeBSD, you can skip the CMake configuration described in [Configuring the system](#Configuring-the-system) above and install Quantum++ with the one-liner

    pkg install quantum++
    
Note that the FreeBSD version of Quantum++ may not be the latest release.

## Building and running a standalone application that uses Quantum++

Below is a minimal `CMakeLists.txt` of a standalone application that uses
Quantum++. For simplicity, we assume that the whole application is located in a
single file `src/main.cpp`, and the `CMakeLists.txt` is located in the project's
root directory.

```cmake
cmake_minimum_required(VERSION 3.12)
project(my_qpp_app)
set(CMAKE_CXX_STANDARD 11) # can set the C++ standard to C++17 or later as well

# Disable spurious Eigen warnings if using C++17 with MSVC (warning STL4007)
if (MSVC)
	add_compile_definitions(_SILENCE_CXX17_ADAPTOR_TYPEDEFS_DEPRECATION_WARNING)
endif ()

# The line below is optional (uncomment if Quantum++ was installed in a
# non-standard location)

# set(CMAKE_PREFIX_PATH "/path/to/installed/qpp")

find_package(qpp REQUIRED)
add_executable(my_qpp_app src/main.cpp)
target_link_libraries(my_qpp_app PUBLIC ${QPP_LINK_DEPS} libqpp)
```

**Do not forget to ALWAYS add `${QPP_LINK_DEPS}` (verbatim) to `target_link_libraries()`!** (last line of the `CMakeLists.txt` file above.

To build the application, first create an out-of-source build directory, e.g., under UNIX/UNIX-like

```bash
mkdir build && cd build
```

followed by configuring and building the application

```bash
cmake .. [optional arguments] && make -j8
```

where the optional arguments were described in the [Configuring the system](#Configuring-the-system) section. The commands above will produce the `my_qpp_app` executable inside the `build` folder.


## Building and running a standalone application that uses Quantum++ without a build system

Quantum++ is a header-only library. So, technically you can build an application that uses Quantum++ without using a building system, by simply using the compiler and specifying the location to all required dependencies, like below (assumes UNIX/UNIX-like, adapt accordingly for Windows)

```bash
c++ -pedantic -std=c++11 -Wall -Wextra -Weffc++ -fopenmp \
    -O3 -DNDEBUG -DEIGEN_NO_DEBUG \
    -isystem $HOME/eigen3 -I $HOME/qpp/include \
     src/main.cpp -o my_qpp_app
```

If you intend to go this route, we assume that you are familiar with how compilers work, so we won't add any more explanations to what the line above does. We still highly recommend to compile and build your applications using a modern build system such as [CMake](http://www.cmake.org/).


## Additional platform-specific instructions

#### macOS/OS X specific instructions

- We highly recommend to install [clang](http://clang.llvm.org/) version 3.8 or later via [Homebrew](https://brew.sh/) or [MacPorts](https://www.macports.org/), as the native AppleClang does not offer OpenMP support.
- In case you get any compiler or linker errors when OpenMP is enabled, you need to install the `libomp` package, e.g.

```bash
brew install libomp
```

#### MATLAB support under Windows
If building under Windows with
[MATLAB](http://www.mathworks.com/products/matlab/)
support, please add the location of `libmx.dll` and `libmat.dll`
(the `.dll` **and not** the `.lib` files) to your
`PATH` environment variable. In our case they are located under
`C:\Program Files\MATLAB\R2021a\bin\win64`.
