# Installation instructions

Quantum++ is a header-only library that uses [CMake](https://cmake.org/) as
its build/install system. Quantum++ is platform-independent,
supporting [UNIX](https://www.opengroup.org/membership/forums/platform/unix)
(including
[macOS](https://www.apple.com/macos/)) and UNIX-like operating systems
(e.g., [Linux](https://www.linux.org)), as well
as [Windows](https://www.microsoft.com/en-us/windows).

## Pre-requisites

- [CMake](https://cmake.org/)
- [Eigen 3](https://eigen.tuxfamily.org) linear algebra library
    - Preferably install Eigen3 with via CMake or with a package manager,
      e.g. `sudo apt install libeigen3-dev` to install on Ubuntu/Debian Linux,
      so it is visible system-wide
    - **Important**: If, when building with Quantum++, your system is unable to
      detect the location of the Eigen3 matrix library, set the environment
      variable `EIGEN3_INSTALL_DIR` to point to the location of the Eigen3
      library (include the `include/eigen3` part of the path), or pass the
      argument `-DEIGEN3_INSTALL_DIR=/path/to/eigen3` to CMake
- C++17 compliant compiler, e.g., [gcc](https://gcc.gnu.org/)
  , [clang](https://clang.llvm.org)
  , [MSVC](https://visualstudio.microsoft.com/vs/) etc.

## Optional

- [Python 3](https://www.python.org/) for running the `pyqpp` Python 3 wrapper
- [MATLAB](https://www.mathworks.com/products/matlab/) compiler shared libraries
  and include header files, in case you want to enable interoperability with
  MATLAB. If enabled, allows applications build with Quantum++ to save/load
  Quantum++ matrices and vectors to/from MATLAB. The locations of the MATLAB
  compiler shared libraries and header files are platform-specific, e.g., under
  Linux, they may be located under `/usr/local/MATLAB/R2021a/bin/glnxa64` and
  `/usr/local/MATLAB/R2021a/extern/include`, respectively. On your platform the
  locations may of course differ.

## Configuring the system

First configure the system via CMake to use an out-of-source build directory 
(e.g., `./build`) by executing in a terminal/console/command prompt

```bash
cmake -B build
```

## Building the examples and/or unit_tests

To build
the [examples](https://github.com/softwareQinc/qpp/tree/main/examples) and/or 
the [unit tests](https://github.com/softwareQinc/qpp/tree/main/unit_tests), you 
need to pass the additional 
[optional flags](#cmake-optional-arguments-and-flags)  
`WITH_EXAMPLES=ON` and/or `WITH_UNIT_TESTS=ON` to CMake, e.g.,

```bash
cmake -B build -DWITH_EXAMPLES=ON -DWITH_UNIT_TESTS=ON
```

followed by the build command 

```bash
cmake --build build --target=examples --target=unit_tests --parallel 8
```

The above command builds all examples as executables in `./build`, 
followed by building the unit tests executable in 
`./build/unit_tests/unit_tests`. The `--parallel 8` instructs CMake to build
in parallel using 8 threads. 

Tu run the unit tests, execute

```bash
ctest --test-dir build
```

To build **only** a specific target, execute, e.g.,

```bash
cmake --build build --target=bb84
```

The command above only builds the example 
[examples/bb84](https://github.com/softwareQinc/qpp/tree/main/examples/bb84) 
and output the executable in `./build/bb84`. Reminder: for this to work, 
do not forget to configure the system with the `-DWITH_EXAMPLES=ON` flag.

## CMake optional arguments and flags

| Optional argument      | Value                                   | Description                                                                                                                                                                                                                      |
|------------------------|-----------------------------------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `CMAKE_INSTALL_PREFIX` | `/path/to/install`                      | Installs Quantum++ header files in a non-standard location (e.g., due to lack of admin. rights)                                                                                                                                  |
| `EIGEN3_INSTALL_DIR`   | `/path/to/eigen3`                       | Path to Eigen3 installation, if not automatically detected. This path can alternatively be enforced by setting the environment variable with the same name, e.g., via `export EIGEN3_INSTALL_DIR=/path/to/eigen3` in UNIX/Linux. |
| `SANITIZE`             | `ON/OFF` [`OFF` by default]             | Enable code sanitizing (only for gcc/clang)                                                                                                                                                                                      |
| `USE_OPENQASM2_SPECS`  | `ON/OFF` [`OFF` by default]             | Enables/disables using the OpenQASM 2.0 standard instead of Qiskit specifications; see [`DISCREPANCIES.md`](https://github.com/softwareQinc/qpp/blob/main/DISCREPANCIES.md)                                                      |
| `WITH_EXAMPLES`        | `ON/OFF` [`OFF` by default]             | Enables/disables examples as a CMake build target                                                                                                                                                                                |
| `WITH_UNIT_TESTS`      | `ON/OFF` [`OFF` by default]             | Enables/disables unit tests as a CMake build target                                                                                                                                                                              |
| `WITH_OPENMP`          | `ON/OFF` [`ON` by default]              | Enables (if available)/disables OpenMP multi-processing library                                                                                                                                                                  |
| `WITH_MATLAB=ON/OFF`   | `ON/OFF` [`OFF` by default]             | Enables (if available)/disables interoperability with MATLAB, allowing to detect MATLAB installation automatically. If enabled, allows applications to save/load Quantum++ matrices and vectors to/from MATLAB.                  |
| `TYPE_IDX`             | `default`, etc. [`default` by default]  | Integer index type (`qpp::idx`)                                                                                                                                                                                                  |
| `TPYE_BIGINT`          | `default`, etc. [`default` by default]  | Signed big integer type (`qpp::bigint`)                                                                                                                                                                                          |
| `TPYE_FP`              | `default`, etc. [`default` by default]  | Floating-point type (`qpp::realT`)                                                                                                                                                                                               |

If `WITH_MATLAB=ON` and the system could not detect your MATLAB installation,
you can manually specify the path to MATLAB's installation directory via the
additional argument

    MATLAB_INSTALL_DIR=/path/to/MATLAB

If you are still receiving errors, you can manually specify the path to MATLAB's
required libraries and header files via the additional arguments

    MATLAB_LIB_DIR=/path/to/MATLAB/libs
    MATLAB_INCLUDE_DIR=/path/to/MATLAB/headers

## Installing Quantum++

To install Quantum++ (after [Configuring the system](#configuring-the-system)), 
execute in a terminal/console (UNIX/UNIX-like systems)

```bash
sudo cmake --build build --target install
```

which installs Quantum++ in `/usr/local/include/qpp`.

### Windows

To install Quantum++ on a Windows platform, execute in an
Administrator Developer Command Prompt

```shell
cmake --build build --target INSTALL
```

which installs Quantum++ in `C:\Program Files (x86)\qpp`.

### FreeBSD

We are proud to be part of the [FreeBSD](https://www.freebsd.org/) operating
system as an official package. If you are running FreeBSD, you can install 
Quantum++ with

    pkg install quantum++

### macOS

If you are running macOS, you can install Quantum++ via 
[Homebrew](https://brew.sh) with

    brew install quantum++

## Building and running a standalone application that uses Quantum++

Below is a minimal `CMakeLists.txt` of a standalone application that uses
Quantum++. For simplicity, we assume that the whole application is located in a
single file `src/main.cpp`, and the `CMakeLists.txt` is located in the project's
root directory.

```cmake
cmake_minimum_required(VERSION 3.15)
project(standalone)
set(CMAKE_CXX_STANDARD 17)

# If the Quantum++ installation path was non-standard, i.e., specified by
#
# cmake -B build -DCMAKE_INSTALL_PREFIX=/path/to/installed/qpp
#
# then uncomment the following line and replace the installation path with yours

# set(CMAKE_PREFIX_PATH "/path/to/installed/qpp")

find_package(qpp REQUIRED)
add_executable(standalone src/main.cpp)
target_link_libraries(standalone PUBLIC ${QPP_LINK_DEPS} libqpp)

```

**Do not forget to ALWAYS add `${QPP_LINK_DEPS}` (verbatim)
to `target_link_libraries()`!** (last line of the
`CMakeLists.txt` file above).

Configure the application in an out-of-source directory by executing

```bash
cmake -B build 
```

followed by building the application with

```bash
cmake --build build
```

The commands above builds the `standalone` executable inside the `build` 
directory.

## Building and running a standalone application that uses Quantum++ without a build system

Quantum++ is a header-only library. Hence, you can technically build an 
application that uses Quantum++ without using a building system, by simply using 
the compiler and specifying the location to all required dependencies, like 
below (assumes UNIX/UNIX-like, adapt accordingly for Windows)

```bash
c++ -pedantic -std=c++17 -Wall -Wextra -Weffc++ -fopenmp \
    -O3 -DNDEBUG -DEIGEN_NO_DEBUG \
    -isystem $HOME/eigen3 -I $HOME/qpp/include -I $HOME/qpp/qasmtools/include \
     src/main.cpp -o my_qpp_app
```

If you intend to go via this route, we assume that you are familiar with how
compilers work, so we won't add any more explanations to what the line above
does.

## Additional platform-specific instructions

### Eigen 3 installation under Windows

- We **strongly** recommend installing [Eigen3](https://eigen.tuxfamily.org)
  using the [CMake](https://cmake.org) system, according to the installation
  instructions file INSTALL from the [Eigen3](https://eigen.tuxfamily.org) root
  directory (which you obtain after un-zipping the Eigen distribution archive).
  For MSVC, this translates into downloading the Eigen3 archive
  form [https://eigen.tuxfamily.org](https://eigen.tuxfamily.org), unziping it
  to e.g. `C:\path\to\eigen-3.x.x\`, followed by executing the following in an
  Administrator Developer Command Prompt

```shell
cd C:\path\to\eigen-3.x.x\
cmake -B build
cmake --build build --target INSTALL
```

Running the last line may require Administrator privileges. For other compilers,
you may need to replace the last line `cmake --build . --target INSTALL`
with `make install`.

### macOS/OS X specific instructions

- We highly recommend installing [clang](https://clang.llvm.org/) 
  via [Homebrew](https://brew.sh/), since the native AppleClang does not offer 
  OpenMP support.
- In case you get any compiler or linker errors when OpenMP is enabled, you need
  to install the `libomp` package, e.g., execute

```bash
brew install libomp
```

### MATLAB support under Windows

If building under Windows
with [MATLAB](https://www.mathworks.com/products/matlab/) support, please add the
location of
`libmx.dll` and `libmat.dll` (the `.dll` **and not** the `.lib` files) to
your `PATH` environment variable. In our case they are located
under `C:\Program Files\MATLAB\R2021a\bin\win64`.

## Python 3 wrapper

[pyqpp](https://github.com/softwareQinc/qpp/blob/main/pyqpp) is a Python 3
wrapper for Quantum++. pyqpp requires the same dependencies as Quantum++, and
can be installed using `pip`

```
pip install git+https://github.com/softwareQinc/qpp
```

**Important**: If the installation fails due to your system being unable to
detect the location of the Eigen3 matrix library, set the environment variable
`EIGEN3_INSTALL_DIR` to point to the location of the Eigen3 library
(include the `include/eigen3` part of the path).

For more details, please see
[pyqpp/README.md](https://github.com/softwareQinc/qpp/blob/main/pyqpp/README.md).
