# Installation instructions

**Quantum++** is a header-only library that uses [CMake](https://cmake.org/) as
its build/install system. **Quantum++** is platform-independent, supporting
[UNIX](https://www.opengroup.org/membership/forums/platform/unix) (including
[macOS](https://www.apple.com/macos/)) and UNIX-like operating systems (e.g.,
[Linux](https://www.linux.org)), as well as
[Windows](https://www.microsoft.com/en-us/windows).

---

## Pre-requisites

- C++17 compliant compiler, e.g., [GCC](https://gcc.gnu.org/)
  , [Clang](https://clang.llvm.org)
  , [MSVC](https://visualstudio.microsoft.com/vs/) etc.
- [CMake](https://cmake.org/)
- [Eigen 3](https://eigen.tuxfamily.org) linear algebra library. If missing, it
  will be installed automatically by CMake as a build dependency.

### Optional

- [Python 3](https://www.python.org/) for building and running the **pyqpp**
  Python 3 wrapper
- [MATLAB](https://www.mathworks.com/products/matlab/) compiler shared
  libraries and include header files, in case you want to enable
  interoperability with MATLAB. If enabled, allows applications built with
  **Quantum++** to save/load **Quantum++** matrices and vectors to/from MATLAB.

---

## Configuring the system

First configure the system via CMake to use an out-of-source build directory
(e.g., `./build`) by executing (in a terminal/console/command prompt) under the
project's root directory

```shell
cmake -B build
```

---

## Building the examples and/or unit tests

To build the
[examples](https://github.com/softwareQinc/qpp/tree/main/examples), execute

```shell
cmake --build build --target examples --parallel 8
```

The above command builds all examples as executables in `./build`. The
`--parallel 8` flag instructs CMake to build in parallel using 8 threads,
modify accordingly.

To build the
[unit tests](https://github.com/softwareQinc/qpp/tree/main/unit_tests), execute

```shell
cmake --build build/unit_tests --target unit_tests --parallel 8
```

Tu run the unit tests, execute

```shell
ctest --test-dir build
```

To build **only** a specific target, execute, e.g.,

```shell
cmake --build build --target bb84
```

The command above builds only the example
[examples/bb84.cpp](https://github.com/softwareQinc/qpp/tree/main/examples/bb84.cpp)
and outputs the executable `./build/bb84[.exe]`.

---

## CMake optional arguments and flags

> **Note:**
> All CMake flags below **do not propagate** to projects that use **Quantum++**
> via `find_package(qpp ...)`.
>
> Each consumer project must explicitly define these flags in its own
> `CMakeLists.txt` if needed.

| Optional argument       | Value                                  | Description                                                                                                                                                                                                         |
| ----------------------- | -------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `CMAKE_INSTALL_PREFIX`  | `/path/to/install`                     | Specifies a custom installation directory for **Quantum++** header files -- useful when you lack administrative privileges or want a non-default install location                                                   |
|                         |                                        |                                                                                                                                                                                                                     |
| `QPP_MATLAB`            | `ON/OFF` [`OFF` by default]            | Enables (if available)/disables interoperability with MATLAB, allowing to detect MATLAB installation automatically. If enabled, allows applications to save/load **Quantum++** matrices and vectors to/from MATLAB. |
| `QPP_OPENMP`            | `ON/OFF` [`ON` by default]             | Enables (if available)/disables OpenMP multi-processing library                                                                                                                                                     |
| `QASMTOOLS_QASM2_SPECS` | `ON/OFF` [`OFF` by default]            | Enables/disables using the OpenQASM 2.0 standard instead of Qiskit specifications -- see [`DISCREPANCIES.md`](https://github.com/softwareQinc/qpp/blob/main/DISCREPANCIES.md)                                       |
|                         |                                        |                                                                                                                                                                                                                     |
| `QPP_SANITIZE`          | `ON/OFF` [`OFF` by default]            | Enable code sanitizing                                                                                                                                                                                              |
|                         |                                        |                                                                                                                                                                                                                     |
| `QPP_BIGINT`            | `default`, etc. [`default` by default] | Signed big integer type (`qpp::bigint`)                                                                                                                                                                             |
| `QPP_FP`                | `default`, etc. [`default` by default] | Floating-point type (`qpp::realT`)                                                                                                                                                                                  |
| `QPP_IDX`               | `default`, etc. [`default` by default] | Integer index type (`qpp::idx`)                                                                                                                                                                                     |

---

## Installing Quantum++

### UNIX/UNIX-like/Windows

To install **Quantum++** (after
[Configuring the system](#configuring-the-system)), execute in a
terminal/console (UNIX/UNIX-like systems)

```shell
sudo cmake --build build --target install
```

or in an Administrator Command Prompt (Windows)

```shell
cmake --build build --target install
```

The above commands install **Quantum++** in `/usr/local/include/qpp` on
UNIX/UNIX-like platforms, and in `C:\Program Files (x86)\qpp` on Windows
platforms.

To uninstall, execute in a terminal/console (UNIX/UNIX-like systems)

```shell
sudo cmake --build build --target uninstall
```

or in an Administrator Command Prompt (Windows)

```shell
cmake --build build --target uninstall
```

### FreeBSD

We are proud to be part of the [FreeBSD](https://www.freebsd.org/) operating
system as an official package. If you are running FreeBSD, you can install
**Quantum++** with

```shell
sudo pkg install quantum++
```

and uninstall it with

```shell
sudo pkg remove quantum++
```

### macOS/Linux

If you are running macOS or Linux, you can install **Quantum++** via
[Homebrew](https://brew.sh) with

```shell
brew install quantum++
```

and uninstall it with

```shell
brew uninstall quantum++
```

---

## Building and running a standalone application that uses Quantum++

Below is a minimal `CMakeLists.txt` of a standalone application that uses
**Quantum++**. For simplicity, we assume that the whole application is located
in a single file `src/main.cpp`, and the `CMakeLists.txt` is located in the
project's root directory.

```cmake
cmake_minimum_required(VERSION 3.20)
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
target_link_libraries(standalone PRIVATE libqpp)

```

Configure the application in an out-of-source directory by executing

```shell
cmake -B build
```

followed by building the application with

```shell
cmake --build build
```

The commands above builds the `standalone` executable inside the `build`
directory.

---

## Building and running a standalone application that uses Quantum++ without a build system

**Quantum++** is a header-only library. Hence, you can technically build an
application that uses **Quantum++** without using a building system, by simply
using the compiler and specifying the location to all required dependencies,
like below (assumes UNIX/UNIX-like, adapt accordingly for Windows)

```shell
g++ -pedantic -std=c++17 -Wall -Wextra -Weffc++ -fopenmp \
    -O3 -DNDEBUG -DEIGEN_NO_DEBUG
    -isystem $HOME/eigen3 -I $HOME/qpp/include -I $HOME/qpp/qasmtools/include \
     src/main.cpp -o my_qpp_app
```

If you intend to go via this route, we assume that you are familiar with how
compilers work, so we won't add any more explanations to what the line above
does.

---

## Additional platform-specific instructions

### macOS/OS X specific instructions

- We highly recommend installing [clang](https://clang.llvm.org/)
  via [Homebrew](https://brew.sh/), since the native AppleClang does not offer
  OpenMP support.
- In case you get any compiler or linker errors when OpenMP is enabled, you
  need to install the `libomp` package, e.g., execute

```shell
brew install libomp
```

### MATLAB support under Windows

If building under Windows
with [MATLAB](https://www.mathworks.com/products/matlab/) support, please add
the
location of
`libmx.dll` and `libmat.dll` (the `.dll` **and not** the `.lib` files) to your
`PATH` environment variable. For example, on our platform they are located
under `C:\Program Files\MATLAB\R2021a\bin\win64`.

---

## Python 3 wrapper

[**pyqpp**](https://github.com/softwareQinc/qpp/blob/main/pyqpp) is a Python 3
wrapper for **Quantum++**. **pyqpp** requires the same dependencies as
**Quantum++**, and can be installed using `pip`

```shell
pip install git+https://github.com/softwareQinc/qpp
```

For more details, please see
[pyqpp/README.md](https://github.com/softwareQinc/qpp/blob/main/pyqpp/README.md).

---

## Platform-dependent issues

### SunOS/OpenIndiana

The Python3 wrapper
[**pyqpp**](https://github.com/softwareQinc/qpp/blob/main/pyqpp) doesn't
compile under SunOS/OpenIndiana due to errors in `<cmath>` such as

```
no member named 'llround' in the global namespace; did you mean 'lround'?
```
