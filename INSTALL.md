# Installation instructions

Quantum++ is a header-only library that uses [CMake](https://cmake.org/) as
its build/install system. Quantum++ is platform-independent,
supporting [UNIX](https://www.opengroup.org/membership/forums/platform/unix)
(including
[macOS](https://www.apple.com/macos/)) and UNIX-like operating systems
(e.g., [Linux](https://www.linux.org)), as well
as [Windows](https://www.microsoft.com/en-us/windows).

## Pre-requisites

- [CMake](https://cmake.org/) version 3.12 or later
- [Eigen 3](https://eigen.tuxfamily.org) linear algebra library
    - Preferably install Eigen3 with via CMake or with a package manager,
      e.g. `sudo apt install libeigen3-dev` to install on Ubuntu/Debian Linux,
      so it is visible system-wide
    - **Important**: If, when building with Quantum++, your system is unable to
      detect the location of the Eigen3 matrix library, set the environment
      variable `EIGEN3_INSTALL_DIR` to point to the location of the Eigen3
      library (include the `include/eigen3` part of the path), or pass the
      flag `-DEIGEN3_INSTALL_DIR=/path/to/eigen3` to CMake
- A C++17 compliant compiler, e.g., [gcc](https://gcc.gnu.org/)
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

For UNIX/UNIX-like/Windows, first create an out-of-source build directory, e.g.,
from the project's root directory type in a terminal/console/command prompt

```bash
mkdir build && cd build
```

Then, configure the system by typing (from inside the newly created `build`
directory)

```bash
cmake ..
```

If you plan to compile
the [examples](https://github.com/softwareQinc/qpp/tree/main/examples) or the
[unit tests](https://github.com/softwareQinc/qpp/tree/main/unit_tests) you may
need to pass additional arguments to CMake,

```bash
cmake .. [optional arguments]
```

where [optional arguments] are passed as `-DOPTIONAL_ARGUMENT=VALUE`. The
Quantum++-specific optional arguments are:

| Optional argument      | Value                       | Description                                                                                                                                                                                                                      |
|------------------------|-----------------------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `CMAKE_INSTALL_PREFIX` | `/path/to/install`          | Installs Quantum++ header files in a non-standard location (e.g., due to lack of admin. rights)                                                                                                                                  |
| `EIGEN3_INSTALL_DIR`   | `/path/to/eigen3`           | Path to Eigen3 installation, if not automatically detected. This path can alternatively be enforced by setting the environment variable with the same name, e.g., via `export EIGEN3_INSTALL_DIR=/path/to/eigen3` in UNIX/Linux. |
| `SANITIZE`             | `ON/OFF` [`OFF` by default] | Enable code sanitizing (only for gcc/clang)                                                                                                                                                                                      |
| `USE_OPENQASM2_SPECS`  | `ON/OFF` [`OFF` by default] | Enables/disables using the OpenQASM 2.0 standard instead of Qiskit specifications; see [`DISCREPANCIES.md`](https://github.com/softwareQinc/qpp/blob/main/DISCREPANCIES.md)                                                      |
| `WITH_EXAMPLES`        | `ON/OFF` [`OFF` by default] | Enables/disables examples as a CMake build target                                                                                                                                                                                |
| `WITH_UNIT_TESTS`      | `ON/OFF` [`OFF` by default] | Enables/disables unit tests as a CMake build target                                                                                                                                                                              |
| `WITH_OPENMP`          | `ON/OFF` [`ON` by default]  | Enables (if available)/disables OpenMP multi-processing library                                                                                                                                                                  |
| `WITH_MATLAB=ON/OFF`   | `ON/OFF` [`OFF` by default] | Enables (if available)/disables interoperability with MATLAB, allowing to detect MATLAB installation automatically. If enabled, allows applications to save/load Quantum++ matrices and vectors to/from MATLAB.                  |

If `WITH_MATLAB=ON` and the system could not detect your MATLAB installation,
you can manually specify the path to MATLAB's installation directory via the
additional argument

    MATLAB_INSTALL_DIR=/path/to/MATLAB

If you still receive errors, you can manually specify the path to MATLAB's
required libraries and header files via the additional arguments

    MATLAB_LIB_DIR=/path/to/MATLAB/libs
    MATLAB_INCLUDE_DIR=/path/to/MATLAB/headers

## Installing Quantum++ into the default/custom target directory

To install Quantum++ into a target directory (
after [Configuring the system](#Configuring-the-system)) follow the
platform-specific instructions below.

### UNIX/UNIX-like

From the same out-of-source build directory execute

```bash
sudo make install 
```

or

```bash
sudo cmake --build . --target install
```

The commands above will install Quantum++ in: `/usr/local/include/qpp`
and `/usr/local/lib/cmake/qpp`.

### Windows

From the same out-of-source build directory execute under an Administrator
Command Prompt

```shell
cmake --build . --target INSTALL
```

The commands above will install Quantum++ in: `C:\Program Files (x86)\qpp`

### FreeBSD

We are proud to be part of the [FreeBSD](https://www.freebsd.org/) operating
system as an official package. If you are running FreeBSD, you can skip the
CMake configuration described
in [Configuring the system](#Configuring-the-system)
above and install Quantum++ with the one-liner

    pkg install quantum++

Note that the FreeBSD version of Quantum++ may not correspond to the latest
release.

## Building and running a standalone application that uses Quantum++

Below is a minimal `CMakeLists.txt` of a standalone application that uses
Quantum++. For simplicity, we assume that the whole application is located in a
single file `src/main.cpp`, and the `CMakeLists.txt` is located in the project's
root directory.

```cmake
cmake_minimum_required(VERSION 3.12)
project(my_qpp_app)
set(CMAKE_CXX_STANDARD 17)

# If the Quantum++ installation path was non-standard, i.e., specified by
#
# cmake .. -DCMAKE_INSTALL_PREFIX=/path/to/installed/qpp
#
# then uncomment the following line and replace the installation path with yours

# set(CMAKE_PREFIX_PATH "/path/to/installed/qpp")

find_package(qpp REQUIRED)
add_executable(my_qpp_app src/main.cpp)
target_link_libraries(my_qpp_app PUBLIC ${QPP_LINK_DEPS} libqpp)
```

**Do not forget to ALWAYS add `${QPP_LINK_DEPS}` (verbatim)
to `target_link_libraries()`!** (last line of the
`CMakeLists.txt` file above).

To build the application, first create an out-of-source build directory, e.g.,
under UNIX/UNIX-like execute

```bash
mkdir build && cd build
```

followed by configuring and building the application

```bash
cmake .. [optional arguments] && make -j8
```

where the optional arguments were described in
the [Configuring the system](#Configuring-the-system) section. The commands
above will produce the `my_qpp_app` executable inside the `build` folder.

## Building and running a standalone application that uses Quantum++ without a build system

Quantum++ is a header-only library. So, technically you can build an application
that uses Quantum++ without using a building system, by simply using the
compiler and specifying the location to all required dependencies, like below
(assumes UNIX/UNIX-like, adapt accordingly for Windows)

```bash
c++ -pedantic -std=c++17 -Wall -Wextra -Weffc++ -fopenmp \
    -O3 -DNDEBUG -DEIGEN_NO_DEBUG \
    -isystem $HOME/eigen3 -I $HOME/qpp/include -I $HOME/qpp/qasmtools/include \
     src/main.cpp -o my_qpp_app
```

If you intend to go this route, we assume that you are familiar with how
compilers work, so we won't add any more explanations to what the line above
does. We still highly recommend compiling and building your applications using a
modern build system such as [CMake](https://cmake.org/).

## Building and running examples

We assume that you have already installed Quantum++, if not please go back and
read the instructions above. We also assume that you are now familiar with the
CMake build system, so in the following we only showcase the commands, with no
further explanations. If the lines below are not clear, please read the
instructions above, or look into the
[AppVeyor](https://github.com/softwareQinc/qpp/tree/main/.appveyor) or
[CircleCI](https://github.com/softwareQinc/qpp/tree/main/.circleci) continuous
integration scripts to understand better how to compile and run under various
platforms (and using alternatives build systems such as,
e.g., [Ninja](https://ninja-build.org/)).

Finally, we assume that you type the commands below in a
terminal/console/command prompt, inside an out-of-source build directory,
e.g., `./build`.

### UNIX/UNIX-like

```bash
cmake .. -DWITH_EXAMPLES=ON [optional arguments]
make -j8 examples
./bb84 # runs, e.g., the BB84 example
```

If you want to compile a single example/subset of examples,
replace `make -j8 examples` with
`make -j8 example_target[s]`, e.g.,

```bash
make bb84 # builds the BB84 example
```

or

```bash
make -j8 bb84 grover # builds both BB84 and Grover examples
```

### Windows

```shell
cmake .. -DWITH_EXAMPLES=ON [optional arguments]
msbuild -verbosity:minimal -m:8 examples.vcxproj
.\Debug\bb84 # runs, e.g., the BB84 example
```

The examples will be built under the `./Debug` directory (default). You can also
compile and build Release versions, please go ahead and read more about CMake
under Windows from your favourite source.

If you want to compile a single example/subset of examples,
replace `msbuild -verbosity:minimal -m:8 examples.vcxproj`
with
`msbuild -verbosity:minimal -m:8 example_target[s].vcxproj`, e.g.,

```shell
msbuild -verbosity:minimal bb84.vcxproj # builds the BB84 example
```

or

```shell
msbuild -verbosity:minimal -m:8 bb84.vcxproj grover.vcxproj # builds both BB84 and Grover examples
```    

## Building and running unit tests

### UNIX/UNIX-like

```bash
cmake .. -DWITH_UNIT_TESTS=ON [optional arguments]
make -j8 unit_tests
ctest # or ./unit_tests/unit_tests
```

### Windows

```shell
msbuild -verbosity:minimal -m:8 .\unit_tests\unit_tests.vcxproj
ctest # or .\unit_tests\Debug\unit_tests.exe
```    

## Additional platform-specific instructions

#### Eigen 3 installation under Windows

- We **strongly** recommend installing [Eigen3](https://eigen.tuxfamily.org)
  using the [CMake](https://cmake.org) system, according to the installation
  instructions file INSTALL from the [Eigen3](https://eigen.tuxfamily.org) root
  directory (which you obtain after un-zipping the Eigen distribution archive).
  For MSVC, this translates into downloading the Eigen3 archive
  form [https://eigen.tuxfamily.org](https://eigen.tuxfamily.org), unziping it
  to e.g. `C:\path\to\eigen-3.x.x\`, followed by typing the following in a
  Command Prompt

```shell
cd C:\path\to\eigen-3.x.x\
mkdir build && cd build
cmake ..
cmake --build . --target INSTALL
```

Running the last line may require Administrator privileges. For other compilers,
you may need to replace the last line `cmake --build . --target INSTALL`
with `make install`.

#### macOS/OS X specific instructions

- We highly recommend installing [clang](https://clang.llvm.org/) version 3.8 or
  later via [Homebrew](https://brew.sh/)
  or [MacPorts](https://www.macports.org/), as the native AppleClang does not
  offer OpenMP support.
- In case you get any compiler or linker errors when OpenMP is enabled, you need
  to install the `libomp` package, e.g., execute

```bash
brew install libomp
```

#### MATLAB support under Windows

If building under Windows
with [MATLAB](https://www.mathworks.com/products/matlab/) support, please add the
location of
`libmx.dll` and `libmat.dll` (the `.dll` **and not** the `.lib` files) to
your `PATH` environment variable. In our case they are located
under `C:\Program Files\MATLAB\R2021a\bin\win64`.

## Python 3 wrapper

pyqpp is a Python 3 wrapper for Quantum++. pyqpp requires the same dependencies
as Quantum++, and can be installed using `pip`:

```
pip install git+https://github.com/softwareQinc/qpp
```

**Important**: If the installation fails due to your system being unable to
detect the location of the Eigen3 matrix library, set the environment variable
`EIGEN3_INSTALL_DIR` to point to the location of the Eigen3 library
(include the `include/eigen3` part of the path).

For more details, please see
[pyqpp/README.md](https://github.com/softwareQinc/qpp/blob/main/pyqpp/README.md)
.
