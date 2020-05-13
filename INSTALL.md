# Building instructions for the truly impatient

If you really cannot wait anymore, continue reading this section. If you afford to be a bit patient, please skip the reminder of this section and read the detailed 
[Building instructions for POSIX compliant systems](https://github.com/softwareqinc/qpp/wiki/2.-Building-instructions-for-POSIX-compliant-platforms) 
or 
[Building instructions for Windows](https://github.com/softwareqinc/qpp/wiki/3.-Building-instructions-for-Windows-platforms) 
in the Wiki or in the reminder of this document.

## Quick instructions for building on POSIX compliant systems

In this section I assume that you use a POSIX compliant system (e.g. UNIX-like).
To get started with [Quantum++](https://github.com/softwareqinc/qpp), first
install the [Eigen 3](http://eigen.tuxfamily.org/) library into your home directory, 
as `$HOME/eigen`. You can change the name of the directory, but in the
current document I will use `$HOME/eigen` as the location of the
[Eigen 3](http://eigen.tuxfamily.org/) library. Next, download or clone
[Quantum++](https://github.com/softwareqinc/qpp) library into the home
directory as `$HOME/qpp`. Finally, make sure that your compiler supports
C++11 and preferably [OpenMP](http://openmp.org/). As a compiler I
recommend [g++](https://gcc.gnu.org/) version 5.0 or later or
[clang++](http://clang.llvm.org) version 3.8 or later. You are now ready to go!

Next, we will build a simple minimal example to test that the installation was
successful. Create a directory called `$HOME/testing`, and inside it
create the file `minimal.cpp`, available for download in 
[`examples/minimal.cpp`](https://github.com/softwareqinc/qpp/blob/master/examples/minimal.cpp) and with 
the content listed below.

```CPP
// Minimal example
// Source: ./examples/minimal.cpp
#include <iostream>

#include "qpp.h"

int main() {
    using namespace qpp;
    std::cout << "Hello Quantum++!\nThis is the |0> state:\n";
    std::cout << disp(st.z0) << '\n';
}
```

Next, compile the file using the C++11 compliant compiler. In the
following I assume you use [g++](https://gcc.gnu.org/), but the building
instructions are similar for other compilers. From the directory
`$HOME/testing` type

```bash
g++ -std=c++11 -isystem $HOME/eigen -I $HOME/qpp/include minimal.cpp -o minimal
```

Your compile command may differ from the above, depending on the C++
compiler and operating system. If everything went fine, the above
command should build an executable `minimal` in `$HOME/testing`, which
can be run by typing `./minimal`. The output should be similar to the
following:

```
Hello Quantum++!
This is the |0> state:
1
0
```

The line

```CPP
#include "qpp.h"
```

includes the main header file of the
library `qpp.h`. This header file includes all other necessary internal
[Quantum++](https://github.com/softwareqinc/qpp) header files. The line
```CPP
using namespace qpp;
```
injects the [Quantum++](https://github.com/softwareqinc/qpp)'s namespace into
the working namespace, so we don't need to prefix the library functions by
`qpp::`. Finally the line 

```CPP
std::cout << disp(st.z0) << '\n';
```

displays the state |0> represented by the Singleton `qpp::st.z0` in a
nice format using the display manipulator `qpp::disp()`.

## Quick instructions for building on [Windows](https://www.microsoft.com/en-us/windows)

We provide [CMake](http://www.cmake.org/) support for [Visual Studio](https://www.visualstudio.com). We recommend using 
Visual Studio 2017 or later. For comprehensive details about using [CMake](http://www.cmake.org/) with 
[Visual Studio](https://www.visualstudio.com) please read 
[this page](https://docs.microsoft.com/en-us/cpp/build/cmake-projects-in-visual-studio?view=vs-2019) and see the [Wiki](https://github.com/softwareQinc/qpp/wiki/3.-Building-instructions-for-Windows-platforms).

---
# Building instructions for POSIX compliant systems

## Pre-requisites

- Recommended compiler: [g++](https://gcc.gnu.org/) version 5.0 or later
(for good C++11 support). You may also use [clang++](http://clang.llvm.org/) version 3.8 or later, please read the ["Additional remarks/Building with clang++"](https://github.com/softwareqinc/qpp/blob/master/INSTALL.md#clang) subsection near the end of this document for more platform-dependent details.
- [Eigen 3](http://eigen.tuxfamily.org) linear algebra library. I assume here that 
the library is installed in `$HOME/eigen`, although the location may vary, e.g. if 
the library was installed using a package manager.
- [Quantum++](https://github.com/softwareqinc/qpp) library located in `$HOME/qpp`.

## Optional

- [CMake](http://www.cmake.org/) version 3.0 or later, highly recommended
- [MATLAB](http://www.mathworks.com/products/matlab/) compiler include header 
files, e.g.:
`/Applications/MATLAB_R2017b.app/extern/include`
- [MATLAB](http://www.mathworks.com/products/matlab/) compiler shared library 
files, e.g.:
`/Applications/MATLAB_R2017b.app/bin/maci64`

## Building using [CMake](http://www.cmake.org/) (version 3.0 or later)

The current version of the repository has a 
[`CMakeLists.txt`](https://github.com/softwareqinc/qpp/blob/master/CMakeLists.txt)
configuration file for building examples using [CMake](http://www.cmake.org/). 
To build an(the) example(s) using [CMake](http://www.cmake.org/), 
I recommend an out-of-source build, i.e., from the root of the project 
(where the `include` folder is located), type

```bash
mkdir build
cd build
cmake ..
make [<target>]
```

If the optional `<target>` is specified, the commands above build the release version (default) executable corresponding to `<target>`, 
from the corresponding example source file `examples/<target>.cpp`,
without [MATLAB](http://www.mathworks.com/products/matlab/) support (default), 
inside the directory `build`. If `<target>` is not specified, all target examples from the directory `examples` are built.

If the location of [Eigen 3](http://eigen.tuxfamily.org) is not detected
automatically by the [CMake](http://www.cmake.org/) build script, 
then the build script will fail (with an error message). In this case 
the location of [Eigen 3](http://eigen.tuxfamily.org) needs to be specified
manually in the [CMake](http://www.cmake.org/) build command line by passing
the `-DEIGEN3_INSTALL_DIR=/path/to/eigen3` flag, e.g.

```bash
cmake .. -DEIGEN3_INSTALL_DIR=/usr/local/eigen3
```

To build a different configuration, 
e.g. the debug version with [MATLAB](http://www.mathworks.com/products/matlab/)
support, type from the root of the project

```bash
cd build
rm -rf *
cmake .. -DCMAKE_BUILD_TYPE=Debug -DWITH_MATLAB=ON
make [<target>]
```

The flag `-DWITH_MATLAB=ON` instructs [CMake](http://www.cmake.org/)  to detect the [MATLAB](http://www.mathworks.com/products/matlab/) installation automatically. If successful, it builds and links with the [MATLAB](http://www.mathworks.com/products/matlab/) compiler libraries and injects the macro `HAS_MATLAB_ENABLED` in the source code. If [MATLAB](http://www.mathworks.com/products/matlab/) was not detected automatically, specify its location manually via the additional flag `-DMATLAB_INSTALL_DIR=/path/to/MATLAB`, e.g. `-DMATLAB_INSTALL_DIR=/Applications/MATLAB_R2017b.app`. If you still get some [CMake](http://www.cmake.org/) errors, try adding the flags `-DMATLAB_INCLUDE_DIR=/path/to/MATLAB/headers` and `-DMATLAB_LIB_DIR=/path/to/MATLAB/libs` to specify the exact location of the [MATLAB](http://www.mathworks.com/products/matlab/) compiler include header files and compiler shared library files, respectively.
    
To disable [OpenMP](http://openmp.org/) support (enabled by default), type
   
```bash
cd build
rm -rf *
cmake .. -DWITH_OPENMP=OFF
make [<target>]
```

Inspect also [`CMakeLists.txt`](https://github.com/softwareqinc/qpp/blob/master/CMakeLists.txt) 
for additional fine-tuning options. Do not forget to clean the `build` 
directory before a fresh build!

## Building without an automatic build system

- Example file: `$HOME/qpp/examples/minimal.cpp`
- Output executable: `$HOME/qpp/examples/minimal`
- You must run the commands below from inside the directory 
`$HOME/qpp/examples` 

### Release version (without [MATLAB](http://www.mathworks.com/products/matlab/) support) 

```bash
g++ -pedantic -std=c++11 -Wall -Wextra -Weffc++ -fopenmp \
    -O3 -DNDEBUG -DEIGEN_NO_DEBUG \
    -isystem $HOME/eigen -I $HOME/qpp/include \
     minimal.cpp -o minimal
```

### Debug version (without [MATLAB](http://www.mathworks.com/products/matlab/) support)

```bash
g++ -pedantic -std=c++11 -Wall -Wextra -Weffc++ -fopenmp \
    -g3 -DDEBUG \
    -isystem $HOME/eigen -I $HOME/qpp/include \
     minimal.cpp -o minimal
```

### Release version (with [MATLAB](http://www.mathworks.com/products/matlab/) support)

```bash
g++ -pedantic -std=c++11 -Wall -Wextra -Weffc++ -fopenmp \
    -O3 -DNDEBUG -DEIGEN_NO_DEBUG \
    -isystem $HOME/eigen -I $HOME/qpp/include \
    -I/Applications/MATLAB_R2017b.app/extern/include \
    -L/Applications/MATLAB_R2017b.app/bin/maci64 \
    -lmx -lmat minimal.cpp -o minimal
```

### Debug version (with [MATLAB](http://www.mathworks.com/products/matlab/) support)

```bash
g++ -pedantic -std=c++11 -Wall -Wextra -Weffc++ -fopenmp \
    -g3 -DDEBUG \
    -isystem $HOME/eigen -I $HOME/qpp/include \
    -I /Applications/MATLAB_R2017b.app/extern/include \
    -L /Applications/MATLAB_R2017b.app/bin/maci64 \
    -lmx -lmat minimal.cpp -o minimal
```

## Additional remarks

### <a name="clang"></a>Building with [clang++](http://clang.llvm.org/)
- Make sure to use [clang++](http://clang.llvm.org/) version 3.8 or later (previous versions lack or do not have good [OpenMP](http://openmp.org/) support).

#### Linux specific instructions
- On Linux, you can install [clang++](http://clang.llvm.org/) version 3.8 or later via the package manager, e.g. in Debian 9 or Ubuntu 16.10 use

    ```bash
    sudo apt-get install clang-3.8
    ```

- In case you get any compiler or linker errors when 
[OpenMP](http://openmp.org/) is enabled, you need to install the `libomp-dev` package, e.g in Debian 9 or Ubuntu 16.10 use

    ```bash
    sudo apt-get install libomp-dev
    ```

#### [OS X/macOS](http://www.apple.com/osx) specific instructions

- I highly recommend to install [clang++](http://clang.llvm.org/) version 3.8 or later via [MacPorts](https://www.macports.org/).
- In case you get any compiler or linker errors when [OpenMP](http://openmp.org/) is enabled, you need to install the `libomp` package, e.g using [MacPorts](https://www.macports.org/)

    ```bash
    sudo port install libomp
    ```

### Building debug versions with [g++](https://gcc.gnu.org/) on [OS X/macOS](http://www.apple.com/osx)

- If you build a debug version with [g++](https://gcc.gnu.org/) and use 
[gdb](http://www.gnu.org/software/gdb/) to step inside template functions 
you may want to add `-fno-weak` compiler flag. See 
<http://stackoverflow.com/questions/23330641/gnu-gdb-can-not-step-into-template-functions-os-x-mavericks>
for more details about this problem.

---
# Building instructions for Windows

## Via [Visual Studio](https://www.visualstudio.com)

We provide [CMake](http://www.cmake.org/) support for [Visual Studio](https://www.visualstudio.com). We recommend using 
Visual Studio 2017 or later. For comprehensive details about using [CMake](http://www.cmake.org/) with 
[Visual Studio](https://www.visualstudio.com) please read 
[this page](https://docs.microsoft.com/en-us/cpp/build/cmake-projects-in-visual-studio?view=vs-2019) and see the [Wiki](https://github.com/softwareQinc/qpp/wiki/3.-Building-instructions-for-Windows-platforms).

## Via [Cygwin](https://www.cygwin.com) or [MSYS2](https://www.msys2.org)

- Use the [same building instructions as for POSIX systems](https://github.com/softwareqinc/qpp/wiki/2.-Building-instructions-for-POSIX-compliant-platforms).
- **Note:** some earlier versions of
[Cygwin](https://www.cygwin.com) had a bug related to lack of support for some
C++11 math functions, see
<http://stackoverflow.com/questions/28997206/cygwin-support-for-c11-in-g4-9-2>
for more details. Quick fix: patch the standard library header file `<cmath>`
using the provided patch [`cmath_cygwin.patch`](https://github.com/softwareqinc/qpp/blob/master/cmath_cygwin.patch).
Later [Cygwin](https://www.cygwin.com) versions seem to have fixed the issue (as of Nov. 2016).

---
**Important**: If building under Windows with 
[MATLAB](http://www.mathworks.com/products/matlab/)
support, please add the location of `libmx.dll` and `libmat.dll` 
(the `.dll` **and not** the `.lib` files) to your
`PATH` environment variable. In my case they are located under 
`C:\Program Files\MATLAB\R2020a\bin\win64`.
