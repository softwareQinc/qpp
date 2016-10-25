# Quantum++ 
## Version 1.0.0-devel - development

**Build status:**
Master
[![Build Status](https://api.travis-ci.org/vsoftco/qpp.svg?branch=master)](https://travis-ci.org/vsoftco/qpp)
Devel
[![Build Status](https://api.travis-ci.org/vsoftco/qpp.svg?branch=v1.0.0-devel)](https://travis-ci.org/vsoftco/qpp)

Quantum++ is a modern C++11 general purpose quantum computing library, composed 
solely of template header files. Quantum++ is written in standard C++11 and 
has very low external dependencies, using only the 
[Eigen 3](http://eigen.tuxfamily.org) linear algebra header-only template 
library and, if available, the [OpenMP](http://openmp.org/) multi-processing 
library. 

Quantum++ is not restricted to qubit systems or specific quantum 
information processing tasks, being capable of simulating arbitrary quantum 
processes. The main design factors taken in consideration were the ease of 
use, high portability, and high performance. The library's simulation
capabilities are only restricted by the amount of available physical memory. 
On a typical machine (Intel i5 8Gb RAM) Quantum++ can successfully simulate 
the evolution of 25 qubits in a pure state or of 12 qubits in a mixed state 
reasonably fast.

To report any bugs or ask for additional features/enhancements, please 
[submit an issue](https://github.com/vsoftco/qpp/issues) with an appropriate 
label.

If you are interesting in contributing to this project, please contact me. 
To contribute, you need to have a solid knowledge of C++ (preferably C++11), 
including templates and the standard library, a basic knowledge of 
quantum computing and linear algebra, and working experience with 
[Eigen 3](http://eigen.tuxfamily.org).

For additional [Eigen 3](http://eigen.tuxfamily.org) documentation 
see <http://eigen.tuxfamily.org/dox/>. For a simple 
[Eigen 3](http://eigen.tuxfamily.org) quick ASCII reference see
<http://eigen.tuxfamily.org/dox/AsciiQuickReference.txt>.

Copyright (c) 2013 - 2017 Vlad Gheorghiu, vgheorgh AT gmail DOT com.

---
Quantum++ is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Quantum++ is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Quantum++.  If not, see <http://www.gnu.org/licenses/>.

---
## Building instructions

### Configuration

- Compiler: [g++](https://gcc.gnu.org/) version 4.8.2 or later 
(for good C++11 support)
- [Eigen 3](http://eigen.tuxfamily.org) library located in `$HOME/eigen`
- Quantum++ library located in `$HOME/qpp`

##### Optional

- [MATLAB](http://www.mathworks.com/products/matlab/) compiler 
include header files:
`/Applications/MATLAB_R2016a.app/extern/include`
- [MATLAB](http://www.mathworks.com/products/matlab/) compiler 
shared library files:
`/Applications/MATLAB_R2016a.app/bin/maci64`

### Building using [CMake](http://www.cmake.org/) (version 3.0.0 or later)

The current version of the repository has a `./CMakeLists.txt` configuration 
file for building examples using [CMake](http://www.cmake.org/). 
To build an example using [CMake](http://www.cmake.org/), 
I recommend an out-of-source build, i.e., from the root of the project 
(where `./include` is located), type

    mkdir ./build
    cd ./build
    cmake ..
    make

The commands above build the relase version (default) executable `qpp`, 
from the source file `./examples/minimal.cpp`,
without [MATLAB](http://www.mathworks.com/products/matlab/) support (default), 
inside the directory `./build`. To build a different configuration, 
e.g. debug version with [MATLAB](http://www.mathworks.com/products/matlab/) 
support, type from the root of the project

    cd ./build
    rm -rf *
    cmake -DCMAKE_BUILD_TYPE=Debug -DWITH_MATLAB=ON ..
    make
    
Or, to disable [OpenMP](http://openmp.org/) support (enabled by default), type
   
    cd ./build
    rm -rf *
    cmake -DWITH_OPENMP=OFF ..
    make

To change the name of the example file, the location of the
[Eigen 3](http://eigen.tuxfamily.org)
library or the location of [MATLAB](http://www.mathworks.com/products/matlab/) 
installation, edit the `./CMakeLists.txt` file. See also `./CMakeLists.txt` 
for additional options. Do not forget to clean the `./build` directory before 
a fresh build!

### Building without an automatic build system

- Example file: `$HOME/qpp/examples/minimal.cpp`
- Output executable: `$HOME/qpp/examples/minimal`
- You must run the commands below from inside the directory 
`$HOME/qpp/examples` 

#### Release version (without [MATLAB](http://www.mathworks.com/products/matlab/) support) 

	g++ -pedantic -std=c++11 -Wall -Wextra -Weffc++ -fopenmp \
         -O3 -DNDEBUG -DEIGEN_NO_DEBUG \
         -isystem $HOME/eigen -I $HOME/qpp/include \
         minimal.cpp -o minimal

#### Debug version (without [MATLAB](http://www.mathworks.com/products/matlab/) support)

	g++ -pedantic -std=c++11 -Wall -Wextra -Weffc++ -fopenmp \
         -g3 -DDEBUG \
         -isystem $HOME/eigen -I $HOME/qpp/include \
          minimal.cpp -o minimal

#### Release version (with [MATLAB](http://www.mathworks.com/products/matlab/) support)

	g++ -pedantic -std=c++11 -Wall -Wextra -Weffc++ -fopenmp \
         -O3 -DNDEBUG -DEIGEN_NO_DEBUG \
         -isystem $HOME/eigen -I $HOME/qpp/include \
         -I/Applications/MATLAB_R2016a.app/extern/include \
         -L/Applications/MATLAB_R2016a.app/bin/maci64 \
         -lmx -lmat minimal.cpp -o minimal

#### Debug version (with [MATLAB](http://www.mathworks.com/products/matlab/) support)

	g++ -pedantic -std=c++11 -Wall -Wextra -Weffc++ -fopenmp \
         -g3 -DDEBUG \
         -isystem $HOME/eigen -I $HOME/qpp/include \
         -I /Applications/MATLAB_R2016a.app/extern/include \
         -L /Applications/MATLAB_R2016a.app/bin/maci64 \
         -lmx -lmat minimal.cpp -o minimal

### Unit testing

Quantum++ was extensively tested via a suite of unit tests constructed with
[Google Test 1.8.0](https://github.com/google/googletest) (included with the 
project in `./unit_tests/lib/gtest-1.8.0`). The source code of the unit tests 
is provided under `./unit_tests/tests`. To build and run the unit tests, I 
strongly recommend to use [CMake](http://www.cmake.org/) version 3.0.0 or 
later. Assuming you do use [CMake](http://www.cmake.org/), switch to the  
`./unit_tests` directory, create a `build` directory inside it, then from the 
newly created `./unit_tests/build` type

    cmake ..
    make
    
The commands above build `./unit_tests/build/tests/qpp_testing`, which you 
then may run.

##### Note

The [CMake](http://www.cmake.org/) configuration file 
`./unit_tests/CMakeLists.txt` defines the same building options and default 
choices as the main `./CMakeLists.txt` of Quantum++.  Therefore you can use the 
same flags as the ones mentioned at the beginning of this document when 
customizing the build. You should modify `./unit_tests/CMakeLists.txt` 
accordingly in case your [Eigen 3](http://eigen.tuxfamily.org) library or 
[MATLAB](http://www.mathworks.com/products/matlab/) include/library files are 
in a different location than the one assumed in this document.

### Additional remarks

- The C++ compiler must be fully standard-C++11 compliant.

- If using [Windows](http://windows.microsoft.com/), I recommend compiling 
under [cygwin](https://www.cygwin.com) via [CMake](http://www.cmake.org/)
and [g++](https://gcc.gnu.org/). See also 
<http://stackoverflow.com/questions/28997206/cygwin-support-for-c11-in-g4-9-2>
for a bug related to lack of support for some C++11 math functions, and
how to fix it. Quick fix: patch the standard library header file `<cmath>` 
using the provided patch `./cmath_cygwin.patch`.

- In case you use [OS X](http://www.apple.com/osx) and want to install
[clang++](http://clang.llvm.org/) version 3.7 or later, I highly recommend 
to install it via [macports](https://www.macports.org/). 

- If you use [clang++](http://clang.llvm.org/) version 3.7 or later and want 
to use [OpenMP](http://openmp.org/) (enabled by default), make sure to modify 
`CLANG_LIBOMP` and `CLANG_LIBOMP_INCLUDE` in `CMakeLists.txt` so they point to 
the correct location of the [OpenMP](http://openmp.org/) library, as otherwise 
[clang++](http://clang.llvm.org/) will not find `<omp.h>` and the `libomp` 
shared library. 

- If you run the program on [OS X](http://www.apple.com/osx) with 
[MATLAB](http://www.mathworks.com/products/matlab/) support, make sure that 
the environment variable `DYLD_LIBRARY_PATH` is set to point to the 
[MATLAB](http://www.mathworks.com/products/matlab/) 
compiler library location, see the `run_OSX_MATLAB` script. 
Otherwise, you get a runtime error similar to  

    > dyld: Library not loaded: @rpath/libmat.dylib.
    
   * I recommend running via a script, as otherwise setting the 
    `DYLD_LIBRARY_PATH` globally may interfere with 
    [macports](https://www.macports.org/)' [CMake](http://www.cmake.org/) 
    installation (in case you use [CMake](http://www.cmake.org/) from 
    [macports](https://www.macports.org/)). If you use a script, 
    then the environment variable is local to the script and 
    does not interfere with the rest of the system.

   * Example of script, assumed to be located in the root directory 
    of Quantum++
        
            #!/bin/sh
            
            MATLAB=/Applications/MATLAB_R2016a.app
            export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:$MATLAB/bin/maci64
            
            ./build/qpp

- If you build a debug version with [g++](https://gcc.gnu.org/) under 
[OS X](http://www.apple.com/osx) and use 
[gdb](http://www.gnu.org/software/gdb/) to step inside template functions 
you may want to add `-fno-weak` compiler flag. See 
<http://stackoverflow.com/questions/23330641/gnu-gdb-can-not-step-into-template-functions-os-x-mavericks>
for more details about this problem.
