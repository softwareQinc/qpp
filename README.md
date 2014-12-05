# Quantum++ 

Quantum++ is a template-based header-only C++11 quantum computing library, 
developed using [Eigen 3](http://eigen.tuxfamily.org) linear algebra library.
Copyright (c) 2013 - 2014 Vlad Gheorghiu, vgheorgh AT gmail DOT com.

If you are interesting in contributing please let me know. 
There is still work left to be done, and I can provide you with more details 
about what I have in mind. To contribute, you need to have a decent knowledge 
of C++ (preferably C++11), including templates and the standard library, 
a basic knowledge of quantum computing and linear algebra, 
and some working experience with [Eigen 3](http://eigen.tuxfamily.org).

The ultimate goal of this project is to build a universal quantum simulator, 
applicable to a vast majority of problems in quantum information/computation.
The simulator should be fast but nevertheless user-friendly 
for anyone with a basic knowledge of C/C++. 

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

### Configuration:

- Compiler: `g++` >= 4.8 (for good C++11 support)
- [Eigen 3](http://eigen.tuxfamily.org) library located in `$HOME/eigen`
- Quantum++ library located in `$HOME/qpp`
- MATLAB compiler include header files:
`/Applications/MATLAB_R2014b.app/extern/include`
- MATLAB compiler shared library files:
`/Applications/MATLAB_R2014b.app/bin/maci64`


### Building without a build system

- Example file: `$HOME/qpp/examples/example.cpp`
- Output executable: `$HOME/qpp/examples/example`
- Must run the commands below from inside the directory `$HOME/qpp/examples` 

#### Release version (without MATLAB support): 

	g++ -pedantic -std=c++11 -Wall -Wextra -Weffc++ -fopenmp \
         -O3 -DNDEBUG -DEIGEN_NO_DEBUG \
         -isystem $HOME/eigen -I $HOME/qpp/include \
         example.cpp -o example

#### Debug version (without MATLAB support): 

	g++ -pedantic -std=c++11 -Wall -Wextra -Weffc++ -fopenmp \
         -g3 -DDEBUG \
         -isystem $HOME/eigen -I $HOME/qpp/include \
          example.cpp -o example

#### Release version (with MATLAB support): 

	g++ -pedantic -std=c++11 -Wall -Wextra -Weffc++ -fopenmp \
         -O3 -DNDEBUG -DEIGEN_NO_DEBUG \
         -isystem $HOME/eigen -I $HOME/qpp/include \
         -I/Applications/MATLAB_R2014b.app/extern/include \
         -L/Applications/MATLAB_R2014b.app/bin/maci64 \
         -lmx -lmat example.cpp -o example

#### Debug version (with MATLAB support): 

	g++ -pedantic -std=c++11 -Wall -Wextra -Weffc++ -fopenmp \
         -g3 -DDEBUG \
         -isystem $HOME/eigen -I $HOME/qpp/include \
         -I /Applications/MATLAB_R2014b.app/extern/include \
         -L /Applications/MATLAB_R2014b.app/bin/maci64 \
         -lmx -lmat example.cpp -o example


### Building using `cmake`

The current version of the repository has a `CMakeLists.txt` configuration file 
for building examples using `cmake` (`cmake` needs to be installed). To build an 
example using `cmake`, I recommend an out-of-source build, 
 i.e., from the root of the project (where `./include` is located), type

    mkdir ./build
    cd ./build
    cmake ..
    make

The above commands build the relase version (default) executable `qpp`, 
from the source file `./examples/example.cpp`,
without MATLAB support (default), inside the directory `./build`. 
To build a different configuration, e.g. debug version with MATLAB support, 
type from the root of the project

    cd ./build
    rm -rf *
    cmake -DCMAKE_BUILD_TYPE=Debug -DWITH_MATLAB=ON ..
    make
    
Or, to disable OpenMP support (enabled by default), type
   
    cd ./build
    rm -rf *
    cmake -DWITH_OPENMP=OFF ..
    make

To change the name of the example file, the location of the Eigen 3
library or the location of MATLAB installation, 
edit the `CMakeLists.txt` file. See also `CMakeLists.txt` 
for additional options. 
Do not forget to remove everything from the `./build` directory 
before a fresh build!


### Additional remarks

- The C++ compiler must be C++11 compliant.

- If your compiler does not support OpenMP 
(as it is the case e.g with `clang++`), disable OpenMP it in your build, 
as otherwise the linker may not find the `gomp` library.

- If you run the program on OS X with MATLAB support, make sure that 
the environment variable `DYLD_LIBRARY_PATH` is set to point to the MATLAB 
compiler library location, see the `run_OSX_MATLAB` script. 
Otherwise, you will get a runtime error like 
`dyld: Library not loaded: @rpath/libmat.dylib`.

    * I recommend running via a script, as otherwise setting the 
    `DYLD_LIBRARY_PATH` globally may interfere with Macports' `cmake` 
    installation (in case you use `cmake` from `macports`). If you use a 
    script, then the environment variable is local to the script and does not
    interfere with the rest of the system.

    * Example of running script, run from inside the directory where 
    the executable `qpp` is located:
	    
            #!/bin/sh # Run Quantum++ under OS X with MATLAB support
            
            export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:"/Applications/MATLAB_R2014b.app/bin/maci64"
            ./qpp

- If you build a debug version with `g++` under OS X and use `gdb` to step 
inside template functions you may want to add `-fno-weak` compiler flag. See 
<http://stackoverflow.com/questions/23330641/gnu-gdb-can-not-step-into-template-functions-os-x-mavericks>
for more details about this problem.

