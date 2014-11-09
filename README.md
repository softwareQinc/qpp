Quantum++
===

Quantum++ is a template-based header-only C++11 quantum computing library, 
developed using Eigen3 linear algebra library http://eigen.tuxfamily.org/. 
Copyright (c) 2013 - 2014 Vlad Gheorghiu, vgheorgh AT gmail DOT com.

If anyone else is interesting in contributing please let me known. 
There is still work left to be done, and I can provide you with more details about what I have in mind. 
To contribute, you need to have a decent knowledge of C++ (preferably C++11), 
including templates and STL, a basic knowledge of quantum computing and linear algebra, 
and some working experience with Eigen3.

The ultimate goal of this project is to build a universal quantum simulator, 
applicable to a vast majority of problems in quantum information/computation. 
The simulator should be fast but nevertheless user-friendly for anyone with a basic knowledge of C/C++. 

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

## Compiling/linking instructions for a simple example using g++

- Example file: `example.cpp` (current folder, from where g++ is run)
- Output binary: `qpp` (current folder, from where g++ is run)

### Configuration:

- Eigen3 library located in `$HOME/eigen_3.2.2`
- Quantum++ library located in `$HOME/qpp`
- MATLAB compiler include header files: `/Applications/MATLAB_R2014b.app/extern/include`
- MATLAB compiler shared library files: `/Applications/MATLAB_R2014b.app/bin/maci64`

### Release mode (without MATLAB support): 

	g++ -pedantic -std=c++11 -Wall -Wextra -Weffc++ -fopenmp -mtune=native -msse3 -O3 -DNDEBUG -DEIGEN_NO_DEBUG -isystem $HOME/eigen_3.2.2 -I$HOME/qpp/include example.cpp -o qpp

### Debug mode (without MATLAB support): 

	g++ -pedantic -std=c++11 -Wall -Wextra -Weffc++ -fopenmp -mtune=native -msse3 -g3 -DDEBUG -isystem $HOME/eigen_3.2.2 -I $HOME/qpp/include example.cpp -o qpp

### Release mode (with MATLAB support): 

	g++ -pedantic -std=c++11 -Wall -Wextra -Weffc++ -fopenmp -mtune=native -msse3 -O3 -DNDEBUG -DEIGEN_NO_DEBUG -isystem $HOME/eigen_3.2.2 -I$HOME/qpp/include -I/Applications/MATLAB_R2014b.app/extern/include -L/Applications/MATLAB_R2014b.app/bin/maci64 -lmx -lmat example.cpp -o qpp

### Debug mode (with MATLAB support): 

	g++ -pedantic -std=c++11 -Wall -Wextra -Weffc++ -fopenmp -mtune=native -msse3 -g3 -DDEBUG -isystem $HOME/eigen_3.2.2 -I $HOME/qpp/include -I /Applications/MATLAB_R2014b.app/extern/include -L /Applications/MATLAB_R2014b.app/bin/maci64 -lmx -lmat example.cpp -o qpp

Note: if you want to run the program on OS X with MATLAB support, make sure that the environment variable 
`DYLD_LIBRARY_PATH` is set to point to the MATLAB compiler library location, see the `run_OSX_MATLAB` script. 
Example of run command:
	
	export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:"/Applications/MATLAB_R2014b.app/bin/maci64"; ./qpp

Otherwise, you will get a runtime error like `dyld: Library not loaded: @rpath/libmat.dylib`

The current version of the repository has a `CMakeLists.txt` configuration file for building using `cmake`
 (`cmake` needs to be installed). To build the project using `cmake`, I recommend an out-of-source build, 
 i.e., from the root of the project (where `./include` is located), type

    mkdir ./build
    cd ./build
    cmake ..
    make

The above commands create the relase-mode (default) executable `qpp`, 
without MATLAB support (default), inside the directory `./build`. 
To build a different configuration, e.g. debug-mode with MATLAB support, 
type from inside the `./build` directory

    rm -rf *
    cmake -DCMAKE_BUILD_TYPE=Debug -DWITH_MATLAB=ON
    make

To change the location of `Eigen3` library or the location of MATLAB installation, 
edit the `CMakeLists.txt` file. 
Do not forget to remove everything from the `./build` directory before a fresh build!


For convenience, the current version of the repository has also a set of Makefiles available 
under the folder`./Makefile.examples`. To build the executable, copy into the root of the project 
an appropriate Makefile from `./Makefile.examples`, and name it `Makefile` , 
then type `make` (release mode) or `make debug` (debug mode) to produce the executable `qpp`.


`Makefiles.examples` and `CMakeLists.txt` are provided for convenience.
The final version of the library will consist only of header files, 
and it is the user's responsability to create an appropriate Makefile system that suits her/his needs.
