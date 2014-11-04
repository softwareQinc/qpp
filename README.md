quantum++
===

This is a header-only (template-based) C++11 quantum computing library, and it is work in progress.

If anyone else is interesting in contributing please let me known at vgheorgh@gmail.com. There is still a lot of work left to be done, and I can provide you with more details about what I have in mind. If you are interested in contributing, you need to have a decent knowledge of C++ (preferably C++11), including templates and STL, a basic knowledge of quantum computing and linear algebra, and some experience with Eigen linear algebra library http://eigen.tuxfamily.org/ (although this is not mandatory, but highly desirable).

The ultimate goal of this project is to have a universal quantum simulator, applicable to a vast majority of problems in quantum information/computation, and nevertheless to be as fast as possible. Right now the parallelization is done using OpenMP, but probably one can create a new branch using GPU-type multi-processing like CUDA, although this may require a significant amount of code re-writing.

And finally, the simulator should be user-friendly, easy to use for anyone with a basic knowledge of C/C++. My intention is to have it publicly available, so do not expect to make money out of it :)

---

Compiling instructions
===

== Library/header files: ==

Eigen3 library located in $HOME/eigen_3.2.2
Quantum++ library located in $HOME/qpp
MATLAB compiler include header files: /Applications/MATLAB_R2014b.app/extern/include 
MATLAB compiler shared library files: /Applications/MATLAB_R2014b.app/bin/maci64


== Release mode (with MATLAB support): ==

g++ -pedantic -std=c++11 -Wall -Wextra -Weffc++ -fopenmp -mtune=native -msse3 -O3 -DNDEBUG -DEIGEN_NO_DEBUG -isystem $HOME/eigen_3.2.2 -I$HOME/qpp/include -I/Applications/MATLAB_R2014b.app/extern/include -L/Applications/MATLAB_R2014b.app/bin/maci64 -lmx -lmat example.cpp -o qpp

== Release mode (without MATLAB support): ==

g++ -pedantic -std=c++11 -Wall -Wextra -Weffc++ -fopenmp -mtune=native -msse3 -O3 -DNDEBUG -DEIGEN_NO_DEBUG -isystem $HOME/eigen_3.2.2 -I$HOME/qpp/include example.cpp -o qpp


== Debug mode (with MATLAB support): ==

g++ -pedantic -std=c++11 -Wall -Wextra -Weffc++ -fopenmp -mtune=native -msse3 -g3 -DDEBUG -isystem $HOME/eigen_3.2.2 -I $HOME/qpp/include -I /Applications/MATLAB_R2014b.app/extern/include -L /Applications/MATLAB_R2014b.app/bin/maci64 -lmx -lmat example.cpp -o qpp

== Debug mode (without MATLAB support): ==

g++ -pedantic -std=c++11 -Wall -Wextra -Weffc++ -fopenmp -mtune=native -msse3 -g3 -DDEBUG -isystem $HOME/eigen_3.2.2 -I $HOME/qpp/include example.cpp -o qpp