# Quantum++ 
## Version 2.5 - 26 November 2020

**Build status:**

[![Build Status](https://travis-ci.org/softwareQinc/qpp.svg?branch=master)](https://travis-ci.org/softwareQinc/qpp)
[![Build status](https://ci.appveyor.com/api/projects/status/1k2866yffaiaapmw?svg=true)](https://ci.appveyor.com/project/vsoftco/qpp)

**Chat (questions/issues)**

[![Join the chat at https://gitter.im/vsoftco_qpp](https://badges.gitter.im/Join%20Chat.svg)](https://gitter.im/vsoftco_qpp/Lobby/?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)

---
## About

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
[submit an issue](https://github.com/softwareqinc/qpp/issues) with an appropriate 
label.

If you are interesting in contributing to this project, feel free to contact
me. Alternatively, fork the repository, add your contribution, then finally
create a pull request. If I accept the pull request, I will merge your custom
branch with the latest development branch. The latter will eventually be
merged into a future release version.  To contribute, you need to have a solid
knowledge of C++ (preferably C++11), including templates and the standard
library, a basic knowledge of quantum computing and linear algebra, and
working experience with [Eigen 3](http://eigen.tuxfamily.org).

For additional [Eigen 3](http://eigen.tuxfamily.org) documentation 
see <http://eigen.tuxfamily.org/dox/>. For a simple 
[Eigen 3](http://eigen.tuxfamily.org) quick ASCII reference see
<http://eigen.tuxfamily.org/dox/AsciiQuickReference.txt>.

Copyright (c) 2013 - 2020 Vlad Gheorghiu.

---
## License

[Quantum++](https://github.com/softwareqinc/qpp) is distributed under the MIT 
license. Please see the 
[`LICENSE`](https://github.com/softwareqinc/qpp/blob/master/LICENSE) file for more 
details.

---
## Installation instructions and further documentation

Please see the installation guide 
[`INSTALL.md`](https://github.com/softwareqinc/qpp/blob/master/INSTALL.md) 
and the comprehensive [Wiki](https://github.com/softwareqinc/qpp/wiki) for further 
documentation and detailed examples. 

To generate the full official API documentation in both PDF and HTML formats run 
[`doxygen`](http://www.doxygen.nl) on the [`Doxyfile`](https://github.com/softwareqinc/qpp/blob/master/Doxyfile) file. The tool `dot` from the [`Graphviz`](https://www.graphviz.org) package must be installed (`sudo apt-get install graphviz` in Ubuntu/Debian). Running `doxygen` will generate the 
documentation directory `doc` containing both the HTML and LaTeX documentation.

The HTML documentation file will be accessible by opening `doc/html/index.html` with the browser of your choice.

To generate a PDF file of the documentation, run 

```bash
latexmk -pdf refman.tex
```

from the `doc/latex` directory or compile the file `doc/latex/refman.tex` with your LaTeX compiler. This will create the `doc/latex/refman.pdf` documentation file. Consult your favourite LaTeX manual for how to compile/build LaTeX files under your specific operating system.
