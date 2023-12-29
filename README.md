# Quantum++

## Version 4.3.4 - 14 August 2023

[![GitHub actions](https://github.com/softwareqinc/qpp/actions/workflows/cmake.yml/badge.svg)](https://github.com/softwareQinc/qpp/actions)

---

## About

Quantum++ is a modern C++ general purpose quantum computing library, composed
solely of template header files. Quantum++ is written in standard C++17 and has
very low external dependencies, using only
the [Eigen 3](https://eigen.tuxfamily.org) linear algebra header-only template
library and, if available, the [OpenMP](https://www.openmp.org/) multiprocessing
library.

Quantum++ is not restricted to qubit systems or specific quantum information
processing tasks, being capable of simulating arbitrary quantum processes. The
main design factors taken in consideration were the ease of use, high
portability, and high performance. The library's simulation capabilities are
only restricted by the amount of available physical memory. On a typical
machine (Intel i5 8Gb RAM) Quantum++ can successfully simulate the evolution of
25 qubits in a pure state or of 12 qubits in a mixed state reasonably fast.

To report any bugs or ask for additional features/enhancements, please
[submit an issue](https://github.com/softwareQinc/qpp/issues) with an
appropriate label.

If you are interested in contributing to this project, feel free to contact us.
Alternatively, fork the repository, create a custom branch, add your
contribution, then finally create a pull request. If we accept the pull request,
we will merge your custom branch with the latest main/development branch. The
latter will eventually be merged into a future release version. To contribute,
it is preferable to have a solid knowledge of modern C++ (preferably C++17 or
later), including templates and the standard library, a basic knowledge of
quantum computing and linear algebra, and working experience
with [Eigen 3](https://eigen.tuxfamily.org).

For additional [Eigen 3](https://eigen.tuxfamily.org) documentation
see <https://eigen.tuxfamily.org/dox/>. For a simple
[Eigen 3](https://eigen.tuxfamily.org) quick ASCII reference see
<https://eigen.tuxfamily.org/dox/AsciiQuickReference.txt>.

Copyright (c) 2017 - 2024 softwareQ Inc. All rights reserved.

---

## License

[Quantum++](https://github.com/softwareQinc/qpp) is distributed under the MIT
license. Please see the
[`LICENSE`](https://github.com/softwareQinc/qpp/blob/main/LICENSE) file for more
details.

---

## Installation instructions and further documentation

Please see the installation guide
[`INSTALL.md`](https://github.com/softwareQinc/qpp/blob/main/INSTALL.md) and the
comprehensive [Wiki](https://github.com/softwareQinc/qpp/wiki) for further
documentation and detailed examples.

To generate the full official API documentation in both LaTeX and HTML formats
run
[`doxygen`](https://www.doxygen.nl) on
the [`Doxyfile`](https://github.com/softwareQinc/qpp/blob/main/Doxyfile) file.
The tool `dot` from the [`Graphviz`](https://www.graphviz.org) package must be
installed (`sudo apt-get install graphviz` on Ubuntu/Debian Linux,
or `brew install graphviz` on macOS). Running `doxygen` will generate the
documentation directory `doc` containing both the HTML and LaTeX documentation.

The HTML documentation file will be accessible by opening `doc/html/index.html`
with the browser of your choice. To generate a PDF file of the documentation,
run

```shell
latexmk -pdf refman.tex
```

from the `doc/latex` directory or compile the file `doc/latex/refman.tex` with
your LaTeX compiler. This will create the `doc/latex/refman.pdf` documentation
file. Consult your favourite LaTeX manual for how to compile/build LaTeX files
under your specific operating system.

---

## Python 3 wrapper

[pyqpp](https://github.com/softwareQinc/qpp/blob/main/pyqpp) is a Python 3
wrapper for Quantum++. pyqpp requires the same dependencies as Quantum++, and
can be installed using `pip`

```shell
pip install git+https://github.com/softwareQinc/qpp
```

**Important**: If the installation fails due to your system being unable to
detect the location of the Eigen3 matrix library, set the environment variable
`EIGEN3_INSTALL_DIR` to point to the location of the Eigen3 library
(include the `include/eigen3` part of the path).

For more details, please see
[pyqpp/README.md](https://github.com/softwareQinc/qpp/blob/main/pyqpp/README.md).
