// Various physical layouts
// Source: ./examples/layouts.cpp

#include <iostream>

#include "qpp/qpp.hpp"

int main() {
    using namespace qpp;

    std::cout << ">> Various physical qudit layouts\n";

    // a 2 x 2 x 2 orthogonal lattice
    auto lattice = Lattice{2, 2, 2};
    // a 2 x 2 x 2 periodic boundary orthogonal lattice
    auto periodic_boundary_lattice = PeriodicBoundaryLattice{2, 2, 2};

    std::cout << ">> The (1, 1, 1) coordinate in an orthogonal lattice of "
                 "dimensions ";
    std::cout << disp(lattice.get_dims(), IOManipContainerOpts{}.set_sep(", "))
              << " is:\n";
    std::cout << lattice(1, 1, 1) << '\n';

    std::cout << ">> The index " << lattice(1, 1, 1)
              << " in an orthogonal lattice of dimensions ";
    std::cout << disp(lattice.get_dims(), IOManipContainerOpts{}.set_sep(", "))
              << " maps to:\n";
    std::cout
        << disp(lattice.to_coordinates(lattice(1, 1, 1)),
                IOManipContainerOpts{}.set_sep(", ").set_left("(").set_right(
                    ")"))
        << '\n';

    std::cout << ">> The (2, 2, 3) coordinate in a periodic boundary "
                 "orthogonal lattice of dimensions ";
    std::cout << disp(periodic_boundary_lattice.get_dims(),
                      IOManipContainerOpts{}.set_sep(", "))
              << " is:\n";
    std::cout << periodic_boundary_lattice(3, 3, 4) << '\n';
}
