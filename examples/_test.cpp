// Used for testing, do not use it as an example
#include <iostream>

#include "qpp.h"

int main() {
    /////////// testing ///////////
    using namespace qpp;

    QCircuit qc{1, 3};
    qc.cCTRL(gt.X, {0, 1}, 0, {1, 0});
    qc.measureZ(0, 2);
    std::cout << qc << "\n\n";
    std::cout << qc.to_JSON() << "\n\n";

    QEngine q_engine{qc};
    q_engine.set_dit(1, 1);
    q_engine.execute();
    std::cout << q_engine << "\n\n";

    std::vector<idx> dims{31, 42, 150, 12};
    std::vector<idx> coordinates{10, 12, 30, 2};

    // Layouts

    // using vectors
    Lattice l(dims);
    std::cout << "Dims: " << disp(l.get_dims(), " ") << "\n";

    std::cout << "Coords: " << disp(coordinates, " ") << "\n";
    std::cout << "Index: " << l(coordinates) << "\n";
    std::cout << "Back to coords: "
              << disp(l.to_coordinates(l(coordinates)), " ") << "\n\n";

    // using variadics
    Lattice l2(32, 43, 151, 14);
    std::cout << "Dims: " << disp(l.get_dims(), " ") << "\n";

    std::cout << "Coords: " << disp(std ::vector<idx>{11, 20, 31, 9}, " ")
              << "\n";
    std::cout << "Index: " << l(11, 20, 31, 9) << "\n";
    std::cout << "Back to coords: "
              << disp(l.to_coordinates(l(11, 20, 31, 9)), " ") << "\n\n";
}
