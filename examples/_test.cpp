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

    Lattice l(3, 4, 5); // Dx = 3, Dy = 4, Dz = 5
    std::cout << l(1, 2, 3) << "\n";
}
