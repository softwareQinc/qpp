// Used for testing, do not use it as an example
#include <iostream>

#include "qpp.h"

int main() {
    using namespace qpp;
    /////////// testing ///////////

    QCircuitDescription qc{4, 2, 2, "test_circuit"};
    qc.measureZ(3, 0);
    qc.gate(gt.X, 0, "named_X");
    qc.gate(gt.Z, 1);
    qc.CTRL(gt.X, 0, 1, "ctrl_X");
    qc.measureZ(0, 0);
    qc.measureV(gt.H, 1, 1);
    qc.measureZ(2, 0);

    std::cout << qc << '\n';
    std::cout << qc.to_JSON() << "\n\n";

    auto n = QuditDepolarizingNoise(0.2, 3);
    std::cout << n.get_d() << '\n';
    for (auto&& elem : n.get_Ks())
        std::cout << disp(elem) << "\n\n";
}
