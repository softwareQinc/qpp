// Used for testing, do not use it as an example
#include <iostream>

#include "qpp.h"

int main() {
    using namespace qpp;

    /////////// testing ///////////
    // we simulate the action fully depolarizing noise applied multiple times
    QCircuitDescription qc{3,2,2, "test_circuit"};
    qc.gate(gt.X, 0, "named_X");
    qc.gate(gt.Z, 1);
    qc.CTRL(gt.X, 0, 1, "ctrl_X");
    qc.measureV(gt.H, 0, 0);
    qc.measureV(gt.H, 1, 1);

    std::cout << qc << '\n';
    std::cout << qc.to_JSON() << '\n';

}