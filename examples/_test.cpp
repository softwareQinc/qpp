// Used for testing, do not use it as an example
#include <iostream>

#include "qpp.h"

int main() {
    using namespace qpp;
    /////////// testing ///////////

    QCircuitDescription qcd{4, 2, 2, "test_circuit"};
    qcd.gate(gt.X, 0, "named_X");
    qcd.measureZ(3, 0);
    qcd.gate(gt.X, 0, "named_X");
    qcd.gate(gt.Z, 1);
    qcd.CTRL(gt.X, 0, 1, "ctrl_X");
    qcd.gate_fan(gt.H);
    qcd.gate_fan(gt.H, {0, 2, 3});
    qcd.measureZ(0, 0);
    qcd.measureV(gt.H, 1, 1);
    qcd.measureZ(2, 0);

    std::cout << qcd << '\n';
    std::cout << qcd.to_JSON() << "\n\n";

    QCircuit qc{qcd};
    auto it = std::begin(qcd);
    for (; it != std::end(qcd); ++it) {
        std::cout << *it << '\n';
        qc.execute(*it);
    }
    std::cout << qc << '\n';
    std::cout << qc.to_JSON() << "\n\n";

    auto n = QuditDepolarizingNoise(0.2, 3);
    std::cout << n.get_d() << '\n';
    for (auto&& elem : n.get_Ks())
        std::cout << disp(elem) << "\n\n";
}
