// Used for testing, do not use it as an example
#include <iostream>

#include "qpp.h"
#include "experimental/experimental.h"

int main() {
    using namespace qpp;
    using namespace qpp::experimental;

    QCircuitDescription circ{10, 10};
    // circ.gate(gt.X, 7, "X gate on qubit 7");
    circ.measureV(gt.X, 4, 0, "Initial measurement on qubit 4");
    circ.gate(gt.X, 2, "X gate on qubit 2");
    circ.gate(gt.CNOT, 1, 3, "CNOT on 1 3");
    circ.measureZ(2, 0, "Z meas on qubit 2");
    circ.gate(gt.TOF, 5, 3, 6, "TOF on 5 3 6");
    circ.measureZ(6, 0);
    circ.measureZ(9, 0);
    circ.CTRL(gt.X, 6, 7, "Control X 6->7");
    circ.measureV(gt.H, 0, 0, "Measure in X basis");
    circ.gate(gt.Y, 8);
    circ.measureKs({std::sqrt(2) * gt.X, std::sqrt(2) * gt.Z}, 7, 0, "Kraus");
    circ.measureZ(8, 0);
    circ.gate(gt.Z, {1, 2, 3}, "Z on 1 2 3");
    circ.measureZ(5, 9);
    circ.QFT({1, 2, 3});
    circ.TFQ({1, 2, 3});
    // circ.gate(gt.TOF, 4, 6, 7, "TOF 4 6 7");

    std::cout << circ << std::endl;

    QCircuitDescription test{10, 1, 2, "some circuit"};

    test.measureZ(1, 0);
    test.measureZ(3, 0);

    std::cout << test << '\n';

    QCircuitDescription teleport{2, 2, 2, "teleportation"};
    teleport.gate(gt.H, 0);
    teleport.gate(gt.CNOT, 0, 1);

    teleport.measureZ(0, 0);
    teleport.measureZ(1, 1);

    // teleport.run();

    // std::cout << disp(teleport.get_psi()) << '\n';

    std::cout << teleport << '\n';

    QCircuit qc_teleport{teleport};
    std::cout << qc_teleport << '\n';
}
