// Used for testing, do not use it as an example
#include <iostream>

#include "qpp.h"
#include "experimental/experimental.h"

int main() {
    using namespace qpp;
    using namespace qpp::experimental;

    QCircuit circ{10, 12};
    // circ.apply(gt.X, 7, "X gate on qubit 7");
    circ.measureV(gt.X, 4, "Initial measurement on qubit 4");
    circ.apply(gt.X, 2, "X gate on qubit 2");
    circ.apply(gt.CNOT, 1, 3, "CNOT on 1 3");
    circ.measureZ(2, "Z meas on qubit 2");
    circ.apply(gt.TOF, 5, 3, 6, "TOF on 5 3 6");
    circ.measureZ(6);
    circ.measureZ(9);
    circ.CTRL(gt.X, 6, 7, "Control X 6->7");
    circ.measureV(gt.H, 7, "Measure in X basis");
    circ.apply(gt.Y, 8);
    circ.measureKs({std::sqrt(2) * gt.X, std::sqrt(2) * gt.Z}, 9, "Kraus");
    circ.measureZ(9);
    // circ.apply(gt.TOF, 4, 6, 7, "TOF 4 6 7");

    std::cout << circ;

    QCircuit test{4};
    std::cout << test;
}
