// Used for testing, do not use it as an example
#include <iostream>

#include "qpp.h"
#include "experimental/experimental.h"

int main() {
    using namespace qpp;
    using namespace qpp::experimental;

    QCircuit circ{10};
    circ.apply(gt.X, 2, "X gate on qubit 2");
    circ.apply(gt.CNOT, 1, 3, "CNOT on 1 3");
    circ.measureZ(2, "Z meas on qubit 2");
    circ.apply(gt.TOF, 5, 3, 6, "TOF on 5 3 6");
    circ.measureZ(6, "");

    std::cout << circ;
}
