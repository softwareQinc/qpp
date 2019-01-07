// Used for testing, do not use it as an example
#include <iostream>

#include "qpp.h"
#include "experimental/experimental.h"

int main() {
    using namespace qpp;
    using namespace qpp::experimental;

    QCircuit circ{10};
    circ.apply(gt.X, 2, "X gate on qubit 2");
    circ.apply(gt.CNOT, 1, 2, "CNOT on 12");

    for (auto elem : circ.gates_) {
        std::cout << elem.gate_type_ << " " << elem.name_ << '\n';
    }

    auto x = QCircuit::MeasureType::MEASURE_V;
    std::cout << x << " ";
}
