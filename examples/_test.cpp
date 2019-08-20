// Used for testing, do not use it as an example
#include <iostream>

#include "qpp.h"

int main() {
    /////////// testing ///////////
    using namespace qpp;

    QCircuit qc = *qasm::read_from_file(
        PATH "/unit_tests/tests/qasm/circuits/generic/teleport.qasm");
    std::cout << qc << "\n\n";

    QEngine q_engine{qc};
    q_engine.execute();
    std::cout << q_engine << '\n';
}