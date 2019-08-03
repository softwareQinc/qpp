// Used for testing, do not use it as an example
#include <iostream>

#include "qpp.h"

int main() {
    /////////// testing ///////////
    using namespace qpp;

    QCircuit qc{4, 3, 3};
    qc.add_circuit(adjoint(QCircuit{4, 2, 3}.QFT()),0);
    std::cout << qc << "\n\n";
    std::cout << qc.to_JSON() << "\n\n";
}
