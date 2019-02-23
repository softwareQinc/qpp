// Used for testing, do not use it as an example
#include <iostream>

#include "qpp.h"

int main() {
    /////////// testing ///////////
    using namespace qpp;

    // Quantum coin flip
    QCircuit qc{1, 1};
    qc.gate(gt.H, 0).measureZ(0, 0);

    QEngine q_engine{qc};
    q_engine.execute();
    std::cout << q_engine.get_dit(0) << '\n';
}
