// Used for testing, do not use it as an example
#include <iostream>

#include "qpp.h"

int main() {
    /////////// testing ///////////
    using namespace qpp;

    QCircuit qc{1, 2};
    qc.cCTRL(gt.X, 0, 0, 1);
    qc.measureZ(0, 1);
    std::cout << qc << "\n\n";

    QEngine q_engine{qc};
    q_engine.execute();
    std::cout << q_engine << "\n\n";
}
