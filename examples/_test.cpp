// Used for testing, do not use it as an example
#include <iostream>

#include "qpp.h"

int main() {
    /////////// testing ///////////
    using namespace qpp;

    QCircuit qc{3, 3};
    qc.reset(0);
    qc.reset({1, 2});
    qc.cCTRL(gt.X, {0, 1}, 0, {1, 0});
    qc.measureZ(0, 2);
    std::cout << qc << "\n\n";
    std::cout << qc.to_JSON() << "\n\n";

    QEngine q_engine{qc};
    q_engine.set_dit(1, 1);
    q_engine.execute();
    std::cout << q_engine << "\n\n";
}