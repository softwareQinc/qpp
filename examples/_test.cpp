// Used for testing, do not use it as an example
#include <iostream>

#include "qpp.h"

int main() {
    /////////// testing ///////////
    using namespace qpp;

    QCircuit qc{4, 4};
    qc.QFT();
    qc.measureZ(0,0);
    qc.measureZ(1,1);

    std::cout << qc << std::endl << std::endl;

    QEngine q_engine{qc};
    std::cout << q_engine << std::endl;

    std::cout << qc.get_gate_depth("CTRL-R2") << std::endl;
    std::cout << qc.get_gate_depth() << std::endl;
    std::cout << qc.get_gate_count("H") << std::endl;
}
