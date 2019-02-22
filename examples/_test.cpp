// Used for testing, do not use it as an example
#include <iostream>

#include "qpp.h"

int main() {
    /////////// testing ///////////
    using namespace qpp;

    idx N = 4;
    QCircuit qc{N};

    for (idx i = 0; i < N; ++i) {
        idx ctrl = randidx(0, N - 1);
        idx target = randidx(0, N - 1);
        while (ctrl == target)
            target = randidx(0, N - 1);

        qc.gate(gt.H, ctrl);
        qc.gate(gt.CNOT, ctrl, target);
    }

    std::cout << qc << std::endl;
    std::cout << qc.get_gate_count() << std::endl;
    std::cout << qc.get_gate_depth("H") << std::endl;
    std::cout << qc.get_gate_depth("CNOT") << std::endl;
    std::cout << qc.get_gate_depth() << std::endl;
}
