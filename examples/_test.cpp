// Used for testing, do not use it as an example
#include <iostream>

#include "qpp.h"

int main() {
    /////////// testing ///////////
    using namespace qpp;

    const idx N = 10000; // qubits
    const idx w = 100;   // gates per qubit
    QCircuit qc{N};

    std::cout << "filling the circuit...\n";
    for (idx i = 0; i < w * N; ++i) {
        idx ctrl = randidx(0, N - 1);
        idx target = randidx(0, N - 1);
        while (target == ctrl)
            target = randidx(0, N - 1);
        qc.gate(gt.CNOT, ctrl, target);
    }

    std::cout << "gate count: " << qc.get_gate_count() << std::endl;
    std::cout << "gate depth: " << qc.get_gate_depth() << std::endl;
}
