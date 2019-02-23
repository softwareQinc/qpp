// Used for testing, do not use it as an example
#include <iostream>

#include "qpp.h"

int main() {
    /////////// testing ///////////
    using namespace qpp;

    idx nq = 10000;
    QCircuit qc{nq};

    for (idx i = 0; i < nq; ++i) {
        idx ctrl = randidx(0, nq - 1);
        idx target = randidx(0, nq - 1);
        while (ctrl == target)
            target = randidx(0, nq - 1);

        // qc.gate(gt.H, ctrl);
        qc.gate(gt.H, randidx(0, nq - 1));
        qc.gate(gt.CNOT, ctrl, target);
    }

    // std::cout << qc << '\n';
    std::cout << qc.get_gate_count() << '\n';
    std::cout << qc.get_gate_depth("H") << '\n';
    std::cout << qc.get_gate_depth("CNOT") << '\n';
    std::cout << qc.get_gate_depth() << "\n\n";

    idx nc = nq;
    Bit_circuit bit_circuit{nc};
    for (idx i = 0; i < nc; ++i) {
        idx a = randidx(0, nc - 1);
        idx b = randidx(0, nc - 1);
        while (a == b)
            b = randidx(0, nc - 1);

        bit_circuit.NOT(randidx(0, nc - 1));
        bit_circuit.CNOT(a, b);
    }

    std::cout << bit_circuit.get_gate_count() << '\n';
    std::cout << bit_circuit.get_gate_count("X") << '\n';
    std::cout << bit_circuit.get_gate_count("CNOT") << "\n\n";

    std::cout << bit_circuit.get_gate_depth("X") << '\n';
    std::cout << bit_circuit.get_gate_depth("CNOT") << '\n';
    std::cout << bit_circuit.get_gate_depth() << "\n\n";

    bit_circuit.reset();

    std::cout << bit_circuit.get_gate_count() << '\n';
    std::cout << bit_circuit.get_gate_count("X") << '\n';
    std::cout << bit_circuit.get_gate_count("CNOT") << "\n\n";

    std::cout << bit_circuit.get_gate_depth("X") << '\n';
    std::cout << bit_circuit.get_gate_depth("CNOT") << '\n';
    std::cout << bit_circuit.get_gate_depth() << "\n\n";
}
