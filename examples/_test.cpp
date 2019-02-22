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

    // std::cout << qc << std::endl;
    std::cout << qc.get_gate_count() << std::endl;
    std::cout << qc.get_gate_depth("H") << std::endl;
    std::cout << qc.get_gate_depth("CNOT") << std::endl;
    std::cout << qc.get_gate_depth() << std::endl << std::endl;

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

    std::cout << bit_circuit.gate_count.total << std::endl;
    std::cout << bit_circuit.gate_count.X << std::endl;
    std::cout << bit_circuit.gate_count.CNOT << std::endl << std::endl;

    std::cout << bit_circuit.gate_depth.X << std::endl;
    std::cout << bit_circuit.gate_depth.CNOT << std::endl;
    std::cout << bit_circuit.gate_depth.total << std::endl;
}
