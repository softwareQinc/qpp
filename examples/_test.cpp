// Used for testing, do not use it as an example
#include <iostream>

#include "qpp.h"

int main() {
    /////////// testing ///////////
    using namespace qpp;

    idx nq = 2; //, nc = 2;
    QCircuit qc{3, 3};
    qc.CTRL(gt.X, 0, 1);
    qc.CTRL(gt.X, 1, 2);

    qc.gate_fan(gt.H);
    qc.gate_fan(gt.H);
    qc.gate_fan(gt.H);

    // qc.replicate(9);
    for (idx i = 0; i < nq; ++i)
        qc.measureZ(i, i);

    QCircuit qc_add{3, 2};
    qc_add.CTRL(gt.Y, 0, 1);
    qc_add.gate(gt.H, 2);

    std::cout << qc << "\n\n";
    std::cout << qc_add << "\n\n";

    qc.add_circuit(qc_add, 4);
    qc.measureZ(6, 4);

    std::cout << qc << "\n\n";
    //std::cout << qc.to_JSON() << "\n\n";

    QEngine q_engine{qc};
    q_engine.execute(360);
    std::cout << q_engine << "\n\n";

    nq = 3;
    QCircuit q_circuit{nq, nq};
    q_circuit.QFT();

    q_circuit.replicate(24);

    for(idx i = 0 ; i < nq; ++i)
        q_circuit.measureZ(i, i);

    std::cout << q_circuit << "\n\n";

    QEngine q_engine1{q_circuit};
    q_engine1.execute(1024);

    std::cout << q_engine1 << "\n";
}
