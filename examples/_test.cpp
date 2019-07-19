// Used for testing, do not use it as an example
#include <iostream>

#include "qpp.h"

int main() {
    /////////// testing ///////////
    using namespace qpp;

    idx nq = 2, nc = 3;
    QCircuit qc{nq, nc};
    qc.gate_fan(gt.H);
    for (idx i = 0; i < qc.get_nq(); ++i)
        qc.measureZ(i, i);

    auto qc1 = qc;

    std::cout << qc1 << "\n";

    std::cout << '\n' << "engine " << '\n';
    QEngine q_engine{qc1};
    q_engine.execute();
    q_engine.execute(1, false);

    std::cout << q_engine << '\n';
    std::cout << q_engine.to_JSON() << '\n';
}
