// Used for testing, do not use it as an example
#include <iostream>

#include "qpp.h"

int main() {
    /////////// testing ///////////
    using namespace qpp;

    GridView gv{2, 2};

    QCircuit qc{4};
    qc.gate_fan(gt.H);
    qc.CTRL(gt.Z, gv(0, 0), gv(0,1));
    qc.CTRL(gt.Z, gv(0, 1), gv(1,1));
    qc.CTRL(gt.Z, gv(1, 1), gv(1,0));
    qc.CTRL(gt.Z, gv(1, 0), gv(0,0));
    qc.gate(gt.X, gv(0, 0));

    QCircuit qc1{4};
    qc1.gate_fan(gt.H);
    qc1.CTRL(gt.Z, gv(0, 0), gv(0,1));
    qc1.CTRL(gt.Z, gv(0, 1), gv(1,1));
    qc1.CTRL(gt.Z, gv(1, 1), gv(1,0));
    qc1.CTRL(gt.Z, gv(1, 0), gv(0,0));
    qc1.gate(gt.Z, gv(0, 1));
    qc1.gate(gt.Z, gv(1, 0));
    //qc1.gate(gt.Z, gv(1, 1));

    QEngine q_engine{qc};
    QEngine q_engine1{qc1};

    q_engine.execute();
    auto psi = q_engine.get_psi();

    q_engine1.execute();
    auto psi1 = q_engine1.get_psi();
    std::cout << q_engine1 << "\n\n";

    std::cout << norm(psi) << "\n";
    std::cout << disp(transpose(psi)) << "\n";
    std::cout << norm(psi - psi1) << "\n";

    QCircuit qc2{2,2};
    qc2.gate(gt.X, 0);
    qc2.nop();
    qc2.gate(gt.Z, 1);
    qc2.gate(gt.Y, 1);
    std::cout << qc2 << "\n\n";

    auto qc3 = adjoint(qc2);
    std::cout << qc3 << "\n\n";
}
