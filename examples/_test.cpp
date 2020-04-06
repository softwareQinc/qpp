// Used for testing, do not use it as an example
#include <iostream>

#include "qpp.h"
#include "experimental/experimental.h"

int main() {
    /////////// testing ///////////
    using namespace qpp;
    //
    // ket in = st.plus();
    // qram data{1, 2, 3};
    // ket out = qRAM(in, data);
    // std::cout << disp(out) << '\n';

    idx nq = 20, nc = 2;
    QCircuit q_circuit{nq, nc};
    q_circuit.gate(gt.H, 0);
    q_circuit.CTRL(gt.X, 1, 5);
    q_circuit.CTRL(gt.X, 10, 12);

    std::cout << "Before compression\n" << q_circuit << '\n';
    std::cout << disp(q_circuit.get_clean(), ",") << '\n';

    std::cout << "\nAfter compression\n";

    q_circuit.compress();
    std::cout << q_circuit << '\n';
    std::cout << disp(q_circuit.get_clean(), ",") << '\n';

    QEngine q_engine{q_circuit};
    q_engine.execute();
}
