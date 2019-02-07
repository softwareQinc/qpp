// Used for testing, do not use it as an example
#include <iostream>

#include "qpp.h"

int main() {
    using namespace qpp;
    /////////// testing ///////////

    QCircuit qc{4, 4, 2, "test_circuit"};
    qc.gate_fan(gt.H);
    qc.gate(gt.X, 0, "named_X");
    qc.measureZ(3, 0);
    qc.gate(gt.X, 0, "named_X");
    qc.gate(gt.Z, 1);
    qc.cCTRL(gt.X, 0, 1);
    qc.CTRL(gt.X, 0, 1);
    qc.gate_fan(gt.H);
    qc.gate_fan(gt.H, {0, 2});
    qc.measureZ(0, 1);
    qc.measureV(gt.Z, 1, 2);
    qc.measureZ(2, 3);

    std::cout << qc << '\n';
    std::cout << qc.to_JSON() << "\n\n";

    QEngine engine{qc};
    for (auto&& step : qc) {
        engine.execute(step);
    }
    std::cout << engine << '\n';
    std::cout << engine.to_JSON() << "\n\n";

    engine.reset();
    for (auto&& step : qc) {
        engine.execute(step);
    }
    std::cout << engine << '\n';
    std::cout << engine.to_JSON() << "\n\n";

    std::cout << qc.get_gate_count("H") << "\n";
    std::cout << qc.get_gate_count() << "\n";
    std::cout << qc.get_measurement_count("Z") << "\n";
    std::cout << qc.get_measurement_count() << "\n";

    std::cout << hash_eigen_expression(gt.X) << "\n";
    std::cout << std::hash<cmat>{}(gt.X) << "\n";

    dmat a(2, 2);
    a << 1, 2, 3, 4.000000000000001;
    dmat b(2, 2);
    b << 1, 2, 3, 4.000000000000002;
    cmat c(1, 1);
    c << 1.1 + 2.2_i;
    std::cout << std::hash<dmat>{}(a) << "\n";
    std::cout << std::hash<dmat>{}(b) << "\n";
    std::cout << std::hash<cmat>{}(c) << "\n";
}
