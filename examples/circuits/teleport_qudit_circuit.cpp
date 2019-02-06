// Qudit teleporation circuit simulator
// Source: ./examples/circuits/teleport_qudit_circuit.cpp
#include <iostream>

#include "qpp.h"

int main() {
    using namespace qpp;
    idx d = 5; // qudit dimension
    std::cout << ">> Qudit teleportation (d = " << d << ") ";
    std::cout << "quantum circuit simulation\n\n";

    // quantum circuit with 3 qudits, 2 classical bits, (optional) dimension
    // d = 5 and (optional) name = "qudit teleportation"
    QCircuit qc{3, 2, d, "qudit teleportation"};
    // set the qubit 0 to a random state
    cmat U = randU(d);
    // apply the gate U to qubit 0
    qc.gate(U, 0);

    // set the MES between qudits 1 and 2
    qc.gate(gt.Fd(d), 1);
    qc.CTRL(gt.Xd(d), 1, 2);

    // perform the Bell measurement between qudits 0 and 1
    qc.CTRL(adjoint(gt.Xd(d)), 0, 1);
    qc.gate(adjoint(gt.Fd(d)), 0);

    // perform the measurements
    qc.measureZ(0, 0);
    qc.measureZ(1, 1);

    // apply the classical controls
    qc.cCTRL(adjoint(gt.Xd(d)), 1, 2);
    qc.cCTRL(gt.Zd(d), 0, 2);

    // initialize the quantum engine with a circuit
    QEngine engine{qc};

    // display the quantum circuit
    std::cout << ">> BEGIN CIRCUIT\n";
    std::cout << engine.get_circuit() << '\n';
    std::cout << ">> END CIRCUIT\n\n";

    // execute the circuit
    for (auto&& step : qc)
        engine.execute(step);

    // display the measurement statistics
    std::cout << ">> BEGIN AFTER RUNNING\n";
    std::cout << engine << '\n';
    std::cout << ">> END AFTER RUNNING\n\n";

    // verify that the teleportation was successful
    ket psi_initial = U * mket({0}, d);
    ket psi_final = engine.get_psi();
    std::cout << ">> Teleported state:\n";
    std::cout << disp(psi_final) << '\n';
    std::cout << ">> Norm difference: " << norm(psi_final - psi_initial)
              << '\n';
}
