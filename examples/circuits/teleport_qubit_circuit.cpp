// Qubit teleporation circuit simulator
// Source: ./examples/circuits/teleport_qubit_circuit.cpp
#include <iostream>

#include "qpp.h"

int main() {
    using namespace qpp;
    std::cout << ">> Qubit teleportation quantum circuit simulation\n\n";

    // quantum circuit with 3 qubits and 2 classical bits
    QCircuit qc{3, 2};
    // set the qubit 0 to a random state
    cmat U = randU(2);
    // apply the gate U with name randU to qubit 0
    qc.gate(U, 0, "randU");

    // set the MES between qubits 1 and 2
    qc.gate(gt.H, 1);
    qc.CTRL(gt.X, 1, 2);

    // perform the Bell measurement between qubits 0 and 1
    qc.CTRL(gt.X, 0, 1);
    qc.gate(gt.H, 0);
    qc.measureZ(0, 0);
    qc.measureZ(1, 1);

    // apply the classical controls
    qc.cCTRL(gt.X, 1, 2);
    qc.cCTRL(gt.Z, 0, 2);

    // initialize the quantum engine with a circuit
    QEngine engine{qc};

    // display the quantum circuit
    std::cout << ">> BEGIN CIRCUIT\n";
    std::cout << engine.get_circuit() << '\n';
    std::cout << ">> END CIRCUIT\n\n";

    // execute the entire circuit step by step, same effect as engine.execute();
    for (auto&& step : qc)
        engine.execute(step);

    // display the measurement statistics
    std::cout << ">> BEGIN ENGINE STATISTICS\n";
    std::cout << engine << '\n';
    std::cout << ">> END ENGINE STATISTICS\n\n";

    // verify that the teleportation was successful
    ket psi_initial = U * 0_ket;
    ket psi_final = engine.get_psi();
    std::cout << ">> Teleported state:\n";
    std::cout << disp(psi_final) << '\n';
    std::cout << ">> Norm difference: " << norm(psi_final - psi_initial)
              << '\n';
}
