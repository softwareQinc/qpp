// Qubit teleporation circuit simulator
// Source: ./examples/circuits/teleport_qubit_circuit.cpp
#include <iostream>

#include "qpp.h"

int main() {
    using namespace qpp;
    std::cout << ">> Qubit teleportation quantum circuit simulation\n\n";

    // quantum circuit description with 3 qubits and 2 classical bits
    QCircuitDescription teleport_qubit{3, 2};
    // set the qubit 0 to a random state
    cmat U = randU(2);
    teleport_qubit.gate(U, 0, "randU");

    // set the MES between qudits 1 and 2
    teleport_qubit.gate(gt.H, 1);
    teleport_qubit.CTRL(gt.X, 1, 2);

    // perform the Bell measurement between qudits 0 and 1
    teleport_qubit.CTRL(gt.X, 0, 1);
    teleport_qubit.gate(gt.H, 0);
    teleport_qubit.measureZ(0, 0);
    teleport_qubit.measureZ(1, 1);

    // apply the classical controls
    teleport_qubit.cCTRL(gt.X, 1, 2);
    teleport_qubit.cCTRL(gt.Z, 0, 2);

    QCircuit qc_teleport_qubit{teleport_qubit};

    std::cout << ">> BEGIN CIRCUIT DESCRIPTION\n";
    std::cout << qc_teleport_qubit.get_circuit_description() << '\n';
    std::cout << ">> END CIRCUIT DESCRIPTION\n\n";

    // non-verbose run, use qc_teleport_qubit.run(true) for verbose running
    qc_teleport_qubit.run();

    // displays the  measurement statistics
    std::cout << ">> BEGIN AFTER RUNNING\n";
    std::cout << qc_teleport_qubit << '\n';
    std::cout << ">> END AFTER RUNNING\n\n";

    // verify that the teleportation was successful
    ket psi_initial = U * 0_ket;
    ket psi_final = qc_teleport_qubit.get_psi();
    std::cout << ">> Teleported state:\n";
    std::cout << disp(psi_final) << '\n';
    std::cout << ">> Norm difference: " << norm(psi_final - psi_initial)
              << '\n';
}
