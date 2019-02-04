// Qubit teleporation circuit simulator
// Source: ./examples/circuits/teleport_qubit_circuit.cpp
#include <iostream>

#include "qpp.h"

int main() {
    using namespace qpp;
    std::cout << ">> Qubit teleportation quantum circuit simulation\n\n";

    // quantum circuit description with 3 qubits and 2 classical bits
    QCircuitDescription qcd{3, 2};
    // set the qubit 0 to a random state
    cmat U = randU(2);
    // apply the gate U with name randU to qubit 0
    qcd.gate(U, 0, "randU");

    // set the MES between qubits 1 and 2
    qcd.gate(gt.H, 1);
    qcd.CTRL(gt.X, 1, 2);

    // perform the Bell measurement between qubits 0 and 1
    qcd.CTRL(gt.X, 0, 1);
    qcd.gate(gt.H, 0);
    qcd.measureZ(0, 0);
    qcd.measureZ(1, 1);

    // apply the classical controls
    qcd.cCTRL(gt.X, 1, 2);
    qcd.cCTRL(gt.Z, 0, 2);

    QCircuit qc{qcd};

    std::cout << ">> BEGIN CIRCUIT DESCRIPTION\n";
    std::cout << qc.get_circuit_description() << '\n';
    std::cout << ">> END CIRCUIT DESCRIPTION\n\n";

    // execute the circuit
    for (auto&& elem : qcd)
        qc.execute(elem);

    // display the  measurement statistics
    std::cout << ">> BEGIN AFTER RUNNING\n";
    std::cout << qc << '\n';
    std::cout << ">> END AFTER RUNNING\n\n";

    // verify that the teleportation was successful
    ket psi_initial = U * 0_ket;
    ket psi_final = qc.get_psi();
    std::cout << ">> Teleported state:\n";
    std::cout << disp(psi_final) << '\n';
    std::cout << ">> Norm difference: " << norm(psi_final - psi_initial)
              << '\n';
}
