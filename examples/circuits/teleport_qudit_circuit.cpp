// Qudit teleporation circuit simulator
// Source: ./examples/circuits/teleport_qudit_circuit.cpp
#include <iostream>

#include "qpp.h"

int main() {
    using namespace qpp;
    idx d = 5; // qudit dimension
    std::cout << ">> Qudit teleportation (d = " << d << ") ";
    std::cout << "quantum circuit simulation\n\n";

    // quantum circuit description with 3 qudits, 2 classical bits, (optional)
    // dimension d = 5 and (optional) name = "qudit teleportation"
    QCircuitDescription qcd{3, 2, d, "qudit teleportation"};
    // set the qubit 0 to a random state
    cmat U = randU(d);
    // apply the gate U to qubit 0
    qcd.gate(U, 0);

    // set the MES between qudits 1 and 2
    qcd.gate(gt.Fd(d), 1);
    qcd.CTRL(gt.Xd(d), 1, 2);

    // perform the Bell measurement between qudits 0 and 1
    qcd.CTRL(adjoint(gt.Xd(d)), 0, 1);
    qcd.gate(adjoint(gt.Fd(d)), 0);

    // perform the measurements
    qcd.measureZ(0, 0);
    qcd.measureZ(1, 1);

    // apply the classical controls
    qcd.cCTRL(adjoint(gt.Xd(d)), 1, 2);
    qcd.cCTRL(gt.Zd(d), 0, 2);

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
    ket psi_initial = U * mket({0}, d);
    ket psi_final = qc.get_psi();
    std::cout << ">> Teleported state:\n";
    std::cout << disp(psi_final) << '\n';
    std::cout << ">> Norm difference: " << norm(psi_final - psi_initial)
              << '\n';
}
