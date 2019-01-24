// Qubit teleporation circuit simulator
// Source: ./examples/teleport_qubit.cpp
#include <iostream>

#include "qpp.h"

int main() {
    using namespace qpp;
    idx d = 5; // qudit dimension
    std::cout << ">> Qudit teleportation (d = " << d << ") ";
    std::cout << "quantum circuit simulation\n\n";

    // quantum circuit description with 3 qudits, 2 classical bits, (optional)
    // dimension d = 5 and (optional) name = "qudit teleportation"
    QCircuitDescription tele_qudit{3, 2, d, "qudit teleportation"};
    // set the qubit 0 to a random state
    cmat U = randU(d);
    tele_qudit.gate(U, 0);

    // set the MES between qudits 1 and 2
    tele_qudit.gate(gt.Fd(d), 1);
    tele_qudit.CTRL(gt.Xd(d), 1, 2);

    // perform the Bell measurement between qudits 0 and 1
    tele_qudit.CTRL(adjoint(gt.Xd(d)), 0, 1);
    tele_qudit.gate(adjoint(gt.Fd(d)), 0);

    // perform the measurements
    tele_qudit.measureZ(0, 0);
    tele_qudit.measureZ(1, 1);

    // apply the classical controls
    tele_qudit.cCTRL(adjoint(gt.Xd(d)), 1, 2);
    tele_qudit.cCTRL(gt.Zd(d), 0, 2);

    QCircuit qc_tele_qudit{tele_qudit};
    std::cout << ">> BEGIN CIRCUIT DESCRIPTION\n";
    std::cout << qc_tele_qudit.get_circuit_description() << '\n';
    std::cout << ">> END CIRCUIT DESCRIPTION\n\n";

    // non-verbose run, use qc_teleport_qubit.run(true) for verbose running
    qc_tele_qudit.run();

    // displays the  measurement statistics
    std::cout << ">> BEGIN AFTER RUNNING\n";
    std::cout << qc_tele_qudit << '\n';
    std::cout << ">> END AFTER RUNNING\n\n";

    // verify that the teleportation was successful
    ket psi_initial = U * mket({0}, d);
    ket psi_final = qc_tele_qudit.get_psi();
    std::cout << ">> Teleported state:\n";
    std::cout << disp(psi_final) << '\n';
    std::cout << ">> Norm difference: " << norm(psi_final - psi_initial)
              << '\n';
}