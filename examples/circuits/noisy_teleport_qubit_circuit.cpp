// Qubit noisy teleporation circuit simulator
// Source: ./examples/circuits/noisy_teleport_qubit_circuit.cpp
#include <iostream>

#include "qpp.h"

int main() {
    using namespace qpp;
    std::cout << ">> Qubit noisy teleportation quantum circuit simulation\n\n";

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

    // initialize the noisy quantum engine with an amplitude damping noise model
    // and a quantum circuit; in C++17 you can make use of the class template
    // argument deduction rules to simply write
    // QNoisyEngine noisy_engine{qc, QubitAmplitudeDampingNoise{0.99}};
    QNoisyEngine<QubitAmplitudeDampingNoise> noisy_engine{
        qc, QubitAmplitudeDampingNoise{0.99}};

    // display the quantum circuit
    std::cout << ">> BEGIN CIRCUIT\n";
    std::cout << noisy_engine.get_circuit() << '\n';
    std::cout << ">> END CIRCUIT\n\n";

    // execute the entire circuit
    noisy_engine.execute();

    // display the measurement statistics
    std::cout << ">> BEGIN NOISY ENGINE STATISTICS\n";
    std::cout << noisy_engine << '\n';
    std::cout << ">> END NOISY ENGINE STATISTICS\n\n";

    // verify how successful the teleportation was
    ket psi_initial = U * 0_ket;
    ket psi_final = noisy_engine.get_psi();
    std::cout << ">> Initial state:\n";
    std::cout << disp(psi_initial) << '\n';
    std::cout << ">> Teleported state:\n";
    std::cout << disp(psi_final) << '\n';
    std::cout << ">> Norm difference: " << norm(psi_final - psi_initial)
              << '\n';
}
