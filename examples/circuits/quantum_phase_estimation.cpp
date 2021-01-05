// Quantum Phase Estimation circuit simulator
// Source: ./examples/circuits/qpe.cpp
// See also ./examples/qpe.cpp for a low-level API example
#include <cmath>
#include <iostream>
#include <numeric>
#include <vector>

#include "qpp.h"

int main() {
    using namespace qpp;
    idx nq = 4; // number of qubits; qubits 0, 1, ... , n - 2 are counting
                // qubits, qubit n - 1 is the ancilla
    idx nc = 1; // nc stores a 'dit'; increase nq for more precision.
    std::cout << ">> Quantum phase estimation quantum circuit simulation on ";
    std::cout << "n = " << nq << " qubits\n\n";

    cmat U(2, 2); // initialize a unitary operator
    // we use the T gate as an example; we expect estimated theta = 1/8 (0.125).
    double theta = 0.125; // change if you want, increase nq for more precision
    U << 1, 0, 0, std::exp(2 * pi * 1_i * theta);

    QCircuit qc{nq, nc};
    std::vector<idx> counting_qubits(nq - 1);
    std::iota(counting_qubits.begin(), counting_qubits.end(), 0);
    idx ancilla = nq - 1;

    qc.gate_fan(gt.H, counting_qubits).gate(gt.X, ancilla);
    for (idx i = ancilla; i-- > 0;) {
        qc.CTRL(U, i, ancilla);
        U = powm(U, 2);
    }
    qc.TFQ(counting_qubits); // inverse Fourier transform
    // measure many qubits at once, result is a dit
    qc.measureZ(counting_qubits, 0);

    // display the quantum circuit
    std::cout << ">> BEGIN CIRCUIT\n";
    std::cout << qc << '\n';
    std::cout << ">> END CIRCUIT\n\n";

    QEngine engine{qc};
    engine.execute();
    // decimal representation of the measurement result
    auto decimal = engine.get_dit(0);
    double theta_e = decimal / std::pow(2, counting_qubits.size());

    std::cout << ">> Input theta = " << theta << '\n';
    std::cout << ">> Estimated theta = " << theta_e << '\n';
    std::cout << ">> Norm difference: " << std::abs(theta_e - theta) << '\n';
}
