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
    idx nq_c = 4;         // number of counting qubits
    idx nq_a = 2;         // number of ancilla qubits
    idx nq = nq_c + nq_a; // total number of qubits
    idx nc = 1;           // nc stores a 'dit'; increase nq for more precision.

    std::cout << ">> Quantum phase estimation quantum circuit simulation\n";
    std::cout << ">> nq_c = " << nq_c << " counting qubits, nq_a = " << nq_a
              << " ancilla qubits\n\n";

    // we use the T\otimes T gate as an example; we want to estimate its last
    // (4-th) eigenvalue; we expect estimated theta = 1/4 (0.25).
    cmat U = kron(gt.T, gt.T); // diag(1,e^{\pi i/4},e^{\pi i/4},e^{2\pi i/4})
    double theta = 0.25;

    QCircuit qc{nq, nc};
    std::vector<idx> counting_qubits(nq_c);
    std::iota(counting_qubits.begin(), counting_qubits.end(), 0);
    std::vector<idx> ancilla(nq_a);
    std::iota(ancilla.begin(), ancilla.end(), nq_c);

    qc.gate_fan(gt.H, counting_qubits);
    qc.gate_fan(gt.X, ancilla); // prepare |11>, the fourth eigenvector of U
    for (idx i = nq_c; i-- > 0;) {
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
