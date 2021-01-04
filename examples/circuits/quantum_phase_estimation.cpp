// Quantum Phase Estimation circuit simulator
// Source: ./examples/circuits/qpe.cpp
// See also ./examples/qpe.cpp for a low-level API example
#include <cmath>
#include <iostream>
#include <vector>

#include "qpp.h"

int main() {
    using namespace qpp;
    idx nq = 4, nc = 1; // nc stores a 'dit'; increase nq for more precision.
    std::cout << ">> Quantum phase estimation quantum circuit simulation on ";
    std::cout << "n = " << nq << " qubits\n\n";

    cmat U(2, 2); // initialize a unitary operator
    // use T-Gate as example; we expect estimated theta = 1/8.
    double theta = 0.125; // change if you want, increase nq for more precision
    U << 1, 0, 0, std::exp(2 * pi * 1_i * theta);

    QCircuit qc{nq, nc};
    std::vector<idx> first_qubits(nq - 1);
    idx last_qubit = nq - 1;
    std::iota(first_qubits.begin(), first_qubits.end(), 0);

    qc.gate_fan(gt.H, first_qubits).gate(gt.X, last_qubit);
    for (idx i = last_qubit; i-- > 0;) {
        qc.CTRL(U, i, last_qubit);
        U = powm(U, 2);
    }
    qc.TFQ(first_qubits);         // inverse Fourier transform
    qc.measureZ(first_qubits, 0); // measure many qubits, result is a dit

    // display the quantum circuit
    std::cout << ">> BEGIN CIRCUIT\n";
    std::cout << qc << '\n';
    std::cout << ">> END CIRCUIT\n\n";

    QEngine engine{qc};
    engine.execute();
    auto measured_result = engine.get_dit(0);
    double theta_e = measured_result / std::pow(2, first_qubits.size());

    std::cout << ">> Input theta = " << theta << '\n';
    std::cout << ">> Estimated theta = " << theta_e << '\n';
    std::cout << ">> Norm difference: " << std::abs(theta_e - theta) << '\n';
}
