// Quantum Phase Estimation
// Source: ./examples/qpe.cpp
// See also ./examples/circuits/quantum_phase_estimation.cpp for a high-level
// API example

/*
A program to construct the following quantum phase estimator circuit and
execute simulation on the phase of U = diag(1, e^{2*pi*i*theta}).

|0> ----H--------------@----x-------------------R3+--R2+--H----|D---- [q0]
                       |    |                   |    |
|0> ----H---------@----+----+---------R2+--H----+----@---------|D---- [q1]
                  |    |    |         |         |
|0> ----H----@----+----+----x----H----@---------@--------------|D---- [q2]
             |    |    |
|0> ----X----U----U^2--U^4------------------------------------------- [q3]

*/

#include <cmath>
#include <iostream>
#include <numeric>
#include <vector>

#include "qpp.h"

int main() {
    using namespace qpp;
    idx n = 4; // number of qubits; qubits 0, 1, ... , n - 2 are counting
               // qubits, qubit n - 1 is the ancilla
    ket psi = mket(std::vector<idx>(n, 0)); // |0>^\otimes n

    cmat U(2, 2); // initialize a unitary operator
    // we use the T gate as an example; we expect estimated theta = 1/8 (0.125).
    double theta = 0.125; // change if you want, increase n for more precision
    U << 1, 0, 0, std::exp(2 * pi * 1_i * theta);

    ket result = psi;
    std::vector<idx> counting_qubits(n - 1);
    std::iota(counting_qubits.begin(), counting_qubits.end(), 0);
    idx ancilla = n - 1;

    std::cout << ">> QPE on n = " << n << " qubits. ";
    std::cout << "The sequence of applied gates is:\n";

    for (idx i = 0; i < counting_qubits.size(); ++i) {
        // apply Hadamard on counting qubits
        result = apply(result, gt.H, {i});
        std::cout << "H" << i << " ";
    }
    result = apply(result, gt.X, {ancilla}); // apply X on the ancilla
    std::cout << "X" << ancilla << '\n';

    // apply controlled unitary operations
    idx reps = 1;
    for (idx i = 0; i < n - 1; ++i) {
        std::cout << "CU(" << n - i - 2 << ", " << n - 1 << ")^" << reps;
        for (idx r = 0; r < reps; ++r) {
            result = applyCTRL(result, U, {n - i - 2}, {n - 1});
        }
        reps *= 2;
        std::cout << '\n';
    }

    // apply inverse quantum Fourier transform to convert state of the counting
    // register
    result = applyTFQ(result, counting_qubits);
    std::cout << "QFT^{-1}" << disp(counting_qubits, ", ", "(", ")") << "\n";

    // measure the counting register and readout probabilities
    auto measured = measure_seq(result, {counting_qubits});
    auto res = std::get<RES>(measured);
    std::cout << ">> Measurement result [q0 q1 ... ]: " << disp(res, " ");
    std::cout << '\n';

    // decimal representation of the measurement result
    idx decimal = multiidx2n(res, std::vector<idx>(counting_qubits.size(), 2));

    // readout phase estimate
    double theta_e = decimal / std::pow(2, n - 1);
    std::cout << ">> Input theta = " << theta << '\n';
    std::cout << ">> Estimated theta = " << theta_e << '\n';
    std::cout << ">> Norm difference: " << std::abs(theta_e - theta) << '\n';
}
