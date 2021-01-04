// Quantum Phase Estimation
// Source: ./examples/qpe.cpp
// See also ./examples/circuits/quantum_phase_estimation.cpp for a high-level
// API example

/*
A program to construct the following quantum phase estimator circuit and
execute simulation on the phase of U = diag(1, e^{2*pi*i*theta}).

|0> ---H----@--------------x----H----@----------@----------------|D--- [q0]
            |              |         |          |
|0> ---H----+----@---------+---------R2----H----+-----@----------|D--- [q1]
            |    |         |                    |     |
|0> ---H----+----+----@----x--------------------R3----R2----H----|D--- [q2]
            |    |    |
|0> ---X----U----U^2--U^4--------------------------------------------- [q3]

*/

#include <cmath>
#include <iostream>
#include <vector>

#include "qpp.h"

int main() {
    using namespace qpp;
    // initial four qubit state
    // qubits 0, 1, 2 are counting qubits, qubit 3 is the ancilla
    std::vector<idx> qubits{0, 0, 0, 0};
    ket psi = mket(qubits);

    cmat U(2, 2); // initialize a unitary operator
    // use T-Gate as example; we expect estimated theta = 1/8.
    double theta = 1 / 8.;
    U << 1, 0, 0, std::exp(2 * pi * 1_i * theta);

    ket result = psi;
    idx n = qubits.size(); // number of qubits

    std::cout << ">> QPE on n = " << n << " qubits. ";
    std::cout << "The sequence of applied gates is:\n";

    result = apply(result, gt.X, {n - 1}); // apply X on the ancilla
    std::cout << "X" << (n - 1) << " ";

    for (idx i = 0; i < n - 1; ++i) {
        // apply Hadamard on counting qubits
        result = apply(result, gt.H, {i});
        std::cout << "H" << i << " ";
    }
    std::cout << '\n';

    // apply controlled unitary operations
    idx repetitions = 1;
    for (idx i = 0; i < n - 1; ++i) {
        for (idx r = 0; r < repetitions; ++r) {
            result = applyCTRL(result, U, {i}, {n - 1});
            std::cout << "CU"
                      << "(" << i << ", " << n - 1 << ") ";
        }
        repetitions *= 2;
        std::cout << '\n';
    }

    // apply inverse quantum Fourier transform to convert state of the counting
    // register
    for (idx i = 0; i < (n - 1) / 2; ++i) {
        std::cout << "SWAP(" << i << ", " << (n - 1) - i - 1 << ")\n";
        result = apply(result, gt.SWAP, {i, (n - 1) - i - 1});
    }

    for (idx j = 0; j < n - 1; ++j) {
        for (idx m = 0; m < j; ++m) {
            cmat Rj(2, 2);
            Rj << 1, 0, 0, std::exp(-1_i * pi / std::pow(2, j - m));
            result = applyCTRL(result, Rj, {m}, {j});
            std::cout << "R" << j - m + 1 << "(" << m << ", " << j << ") ";
        }
        result = apply(result, gt.H, {j});
        std::cout << "H" << j << "\n";
    }
    std::cout << '\n';

    // measure the counting register and readout probabilities
    double decimal = 0.0;
    for (idx i = 0; i < n - 1; ++i) {
        auto measured = measure(result, gt.Z, {i});
        auto res = std::get<RES>(measured);
        auto prob = std::get<PROB>(measured);
        std::cout << ">> Measurement result q" << i << ": " << res << '\n';
        std::cout << ">> Probabilities: ";
        std::cout << disp(prob, ", ") << '\n';
        if (prob[0] < prob[1]) {
            decimal += std::pow(2, i);
        }
    }
    std::cout << '\n';

    // readout phase estimate
    double theta_e = decimal / std::pow(2, n - 1);
    std::cout << ">> Input theta = " << theta << '\n';
    std::cout << ">> Estimated theta = " << theta_e << '\n';
    std::cout << ">> Norm difference: " << std::abs(theta_e - theta) << '\n';
}
