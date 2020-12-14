/*
 **************************************************************************
 Quantum Phase Estimation
 Source: ./examples/qpe.cpp
 
 A program to construct the following quantum phase estimator circuit and
 execute simulation on the T-gate.
                                 
|0> ---H----@--------------x----H----@----------@----------------M--- [q0]
            |              |         |          |
|0> ---H----+----@---------+---------R2----H----+-----@----------M--- [q1]
            |    |         |                    |     |
|0> ---H----+----+----@----x--------------------R3----R2----H----M--- [q2]
            |    |    |
|0> ---X----U----U^2--U^4-------------------------------------------- [q3]
 
 **************************************************************************
 */

#include <cmath>
#include <iostream>
#include <vector>

#include "qpp.h"

int main() {

    using namespace qpp;

    // Initial four qubit state: qubits 0,1,2 are counting qubits, qubit 3 is the ancilla
    std::vector<idx> qubits{0, 0, 0, 0};

    cmat U(2,2); // initialize a unitary operator
    U << 1, 0, 0, std::exp(1_i*pi/4.); // use T-Gate as example. We expect theta = 1/8.

    ket psi = mket(qubits);
    ket result = psi;
    idx n = qubits.size(); // number of qubits

    std::cout << ">> QPE on N = " << n << " qubits. ";
    std::cout << "The sequence of applied gates is:\n";

    result = apply(result, gt.X, {n-1}); // apply X on the ancilla
    std::cout << "X" << (n-1) << " ";

    for (idx qubit = 0; qubit < n-1; ++qubit) {
        result = apply(result, gt.H, {qubit}); // apply Hadamard on counting qubits
        std::cout << "H" << qubit << " ";
    }
    std::cout << '\n';

    // Apply controlled unitary operations
    idx repetitions = 1;
    for (idx qubit = 0; qubit < n-1; ++qubit) {
        for (idx r = 0; r < repetitions; ++r) {
            result = applyCTRL(result, U, {qubit}, {n-1});
            std::cout << "CU" << "(" << qubit << ", " << n-1 << ") ";
        }
        repetitions *= 2;
        std::cout << '\n';
    }

    // Apply inverse quantum Fourier transform to convert state of the counting register
    for (idx i = 0; i < (n-1) / 2; ++i) {
        std::cout << "SWAP(" << i << ", " << (n-1) - i - 1 << ")\n";
        result = apply(result, gt.SWAP, {i, (n-1) - i - 1});
    }

    for (idx j = 0; j < n-1; ++j) {
        for (idx m = 0; m < j; ++m) {
            cmat Rj(2,2);
            Rj << 1, 0, 0, std::exp(-1_i*pi/pow(2,j-m));
            result = applyCTRL(result, Rj, {m}, {j});
            std::cout << "R" << j << "(" << m << ", " << j << ") ";
        }
        result = apply(result, gt.H, {j});
        std::cout << "H" << j << " ";
        std::cout << '\n';
    }
    std::cout << '\n';

    // Measure the counting register and readout probablilites
    double decimal = 0.0;
    for (idx i = 0; i < n-1; ++i) {
        auto measured = measure(result, gt.Z, {i});
        auto res = std::get<RES>(measured);
        auto prob = std::get<PROB>(measured);
        std::cout << ">> Measurement result q" << i << ": " << res << '\n';
        std::cout << ">> Probabilities: ";
        std::cout << disp(prob, ", ") << '\n';
        if (prob[0] < prob[1]) { decimal += pow(2,i); }
    }
    std::cout << '\n';
    
    // Readout phase estimate
    double theta = decimal/pow(2,n-1);
    std::cout << "theta = " << theta << '\n';
    
}
    


