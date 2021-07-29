// Quantum Fourier transform
// Source: ./examples/qft.cpp
#include <cmath>
#include <iostream>
#include <vector>

#include "qpp.h"

int main() {
    using namespace qpp;

    std::vector<idx> qubits{1, 0, 1, 1, 0}; // initial state
    ket psi = mket(qubits);
    ket result = psi;

    idx n = qubits.size();                                   // number of qubits
    auto D = static_cast<idx>(std::llround(std::pow(2, n))); // dimension 2^n
    std::cout << ">> QFT on n = " << n << " qubits. ";

    std::cout << "The sequence of applied gates is:\n";
    for (idx i = 0; i < n; ++i) {
        std::cout << "H" << i << " ";
        result = apply(result, gt.H, {i}); // apply Hadamard on qubit i
        // apply controlled rotations
        for (idx j = 2; j <= n - i; ++j) {
            cmat Rj(2, 2);
            auto pow_j = static_cast<idx>(std::llround(std::pow(2, j)));
            Rj << 1, 0, 0, omega(pow_j);
            result = applyCTRL(result, Rj, {i + j - 1}, {i});
            std::cout << "R" << j << "(" << i + j - 1 << ", " << i << ") ";
        }
        std::cout << '\n';
    }

    // we have the qubits in reversed order, we must swap them
    for (idx i = 0; i < n / 2; ++i) {
        std::cout << "SWAP(" << i << ", " << n - i - 1 << ")\n";
        result = apply(result, gt.SWAP, {i, n - i - 1});
    }

    // check that we got the Fourier transform, compute the norm difference
    if (n < 14) { // otherwise not enough memory for computing gt.Fd(D) * psi
        std::cout << ">> Norm difference: " << norm(result - gt.Fd(D) * psi)
                  << '\n';
    }
}
