// Source: ./examples/qft_inplace.cpp
//
// Quantum Fourier transform, in place

#include <cmath>
#include <iostream>
#include <vector>

#include <qpp/qpp.hpp>

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
        // apply Hadamard on qubit 'i', in place
        apply_inplace(result, gt.H, {i});
        // apply controlled rotations, in place
        for (idx j = 2; j <= n - i; ++j) {
            cmat Rj(2, 1);
            auto pow_j = static_cast<idx>(std::llround(std::pow(2, j)));
            Rj << 1, omega(pow_j);
            applyCTRL_diag_inplace(result, Rj, {i + j - 1}, {i});
            std::cout << "R" << j << "(" << i + j - 1 << ", " << i << ") ";
        }
        std::cout << '\n';
    }

    // we have the qubits in reversed order, we must swap them, in place
    for (idx i = 0; i < n / 2; ++i) {
        std::cout << "SWAP(" << i << ", " << n - i - 1 << ")\n";
        apply_inplace(result, gt.SWAP, {i, n - i - 1});
    }

    // check that we got the Fourier transform, compute the norm difference
    if (n < 14) { // otherwise, not enough memory for computing gt.Fd(D) * psi
        std::cout << ">> Norm difference: " << norm(result - gt.Fd(D) * psi)
                  << '\n';
    }
}
