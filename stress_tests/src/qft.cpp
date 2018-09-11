// Quantum Fourier transform stress test on a pure state of n qubits
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

#include <omp.h>

#include "qpp.h"

int main(int argc, char **argv) {
    using namespace qpp;
    if (argc != 3) {
        std::cerr << "Please specify the number of cores and qubits!\n";
        exit(EXIT_FAILURE);
    }

    idx num_cores = std::stoi(argv[1]); // number of cores
    idx n = std::stoi(argv[2]);         // number of qubits
    omp_set_num_threads(num_cores);     // number of cores

    std::vector<idx> qubits(n); // initial state
    ket psi = mket(qubits);
    ket result = psi;

    Timer<> t; // start timing
    for (idx i = 0; i < n; ++i) {
        result = apply(result, gt.H, {i}); // apply Hadamard on qubit i
        // apply controlled rotations
        for (idx j = 2; j <= n - i; ++j) {
            cmat Rj(2, 2);
            Rj << 1, 0, 0, omega(std::pow(2, j));
            result = applyCTRL(result, Rj, {i + j - 1}, {i});
        }
    }
    // we have the qubits in reversed order, we must swap them
    for (idx i = 0; i < n / 2; ++i) {
        result = apply(result, gt.SWAP, {i, n - i - 1});
    }
    std::cout << num_cores << ", " << n << ", " << t.toc() << '\n';
}
