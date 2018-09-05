// Quantum Fourier transform stress test
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <vector>

#include <omp.h>

#include "qpp.h"

int main(int argc, char **argv) {
    using namespace qpp;
    if (argc != 3) {
        std::cerr << "Please specify the number of cores and qubits!\n";
        exit(EXIT_FAILURE);
    }

    //std::cout << "Max cores: " << omp_get_num_procs() << "\n";
    idx num_cores = std::stoi(argv[1]); // number of cores
    idx N = std::stoi(argv[2]);         // number of qubits
    omp_set_num_threads(num_cores);     // number of cores

    std::vector<idx> qubits(N); // initial state
    ket psi = mket(qubits);
    ket result = psi;
    Timer<> t;
    for (idx i = 0; i < N; ++i) {
        result = apply(result, gt.H, {i}); // apply Hadamard on qubit i
        // apply controlled rotations
        for (idx j = 2; j <= N - i; ++j) {
            cmat Rj(2, 2);
            Rj << 1, 0, 0, omega(std::pow(2, j));
            result = applyCTRL(result, Rj, {i + j - 1}, {i});
        }
    }
    // we have the qubits in reversed order, we must swap them
    for (idx i = 0; i < N / 2; ++i) {
        result = apply(result, gt.SWAP, {i, N - i - 1});
    }
    std::cout << num_cores << ", " << N << ", " << t.toc() << '\n';
}
