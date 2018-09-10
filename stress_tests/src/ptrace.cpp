// Partial trace stress test on a density matrix of n qubits
// Reduce the number of NQ_END in run.sh to accommodate for density matrices
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

    idx num_cores = std::stoi(argv[1]);   // number of cores
    idx n = std::stoi(argv[2]);           // number of qubits
    idx D = std::round(std::pow(2, n));   // dimension
    omp_set_num_threads(num_cores);       // number of cores

    cmat randcmat = cmat::Random(D, D);   // random matrix
    std::vector<idx> subsys_ptrace = {0}; // partial trace over first qubit

    Timer<> t; // start timing
    ptrace(randcmat, subsys_ptrace);
    std::cout << num_cores << ", " << n << ", " << t.toc() << '\n';
}
