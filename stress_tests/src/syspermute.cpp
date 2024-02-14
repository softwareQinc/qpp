// System permutation stress test on a pure state of n qubits

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

#include <omp.h>

#include "qpp/qpp.h"

int main(int argc, char** argv) {
    using namespace qpp;
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " <n_cores> <n_qubits>\n";
        exit(EXIT_FAILURE);
    }

    int num_cores = std::stoi(argv[1]);                     // number of cores
    idx n = std::stoi(argv[2]);                             // number of qubits
    idx D = static_cast<idx>(std::llround(std::pow(2, n))); // dimension 2^n
    omp_set_num_threads(num_cores);                         // number of cores

    ket randket = ket::Random(D);       // random ket
    std::vector<idx> subsys_syspermute; // partially syspermute reversed order
    for (idx i = 0; i < n; ++i) {
        subsys_syspermute.push_back(i);
    }

    Timer<> t; // start timing
    syspermute(randket, subsys_syspermute);
    std::cout << num_cores << ", " << n << ", " << t.toc() << '\n';
}
