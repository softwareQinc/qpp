// System permutation stress test on a pure state of n qubits
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

    ket randket = ket::Random(D);         // random ket
    std::vector<idx> subsys_syspermute;   // partially syspermute reversed order
    for (idx i = 0; i < n; ++i)
        subsys_syspermute.push_back(i);

    Timer<> t;
    syspermute(randket, subsys_syspermute);
    std::cout << num_cores << ", " << n << ", " << t.toc() << '\n';
}
