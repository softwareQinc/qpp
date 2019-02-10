// Used for testing, do not use it as an example
#include <iostream>

#include "qpp.h"

int main() {
    using namespace qpp;
    /////////// testing ///////////

    QCircuit qc(4, 0);
    idx i;
    idx N = 1000000;

    std::cout << "1 qubit gates...\n";
    try {
        for (i = 0; i < N; ++i)
            qc.gate(rand<cmat>(2,2), 0);
    } catch (...) {
        std::cerr << "Collision at step: " << i << '\n';
        throw;
    }
    std::cout << "2 qubit gates...\n";
    try {
        for (i = 0; i < N; ++i)
            qc.gate(randU(4), 0, 1);
    } catch (...) {
        std::cerr << "Collision at step: " << i << '\n';
        throw;
    }
    std::cout << "3 qubit gates...\n";
    try {
        for (i = 0; i < N; ++i)
            qc.gate(randU(8), 0, 1, 2);
    } catch (...) {
        std::cerr << "Collision at step: " << i << '\n';
        throw;
    }
}
