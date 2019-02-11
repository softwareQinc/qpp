// Used for testing, do not use it as an example
#include <iostream>

#include "qpp.h"

int main() {
    /////////// testing ///////////
    using namespace qpp;

    QCircuit qc(10, 0);
    idx i;
    idx N = 1000;

    std::cout << "1 qubit gates...\n";
    try {
        for (i = 0; i < N; ++i)
            qc.gate(rand<cmat>(2, 2), 0);
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

    try {
        std::cout << "Executing...\n";
        QEngine qEngine{qc};
        for (auto&& elem : qc) {
            std::cout << elem << '\n';
            qEngine.execute(elem);
        }
    }
    catch(std::exception& e)
    {
        std::cerr << e.what() << '\n';
    }
}
