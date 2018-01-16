// Exceptions
// Source: ./examples/exceptions.cpp
#include <exception>
#include <iostream>

#include "qpp.h"

int main() {
    using namespace qpp;
    cmat rho = randrho(16); // 4 qubits (subsystems)
    try {
        // the line below throws qpp::exception::SubsysMismatchDims
        double mInfo = qmutualinfo(rho, {0}, {4});
        std::cout << ">> Mutual information between first and last subsystem: ";
        std::cout << mInfo << '\n';
    } catch (const std::exception& e) {
        std::cout << ">> Exception caught: " << e.what() << '\n';
    }
}
