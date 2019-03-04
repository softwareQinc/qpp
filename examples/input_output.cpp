// Input/output
// Source: ./examples/input_output.cpp
#include <iostream>

#include "qpp.h"

int main() {
    using namespace qpp;
    // Quantum++ native input/output
    cmat rho = randrho(256);                 // an 8 qubit density operator
    save(rho, "rho.dat");                    // save it
    cmat loaded_rho = load<cmat>("rho.dat"); // load it back
    // display the difference in norm, should be 0
    std::cout << ">> Norm difference load/save: ";
    std::cout << norm(loaded_rho - rho) << '\n';
}
