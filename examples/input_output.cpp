// Input/output
// Source: ./examples/input_output.cpp
#include <iostream>
#include <fstream>

#include "qpp.h"

int main() {
    using namespace qpp;
    // Quantum++ native input/output
    cmat rho = randrho(256); // an 8 qubit density operator
    {
        std::ofstream fout("rho.dat", std::ios::out | std::ios::binary);
        save(rho, fout); // save it
    }
    std::ifstream fin("rho.dat", std::ios::in | std::ios::binary);
    cmat loaded_rho = load<cmat>(fin); // load it back
    // display the difference in norm, should be 0
    std::cout << ">> Norm difference load/save: ";
    std::cout << norm(loaded_rho - rho) << '\n';
}
