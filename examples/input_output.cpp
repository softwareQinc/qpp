// Source: ./examples/input_output.cpp
//
// Input/output

#include <fstream>
#include <iostream>

#include "qpp/qpp.hpp"

int main() {
    using namespace qpp;

    // Quantum++ input/output
    cmat rho = randrho(256); // an 8 qubit density operator
    {
        std::ofstream fout("rho.txt");
        save(rho, fout); // save it
    }
    std::ifstream fin("rho.txt");
    cmat loaded_rho = load<cmat>(fin); // load it back
    // display the difference in norm, should be 0
    std::cout << ">> Norm difference load/save: ";
    std::cout << norm(loaded_rho - rho) << '\n';
}
