// MATLAB input/output
// Source: ./examples/matlab_io.cpp
#include <iostream>

#include "qpp.h"
#include "MATLAB/matlab.hpp" // must be explicitly included

int main() {
    using namespace qpp;
    // interfacing with MATLAB
    cmat rho = randrho(256); // an 8 qubit density operator
    saveMATLAB(rho, "rho.mat", "rho", "w");
    cmat loaded_rho = loadMATLAB<cmat>("rho.mat", "rho");
    // display the difference in norm, should be 0
    std::cout << ">> Norm difference MATLAB load/save: ";
    std::cout << norm(loaded_rho - rho) << '\n';
}
