// Input/output
// Source: ./examples/input_output.cpp
#include <qpp.h>
#include <MATLAB/matlab.h> // must be explicitly included

using namespace qpp;

int main()
{
    // Quantum++ native input/output
    cmat rho = randrho(256); // an 8 qubit density operator
    save(rho, "rho.dat"); // save it
    cmat loaded_rho = load<cmat>("rho.dat"); // load it back
    // display the difference in norm, should be 0
    std::cout << ">> Norm difference load/save: ";
    std::cout << norm(loaded_rho - rho) << std::endl;

    // interfacing with MATLAB
    saveMATLABmatrix(rho, "rho.mat", "rho", "w");
    loaded_rho = loadMATLABmatrix<cmat>("rho.mat", "rho");
    // display the difference in norm, should be 0
    std::cout << ">> Norm difference MATLAB load/save: ";
    std::cout << norm(loaded_rho - rho) << std::endl;
}
