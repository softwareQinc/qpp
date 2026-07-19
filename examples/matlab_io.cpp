// Source: ./examples/matlab_io.cpp
//
// MATLAB input/output

#include <iostream>

#include <qpp/qpp.hpp>

#include "qpp/MATLAB/matlab.hpp" // must be explicitly included

int main() {
    using namespace qpp;
    namespace fs = std::filesystem;
    fs::path filename = "rho.mat";
    fs::path rho_mat = fs::temp_directory_path() / filename;
    std::string matlab_name = "rho";

    // interfacing with MATLAB
    cmat rho = randrho(256); // an 8 qubit density operator
    save_MATLAB(rho, rho_mat.string(), matlab_name, "w");
    cmat loaded_rho = load_MATLAB<cmat>(rho_mat.string(), matlab_name);
    // display the difference in norm, should be 0
    std::cout << ">> Norm difference MATLAB load/save: ";
    std::cout << norm(loaded_rho - rho) << '\n';
}
