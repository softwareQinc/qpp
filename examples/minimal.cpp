// Source: ./examples/minimal.cpp
//
// Minimal example

#include <iostream>

#include <qpp/qpp.hpp>

int main() {
    using namespace qpp;

    std::cout << "Hello Quantum++!\nThis is the |0> state:\n";
    std::cout << disp(0_ket) << '\n';

    std::cout << "This is some random ket in Dirac notation\n";
    auto dirac_disp_opt = IOManipDiracOpts{}
                              .set_amplitudes_after(false)
                              .set_mul_op("")
                              .set_plus_op(" + ");
    std::cout << disp(dirac(randket(2)), dirac_disp_opt) << '\n';
}
