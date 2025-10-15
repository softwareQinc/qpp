// Source: ./examples/qram.cpp
//
// Quantumly-accessible Random Access Memory over classical data

#include <iostream>

#include "qpp/qpp.hpp"

int main() {
    using namespace qpp;

    std::cout << ">> qRAM over classical data\n";

    ket in = st.x0;  // |+> input state
    qram data{0, 1}; // qRAM data
    std::cout << ">> qRAM input:\n" << disp(in) << '\n';
    std::cout << ">> Classical data:\n"
              << disp(data, IOManipContainerOpts{}.set_sep(", ")) << '\n';
    ket out = qRAM(in, data); // qRAM output, automatically sets the dimension
                              // of the qRAM subsystem to 2
    std::cout << ">> qRAM output:\n" << disp(out) << '\n';
}
