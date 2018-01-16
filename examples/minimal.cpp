// Minimal example
// Source: ./examples/minimal.cpp
#include <iostream>

#include "qpp.h"

int main() {
    using namespace qpp;
    std::cout << "Hello Quantum++!\nThis is the |0> state:\n";
    std::cout << disp(st.z0) << '\n';
}
