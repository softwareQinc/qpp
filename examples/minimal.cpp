// Minimal example
// Source: ./examples/minimal.cpp
#include <iostream>
#include "qpp.h"

using namespace qpp;

int main()
{
    std::cout << "Hello Quantum++!\nThis is the |0> state:\n";
    std::cout << disp(st.z0);
}
