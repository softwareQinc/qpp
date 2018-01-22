// Gates and states
// Source: ./examples/gates_states.cpp
#include <iostream>

#include "qpp.h"

int main() {
    using namespace qpp;
    ket psi = st.z0; // |0> state
    cmat U = gt.X;
    ket result = U * psi;

    std::cout << ">> The result of applying the bit-flip gate X on |0> is:\n";
    std::cout << disp(result) << '\n';

    psi = 10_ket; // |10> state
    U = gt.CNOT;  // Controlled-NOT
    result = U * psi;

    std::cout << ">> The result of applying the gate CNOT on |10> is:\n";
    std::cout << disp(result) << '\n';

    U = randU();
    std::cout << ">> Generating a random one-qubit gate U:\n";
    std::cout << disp(U) << '\n';

    result = applyCTRL(psi, U, {0}, {1}); // Controlled-U
    std::cout << ">> The result of applying the CTRL-U gate on |10> is:\n";
    std::cout << disp(result) << '\n';
}
