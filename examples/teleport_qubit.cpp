// Qubit teleporation
// Source: ./examples/teleport_qubit.cpp
#include <iostream>
#include <tuple>

#include "qpp.h"

int main() {
    using namespace qpp;
    // input state
    ket psi_a = randket();

    std::cout << ">> Qubit teleporation\n";
    std::cout << ">> Initial state:\n" << disp(psi_a) << '\n';

    // the entangled resource
    ket phi_AB = st.b00;

    // global input state
    ket input_aAB = kron(psi_a, phi_AB);

    // apply a CNOT on qubits 'AB' followed by an H on qubit 'a'
    input_aAB = applyCTRL(input_aAB, gt.X, {0}, {1});
    input_aAB = apply(input_aAB, gt.H, {0});

    // measure the aA part
    auto results = measure_seq(input_aAB, {0, 1});

    // measurement results
    idx z = std::get<0>(results)[0];
    idx x = std::get<0>(results)[1];
    std::cout << ">> Alice's measurement results: ";
    std::cout << "x = " << x << " and z = " << z;

    // probability of obtaining the measurement results x and z
    double p = std::get<1>(results);
    std::cout << ", obtained with probability: " << p << '\n';

    // the output state (before correction)
    ket out_B = std::get<2>(results);
    std::cout << ">> Bob's state (before correction):\n";
    std::cout << disp(out_B) << '\n';

    // perform the correction on B
    out_B = powm(gt.Z, z) * powm(gt.X, x) * out_B;
    std::cout << ">> Bob must apply the correction operator Z^" << z << " X^"
              << x << '\n';

    // display the output
    std::cout << ">> Bob's final state (after correction):\n";
    std::cout << disp(out_B) << '\n';

    // verification
    std::cout << ">> Norm difference: " << norm(out_B - psi_a) << '\n';
}
