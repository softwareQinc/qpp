// Qubit teleportation
// Source: ./examples/teleport_qubit.cpp
// See also: ./examples/teleport_qudit.cpp

#include <iostream>
#include <tuple>

#include "qpp/qpp.hpp"

int main() {
    using namespace qpp;

    // input state
    ket psi_a = randket();

    std::cout << ">> Qubit teleportation\n";
    std::cout << ">> Initial state:\n" << disp(psi_a) << '\n';

    // the entangled resource
    ket phi_AB = st.b00;

    // global input state
    ket input_aAB = kron(psi_a, phi_AB);

    // apply a CNOT on qubits 'AB' followed by an H on qubit 'a'
    input_aAB = applyCTRL(input_aAB, gt.X, {0}, {1});
    input_aAB = apply(input_aAB, gt.H, {0});

    // measure the aA part
    auto [ms_aA, probs_aA, psi_B] = measure_seq(input_aAB, {0, 1});

    // measurement results
    idx z = ms_aA[0];
    idx x = ms_aA[1];
    std::cout << ">> Alice's measurement result: ";
    std::cout << "x = " << x << " and z = " << z;

    // probability of obtaining the measurement results x and z
    realT p = prod(probs_aA);
    std::cout << ", obtained with probability: " << p << '\n';

    // the output state (before correction)
    std::cout << ">> Bob's state (before correction):\n";
    std::cout << disp(psi_B) << '\n';

    // perform the correction on B
    psi_B = powm(gt.Z, z) * powm(gt.X, x) * psi_B;
    std::cout << ">> Bob must apply the correction operator Z^" << z << " X^"
              << x << '\n';

    // display the output
    std::cout << ">> Bob's final state (after correction):\n";
    std::cout << disp(psi_B) << '\n';

    // verification
    std::cout << ">> Norm difference: " << norm(psi_B - psi_a) << '\n';
}
