// Qudit teleporation
// Source: ./examples/teleport_qudit.cpp
#include <cmath>
#include <iostream>
#include <tuple>
#include <vector>

#include "qpp.h"

int main() {
    using namespace qpp;
    idx D = 3; // size of the system
    std::cout << ">> Qudit teleportation, D = " << D << '\n';

    ket mes_AB = st.mes(D); // maximally entangled state resource

    // circuit that measures in the qudit Bell basis
    cmat Bell_aA =
        adjoint(gt.CTRL(gt.Xd(D), {0}, {1}, 2, D) * kron(gt.Fd(D), gt.Id(D)));

    ket psi_a = randket(D); // random qudit state
    std::cout << ">> Initial state:\n";
    std::cout << disp(psi_a) << '\n';

    ket input_aAB = kron(psi_a, mes_AB); // joint input state aAB
    // output before measurement
    ket output_aAB = apply(input_aAB, Bell_aA, {0, 1}, D);

    // measure on aA
    auto measured_aA = measure(output_aAB, gt.Id(D * D), {0, 1}, D);
    idx m = std::get<0>(measured_aA); // measurement result

    std::vector<idx> midx = n2multiidx(m, {D, D});
    std::cout << ">> Alice's measurement result: ";
    std::cout << m << " -> " << disp(midx, " ") << '\n';
    std::cout << ">> Alice's measurement probabilities: ";
    std::cout << disp(std::get<1>(measured_aA), ", ") << '\n';

    // conditional result on B before correction
    ket output_m_B = std::get<2>(measured_aA)[m];

    // perform the correction on B
    cmat correction_B =
        powm(gt.Zd(D), midx[0]) * powm(adjoint(gt.Xd(D)), midx[1]);
    std::cout << ">> Bob must apply the correction operator Z^" << midx[0]
              << " X^" << (D - midx[1]) % D << '\n';
    ket psi_B = correction_B * output_m_B;

    // display the output
    std::cout << ">> Bob's final state (after correction):\n";
    std::cout << disp(psi_B) << '\n';

    // verification
    std::cout << ">> Norm difference: " << norm(psi_B - psi_a) << '\n';
}
