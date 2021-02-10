// Qudit teleportation
// Source: ./examples/teleport_qudit.cpp
// See also: ./examples/teleport_qubit.cpp
#include <iostream>
#include <tuple>
#include <vector>

#include "qpp.h"

int main() {
    using namespace qpp;
    idx d = 3; // size of the system
    std::cout << ">> Qudit teleportation, d = " << d << '\n';

    ket mes_AB = st.mes(d); // maximally entangled state resource

    // circuit used to measure in the qudit Bell basis
    cmat Bell_aA = kron(adjoint(gt.Fd(d)), gt.Id(d)) *
                   gt.CTRL(adjoint(gt.Xd(d)), {0}, {1}, 2, d);

    ket psi_a = randket(d); // random qudit state
    std::cout << ">> Initial state:\n";
    std::cout << disp(psi_a) << '\n';

    ket input_aAB = kron(psi_a, mes_AB); // joint input state aAB
    // output before measurement
    ket output_aAB = apply(input_aAB, Bell_aA, {0, 1}, d);

    // measure on aA
    auto measured_aA = measure(output_aAB, gt.Id(d * d), {0, 1}, d);
    idx m = std::get<RES>(measured_aA); // measurement result

    std::vector<idx> midx = n2multiidx(m, {d, d});
    std::cout << ">> Alice's measurement result: ";
    std::cout << m << " -> " << disp(midx, " ") << '\n';
    std::cout << ">> Alice's measurement probabilities: ";
    std::cout << disp(std::get<PROB>(measured_aA), ", ") << '\n';

    // conditional result on B before correction
    ket output_m_B = std::get<ST>(measured_aA)[m];

    // perform the correction on B
    cmat correction_B =
        powm(gt.Zd(d), midx[0]) * powm(adjoint(gt.Xd(d)), midx[1]);
    std::cout << ">> Bob must apply the correction operator Z^" << midx[0]
              << " X^" << (d - midx[1]) % d << '\n';
    ket psi_B = correction_B * output_m_B;

    // display the output
    std::cout << ">> Bob's final state (after correction):\n";
    std::cout << disp(psi_B) << '\n';

    // verification
    std::cout << ">> Norm difference: " << norm(psi_B - psi_a) << '\n';
}
