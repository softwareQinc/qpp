// Qudit dense coding
// Source: ./examples/dense_coding.cpp
#include <cmath>
#include <iostream>
#include <tuple>
#include <vector>
#include "qpp.h"

using namespace qpp;

int main()
{
    idx D = 3; // size of the system
    std::cout << ">> Qudit dense coding, D = " << D << std::endl;

    ket mes_AB = ket::Zero(D * D); // maximally entangled state resource
    for (idx i = 0; i < D; ++i)
        mes_AB += mket({i, i}, D);
    mes_AB /= std::sqrt((double) D);

    // circuit that measures in the qudit Bell basis
    cmat Bell_AB = adjoint(gt.CTRL(gt.Xd(D), {0}, {1}, 2, D)
                           * kron(gt.Fd(D), gt.Id(D)));

    // equal probabilities of choosing a message
    idx m_A = randidx(0, D * D - 1);
    std::vector<idx> midx = n2multiidx(m_A, {D, D});
    std::cout << ">> Alice sent: " << m_A << " -> ";
    std::cout << disp(midx, " ") << std::endl;

    // Alice's operation
    cmat U_A = powm(gt.Zd(D), midx[0]) * powm(adjoint(gt.Xd(D)), midx[1]);

    // Alice encodes the message
    ket psi_AB = apply(mes_AB, U_A, {0}, D);

    // Bob measures the joint system in the qudit Bell basis
    psi_AB = apply(psi_AB, Bell_AB, {0, 1}, D);

    auto measured = measure(psi_AB, gt.Id(D * D));
    std::cout << ">> Bob's measurement probabilities: ";
    std::cout << disp(std::get<1>(measured), ", ") << std::endl;

    // Bob samples according to the measurement probabilities
    idx m_B = std::get<0>(measured);
    std::cout << ">> Bob received: ";
    std::cout << m_B << " -> " << disp(n2multiidx(m_B, {D, D}), " ")
              << std::endl;
}
