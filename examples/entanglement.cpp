// Entanglement
// Source: ./examples/entanglement.cpp
#include <iostream>
#include "qpp.h"

using namespace qpp;

int main()
{
    cmat rho = 0.2 * st.pb00 + 0.8 * st.pb11;
    std::cout << ">> State rho: " << std::endl;
    std::cout << disp(rho) << std::endl;

    std::cout << ">> Concurrence of rho: " << concurrence(rho) << std::endl;
    std::cout << ">> Negativity of rho: " << negativity(rho, {2, 2})
              << std::endl;
    std::cout << ">> Logarithimc negativity of rho: "
              << lognegativity(rho, {2, 2}) << std::endl;

    ket psi = 0.8 * mket({0, 0}) + 0.6 * mket({1, 1});

    // apply some local random unitaries
    psi = kron(randU(2), randU(2)) * psi;

    std::cout << ">> State psi: " << std::endl;
    std::cout << disp(psi) << std::endl;

    std::cout << ">> Entanglement of psi: " << entanglement(psi, {2, 2})
              << std::endl;
    std::cout << ">> Concurrence of psi: " << concurrence(prj(psi))
              << std::endl;
    std::cout << ">> G-Concurrence of psi: " << gconcurrence(psi)
              << std::endl;

    std::cout << ">> Schmidt coefficients of psi: " << std::endl;
    std::cout << disp(schmidtcoeffs(psi, {2, 2})) << std::endl;

    std::cout << ">> Schmidt probabilities of psi: " << std::endl;
    std::cout << disp(schmidtprobs(psi, {2, 2}), ", ") << std::endl;

    cmat UA = schmidtA(psi, {2, 2});
    cmat UB = schmidtB(psi, {2, 2});

    std::cout << ">> Schmidt vectors on Alice's side: " << std::endl;
    std::cout << disp(UA) << std::endl;

    std::cout << ">> Schmidt vectors on Bob's side: " << std::endl;
    std::cout << disp(UB) << std::endl;

    std::cout << ">> State psi in the Schmidt basis: " << std::endl;
    std::cout << disp(adjoint(kron(UA, UB)) * psi) << std::endl;

    // reconstructed state
    ket psi_from_schmidt =
            schmidtcoeffs(psi, {2, 2})(0) * kron(UA.col(0), UB.col(0))
            + schmidtcoeffs(psi, {2, 2})(1)
              * kron(UA.col(1), UB.col(1));
    std::cout << ">> State psi reconstructed from the Schmidt decomposition:\n";
    std::cout << disp(psi_from_schmidt) << std::endl;

    // verification
    std::cout << ">> Norm difference: " << norm(psi - psi_from_schmidt)
              << std::endl;
}
