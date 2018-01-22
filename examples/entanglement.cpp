// Entanglement
// Source: ./examples/entanglement.cpp
#include <iostream>

#include "qpp.h"

int main() {
    using namespace qpp;
    cmat rho = 0.2 * st.pb00 + 0.8 * st.pb11;
    std::cout << ">> State rho:\n";
    std::cout << disp(rho) << '\n';

    std::cout << ">> Concurrence of rho: " << concurrence(rho) << '\n';
    std::cout << ">> Negativity of rho: " << negativity(rho, {2, 2}) << '\n';
    std::cout << ">> Logarithimc negativity of rho: "
              << lognegativity(rho, {2, 2}) << '\n';

    ket psi = 0.8 * 00_ket + 0.6 * 11_ket;

    // apply some local random unitaries
    psi = kron(randU(), randU()) * psi;

    std::cout << ">> State psi:\n";
    std::cout << disp(psi) << '\n';

    std::cout << ">> Entanglement of psi: " << entanglement(psi, {2, 2})
              << '\n';
    std::cout << ">> Concurrence of psi: " << concurrence(prj(psi)) << '\n';
    std::cout << ">> G-Concurrence of psi: " << gconcurrence(psi) << '\n';

    std::cout << ">> Schmidt coefficients of psi:\n";
    std::cout << disp(schmidtcoeffs(psi, {2, 2})) << '\n';

    std::cout << ">> Schmidt probabilities of psi:\n";
    std::cout << disp(schmidtprobs(psi, {2, 2}), ", ") << '\n';

    cmat UA = schmidtA(psi, {2, 2});
    cmat UB = schmidtB(psi, {2, 2});

    std::cout << ">> Schmidt vectors on Alice's side:\n";
    std::cout << disp(UA) << '\n';

    std::cout << ">> Schmidt vectors on Bob's side:\n";
    std::cout << disp(UB) << '\n';

    std::cout << ">> State psi in the Schmidt basis:\n";
    std::cout << disp(adjoint(kron(UA, UB)) * psi) << '\n';

    // reconstructed state
    ket psi_from_schmidt =
        schmidtcoeffs(psi, {2, 2})(0) * kron(UA.col(0), UB.col(0)) +
        schmidtcoeffs(psi, {2, 2})(1) * kron(UA.col(1), UB.col(1));
    std::cout << ">> State psi reconstructed from the Schmidt decomposition:\n";
    std::cout << disp(psi_from_schmidt) << '\n';

    // verification
    std::cout << ">> Norm difference: " << norm(psi - psi_from_schmidt) << '\n';
}
