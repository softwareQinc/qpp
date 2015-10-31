// Entanglement
// Source: ./examples/entanglement.cpp
#include <qpp.h>

using namespace qpp;
using std::cout;
using std::endl;

int main()
{
    cmat rho = 0.2 * st.pb00 + 0.8 * st.pb11;
    cout << ">> State rho: " << endl;
    cout << disp(rho) << endl;

    cout << ">> Concurrence of rho: " << concurrence(rho) << endl;
    cout << ">> Negativity of rho: " << negativity(rho, {2, 2}) << endl;
    cout << ">> Logarithimc negativity of rho: "
            << lognegativity(rho, {2, 2}) << endl;

    ket psi = 0.8 * mket({0, 0}) + 0.6 * mket({1, 1});

    // apply some local random unitaries
    psi = kron(randU(2), randU(2)) * psi;

    cout << ">> State psi: " << endl;
    cout << disp(psi) << endl;

    cout << ">> Entanglement of psi: " << entanglement(psi, {2, 2}) << endl;
    cout << ">> Concurrence of psi: " << concurrence(prj(psi)) << endl;
    cout << ">> G-Concurrence of psi: " << gconcurrence(psi) << endl;

    cout << ">> Schmidt coefficients of psi: " << endl;
    cout << disp(schmidtcoeffs(psi, {2, 2})) << endl;

    cout << ">> Schmidt probabilities of psi: " << endl;
    cout << disp(schmidtprobs(psi, {2, 2}), ", ") << endl;

    cmat UA = schmidtA(psi, {2, 2});
    cmat UB = schmidtB(psi, {2, 2});

    cout << ">> Schmidt vectors on Alice's side: " << endl;
    cout << disp(UA) << endl;

    cout << ">> Schmidt vectors on Bob's side: " << endl;
    cout << disp(UB) << endl;

    cout << ">> State psi in the Schmidt basis: " << endl;
    cout << disp(adjoint(kron(UA, UB)) * psi) << endl;

    // reconstructed state
    ket psi_from_schmidt =
            schmidtcoeffs(psi, {2, 2})(0) * kron(UA.col(0), UB.col(0))
            + schmidtcoeffs(psi, {2, 2})(1)
              * kron(UA.col(1), UB.col(1));
    cout << ">> State psi reconstructed from the Schmidt decomposition:\n";
    cout << disp(psi_from_schmidt) << endl;

    // verification
    cout << ">> Norm difference: " << norm(psi - psi_from_schmidt) << endl;
}
