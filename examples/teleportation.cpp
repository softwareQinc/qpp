// Qudit teleporation
// Source: ./examples/teleportation.cpp
#include <qpp.h>
using namespace qpp;
using std::cout;
using std::endl;

int main()
{
    idx D = 3; // size of the system
    cout << ">> Qudit teleportation, D = " << D << endl;

    ket mes_AB = ket::Zero(D * D); // maximally entangled state resource
    for (idx i = 0; i < D; ++i)
        mes_AB += mket({i, i}, D);
    mes_AB /= std::sqrt((double) D);

    // circuit that measures in the qudit Bell basis
    cmat Bell_aA = adjoint(gt.CTRL(gt.Xd(D), {0}, {1}, 2, D)
                           * kron(gt.Fd(D), gt.Id(D)));

    ket psi_a = randket(D); // random qudit state
    cout << ">> Initial state:" << endl;
    cout << disp(psi_a) << endl;

    ket input_aAB = kron(psi_a, mes_AB); // joint input state aAB
    // output before measurement
    ket output_aAB = apply(input_aAB, Bell_aA, {0, 1}, D);

    // measure on aA
    auto measured_aA = measure(output_aAB, gt.Id(D * D), {0, 1}, D);
    idx m = std::get<0>(measured_aA); // measurement result

    auto midx = n2multiidx(m, {D, D});
    cout << ">> Alice's measurement result: ";
    cout << m << " -> " << disp(midx, " ") << endl;
    cout << ">> Alice's measurement probabilities: ";
    cout << disp(std::get<1>(measured_aA), ", ") << endl;

    // conditional result on B before correction
    ket output_m_B = std::get<2>(measured_aA)[m];
    // correction operator
    cmat correction_B = powm(gt.Zd(D), midx[0]) *
                        powm(adjoint(gt.Xd(D)), midx[1]);
    // apply correction on B
    cout << ">> Bob must apply the correction operator Z^" << midx[0]
        << " X^" << (D - midx[1]) % D << endl;
    ket psi_B = correction_B * output_m_B;

    cout << ">> Bob's final state (after correction): " << endl;
    cout << disp(psi_B) << endl;

    // verification
    cout << ">> Norm difference: " << norm(psi_B - psi_a) << endl;
}
