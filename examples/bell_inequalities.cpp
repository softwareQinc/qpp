// Bell inequalities (CHSH) violation
// Source: ./examples/bell_inequalities.cpp
#include <qpp.h>

using namespace qpp;
using std::cout;
using std::endl;

int main()
{
    ket psi = st.b11; // Bell singlet state (|01> - |10>) / sqrt(2)

    // detector settings (Q and R on Alice's side, S and T on Bob's side)
    cmat Q = gt.Z;
    cmat R = gt.X;
    cmat S = (-gt.Z - gt.X) / sqrt(2);
    cmat T = (gt.Z - gt.X) / sqrt(2);

    // number of "experiments" for each of the 4 detector settings
    idx N = 10000;
    cout << ">> Number N of experiments for each of the 4 measurement";
    cout << " settings = " << N << endl;

    idx statistics[4][4] = {0}; // total statistics
    long E[4] = {0}; // experimental estimate

    idx gate_idx = 0; // gate index (0, 1, 2 or 3)
    for (auto&& gateA: {Q, R}) // measure Alice's side
    {
        // eigenvalues, so we know the order
        dyn_col_vect<double> evalsA = hevals(gateA);
        cmat basisA = hevects(gateA); // eigenvectors, ordered by eigenvalues
        for (auto&& gateB: {S, T}) // measure Bob's side
        {
            // eigenvalues, so we know the order
            dyn_col_vect<double> evalsB = hevals(gateB);
            cmat basisB = hevects(gateB);
            for (idx i = 0; i < N; ++i) // repeat the "experiment" N times
            {
                auto measuredA = measure(psi, basisA, {0});
                idx mA = std::get<0>(measuredA); // result on A
                // the eigenvalues corresponding to the measurement results
                double evalA = evalsA[mA];
                // resulting state on B
                ket psiB = std::get<2>(measuredA)[mA];
                auto measuredB = measure(psiB, basisB);
                idx mB = std::get<0>(measuredB); // measurement result B
                double evalB = evalsB[mB];
                // count the correlations
                if (evalA > 0 && evalB > 0)        // +1 +1 correlation
                {
                    statistics[gate_idx][0]++;
                    E[gate_idx]++;
                }
                else if (evalA > 0 && evalB < 0)  // +1 -1 anti-correlation
                {
                    statistics[gate_idx][1]++;
                    E[gate_idx]--;
                }
                else if (evalA < 0 && evalB > 0)  // -1 +1 anti-correlation
                {
                    statistics[gate_idx][2]++;
                    E[gate_idx]--;
                }
                else if (evalA < 0 && evalB < 0) // -1 -1 correlation
                {
                    statistics[gate_idx][3]++;
                    E[gate_idx]++;
                }
            } // N experiments are done
            ++gate_idx;
        }
    }
    cout << "[N++ | N+- | N-+ | N-- | (N++ + N-- - N+- - N-+)]\n";
    cout << "QS: " << disp(statistics[0], 4, " ");
    cout << "  " << E[0] << endl;
    cout << "QT: " << disp(statistics[1], 4, " ");
    cout << "  " << E[1] << endl;
    cout << "RS: " << disp(statistics[2], 4, " ");
    cout << "  " << E[2] << endl;
    cout << "RT: " << disp(statistics[3], 4, " ");
    cout << "  " << E[3] << endl;

    // Experimental average
    double exp_avg = (E[0] - E[1] + E[2] + E[3]) / static_cast<double>(N);
    cout << ">> Experimental estimate of <QS> + <RS> + <RT> - <QT> = ";
    cout << exp_avg << endl;

    // Theoretical average
    double th_avg = (adjoint(psi) *
                     (kron(Q, S) + kron(R, S) + kron(R, T) - kron(Q, T)) *
                     psi).value().real();
    cout << ">> Theoretical value of <QS> + <RS> + <RT> - <QT> = ";
    cout << th_avg << endl;

    cout << ">> 2 * sqrt(2) = " << 2 * sqrt(2) << endl;
}