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
    cmat S = (-gt.Z - gt.X) / std::sqrt(2);
    cmat T = (gt.Z - gt.X) / std::sqrt(2);

    // number of "experiments" for each of the 4 detector settings
    idx N = 10000;
    std::cout << "Number N of experiments for each of the 4 measurement";
    std::cout << " settings = " << N << std::endl;

    idx statistics[4][4] = {0}; // total statistics
    long E[4] = {0}; // experimental estimate

    idx gate_idx = 0; // gate index (0, 1, 2 or 3)
    for (auto&& gateA: {Q, R}) // measure Alice's side
    {
        auto evalsA = hevals(gateA); // eigenvalues, so we know the order
        auto basisA = hevects(gateA); // eigenvectors, ordered by eigenvalues
        for (auto&& gateB: {S, T}) // measure Bob's side
        {
            // eigenvalues, so we know the order
            auto evalsB = hevals(gateB);
            auto basisB = hevects(gateB);
            for (idx i = 0; i < N; ++i) // repeat the "experiment" N times
            {
                auto measurementA = measure(psi, basisA, {0});
                auto mA = std::get<0>(measurementA); // result on A
                // the eigenvalues corresponding to the measurement results
                auto evalA = evalsA[mA];
                // resulting state on B
                auto rhoB = std::get<2>(measurementA)[mA];
                auto measurementB = measure(rhoB, basisB);
                auto mB = std::get<0>(measurementB); // measurement result B
                auto evalB = evalsB[mB];
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
    std::cout << "[N++ | N+- | N-+ | N-- | (N++ + N-- - N+- - N-+)]\n";
    std::cout << "QS: " << disp(statistics[0], 4, " ");
    std::cout << "  " << E[0] << std::endl;
    std::cout << "QT: " << disp(statistics[1], 4, " ");
    std::cout << "  " << E[1] << std::endl;
    std::cout << "RS: " << disp(statistics[2], 4, " ");
    std::cout << "  " << E[2] << std::endl;
    std::cout << "RT: " << disp(statistics[3], 4, " ");
    std::cout << "  " << E[3] << std::endl;

    // Experimental average
    auto exp_avg = (E[0] - E[1] + E[2] + E[3]) / static_cast<double>(N);
    std::cout << "Experimental estimate of <QS> + <RS> + <RT> - <QT> = ";
    std::cout << exp_avg << std::endl;

    // Theoretical average
    double th_avg = (adjoint(psi) *
                     (kron(Q, S) + kron(R, S) + kron(R, T) - kron(Q, T)) *
                     psi).value().real();
    std::cout << "Theoretical value of <QS> + <RS> + <RT> - <QT> = ";
    std::cout << th_avg << std::endl;

    std::cout << "2 * sqrt(2) = " << 2 * std::sqrt(2) << std::endl;
}