// Bell inequalities (CHSH) violation
// Source: ./examples/bell_inequalities.cpp
#include <qpp.h>
using namespace qpp;
using std::cout;
using std::endl;

//TODO: finish the example, incomplete now

int main()
{
    ket psi = st.b11; // Measure the singlet Bell state (|01>-|10>)/sqrt(2)
    idx N = 1; // number of measurements each party does

    // gates
    cmat Q = gt.Z;
    cmat R = gt.X;
    cmat S = (-gt.Z - gt.X) / std::sqrt(2);
    cmat T = (gt.Z - gt.X) / std::sqrt(2);

    //Q = S = gt.Z;
    //R = T = gt.X;

    idx statistics[4][4] = {0}; // total statistics
    int E[4] = {0}; // experimental estimate

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
            for (idx i = 0; i < N; ++i)
            {
                auto measurementA = measure(psi, basisA, {0});
                auto mA = std::get<0>(measurementA); // result on A
                // the eigenvalues corresponding to the measurement results
                short evalA = static_cast<short>(std::round(evalsA[mA]));
                // resulting state on B
                auto rhoB = std::get<2>(measurementA)[mA];
                auto measurementB = measure(rhoB, basisB);
                auto mB = std::get<0>(measurementB); // measurement result B
                short evalB = static_cast<short>(std::round(evalsB[mB]));
                // count the coincidences
                if (evalA == 1 && evalB == 1)        // +1 +1
                {
                    statistics[gate_idx][0]++;
                    E[gate_idx]++;
                }
                else if (evalA == 1 && evalB == -1)  // +1 -1
                {
                    statistics[gate_idx][1]++;
                    E[gate_idx]--;
                }
                else if (evalA == -1 && evalB == 1)  // -1 +1
                {
                    statistics[gate_idx][2]++;
                    E[gate_idx]--;
                }
                else if (evalA == -1 && evalB == -1) // -1 -1
                {
                    statistics[gate_idx][3]++;
                    E[gate_idx]++;
                }
            }
            ++gate_idx;
        }
    }
    std::cout << "Coincidences N = " << N << std::endl;
    std::cout << "(N++ N+- N-+ N-- E)" << std::endl;
    std::cout << "QS: " << disp(statistics[0], 4, " ");
    std::cout << "  " << E[0] << std::endl;
    std::cout << "QT: " << disp(statistics[1], 4, " ");
    std::cout << "  " << E[1] << std::endl;
    std::cout << "RS: " << disp(statistics[2], 4, " ");
    std::cout << "  " << E[2] << std::endl;
    std::cout << "RT: " << disp(statistics[3], 4, " ");
    std::cout << "  " << E[3] << std::endl;

    double val = (E[0] - E[1] + E[2] + E[3]) / static_cast<double>(N);
    std::cout << "<QS> + <RS> + <RT> - <QT> = " << val << std::endl;
    std::cout << "Theoretical value 2 * sqrt(2) = " << 2 * std::sqrt(2);
}