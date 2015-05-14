// Bell inequality (CHSH) violations
// Source: ./examples/bell_inequalities.cpp
#include <qpp.h>

using namespace qpp;

int main()
{
    ket psi = st.b00; // Measure Bell0 state
    idx N = 20; // number of measurements each party does

    cmat Q = gt.Z;
    cmat R = gt.X;
    cmat S = (-gt.Z - gt.X) / std::sqrt(2);
    cmat T = (gt.Z - gt.X) / std::sqrt(2);

    double qs, rs, rt, qt; // averages
    qs = rs = rt = qt = 0;
    // we compute <Q x S> + <R x S> + <R x T> - <Q x T>
    for (idx i = 0; i < N; ++i)
    {
        // Q x S
        auto result = measure(psi, kron(Q, gt.Id2));

    }
    std::cout << "|psi> = " << std::endl << disp(psi) << std::endl;
    std::cout << "Entanglement(|psi>) = " << entanglement(psi, {2, 2})
    << std::endl;
    std::cout << "<Q x S> + <R x S> + <R x T> - <Q x T> = ";
    std::cout << qs + rs + rt - qt << std::endl;
    std::cout << "2 * sqrt(2) = " << 2 * std::sqrt(2) << std::endl;

    std::vector<double> v{1.1,2.2,3.3};
    std::cout << sum(v) << std::endl;
    std::cout << prod(v) << std::endl;
}