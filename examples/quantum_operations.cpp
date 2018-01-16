// Quantum operations
// Source: ./examples/quantum_operations.cpp
#include <iostream>
#include <vector>

#include "qpp.h"

int main() {
    using namespace qpp;
    cmat rho = st.pb00; // projector onto the Bell state (|00> + |11>) / sqrt(2)
    std::cout << ">> Initial state:\n";
    std::cout << disp(rho) << '\n';

    // partial transpose of first subsystem
    cmat rhoTA = ptranspose(rho, {0});
    std::cout << ">> Eigenvalues of the partial transpose "
              << "of Bell-0 state are:\n";
    std::cout << disp(transpose(hevals(rhoTA))) << '\n';

    std::cout << ">> Measurement channel with 2 Kraus operators:\n";
    std::vector<cmat> Ks{st.pz0, st.pz1}; // 2 Kraus operators
    std::cout << disp(Ks[0]) << "\nand\n" << disp(Ks[1]) << '\n';

    std::cout << ">> Superoperator matrix of the channel:\n";
    std::cout << disp(kraus2super(Ks)) << '\n';

    std::cout << ">> Choi matrix of the channel:\n";
    std::cout << disp(kraus2choi(Ks)) << '\n';

    // apply the channel onto the first subsystem
    cmat rhoOut = apply(rho, Ks, {0});
    std::cout << ">> After applying the measurement channel "
              << "on the first qubit:\n";
    std::cout << disp(rhoOut) << '\n';

    // take the partial trace over the second subsystem
    cmat rhoA = ptrace(rhoOut, {1});
    std::cout << ">> After partially tracing down the second subsystem:\n";
    std::cout << disp(rhoA) << '\n';

    // compute the von-Neumann entropy
    double ent = entropy(rhoA);
    std::cout << ">> Entropy: " << ent << '\n';
}
