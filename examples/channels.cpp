// Quantum channels
// Source: ./examples/channels.cpp
#include <complex>
#include <iostream>
#include <vector>

#include "qpp.h"

int main() {
    using namespace qpp;
    idx nk = 5;
    idx D = 3; // nk Kraus on d-dimensional system
    std::cout << ">> Generating a random channel with " << nk
              << " Kraus operators on a " << D << " dimensional space\n";
    std::vector<cmat> Ks = randkraus(nk, D);

    cmat rho_in = randrho(D);         // random input state
    cmat rho_out = apply(rho_in, Ks); // output state

    std::cout << ">> Computing its Choi matrix...\n";
    cmat choim = kraus2choi(Ks);
    std::cout << ">> Choi matrix:\n" << disp(choim) << '\n';

    std::cout << ">> The eigenvalues of the Choi matrix are:\n"
              << disp(transpose(hevals(choim))) << '\n';

    std::cout << ">> Their sum is: " << sum(hevals(choim)) << '\n';

    std::vector<cmat> Kperps = choi2kraus(choim);
    std::cout << ">> The Kraus rank of the channel is: " << Kperps.size()
              << '\n';

    cmat rho_out1 = apply(rho_in, Kperps);
    // verification
    std::cout << ">> Norm difference on output states: "
              << norm(rho_out1 - rho_out) << '\n';

    std::cout << ">> Superoperator matrix:\n";
    cmat smat = kraus2super(Ks);
    std::cout << disp(smat) << '\n';

    std::cout << ">> The eigenvalues of the superoperator matrix are:\n";
    dyn_col_vect<cplx> evalsupop = evals(smat);
    std::cout << disp(transpose(evalsupop)) << '\n';

    std::cout << ">> Their absolute values are:\n";
    for (idx i = 0; i < (idx) evalsupop.size(); ++i)
        std::cout << std::abs(evalsupop(i)) << " ";

    // verification
    std::cout << "\n>> Norm difference for the superoperator action: ";
    cmat rho_out2 =
        transpose(reshape(smat * reshape(transpose(rho_in), D * D, 1), D, D));
    std::cout << norm(rho_out - rho_out2) << '\n';
}
