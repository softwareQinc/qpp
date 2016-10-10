// Quantum channels
// Source: ./examples/channels.cpp
#include <complex>
#include <iostream>
#include <vector>
#include "qpp.h"

using namespace qpp;

int main()
{
    idx nk = 5;
    idx D = 3; // nk Kraus on d-dimensional system
    std::cout << ">> Generating a random channel with " << nk
              << " Kraus operators on a " << D << " dimensional space"
              << std::endl;
    std::vector<cmat> Ks = randkraus(nk, D);

    cmat rho_in = randrho(D); // random input state
    cmat rho_out = apply(rho_in, Ks); // output state

    std::cout << ">> Computing its Choi matrix..." << std::endl;
    cmat choim = kraus2choi(Ks);
    std::cout << ">> Choi matrix:" << std::endl << disp(choim) << std::endl;

    std::cout << ">> The eigenvalues of the Choi matrix are: "
              << std::endl << disp(transpose(hevals(choim))) << std::endl;

    std::cout << ">> Their sum is: " << sum(hevals(choim)) << std::endl;

    std::vector<cmat> Kperps = choi2kraus(choim);
    std::cout << ">> The Kraus rank of the channel is: "
              << Kperps.size() << std::endl;

    cmat rho_out1 = apply(rho_in, Kperps);
    // verification
    std::cout << ">> Norm difference on output states: "
              << norm(rho_out1 - rho_out) << std::endl;

    std::cout << ">> Superoperator matrix:" << std::endl;
    cmat smat = kraus2super(Ks);
    std::cout << disp(smat) << std::endl;

    std::cout << ">> The eigenvalues of the superoperator matrix are: "
              << std::endl;
    dyn_col_vect<cplx> evalsupop = evals(smat);
    std::cout << disp(transpose(evalsupop)) << std::endl;

    std::cout << ">> Their absolute values are: " << std::endl;
    for (idx i = 0; i < (idx) evalsupop.size(); ++i)
        std::cout << std::abs(evalsupop(i)) << " ";

    // verification
    std::cout << std::endl
              << ">> Norm difference for the superoperator action: ";
    cmat rho_out2 = transpose(
            reshape(smat * reshape(transpose(rho_in), D * D, 1), D, D));
    std::cout << norm(rho_out - rho_out2) << std::endl;
}
