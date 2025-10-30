// Source: ./examples/experimental_average.cpp
//
// Experimental average

#include <iostream>
#include <tuple>

#include "qpp/qpp.hpp"

int main() {
    using namespace qpp;

    ket psi = 0_ket; // same as st.z0;
    cmat X = gt.X;
    dyn_col_vect<realT> evals = hevals(X);
    cmat evects = hevects(X);

    long res = 0;
    idx N = 10000; // number of "measurement experiments"
    for (idx i = 0; i < N; ++i) {
        auto measured = measure(psi, evects);
        idx m = std::get<RES>(measured); // measurement result
        if (evals[m] < 0) {
            --res; // -1
        } else {
            ++res; // +1
        }
    }
    std::cout << ">> N = " << N << " measurements\n";
    std::cout << ">> The experimental average of the observable X\n";
    std::cout << disp(X) << '\n';
    std::cout << "on the state psi\n";
    std::cout << disp(psi) << '\n';
    std::cout << "is: " << static_cast<realT>(res) / static_cast<realT>(N)
              << '\n';
    std::cout << ">> Theoretical average <psi | X | psi> = ";
    std::cout << disp((adjoint(psi) * X * psi).value()) << '\n';
}
