// Experimental average
// Source: ./examples/experimental_average.cpp
#include <iostream>
#include <tuple>

#include "qpp.h"

int main() {
    using namespace qpp;
    ket psi = 0_ket; // same as st.z0;
    cmat A = gt.X;
    dyn_col_vect<double> evals = hevals(A);
    cmat evects = hevects(A);

    long res = 0;
    idx N = 10000; // number of "measurement experiments"
    for (idx i = 0; i < N; ++i) {
        auto measured = measure(psi, evects);
        idx m = std::get<0>(measured); // measurement result
        if (evals[m] < 0)              // -1
            --res;
        else // +1
            ++res;
    }
    std::cout << ">> N = " << N << " measurements\n";
    std::cout << ">> The experimental average of the observable\n";
    std::cout << disp(A) << '\n';
    std::cout << "on the state\n";
    std::cout << disp(psi) << '\n';
    std::cout << "is: " << res / static_cast<double>(N) << '\n';
    std::cout << ">> Theoretical average <psi | A | psi> = ";
    std::cout << disp((adjoint(psi) * A * psi).value()) << '\n';
}
