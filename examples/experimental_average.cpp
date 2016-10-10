// Experimental average
// Source: ./examples/experimental_average.cpp
#include <iostream>
#include <tuple>
#include "qpp.h"

using namespace qpp;

int main()
{
    ket psi = mket({0}); // same as st.z0;
    cmat A = gt.X;
    dyn_col_vect<double> evals = hevals(A);
    cmat evects = hevects(A);

    long res = 0;
    idx N = 10000; // number of "measurement experiments"
    for (idx i = 0; i < N; ++i)
    {
        auto measured = measure(psi, evects);
        idx m = std::get<0>(measured); // measurement result
        if (evals[m] < 0) // -1
            res--;
        else              // +1
            res++;
    }
    std::cout << ">> N = " << N << " measurements" << std::endl;
    std::cout << ">> The experimental average of the observable" << std::endl;
    std::cout << disp(A) << std::endl;
    std::cout << "on the state" << std::endl;
    std::cout << disp(psi) << std::endl;
    std::cout << "is: " << res / static_cast<double>(N) << std::endl;
    std::cout << ">> Theoretical average <psi | A | psi> = ";
    std::cout << disp((adjoint(psi) * A * psi).value()) << std::endl;
}
