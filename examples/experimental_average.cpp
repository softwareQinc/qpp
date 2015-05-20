// Experimental average
// Source: ./examples/experimental_average.cpp
#include <qpp.h>
using namespace qpp;
using std::cout;
using std::endl;

int main()
{
    ket psi = mket({0}); // same as st.z0;
    cmat A = gt.X;
    auto evals = hevals(A);
    auto evects = hevects(A);

    long res = 0;
    idx N = 10000; // number of "measurement experiments"
    for (idx i = 0; i < N; ++i)
    {
        auto measurement = measure(psi, evects);
        idx m = std::get<0>(measurement); // measurement result
        if (evals[m] < 0) // -1
            res--;
        else              // +1
            res++;
    }
    cout << "N = " << N << " measurements" << std::endl;
    cout << "The experimental average of the observable" << endl;
    cout << disp(A) << endl;
    cout << "on the state" << endl;
    cout << disp(psi) << endl;
    cout << "is: " << res / static_cast<double>(N) << endl;
    cout << "Theoretical average <psi | A | psi> = ";
    cout << disp((adjoint(psi) * A * psi).value()) << endl;
}

