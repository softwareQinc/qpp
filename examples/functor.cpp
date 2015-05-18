// Functor
// Source: ./examples/functor.cpp
#include <qpp.h>
using namespace qpp;
using std::cout;
using std::endl;

// test function used by qpp::cwise()
cplx pow3(const cplx& z)
{
    return std::pow(z, 3);
}

int main()
{
    // functor test
    cout << ">> Functor z^3 acting component-wise on:" << endl;
    cmat A(2, 2);
    A << 1, 2, 3, 4;
    cout << disp(A) << endl;

    cout << ">> Result (with lambda):" << endl;
    // functor z^3 componentwise, specify OutputScalar and Derived for lambdas
    cout << disp(cwise<cplx, cmat>(A, [](const cplx& z) -> cplx
    {
        return z * z * z;
    })) << endl;

    cout << ">> Result (with genuine function):" << endl;
    // automatic type deduction for proper functions
    cout << disp(cwise(A, &pow3)) << endl;
}
