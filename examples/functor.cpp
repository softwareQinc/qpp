// Functor
// Source: ./examples/functor.cpp
#include <complex>
#include <iostream>
#include "qpp.h"

using namespace qpp;

// test function used by qpp::cwise()
cplx pow3(const cplx& z)
{
    return std::pow(z, 3);
}

int main()
{
    // functor test
    std::cout << ">> Functor z^3 acting component-wise on:" << std::endl;
    cmat A(2, 2);
    A << 1, 2, 3, 4;
    std::cout << disp(A) << std::endl;

    std::cout << ">> Result (with lambda):" << std::endl;
    // functor z^3 componentwise, specify OutputScalar and Derived for lambdas
    std::cout << disp(cwise<cplx, cmat>(A, [](const cplx& z) -> cplx
    {
        return z * z * z;
    })) << std::endl;

    std::cout << ">> Result (with genuine function):" << std::endl;
    // automatic type deduction for proper functions
    std::cout << disp(cwise(A, &pow3)) << std::endl;
}
