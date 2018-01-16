// Functor
// Source: ./examples/functor.cpp
#include <complex>
#include <iostream>

#include "qpp.h"

// test function used by qpp::cwise()
qpp::cplx pow3(const qpp::cplx& z) { return std::pow(z, 3); }

int main() {
    using namespace qpp;
    // functor test
    std::cout << ">> Functor z^3 acting component-wise on:\n";
    cmat A(2, 2);
    A << 1, 2, 3, 4;
    std::cout << disp(A) << '\n';

    std::cout << ">> Result (with lambda):\n";
    // functor z^3 componentwise, specify OutputScalar and Derived for lambdas
    std::cout << disp(cwise<cplx, cmat>(A, [](const cplx& z) -> cplx {
        return z * z * z;
    })) << '\n';

    std::cout << ">> Result (with genuine function):\n";
    // automatic type deduction for proper functions
    std::cout << disp(cwise(A, &pow3)) << '\n';
}
