// Functor
// Source: ./examples/functor.cpp
#include <complex>
#include <iostream>

#include "qpp.h"

// test function used by qpp::cwise()
qpp::cplx pow3(qpp::cplx z) { return std::pow(z, 3); }

int main() {
    using namespace qpp;

    // functor test
    std::cout << ">> Functor z^3 acting component-wise on:\n";
    cmat A(2, 2);
    A << 1, 2, 3, 4;
    std::cout << disp(A) << '\n';

    std::cout << ">> Result (with lambda):\n";
    auto lambda = [](cplx z) -> cplx { return z * z * z; };
    // functor z^3 component-wise, you must specify both template arguments
    // OutputScalar and Derived (Eigen expression input type) for lambdas passed
    // to qpp::cwise<OutputScalar, Derived>()
    std::cout << disp(cwise<cplx, cmat>(A, lambda)) << '\n';

    std::cout << ">> Result (with genuine function):\n";
    // automatic type deduction for proper functions
    std::cout << disp(cwise(A, &pow3)) << '\n';
}
