// Gram-Schmidt orthogonalization
// Source: ./examples/gram_schmidt.cpp
#include <iostream>

#include "qpp.h"

int main() {
    using namespace qpp;
    cmat A(3, 3);
    A << 1, 1, 0, 0, 2, 0, 0, 0, 0;
    std::cout << ">> Input matrix:\n" << disp(A) << '\n';

    cmat Ags = grams(A);
    std::cout << ">> Result:\n" << disp(Ags) << '\n';

    std::cout << ">> Projector onto G.S. vectors:\n";
    std::cout << disp(Ags * adjoint(Ags)) << '\n';
}
