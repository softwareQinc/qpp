// Gram-Schmidt orthogonalization
// Source: ./examples/gram_schmidt.cpp
#include <iostream>
#include "qpp.h"

using namespace qpp;

int main()
{
    cmat A(3, 3);
    A << 1, 1, 0, 0, 2, 0, 0, 0, 0;
    std::cout << ">> Input matrix:" << std::endl << disp(A) << std::endl;

    cmat Ags = grams(A);
    std::cout << ">> Result:" << std::endl << disp(Ags) << std::endl;

    std::cout << ">> Projector onto G.S. vectors:" << std::endl;
    std::cout << disp(Ags * adjoint(Ags)) << std::endl;
}
