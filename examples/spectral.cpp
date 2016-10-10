// Spectral decomposition
// Source: ./examples/spectral.cpp
#include <iostream>
#include "qpp.h"

using namespace qpp;

int main()
{
    idx D = 4;
    cmat rH = randH(D); // random Hermitian matrix
    std::cout << ">> Original matrix: " << std::endl << disp(rH) << std::endl;

    // spectral decomposition here
    dyn_col_vect<double> evalsH = hevals(rH);
    cmat evectsH = hevects(rH);
    cmat spec = cmat::Zero(D, D);
    // reconstruct the matrix
    for (idx i = 0; i < D; ++i)
        spec += evalsH(i) * prj(evectsH.col(i));

    std::cout << ">> Reconstructed from spectral decomposition: " << std::endl;
    std::cout << disp(spec) << std::endl;

    // verification
    std::cout << ">> Norm difference: " << norm(spec - rH) << std::endl;
}
