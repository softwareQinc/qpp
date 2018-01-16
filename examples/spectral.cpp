// Spectral decomposition
// Source: ./examples/spectral.cpp
#include <iostream>

#include "qpp.h"

int main() {
    using namespace qpp;
    idx D = 4;
    cmat rH = randH(D); // random Hermitian matrix
    std::cout << ">> Original matrix:\n" << disp(rH) << '\n';

    // spectral decomposition here
    dyn_col_vect<double> evalsH = hevals(rH);
    cmat evectsH = hevects(rH);
    cmat spec = cmat::Zero(D, D);
    // reconstruct the matrix
    for (idx i = 0; i < D; ++i)
        spec += evalsH(i) * prj(evectsH.col(i));

    std::cout << ">> Reconstructed from spectral decomposition:\n";
    std::cout << disp(spec) << '\n';

    // verification
    std::cout << ">> Norm difference: " << norm(spec - rH) << '\n';
}
