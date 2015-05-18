// Spectral decomposition
// Source: ./examples/spectral.cpp
#include <qpp.h>
using namespace qpp;
using std::cout;
using std::endl;

int main()
{
    idx D = 4;
    cmat rH = randH(D);
    cout << ">> Original matrix: " << endl << disp(rH) << endl;

    // spectral decomposition here
    dyn_col_vect<double> evalsH = hevals(rH);
    cmat evectsH = hevects(rH);
    cmat spec = cmat::Zero(D, D);
    // reconstruct the matrix
    for (idx i = 0; i < D; ++i)
        spec += evalsH(i) * prj(evectsH.col(i));

    cout << ">> Reconstructed from spectral decomposition: " << endl;
    cout << disp(spec) << endl;

    // verification
    cout << ">> Norm difference: " << norm(spec - rH) << endl;
}
