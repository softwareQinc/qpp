// Gram-Schmidt orthogonalization
// Source: ./examples/gram_schmidt.cpp
#include <qpp.h>
using namespace qpp;
using std::cout;
using std::endl;

int main()
{
    cmat A(3, 3);
    A << 1, 1, 0, 0, 2, 0, 0, 0, 0;
    cout << ">> Input matrix:" << endl << disp(A) << endl;

    cmat Ags = grams(A);
    cout << ">> Result:" << endl << disp(Ags) << endl;

    cout << ">> Projector onto G.S. vectors:" << endl;
    cout << disp(Ags * adjoint(Ags)) << endl;
}
