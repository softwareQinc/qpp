// Bell inequalities (CHSH) violation
// Source: ./examples/bell_inequalities.cpp
#include <qpp.h>
using namespace qpp;
using std::cout;
using std::endl;

int main()
{
    ket psi = st.b00; // Measure Bell0 state
    idx N = 20; // number of measurements each party does

    cmat Q = gt.Z;
    cmat R = gt.X;
    cmat S = (-gt.Z - gt.X) / std::sqrt(2);
    cmat T = (gt.Z - gt.X) / std::sqrt(2);

    // TODO: complete the example
}