// source: ./examples/ex2.cpp

#include <qpp.h>
using namespace qpp;

int main()
{
    ket psi = mket({1, 0});
    cmat U = gt.CNOTab;
    ket result = U * psi;

    std::cout << "The result of applying the Controlled-NOT gate CNOTab on |10> is:\n";
    std::cout << disp(result) << std::endl;

    ket phi = st.z0;
    U = gt.X;
    result = U * phi;
    
    std::cout << "The result of applying the bit-flip gate X on |0> is:\n";
    std::cout << disp(result) << std::endl;
}