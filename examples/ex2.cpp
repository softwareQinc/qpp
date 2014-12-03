// Gates and states
// Source: ./examples/ex2.cpp
#include <qpp.h>
using namespace qpp;

int main()
{
    ket psi = st.z0; // |0> state
    cmat U = gt.X;
    ket result = U * psi;

    std::cout << "The result of applying the bit-flip gate X on |0> is:\n";
    std::cout << disp(result) << std::endl;

    psi = mket({1, 0}); // |10> state
    U = gt.CNOTab; // Controlled-NOT
    result = U * psi;

    std::cout << "The result of applying the gate CNOTab on |10> is:\n";
    std::cout << disp(result) << std::endl;

    U = randU(2);
    std::cout << "Generating a random one-qubit gate U:\n";
    std::cout << disp(U) << std::endl;

    result = applyCTRL(psi, U, {0}, {1});
    std::cout << "The result of applying the Controlled-U gate on |10> is:\n";
    std::cout << disp(result) << std::endl;
}