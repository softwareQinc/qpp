// Quantum error correcting codes
// Source: ./examples/qecc.cpp
#include <iostream>

#include "qpp.h"

int main() {
    using namespace qpp;
    ket a0 = Codes::codeword(Codes::Type::FIVE_QUBIT, 0);
    ket a1 = Codes::codeword(Codes::Type::FIVE_QUBIT, 1);

    ket b0 = Codes::codeword(Codes::Type::STEANE_SEVEN_QUBIT, 0);
    ket b1 = Codes::codeword(Codes::Type::STEANE_SEVEN_QUBIT, 1);

    ket c0 = Codes::codeword(Codes::Type::SHOR_NINE_QUBIT, 0);
    ket c1 = Codes::codeword(Codes::Type::SHOR_NINE_QUBIT, 1);

    std::cout << ">> [[5, 1, 3]] Five qubit code.\n";
    std::cout << ">> Checking codeword orthogonality.\n";
    std::cout << ">> |<0L|1L>| = ";
    std::cout << disp(adjoint(a0) * a1) << '\n';

    std::cout << ">> [[7, 1, 3]] Seven qubit Steane code.\n";
    std::cout << ">> Checking codeword orthogonality.\n";
    std::cout << ">> |<0L|1L>| = ";
    std::cout << disp(adjoint(b0) * b1) << '\n';

    std::cout << ">> [[9, 1, 3]] Nine qubit Shor code.\n";
    std::cout << ">> Checking codeword orthogonality.\n";
    std::cout << ">> |<0L|1L>| = ";
    std::cout << disp(adjoint(c0) * c1) << '\n';
}
