// Quantum error correcting codes
// Source: ./examples/qecc.cpp
#include <iostream>
#include "qpp.h"

using namespace qpp;

int main()
{
    ket a0 = codes.codeword(Codes::Type::FIVE_QUBIT, 0);
    ket a1 = codes.codeword(Codes::Type::FIVE_QUBIT, 1);

    ket b0 = codes.codeword(Codes::Type::SEVEN_QUBIT_STEANE, 0);
    ket b1 = codes.codeword(Codes::Type::SEVEN_QUBIT_STEANE, 1);

    ket c0 = codes.codeword(Codes::Type::NINE_QUBIT_SHOR, 0);
    ket c1 = codes.codeword(Codes::Type::NINE_QUBIT_SHOR, 1);

    std::cout << ">> [[5, 1, 3]] Five qubit code. ";
    std::cout << ">> Checking codeword orthogonality." << std::endl;
    std::cout << ">> |<0L|1L>| = ";
    std::cout << disp(adjoint(a0) * a1) << std::endl;

    std::cout << ">> [[7, 1, 3]] Seven qubit Steane code. ";
    std::cout << ">> Checking codeword orthogonality." << std::endl;
    std::cout << ">> |<0L|1L>| = ";
    std::cout << disp(adjoint(b0) * b1) << std::endl;

    std::cout << ">> [[9, 1, 3]] Nine qubit Shor code. ";
    std::cout << ">> Checking codeword orthogonality." << std::endl;
    std::cout << ">> |<0L|1L>| = ";
    std::cout << disp(adjoint(c0) * c1) << std::endl;
}
